
using Printf
using Statistics
using DelimitedFiles
using Plots

# -------------------------------------------------------
# User configuration
# -------------------------------------------------------

const CASE_DIR = get(ENV, "CASE_DIR",
                     "/home/jfowler/scratch/sim/normalShockVerification")
const OUTPUT_DIR = get(ENV, "OUTPUT_DIR",
                       "/home/jfowler/scratch/sim/normalShockVerification/verificationResults")
const X0 = 0.2
const DOMAIN_LENGTH = 1.0

# Consistent with thermophysicalProperties
const RU = 8314.462618  # J/(kmol K)
const MW = 28.970        # kg/kmol
const CP = 1004.5        # J/(kg K)

const R = RU / MW
const CV = CP - R
const GAMMA = CP / CV

const P1 = 100_000.0
const T1 = 300.0
const U1 = 0.0
const MS = 2.0

# -------------------------------------------------------
# Analytical Rankine-Hugoniot solution
# -------------------------------------------------------

function exact_states(time)
    rho1 = P1 / (R * T1)
    a1 = sqrt(GAMMA * R * T1)
    vs = MS * a1

    pr = 1 + 2GAMMA / (GAMMA + 1) *
             (MS^2 - 1)

    rr = (GAMMA + 1) * MS^2 /
         ((GAMMA - 1) * MS^2 + 2)

    p2 = P1 * pr
    rho2 = rho1 * rr
    T2 = p2 / (rho2 * R)
    u2 = vs * (1 - 1 / rr)

    xs = X0 + vs * time

    state1 = (p=P1, T=T1, u=U1, rho=rho1)
    state2 = (p=p2, T=T2, u=u2, rho=rho2)

    return (; state1, state2, vs, xs)
end

# -------------------------------------------------------
# OpenFOAM ASCII field reader
# Supports uniform and nonuniform internal fields.
# -------------------------------------------------------

function read_foam_field(path, n_cells; vector=false)
    s = read(path, String)

    # Remove comments
    s = replace(s, r"/\*.*?\*/"s => "")
    s = replace(s, r"//[^\n]*" => "")

    pattern = r"internalField\s+(.*?);"s
    m = match(pattern, s)
    isnothing(m) && error("No internalField: $path")

    block = strip(m.captures[1])
    number = raw"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
    values = parse.(Float64,
        [x.match for x in eachmatch(Regex(number), block)])

    if startswith(block, "uniform")
        if vector
            length(values) == 3 ||
                error("Invalid uniform vector field")
            return repeat(reshape(values, 1, 3), n_cells, 1)
        end

        return fill(only(values), n_cells)
    end

    startswith(block, "nonuniform") ||
        error("Unsupported field format: $path")

    # First numeric value is the declared list length.
    count = Int(values[1])
    count == n_cells ||
        error("Cell count mismatch: $path")

    components = vector ? 3 : 1
    data = values[2:end]

    length(data) == components * n_cells ||
        error("Unexpected data length: $path")

    if vector
        return permutedims(reshape(data, 3, :))
    end

    return data
end

function time_directory(case_dir)
    options = filter(readdir(case_dir)) do name
        path = joinpath(case_dir, name)
        time = tryparse(Float64, name)
        isdir(path) && !isnothing(time) && isfinite(time) && time > 0
    end

    isempty(options) && error("No positive solution time directories found")

    times = parse.(Float64, options)
    # Match postProcess -latestTime and retain the actual solver output time.
    i = argmax(times)
    return (; dir=joinpath(case_dir, options[i]), time=times[i])
end

function load_solution(case_dir)
    selected = time_directory(case_dir)
    dir = selected.dir

    xfile = joinpath(dir, "Cx")
    isfile(xfile) ||
        error("Run: postProcess -func writeCellCentres -latestTime")

    # Determine cell count from nonuniform Cx list.
    s = read(xfile, String)
    m = match(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)"s, s)
    isnothing(m) && error("Expected nonuniform Cx field")

    n = parse(Int, m.captures[1])

    x = read_foam_field(xfile, n)
    p = read_foam_field(joinpath(dir, "p"), n)
    T = read_foam_field(joinpath(dir, "T"), n)
    U = read_foam_field(joinpath(dir, "U"), n; vector=true)

    idx = sortperm(x)

    x = x[idx]
    p = p[idx]
    T = T[idx]
    u = U[idx, 1]

    rho = p ./ (R .* T)

    return (; x, p, T, u, rho, time=selected.time)
end


# -------------------------------------------------------
# Numerical shock location
# -------------------------------------------------------

function crossing(x, q, target)
    for i in 1:(length(x)-1)
        qL = q[i] - target
        qR = q[i+1] - target

        if qL * qR <= 0 && q[i] != q[i+1]
            fraction = (target - q[i]) /
                       (q[i+1] - q[i])

            return x[i] + fraction *
                   (x[i+1] - x[i])
        end
    end

    error("Target crossing not found")
end

# -------------------------------------------------------
# Exact step profile
# -------------------------------------------------------

function exact_profile(x, left, right, xs)
    return ifelse.(x .< xs, left, right)
end

# -------------------------------------------------------
# Error norms
# For uniform 1D cells of width dx.
# -------------------------------------------------------

function error_norms(numerical, exact, dx)
    e = numerical .- exact

    L1 = sum(abs.(e)) * dx / DOMAIN_LENGTH
    L2 = sqrt(sum(abs2, e) * dx / DOMAIN_LENGTH)
    Linf = maximum(abs.(e))

    return (; L1, L2, Linf)
end

# -------------------------------------------------------
# Main verification
# -------------------------------------------------------

function verify()
    mkpath(OUTPUT_DIR)
    sol = load_solution(CASE_DIR)
    ref = exact_states(sol.time)

    x = sol.x
    dx = DOMAIN_LENGTH / length(x)

    fields = (:p, :rho, :T, :u)

    println("\nNORMAL SHOCK VERIFICATION")
    println("="^70)

    @printf("Cells:               %d\n", length(x))
    @printf("Time:                %.12e s\n", sol.time)
    @printf("Gamma:               %.10f\n", GAMMA)
    @printf("Exact shock speed:   %.6f m/s\n", ref.vs)

    p_mid = (ref.state1.p + ref.state2.p) / 2
    xs_num = crossing(x, sol.p, p_mid)

    @printf("\nExact shock location: %.9f m\n", ref.xs)
    @printf("CFD shock location:   %.9f m\n", xs_num)
    @printf("Absolute error:       %.9e m\n",
            abs(xs_num - ref.xs))

    @printf("Shock error / dx:     %.6f cells\n",
            abs(xs_num - ref.xs) / dx)

    # Exclude cells near shock for plateau-state tests.
    buffer = 0.03

    mask2 = x .< (ref.xs - buffer)
    mask1 = x .> (ref.xs + buffer)

    println("\nFIELD ERRORS")
    println("-"^70)

    for field in fields
        numerical = getproperty(sol, field)
        q1 = getproperty(ref.state1, field)
        q2 = getproperty(ref.state2, field)

        exact = exact_profile(x, q2, q1, ref.xs)
        norms = error_norms(numerical, exact, dx)

        err1 = mean(numerical[mask1]) - q1
        err2 = mean(numerical[mask2]) - q2

        @printf("\nField: %s\n", string(field))
        @printf("  L1:             %.8e\n", norms.L1)
        @printf("  L2:             %.8e\n", norms.L2)
        @printf("  Linf:           %.8e\n", norms.Linf)
        @printf("  State 1 bias:   %.8e\n", err1)
        @printf("  State 2 bias:   %.8e\n", err2)
    end

    # Pressure transition width: 90% to 10%
    jump = ref.state2.p - ref.state1.p

    x90 = crossing(x, sol.p,
                   ref.state1.p + 0.9jump)
    x10 = crossing(x, sol.p,
                   ref.state1.p + 0.1jump)

    thickness = abs(x10 - x90)

    @printf("\n10-90%% shock thickness: %.8e m\n",
            thickness)
    @printf("Thickness in cells:   %.3f\n",
            thickness / dx)

    # Export numerical and analytical profiles.
    exact_p = exact_profile(x, ref.state2.p,
                              ref.state1.p, ref.xs)
    exact_rho = exact_profile(x, ref.state2.rho,
                                ref.state1.rho, ref.xs)
    exact_T = exact_profile(x, ref.state2.T,
                              ref.state1.T, ref.xs)
    exact_u = exact_profile(x, ref.state2.u,
                              ref.state1.u, ref.xs)

    output = hcat(
        x, sol.p, exact_p,
        sol.rho, exact_rho,
        sol.T, exact_T,
        sol.u, exact_u
    )

    output_file = "$OUTPUT_DIR/verification_profiles.csv"

    open(output_file, "w") do io
        println(io,
            "x,p_cfd,p_exact,rho_cfd,rho_exact," *
            "T_cfd,T_exact,u_cfd,u_exact")

        writedlm(io, output, ',')
    end

    println("\nSaved: $output_file")

    return (; ref, sol, xs_num, thickness)
end

results = verify()

ref = results.ref
sol = results.sol

fields = [
    (:p, "Pressure", "Pa"),
    (:rho, "Density", "kg/m³"),
    (:T, "Temperature", "K"),
    (:u, "Velocity", "m/s")
]

figures = []

for (field, name, units) in fields
    numerical = getproperty(sol, field)

    q1 = getproperty(ref.state1, field)
    q2 = getproperty(ref.state2, field)

    exact = exact_profile(sol.x, q2, q1, ref.xs)

    plt = plot(
        sol.x, numerical;
        label="OpenFOAM",
        linewidth=2,
        xlabel="x (m)",
        ylabel="$name ($units)",
        title="$name — Normal Shock Verification"
    )

    plot!(
        plt, sol.x, exact;
        label="Exact",
        linestyle=:dash,
        linewidth=2
    )

    vline!(
        plt, [ref.xs];
        label="Exact shock position",
        linestyle=:dot
    )

    push!(figures, plt)
end

for (i, plt) in enumerate(figures)
    display(plt)
    savefig(plt, "$OUTPUT_DIR/verification_$(fields[i][1]).png")
end
