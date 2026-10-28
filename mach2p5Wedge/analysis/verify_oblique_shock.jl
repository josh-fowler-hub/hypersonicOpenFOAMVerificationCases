
using Printf
using Statistics
using LinearAlgebra
using DelimitedFiles

# ============================================================
# PHASE 2 — NASA MACH 2.5 OBLIQUE SHOCK VERIFICATION
# ============================================================

# ------------------------------------------------------------
# 1. PATHS AND SETTINGS
# ------------------------------------------------------------

const CASE_DIR = normpath(joinpath(@__DIR__, ".."))
const RESULTS_DIR = joinpath(@__DIR__, "results")

mkpath(RESULTS_DIR)

# Freestream and wedge parameters
const M1 = 2.5
const THETA_DEG = 15.0

# NASA reference freestream
const P1 = 101352.9       # Pa
const T1 = 288.8889       # K

# Match constant/thermophysicalProperties
const RU = 8314.462618    # J/(kmol K)
const MW = 28.970         # kg/kmol
const CP = 1004.5         # J/(kg K)

const R = RU / MW
const CV = CP - R
const GAMMA = CP / CV

# Geometry
const X_TIP = 0.0
const Y_TIP = 0.0
const WEDGE_ANGLE = deg2rad(THETA_DEG)

# Sampling ranges; x positions in meters
const SHOCK_X_MIN = 0.06
const SHOCK_X_MAX = 0.22

const PLATEAU_X_MIN = 0.08
const PLATEAU_X_MAX = 0.22

# Exclude regions near shocks and walls
const PLATEAU_BUFFER = 0.025

# Search half-width about predicted shock line
const SHOCK_SEARCH_WIDTH = 0.04


# ------------------------------------------------------------
# 2. EXACT OBLIQUE-SHOCK SOLUTION
# ------------------------------------------------------------

function theta_from_beta(beta, M, gamma)
    numerator = 2 * cot(beta) *
                (M^2 * sin(beta)^2 - 1)

    denominator =
        M^2 * (gamma + cos(2 * beta)) + 2

    return atan(numerator / denominator)
end


function solve_weak_shock_angle(M, theta, gamma)

    # Minimum allowed angle is the Mach angle.
    mu = asin(1 / M)

    # Search for the first root above the Mach angle.
    betas = range(mu + 1e-8, pi/2 - 1e-8;
                  length=2000)

    residual(beta) =
        theta_from_beta(beta, M, gamma) - theta

    for j in 1:(length(betas) - 1)

        a = betas[j]
        b = betas[j+1]

        fa = residual(a)
        fb = residual(b)

        if fa * fb <= 0

            # Bisection method
            for iteration in 1:100

                mid = (a + b) / 2
                fm = residual(mid)

                if abs(fm) < 1e-13
                    return mid
                end

                if fa * fm <= 0
                    b = mid
                else
                    a = mid
                    fa = fm
                end
            end

            return (a + b) / 2
        end
    end

    error("No attached weak shock found.")
end


function analytical_solution()

    theta = deg2rad(THETA_DEG)

    beta = solve_weak_shock_angle(
        M1, theta, GAMMA
    )

    # Normal Mach number upstream
    Mn1 = M1 * sin(beta)

    # Normal shock relations
    pressure_ratio =
        1 + (2 * GAMMA / (GAMMA + 1)) *
        (Mn1^2 - 1)

    density_ratio =
        ((GAMMA + 1) * Mn1^2) /
        ((GAMMA - 1) * Mn1^2 + 2)

    temperature_ratio =
        pressure_ratio / density_ratio

    Mn2_squared =
        (1 + 0.5 * (GAMMA - 1) * Mn1^2) /
        (GAMMA * Mn1^2 - 0.5 * (GAMMA - 1))

    Mn2 = sqrt(Mn2_squared)

    # Downstream Mach number
    M2 = Mn2 / sin(beta - theta)

    # Upstream thermodynamic state
    rho1 = P1 / (R * T1)
    a1 = sqrt(GAMMA * R * T1)
    U1 = M1 * a1

    # Downstream thermodynamic state
    p2 = P1 * pressure_ratio
    T2 = T1 * temperature_ratio
    rho2 = rho1 * density_ratio

    # Downstream velocity magnitude
    a2 = sqrt(GAMMA * R * T2)
    U2 = M2 * a2

    # Cartesian velocity components
    u2 = U2 * cos(theta)
    v2 = U2 * sin(theta)

    state1 = (
        p=P1,
        T=T1,
        rho=rho1,
        u=U1,
        v=0.0,
        M=M1
    )

    state2 = (
        p=p2,
        T=T2,
        rho=rho2,
        u=u2,
        v=v2,
        M=M2
    )

    return (
        beta=beta,
        theta=theta,
        Mn1=Mn1,
        Mn2=Mn2,
        pressure_ratio=pressure_ratio,
        density_ratio=density_ratio,
        temperature_ratio=temperature_ratio,
        state1=state1,
        state2=state2
    )
end


# ------------------------------------------------------------
# 3. OPENFOAM FIELD READER
# ------------------------------------------------------------

# This parser supports ASCII fields with:
#   internalField uniform ...
#   internalField nonuniform List<scalar> ...
#   internalField nonuniform List<vector> ...
#
# It does not support binary output or decomposed processor data.

function read_foam_internal(path; vector=false)

    text = read(path, String)

    # Strip comments
    text = replace(text, r"/\*.*?\*/"s => "")
    text = replace(text, r"//[^\n]*" => "")

    uniform_match =
        match(r"internalField\s+uniform\s+([^;]+);", text)

    if uniform_match !== nothing

        raw = strip(uniform_match.captures[1])

        if vector
            raw = replace(raw, '(' => ' ', ')' => ' ')
            vals = parse.(Float64, split(raw))
            length(vals) == 3 ||
                error("Invalid vector in $path")
            return vals
        else
            return parse(Float64, raw)
        end
    end

    pattern =
        r"internalField\s+nonuniform\s+List<\w+>\s+(\d+)\s*\((.*?)\)\s*;"s

    m = match(pattern, text)

    m === nothing &&
        error("Cannot parse internalField in $path")

    n = parse(Int, m.captures[1])

    payload = replace(
        m.captures[2],
        '(' => ' ',
        ')' => ' '
    )

    values = parse.(Float64, split(payload))

    if vector
        length(values) == 3n ||
            error("Vector field length mismatch: $path")
        return permutedims(reshape(values, 3, n))
    else
        length(values) == n ||
            error("Scalar field length mismatch: $path")
        return values
    end
end


function expand_field(data, n; vector=false)

    if vector
        if data isa AbstractVector
            return repeat(reshape(data, 1, 3), n, 1)
        end

        size(data, 1) == n ||
            error("Vector cell count mismatch")
        return data
    end

    if data isa Number
        return fill(Float64(data), n)
    end

    length(data) == n ||
        error("Scalar cell count mismatch")

    return data
end


function latest_time_directory(case_dir)

    options = Tuple{Float64,String}[]

    for name in readdir(case_dir)

        t = tryparse(Float64, name)

        if t !== nothing &&
           isdir(joinpath(case_dir, name))
            push!(options, (t, name))
        end
    end

    isempty(options) &&
        error("No OpenFOAM time directories found")

    sort!(options, by=first)

    return joinpath(case_dir, last(options[end]))
end


function load_openfoam_solution(case_dir)

    time_dir = latest_time_directory(case_dir)

    println("Reading: $time_dir")

    xraw = read_foam_internal(joinpath(time_dir, "Cx"))
    yraw = read_foam_internal(joinpath(time_dir, "Cy"))

    xraw isa AbstractVector ||
        error("Cx must be a nonuniform cell field")

    n = length(xraw)

    x = xraw
    y = expand_field(yraw, n)

    p = expand_field(
        read_foam_internal(joinpath(time_dir, "p")), n
    )

    T = expand_field(
        read_foam_internal(joinpath(time_dir, "T")), n
    )

    U = expand_field(
        read_foam_internal(
            joinpath(time_dir, "U"); vector=true
        ), n; vector=true
    )

    u = U[:, 1]
    v = U[:, 2]

    rho = p ./ (R .* T)
    speed = sqrt.(u.^2 .+ v.^2)
    sound = sqrt.(GAMMA .* R .* T)
    mach = speed ./ sound

    all(isfinite, p) || error("Non-finite pressure")
    all(T .> 0) || error("Nonpositive temperature")
    all(p .> 0) || error("Nonpositive pressure")

    return (
        x=x, y=y,
        p=p, T=T, rho=rho,
        u=u, v=v, M=mach,
        n=n
    )
end


# ------------------------------------------------------------
# 4. NUMERICAL SHOCK DETECTION
# ------------------------------------------------------------


function detect_shock(sol, ref)

    columns = Dict{Float64, Vector{Int}}()

    for i in eachindex(sol.x)
        xkey = round(sol.x[i]; digits=9)
        push!(get!(columns, xkey, Int[]), i)
    end

    shock_x = Float64[]
    shock_y = Float64[]

    p1 = ref.state1.p
    p2 = ref.state2.p

    p_mid = 0.5 * (p1 + p2)

    for x in sort(collect(keys(columns)))

        if !(SHOCK_X_MIN <= x <= SHOCK_X_MAX)
            continue
        end

        ids = sort(
            columns[x],
            by=i -> sol.y[i]
        )

        ys = sol.y[ids]
        ps = sol.p[ids]

        y_expected = x * tan(ref.beta)

        crossings = Float64[]

        for j in 1:(length(ps)-1)

            # Search for the expected decreasing
            # pressure transition with increasing y.
            if ps[j] >= p_mid &&
               ps[j+1] < p_mid

                fraction =
                    (p_mid - ps[j]) /
                    (ps[j+1] - ps[j])

                y_cross =
                    ys[j] +
                    fraction * (ys[j+1] - ys[j])

                push!(crossings, y_cross)
            end
        end

        isempty(crossings) && continue

        # If several crossings exist because of
        # oscillations, select the one nearest
        # the analytical shock location.
        distances = abs.(crossings .- y_expected)
        y_detected = crossings[argmin(distances)]

        # Reject detections outside the search band.
        if abs(y_detected - y_expected) >
           SHOCK_SEARCH_WIDTH
            continue
        end

        push!(shock_x, x)
        push!(shock_y, y_detected)
    end

    length(shock_x) >= 3 ||
        error("Insufficient shock detections")

    # Least-squares fit y = mx + b.
    A = hcat(shock_x, ones(length(shock_x)))

    coefficients = A \ shock_y

    slope = coefficients[1]
    intercept = coefficients[2]

    beta_num = atan(slope)

    fitted = slope .* shock_x .+ intercept
    residuals = shock_y .- fitted

    rms = sqrt(mean(residuals.^2))

    return (
        x=shock_x,
        y=shock_y,
        gradient=Float64[],
        slope=slope,
        intercept=intercept,
        beta=beta_num,
        residual_rms=rms
    )
end


# ------------------------------------------------------------
# 5. REGION SAMPLING
# ------------------------------------------------------------

function sample_regions(sol, ref)

    beta = ref.beta
    theta = ref.theta

    # Region 1: undisturbed upstream flow.
    upstream_mask = (
        (sol.x .< -0.05) .&
        (sol.y .> 0.03) .&
        (sol.y .< 0.27)
    )

    # Region 2: post-shock region, away from
    # wedge wall, shock, tip, and outlet.
    y_wedge = max.(sol.x, 0.0) .* tan(theta)
    y_shock = max.(sol.x, 0.0) .* tan(beta)

    downstream_mask = (
        (sol.x .> PLATEAU_X_MIN) .&
        (sol.x .< PLATEAU_X_MAX) .&
        (sol.y .> y_wedge .+ PLATEAU_BUFFER) .&
        (sol.y .< y_shock .- PLATEAU_BUFFER)
    )

    count(upstream_mask) > 20 ||
        error("Too few upstream cells")

    count(downstream_mask) > 20 ||
        error("Too few downstream cells")

    fields = (:p, :T, :rho, :u, :v, :M)

    function region_mean(mask)
        return NamedTuple{
            fields
        }(
            Tuple(
                mean(getproperty(sol, f)[mask])
                for f in fields
            )
        )
    end

    return (
        upstream=region_mean(upstream_mask),
        downstream=region_mean(downstream_mask),
        n_upstream=count(upstream_mask),
        n_downstream=count(downstream_mask)
    )
end


# ------------------------------------------------------------
# 6. VERIFICATION METRICS
# ------------------------------------------------------------

relative_error(qnum, qexact) =
    abs(qnum - qexact) / abs(qexact)

function print_comparison(name, numerical, exact)

    error_percent =
        100 * relative_error(numerical, exact)

    @printf(
        "%-24s %15.7f %15.7f %12.5f %%\n",
        name, numerical, exact, error_percent
    )
end


function verify_solution(sol, ref, shock, samples)

    println()
    println("="^75)
    println("NASA OBLIQUE SHOCK VERIFICATION")
    println("="^75)

    @printf("Total cells: %d\n", sol.n)
    @printf("Specific gas constant: %.8f J/(kg K)\n", R)
    @printf("Heat capacity ratio: %.10f\n", GAMMA)

    println("\nSHOCK GEOMETRY")
    println("-"^75)

    @printf(
        "Analytical angle: %.8f degrees\n",
        rad2deg(ref.beta)
    )

    @printf(
        "Numerical angle:  %.8f degrees\n",
        rad2deg(shock.beta)
    )

    @printf(
        "Angle error:      %.8f degrees\n",
        abs(rad2deg(shock.beta-ref.beta))
    )

    @printf(
        "Fitted intercept: %.8e m\n",
        shock.intercept
    )

    @printf(
        "Shock fit RMS:    %.8e m\n",
        shock.residual_rms
    )

    println("\nPOST-SHOCK STATE")
    println("-"^75)

    println(
        rpad("Quantity", 24),
        lpad("OpenFOAM", 15),
        lpad("Analytical", 15),
        lpad("Error", 15)
    )

    for field in (:p, :rho, :T, :u, :M)

        qnum = getproperty(samples.downstream, field)
        qref = getproperty(ref.state2, field)

        print_comparison(
            string(field), qnum, qref
        )
    end

    println("\nJUMP RATIOS")
    println("-"^75)

    p_ratio =
        samples.downstream.p / samples.upstream.p

    rho_ratio =
        samples.downstream.rho / samples.upstream.rho

    T_ratio =
        samples.downstream.T / samples.upstream.T

    print_comparison(
        "p2/p1", p_ratio, ref.pressure_ratio
    )

    print_comparison(
        "rho2/rho1", rho_ratio, ref.density_ratio
    )

    print_comparison(
        "T2/T1", T_ratio, ref.temperature_ratio
    )

    println("\nSAMPLE COUNTS")
    println("-"^75)

    @printf(
        "Upstream cells: %d\n",
        samples.n_upstream
    )

    @printf(
        "Downstream cells: %d\n",
        samples.n_downstream
    )
end


# ------------------------------------------------------------
# 7. EXPORT NUMERICAL DATA
# ------------------------------------------------------------

function write_csv(path, header, data)

    open(path, "w") do io
        println(io, join(header, ","))
        writedlm(io, data, ',')
    end
end


function export_results(sol, ref, shock, samples)

    # Full cell-center solution
    cell_data = hcat(
        sol.x, sol.y,
        sol.p, sol.T, sol.rho,
        sol.u, sol.v, sol.M
    )

    write_csv(
        joinpath(RESULTS_DIR, "flowfield.csv"),
        ["x", "y", "p", "T", "rho", "u", "v", "M"],
        cell_data
    )

    # Detected shock positions
    exact_y = shock.x .* tan(ref.beta)

    shock_data = hcat(
        shock.x,
        shock.y,
        exact_y,
        shock.y .- exact_y
    )

    write_csv(
        joinpath(RESULTS_DIR, "shock_location.csv"),
        ["x", "y_numerical", "y_exact", "difference"],
        shock_data
    )

    # Small summary file
    summary_data = [
        rad2deg(ref.beta) rad2deg(shock.beta)
    ]

    write_csv(
        joinpath(RESULTS_DIR, "shock_angle.csv"),
        ["beta_exact_deg", "beta_numerical_deg"],
        summary_data
    )
end


# ------------------------------------------------------------
# 8. MAIN PROGRAM
# ------------------------------------------------------------

function main()

    println("Starting oblique-shock verification...")

    ref = analytical_solution()

    sol = load_openfoam_solution(CASE_DIR)

    shock = detect_shock(sol, ref)

    samples = sample_regions(sol, ref)

    verify_solution(sol, ref, shock, samples)

    export_results(sol, ref, shock, samples)

    println()
    println("Verification complete.")
    println("Results saved to: $RESULTS_DIR")

    return (
        ref=ref,
        sol=sol,
        shock=shock,
        samples=samples
    )
end


results = main()
