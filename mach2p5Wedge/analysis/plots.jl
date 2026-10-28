
using Plots
using Printf

# Run and retrieve the verification results.
include("verify_oblique_shock.jl")

ref = results.ref
sol = results.sol
shock = results.shock

# ------------------------------------------------------------
# 1. Pressure field and shock geometry
# ------------------------------------------------------------

p1 = scatter(
    sol.x, sol.y;
    marker_z=sol.p ./ ref.state1.p,
    markersize=1,
    markershape=:rect,
    markerstrokewidth=0,
    markeralpha=1,
    colorbar_title="p/p₁",
    label="",
    xlabel="x (m)",
    ylabel="y (m)",
    title="Numerical pressure field",
    aspect_ratio=:equal,
    size=(1000, 650)
)

xline = collect(range(0.0, 0.30; length=200))

# Wedge geometry
plot!(
    p1,
    xline,
    xline .* tan(ref.theta);
    color=:black,
    linewidth=2,
    label="Wedge"
)

# Exact shock
plot!(
    p1,
    xline,
    xline .* tan(ref.beta);
    color=:red,
    linestyle=:dash,
    linewidth=2,
    label="Exact shock"
)

# Detected shock points
scatter!(
    p1,
    shock.x,
    shock.y;
    markersize=3,
    markercolor=:white,
    markerstrokecolor=:black,
    label="Detected shock"
)

display(p1)

savefig(
    p1,
    joinpath(@__DIR__, "results", "pressure_field.png")
)


# ------------------------------------------------------------
# 2. Shock geometry comparison
# ------------------------------------------------------------

p2 = scatter(
    shock.x,
    shock.y;
    markersize=4,
    label="Detected shock",
    xlabel="x (m)",
    ylabel="y (m)",
    title="Numerical versus analytical shock",
    size=(850, 550)
)

plot!(
    p2,
    shock.x,
    shock.x .* tan(ref.beta);
    linewidth=2,
    linestyle=:dash,
    label="Analytical shock"
)

plot!(
    p2,
    shock.x,
    shock.slope .* shock.x .+ shock.intercept;
    linewidth=2,
    label="Least-squares fit"
)

display(p2)

savefig(
    p2,
    joinpath(@__DIR__, "results", "shock_fit.png")
)


# ------------------------------------------------------------
# 3. Shock-position residual
# ------------------------------------------------------------

residual = shock.y .-
           shock.x .* tan(ref.beta)

p3 = scatter(
    shock.x,
    residual .* 1000;
    xlabel="x (m)",
    ylabel="Position difference (mm)",
    title="Shock-location error",
    label="Numerical - analytical",
    markersize=4,
    size=(850, 550)
)

hline!(
    p3,
    [0.0];
    linestyle=:dash,
    label="Exact"
)

display(p3)

savefig(
    p3,
    joinpath(@__DIR__, "results", "shock_error.png")
)


# ------------------------------------------------------------
# 4. Pressure profiles along selected columns
# ------------------------------------------------------------

p4 = plot(
    xlabel="y (m)",
    ylabel="p/p₁",
    title="Vertical pressure profiles",
    size=(850, 550)
)

sampling_x = [0.08, 0.12, 0.16, 0.20]

for xtarget in sampling_x

    # Find nearest available mesh column.
    distances = abs.(sol.x .- xtarget)
    xactual = sol.x[argmin(distances)]

    ids = findall(
        abs.(sol.x .- xactual) .< 1e-8
    )

    ids = sort(ids, by=i -> sol.y[i])

    plot!(
        p4,
        sol.y[ids],
        sol.p[ids] ./ ref.state1.p;
        linewidth=2,
        label=@sprintf("x = %.3f m", xactual)
    )
end

hline!(
    p4,
    [ref.pressure_ratio];
    linestyle=:dash,
    color=:black,
    label="Exact downstream ratio"
)

hline!(
    p4,
    [1.0];
    linestyle=:dot,
    color=:gray,
    label="Exact upstream ratio"
)

display(p4)

savefig(
    p4,
    joinpath(@__DIR__, "results", "pressure_profiles.png")
)
