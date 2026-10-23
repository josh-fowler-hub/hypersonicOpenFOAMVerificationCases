const γ = 1.4
const Cp = 1004.5
const W = 28.970

const Ru = 8314.462618  # J/(kmol K)
const R = Ru / W

p1 = 100_000.0
T1 = 300.0
u1 = 0.0
Ms = 2.0

ρ1 = p1 / (R * T1)
a1 = sqrt(γ * R * T1)

Vs = Ms * a1

pressure_ratio =
    1 + 2γ / (γ + 1) * (Ms^2 - 1)

density_ratio =
    (γ + 1) * Ms^2 /
    ((γ - 1) * Ms^2 + 2)

p2 = pressure_ratio * p1
ρ2 = density_ratio * ρ1

T2 = p2 / (ρ2 * R)

u2 = Vs * (1 - ρ1 / ρ2)

t = 4e-4
x0 = 0.2
xs = x0 + Vs * t

println("p2 = $p2 Pa")
println("rho2 = $ρ2 kg/m³")
println("T2 = $T2 K")
println("u2 = $u2 m/s")
println("Shock position = $xs m")