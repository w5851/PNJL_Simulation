import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Test

# Load source files
include(joinpath(@__DIR__, "..", "..", "src", "Rotation", "Function_Rotation.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "Rotation", "Advanced_ForwardDiff.jl"))

using .Function_Rotation: get_nodes, pressure_solve_core, hc


# Prepare a small smoke test for dP_dT_rotation
nodes1, nodes2 = get_nodes(8, 4)
omega = 50 / hc

# sample x: phi, Phi1, Phi2 and mu
x = [-2.13, 0.06, 0.12]
mu = [0.0]  # pressure_solve_core expects mu as array-like in some paths
T = 50 / hc

# forwarddiff derivative
d1 = dP_dT_rotation(x, mu, T, nodes1, omega)

# central finite difference
δ = 1e-6
f = T -> pressure_solve_core(x, mu, T, nodes1, omega)
fd = (f(T+δ) - f(T-δ)) / (2δ)

@test isfinite(d1)
@test isfinite(fd)
@test abs(d1 - fd) / (abs(fd) + 1e-12) < 1e-2 # within 1% relative

println("dP_dT_rotation forwarddiff = ", d1)
println("dP_dT_rotation finite diff = ", fd)
