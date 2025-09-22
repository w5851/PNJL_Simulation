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
mu = 0.0  # use scalar mu
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


# ------------------ 新增测试：针对 Advanced_ForwardDiff.jl 中新增函数 ------------------
@testset "Advanced ForwardDiff additional derivatives" begin
	# 测试 dP_dmu_B: 与中心差分比较（使用 pressure_solve_core 作为目标函数）
	d_mu = dP_dmu_B(x, mu, T, nodes1, omega)
	g = mu_val -> pressure_solve_core(x, mu_val, T, nodes1, omega)
	δμ = 1e-6
	fd_mu = (g(mu + δμ) - g(mu - δμ)) / (2δμ)

	@test isfinite(d_mu)
	@test isfinite(fd_mu)
	@test abs(d_mu - fd_mu) / (abs(fd_mu) + 1e-12) < 1e-2

	println("dP_dmu_B forwarddiff = ", d_mu)
	println("dP_dmu_B finite diff = ", fd_mu)

	# 测试 dP_domega: 与中心差分比较（使用 pressure_solve_core 作为目标函数）
	d_w = dP_domega(x, mu, T, nodes1, omega)
	h = w -> pressure_solve_core(x, mu, T, nodes1, w)
	δω = 1e-6
	fd_w = (h(omega + δω) - h(omega - δω)) / (2δω)

	@test isfinite(d_w)
	@test isfinite(fd_w)
	@test abs(d_w - fd_w) / (abs(fd_w) + 1e-12) < 1e-2

	println("dP_domega forwarddiff = ", d_w)
	println("dP_domega finite diff = ", fd_w)

	# 测试 generalized_susceptibility_mu 与直接差分关系（n=0 应等于 P*T^-4）
	P0 = pressure_solve_core(x, mu, T, nodes1, omega)
	gs0 = generalized_susceptibility_mu(0, x, mu, T, nodes1, omega)
	@test isfinite(gs0)
	@test abs(gs0 - P0 * T^(-4)) < 1e-12

	# n=1 与 dP_dmu_B * T^(1-4)
	gs1 = generalized_susceptibility_mu(1, x, mu, T, nodes1, omega)
	@test isfinite(gs1)
	@test abs(gs1 - d_mu * T^(1-4)) / (abs(d_mu * T^(1-4)) + 1e-12) < 1e-2

	# 测试 generalized_susceptibility_omega 的 n=0/1 行为
	gs0w = generalized_susceptibility_omega(0, x, mu, T, nodes1, omega)
	@test isfinite(gs0w)
	@test abs(gs0w - P0 * T^(-4)) < 1e-12

	gs1w = generalized_susceptibility_omega(1, x, mu, T, nodes1, omega)
	@test isfinite(gs1w)
	@test abs(gs1w - d_w * T^(1-4)) / (abs(d_w * T^(1-4)) + 1e-12) < 1e-2
end

