
# 高级 ForwardDiff 工具函数（不作为模块封装）
using ForwardDiff
include("Function_Rotation.jl")
using .Function_Rotation: pressure_solve_core

function dP_dT_rotation(x, mu_B, T, nodes1, omega)
	f = T -> pressure_solve_core(x, mu_B, T, nodes1, omega)
	return ForwardDiff.derivative(f, T)
end

function dP_dT2_rotation(x, mu_B, T, nodes1, omega)
	f = T -> dP_dT_rotation(x, mu_B, T, nodes1, omega)
	return ForwardDiff.derivative(f, T)
end

function dP_dT3_rotation(x, mu_B, T, nodes1, omega)
	f = T -> dP_dT2_rotation(x, mu_B, T, nodes1, omega)
	return ForwardDiff.derivative(f, T)
end

function dP_dT4_rotation(x, mu_B, T, nodes1, omega)
	f = T -> dP_dT3_rotation(x, mu_B, T, nodes1, omega)
	return ForwardDiff.derivative(f, T)
end
