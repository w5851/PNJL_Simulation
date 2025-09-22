
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

# --------------------------- 新增：广义磁化率 -----------------
# P 关于 mu_B 的偏导数（一次）
function dP_dmu_B(x, mu_B, T, nodes1, omega)
    f_mu = mu -> pressure_solve_core(x, mu, T, nodes1, omega)
    return ForwardDiff.derivative(f_mu, mu_B)
end

# P 关于 omega 的偏导数（一次）
function dP_domega(x, mu_B, T, nodes1, omega)
    f_omega = w -> pressure_solve_core(x, mu_B, T, nodes1, w)
    return ForwardDiff.derivative(f_omega, omega)
end

# n 阶广义磁化率：定义为第 n 阶关于 mu_B 的导数乘以 T^(n-4)
function generalized_susceptibility_mu(n::Integer, x, mu_B, T, nodes1, omega)
    if n < 0
        throw(ArgumentError("n must be non-negative"))
    elseif n == 0
        # 0 阶就是压力本身（按定义）
        return pressure_solve_core(x, mu_B, T, nodes1, omega) * T^( -4 )
    end

    # 构造一个函数 f(mu) = P(...)，然后用 ForwardDiff 求 n 阶导数
    f_mu = mu -> pressure_solve_core(x, mu, T, nodes1, omega)

    # 使用 repeated differentiation via ForwardDiff
    deriv = f_mu
    for i in 1:n
        prev = deriv
        deriv = mu -> ForwardDiff.derivative(prev, mu)
    end

    # evaluate at provided mu_B and multiply by T^(n-4)
    return deriv(mu_B) * T^(n - 4)
end

# ---------------------------------------------------------------------------

# n 阶关于 omega 的广义磁化率：定义为第 n 阶关于 omega 的导数乘以 T^(n-4)
function generalized_susceptibility_omega(n::Integer, x, mu_B, T, nodes1, omega)
    if n < 0
        throw(ArgumentError("n must be non-negative"))
    elseif n == 0
        # 0 阶就是压力本身（按定义）
        return pressure_solve_core(x, mu_B, T, nodes1, omega) * T^( -4 )
    end

    # 构造一个函数 f(w) = P(..., omega = w)，然后用 ForwardDiff 求 n 阶导数
    f_w = w -> pressure_solve_core(x, mu_B, T, nodes1, w)

    deriv = f_w
    for i in 1:n
        prev = deriv
        deriv = w -> ForwardDiff.derivative(prev, w)
    end

    return deriv(omega) * T^(n - 4)
end
