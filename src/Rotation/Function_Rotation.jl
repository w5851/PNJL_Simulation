module Function_Rotation

include("Constants_Rotation.jl")
include("../init.jl")
using .Constants_Rotation:hc, π, rho0, a0, a1, a2, a3, b3, b4, T0, Nc, Lambda_f, G_f, K_f, m0_q_f, m0_s_f, r0, coefficients
using .init:gauleg
using SpecialFunctions: log, exp, besselj
using ForwardDiff
using NLsolve

using BenchmarkTools
using StaticArrays
using FiniteDifferences

# Minimal public API: only export symbols used by examples and Advanced_FindCep
# Advanced_FindCep.jl and examples require the core numerical routines and a few constants.
export get_nodes,
       calculate_t_rho,
       calculate_core,
       calculate_thermo,
       calculate_mass,
       pressure_solve_core,
       hc,
       rho0,
       T0,
       Nc,
       Lambda_f,
       G_f,
       m0_q_f,
       m0_s_f

function init_bessel(p,theta,n,w)
    t = sin.(theta)
    ptr2 = (p .* t .* r0).^2

    # 贝塞尔函数项，支持数组广播
    bessel_term = besselj.(n .+ 1, ptr2) .+ besselj.(n, ptr2)

    coefficient = w .* p.^2 .* t .* bessel_term

    return coefficient ./(4.0 .* π^2)
end

function get_nodes(p_num::Int, t_num::Int)
    # 动量节点和权重
    p_nodes, p_weights = gauleg(0.0, Lambda_f, p_num)
    p2_nodes, p2_weights = gauleg(0.0, 20.0, p_num)

    # 角度节点和权重（极角 θ ∈ [0, π]）
    t_nodes, t_weights = gauleg(0.0, π, t_num)

    # 离散节点和权重
    n_nodes = collect(-5:1:5)              # -5 到 5 共 11 个整数
    n_weights = ones(length(n_nodes))      # 权重均为 1

    # 创建三维联合网格 (p, t, n)
    p_mesh  = reshape(p_nodes,  p_num, 1, 1) .* ones(1, t_num, length(n_nodes))
    t_mesh  = reshape(t_nodes,  1, t_num, 1) .* ones(p_num, 1, length(n_nodes))
    n_mesh  = reshape(n_nodes,  1, 1, length(n_nodes)) .* ones(p_num, t_num, 1)
    w_mesh  = reshape(p_weights, p_num, 1, 1) .* reshape(t_weights, 1, t_num, 1) .* reshape(n_weights, 1, 1, length(n_nodes))
    coefficient1 = init_bessel(p_mesh, t_mesh, n_mesh, w_mesh)

    p2_mesh = reshape(p2_nodes, p_num, 1, 1) .* ones(1, t_num, length(n_nodes))
    t2_mesh = t_mesh
    n2_mesh = n_mesh
    w2_mesh = reshape(p2_weights, p_num, 1, 1) .* reshape(t_weights, 1, t_num, 1) .* reshape(n_weights, 1, 1, length(n_nodes))
    coefficient2 = init_bessel(p2_mesh, t2_mesh, n2_mesh, w2_mesh)

    # 返回与 Python np.stack 结构一致的四个三维数组
    nodes1 = [p_mesh, n_mesh, coefficient1]
    nodes2 = [p2_mesh, n2_mesh, coefficient2]
    return nodes1,nodes2
end

@inline function calculate_chiral(phi)
    """计算手征相关量"""
    term1 = G_f * mean(phi.^2)
    return term1
end

@inline function calculate_U(T, Phi1, Phi2)
    T_ratio = T0 / T
    T_b2 = a0 + a1 * T_ratio + a2 * T_ratio^2 + a3 * T_ratio^3
    poly = -0.5 * T_b2 * Phi1 * Phi2 - b3/6 * (Phi1^3 + Phi2^3) + b4/4 * (Phi1^2 * Phi2^2)
    return T^4 * poly
end

@inline function calculate_mass(phi)
    return m0_q_f - 2.0 * G_f * phi
end

@inline function calculate_energy(mass, p, n, omega)
    p2 = p^2
    mass2 = mass^2
    return  sqrt(p2 + mass2) - (0.5 + n) * omega
end

@inline function AA(x,T,Phi1,Phi2)
    exp1 =  exp(-x / T)
    exp2 = exp1 * exp1
    exp3 = exp1 * exp2
    f1 =  1.0 + 3.0 * Phi1 * exp1 + 3.0 * Phi2 * exp2 + exp3
    return f1
end

@inline function AAbar(x,T,Phi1,Phi2)
    exp1 =  exp(-x / T)
    exp2 = exp1 * exp1
    exp3 = exp1 * exp2
    f2 = 1.0 + 3.0 * Phi2 * exp1 + 3.0 * Phi1 * exp2 + exp3
    return f2
end

@inline function calculate_log_term(E_i, mu_i, T, Phi1, Phi2)
    f1 = AA(E_i - mu_i, T, Phi1, Phi2)
    f2 = AA(-E_i - mu_i, T, Phi1, Phi2)
    f3 = AAbar(-E_i + mu_i, T, Phi1, Phi2)
    f4 = AAbar(E_i + mu_i, T, Phi1, Phi2)
    return  log(f1) + log(f2)+ log(f3) + log(f4)
end

@inline function calculate_log_sum(masses, p_nodes, Phi1, Phi2, mu, T, coefficient, n_nodes, omega)
    """逐元素计算对数项"""
    total = 0.0
    
    # 完全展开嵌套循环，每次操作单个元素
    @inbounds for i in eachindex(masses)
        mass_i = masses[i]
        mu_i = mu[i]  # 缓存数组元素
        @inbounds @simd for j in eachindex(p_nodes)
            p = p_nodes[j]
            n = n_nodes[j]
            coefficient_j = coefficient[j]
            E_i = calculate_energy(mass_i, p,n,omega)
            # 直接计算单个元素
            log_term = calculate_log_term(E_i, mu_i, T, Phi1, Phi2)
            total += log_term * coefficient_j
        end
    end
    return total * (-T)
end

# f 函数的系数计算
function calc_factors(T, omega)
    c = Constants.coefficients
    a = c["a"][1] + c["a"][2] * omega^2 + c["a"][3] * omega^4
    b = c["b"][1] + c["b"][2] * omega^2 + c["b"][3] * omega^4
    c_ = c["c"][1] + c["c"][2] * omega^2 + c["c"][3] * omega^4
    d = c["d"][1] + c["d"][2] * omega^2 + c["d"][3] * omega^4
    f = a .* tanh.(b * (T / Constants.T0 - c_)) .+ d
    return f, 1 ./ f
end

# 交互 U
function calc_U(T, Phi1, Phi2, omega)
    f, f_inv = calc_factors(T, omega)
    C = Constants.C
    T0 = Constants.T0
    term = -C * f * (T/T0)^2 * Phi1 * Phi2 - (Phi1^3 + Phi2^3)/3 + (1/C) * f_inv * (T0/T)^2 * (Phi1*Phi2)^2
    return T^4 * term
end

function calculate_pressure(phi,Phi1,Phi2,mu,T,nodes1,omega)
    """计算压力=-omega,T和mu传入前需归一化"""
    # 在函数开始时解包 nodes，并将数组部分转换为视图
    p_nodes2 = @view nodes1[1][:]  # 假设 nodes2[1] 是数组
    n_nodes2 = @view nodes1[2][:]  # 假设 nodes2[2] 是数组
    coef2 = @view nodes1[3][:]  # 假设 nodes2[3] 是数组

    chi = calculate_chiral(phi)
    U = calculate_U(T,Phi1,Phi2)
    
    masses = calculate_mass(phi)
    
    # 计算 log 部分
    log_sum = calculate_log_sum(masses, p_nodes2, Phi1, Phi2, mu, T, coef2,n_nodes2,omega)
    
    return -(chi+U+log_sum)
end

@inline function pressure_wrapper(x, mu, T, nodes1, omega)
    """压力计算的包装函数"""
    phi = x[1]
    Phi1 = x[2]
    Phi2 = x[3]
    return calculate_pressure(phi, Phi1, Phi2, mu, T, nodes1, omega)
end

function calculate_core(x, mu, T, nodes1, omega)
    """核心计算函数"""
    f = x -> pressure_wrapper(x, mu, T, nodes1, omega)
    return ForwardDiff.gradient(f, x)
end

@inline function calculate_rho(x,mu,T,nodes1,omega)
    f_mu = mu -> pressure_wrapper(x, mu, T, nodes1,omega)
    rho = ForwardDiff.derivative(f_mu, mu)
    return rho
end

@inline function calculate_thermo(x , mu,T,nodes1,omega)
    rho = calculate_rho(x, mu, T, nodes1,omega) / rho0

    f_T = T -> pressure_wrapper(x, mu, T, nodes1,omega)
    entropy = ForwardDiff.derivative(f_T, T)

    pressure = pressure_wrapper(x, mu, T, nodes1,omega)
    energy = -pressure + sum(mu .* rho) + T * entropy  # 使用热力学关系计算能量

    return pressure,rho, entropy,energy
end

function calculate_t_rho(x,T,rho,nodes1,omega, fvec=Vector{eltype(x)}(undef, 4))
    x_phi = SVector{3}(x[1:3])
    x_mu = x[4]
    fvec[1:3] .= calculate_core(x_phi, x_mu, T, nodes1,omega)
    fvec[4] = calculate_rho(x_phi, x_mu, T, nodes1,omega) / rho0 - rho
    return fvec
end


function pressure_solve_core(x, mu, T, nodes1, omega)
    """Rotation模型的压力求解核心函数"""
    X0_typed = convert.(promote_type(eltype(x), typeof(T)), x)
    res = nlsolve(x -> calculate_core(x, mu, T, nodes1, omega), X0_typed, autodiff=:forward)
    return pressure_wrapper(res.zero, mu, T, nodes1, omega)
end

# 压力对温度的导数函数（用于分析）
function dP_dT_rotation(x, mu, T, nodes1, omega)
    f = T -> pressure_solve_core(x, mu, T, nodes1, omega)
    return ForwardDiff.derivative(f, T)
end

function dP_dT2_rotation(x, mu, T, nodes1, omega)
    f = T -> dP_dT_rotation(x, mu, T, nodes1, omega)
    return ForwardDiff.derivative(f, T)
end

function dP_dT3_rotation(x, mu, T, nodes1, omega)
    f = T -> dP_dT2_rotation(x, mu, T, nodes1, omega)
    return ForwardDiff.derivative(f, T)
end

function dP_dT4_rotation(x, mu, T, nodes1, omega)
    f = T -> dP_dT3_rotation(x, mu, T, nodes1, omega)
    return ForwardDiff.derivative(f, T)
end


end # module Function_Rotation