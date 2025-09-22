include("Constants_PNJL.jl")
include("init.jl")
using .Constants_PNJL:hc, π, rho0, a0, a1, a2, b3, b4, T0, Nc, Lambda_f, G_f, K_f, m0,m0_q_f,m0_s_f
using .init:gauleg
using SpecialFunctions: log, exp
using ForwardDiff
using NLsolve

using BenchmarkTools
using StaticArrays
using FiniteDifferences

function get_nodes(n_p::Int)
    """获取积分"""
    # p:对应积分节点
    p_nodes, p_weights = gauleg(0.0, Lambda_f, n_p)
    p_nodes2, p_weights2 = gauleg(0.0, 20.0, n_p)
    coefficient1 = p_weights.*p_nodes.^2 ./ π^2  # 转换为球坐标后相同的系数p^2*4π
    coefficient2 = p_weights2.*p_nodes2.^2 ./ π^2  # 转换为球坐标后相同的系数p^2*4π
    return [p_nodes,p_nodes2,coefficient1,coefficient2]
end

@inline function calculate_chiral(phi)
    """计算手征相关量"""
    term1 = 2 * G_f * sum(phi .^ 2) - 4 * K_f * prod(phi)
    return term1
end

@inline function calculate_U(T, Phi1, Phi2)
    """计算极化Polyakov-loop势能"""
    T_ratio = T0 / T
    Ta = a0 + a1 * T_ratio + a2 * T_ratio^2
    Tb = b3 * T_ratio^3
    value = 1 - 6 * Phi2 * Phi1 + 4 * (Phi2^3 + Phi1^3) - 3 * (Phi2 * Phi1)^2
    log_term =  log(value)
    U = T^4 * (-1/2 * Ta * Phi2 * Phi1 + Tb * log_term)  # 对数有效势
    return U
end

@inline function calculate_mass_vec(phi)
    """计算三种夸克的有效质量静态向量，兼容自动微分"""
    phiu, phid, phis = phi
    return SVector{3, eltype(phi)}(
        m0_q_f - 4 * G_f * phiu + 2 * K_f * phid * phis,
        m0_q_f - 4 * G_f * phid + 2 * K_f * phiu * phis,
        m0_s_f - 4 * G_f * phis + 2 * K_f * phiu * phid
    )
end

@inline function calculate_energy(mass_i, p)
    p2 = p^2
    mass_i2 = mass_i^2
    return  sqrt(p2 + mass_i2)
end

@inline function calculate_log_term(E_i, mu_i, T, Phi1, Phi2)
    invT = 1.0 / T  # 预计算倒数
    x_i = (E_i - mu_i) * invT
    x_i_anti = (E_i + mu_i) * invT
    
    # 一次性计算所有指数项
    exp1 =  exp(-x_i)
    exp2 = exp1 * exp1
    exp3 = exp1 * exp2
    exp1_anti =  exp(-x_i_anti)
    exp2_anti = exp1_anti * exp1_anti
    exp3_anti = exp1_anti * exp2_anti
    
    f1_val = 1.0 + 3.0 * Phi1 * exp1 + 3.0 * Phi2 * exp2 + exp3
    f2_val = 1.0 + 3.0 * Phi2 * exp1_anti + 3.0 * Phi1 * exp2_anti + exp3_anti
    
    return  log(f1_val) + log(f2_val)
end

@inline function calculate_energy_sum(masses, p_nodes, coefficient)
    """逐元素计算能量和"""
    total = 0.0
    # 完全展开嵌套循环，每次操作单个元素
    @inbounds for i in eachindex(masses)
        mass_i = masses[i]
        
        @inbounds @simd for j in eachindex(p_nodes)
            E = calculate_energy(mass_i, p_nodes[j])
            total += E * coefficient[j]
        end
    end
    return total * (-Nc)
end

@inline function calculate_log_sum(masses, p_nodes, Phi1, Phi2, mu, T, coefficient)
    """逐元素计算对数项"""
    total = 0.0
    
    # 完全展开嵌套循环，每次操作单个元素
    @inbounds for i in eachindex(masses)
        mass_i = masses[i]
        mu_i = mu[i]  # 缓存数组元素
        @inbounds @simd for j in eachindex(p_nodes)
            p = p_nodes[j]
            E_i = calculate_energy(mass_i, p)
            coefficient_j = coefficient[j]
            # 直接计算单个元素
            log_term = calculate_log_term(E_i, mu_i, T, Phi1, Phi2)
            total += log_term * coefficient_j
        end
    end
    return total * (-T)
end

function calculate_pressure(phi,Phi1,Phi2,mu,T,nodes)
    """计算压力=-omega,T和mu传入前需归一化"""
    # 在函数开始时解包 nodes，并将数组部分转换为视图
    p_nodes1 = @view nodes[1][:]  # 假设 nodes[1] 是数组
    p_nodes2 = @view nodes[2][:]  # 假设 nodes[2] 是数组
    coef1 = @view nodes[3][:]  # 假设 nodes[3] 是数组
    coef2 = @view nodes[4][:]  # 假设 nodes[4] 是数组

    chi = calculate_chiral(phi)
    U = calculate_U(T,Phi1,Phi2)
    
    masses = calculate_mass_vec(phi)
    # 计算能量部分
    energy_sum = calculate_energy_sum(masses, p_nodes1, coef1)
    # 计算 log 部分
    log_sum = calculate_log_sum(masses, p_nodes2, Phi1, Phi2, mu, T, coef2)
    
    return -(chi+U+energy_sum+log_sum)
end

@inline function pressure_wrapper(x, mu, T, nodes)
    phi = SVector{3}(x[1], x[2], x[3])
    Phi1, Phi2 = x[4], x[5]
    return calculate_pressure(phi, Phi1, Phi2, mu, T, nodes)
end

function calculate_core(x, mu, T, nodes)
    # 创建闭包函数，它捕获T、mu和nodes值
    f = x -> pressure_wrapper(x, mu, T, nodes)
    
    return ForwardDiff.gradient(f, x)
end

@inline function calculate_rho(x,mu,T,nodes)
    f_mu = mu -> pressure_wrapper(x, mu, T, nodes)
    rho = ForwardDiff.gradient(f_mu, mu)
    return rho
end

@inline function calculate_thermo(x , mu,T,nodes)
    rho = sum(calculate_rho(x, mu, T, nodes))/(3.0*rho0)

    f_T = T -> pressure_wrapper(x, mu, T, nodes)
    entropy = ForwardDiff.derivative(f_T, T)

    pressure = pressure_wrapper(x, mu, T, nodes)
    energy = -pressure + sum(mu .* rho) + T * entropy  # 使用热力学关系计算能量

    return pressure,rho, entropy,energy
end

function calculate_t_rho(x,T,rho,nodes, fvec=Vector{eltype(x)}(undef, 8))
    x_phi = SVector{5}(x[1:5])
    x_mu = SVector{3}(x[6:8])
    fvec[1:5] .= calculate_core(x_phi, x_mu, T, nodes)
    fvec[6] = x_mu[1] - x_mu[2]  # μ_u - μ_d
    fvec[7] = x_mu[2] - x_mu[3]  # μ_d - μ_s
    fvec[8] = sum(calculate_rho(x_phi, x_mu, T, nodes)) / (3.0*rho0) - rho
    return fvec
end

function Trho(T_start, T_end)
    # 获取节点（p_num=128）
    nodes = get_nodes(128)

    # 输出目录和文件
    outdir = joinpath(@__DIR__, "..", "output")
    mkpath(outdir)
    outfile = joinpath(outdir, "trho.csv")

    # 初始x值（包括 μ_u, μ_d, μ_s）
    x_initial = [-1.8, -1.8, -2.1, 0.8, 0.8, 320 / hc, 320 / hc, 320 / hc]
    x_rho_3 = copy(x_initial)

    open(outfile, "w") do io
        println(io, "T,rho,phi_u,phi_d,phi_s,Phi1,Phi2,mu_u,mu_d,mu_s,pressure,entropy,energy,converged")

        for T in T_start:1/hc:T_end
            # 使用上一个T的rho=3.00的解作为初始值
            x = copy(x_rho_3)

            # 首先单独计算rho=3.00的情况
            rho = 3.00
            converged = false
            try
                res = nlsolve(x -> calculate_t_rho(x, T, rho, nodes), x)
                converged = res.f_converged
                if converged
                    copyto!(x, res.zero)
                    copyto!(x_rho_3, x)
                else
                    @warn "Root finding did not converge for T=$T and rho=$rho"
                end
            catch err
                @warn "Exception in root finding for T=$T and rho=$rho: $err"
                converged = false
            end

            if converged
                x_phi = SVector{5}(x[1:5])
                x_mu = SVector{3}(x[6:8])
                pressure, _, entropy, energy = calculate_thermo(x_phi, x_mu, T, nodes)
            else
                pressure = NaN
                entropy = NaN
                energy = NaN
            end
            println(io, join([T*hc, rho, x..., pressure, entropy, energy, converged], ","))
            flush(io)

            # 然后计算剩余的rho值
            for rho in 2.99:-0.01:0.10
                try
                    res = nlsolve(x -> calculate_t_rho(x, T, rho, nodes), x)
                    converged = res.f_converged
                    if converged
                        copyto!(x, res.zero)
                    else
                        @warn "Root finding did not converge for T=$T and rho=$rho"
                    end
                catch err
                    @warn "Exception in root finding for T=$T and rho=$rho: $err"
                    converged = false
                end

                if converged
                    x_phi = SVector{5}(x[1:5])
                    x_mu = SVector{3}(x[6:8])
                    pressure, _, entropy, energy = calculate_thermo(x_phi, x_mu, T, nodes)
                else
                    pressure = NaN
                    entropy = NaN
                    energy = NaN
                end
                println(io, join([T*hc, rho, x..., pressure, entropy, energy, converged], ","))
                flush(io)
            end
        end
    end

    return nothing
end

function pressure_solve_core(x, mu, T, nodes)
    X0_typed = convert.(promote_type(eltype(x), typeof(T)), x)
    res = nlsolve(x -> calculate_core(x,mu,T,nodes), X0_typed, autodiff=:forward)
    return pressure_wrapper(res.zero,mu,T,nodes)
end
"""
x = [-0.1, -0.1, -1.7, 0.5, 0.5]
mu = [320 / hc, 320 / hc, 320 / hc]
T = 150 / hc
nodes = get_nodes(128)
@show res = pressure_wrapper(x, mu, T, nodes)
#pressure_solve_core(x, mu, T, nodes)
#@show p = pressure_solve_core(x, mu, T, nodes)
#@code_warntype pressure_solve_core(x, mu, T, nodes)
#res = @benchmark pressure_solve_core(x, mu, T, nodes) samples=100 seconds=10
#display(res)
"""
function dP_dT(x, mu, T, nodes)
    f = T -> pressure_solve_core(x, mu, T, nodes)
    return ForwardDiff.derivative(f, T)
end
function dP_dT2(x, mu, T, nodes)
    f = T -> dP_dT(x, mu, T, nodes)
    return ForwardDiff.derivative(f, T)
end
function dP_dT3(x, mu, T, nodes)
    f = T -> dP_dT2(x, mu, T, nodes)
    return ForwardDiff.derivative(f, T)
end
function dP_dT4(x, mu, T, nodes)
    f = T -> dP_dT3(x, mu, T, nodes)
    return ForwardDiff.derivative(f, T)
end
#res = @benchmark dP_dT(x,mu,T,nodes) samples=100 seconds=10
#res = @benchmark dP_dT2(x,mu,T,nodes) samples=100 seconds=10
#res = @benchmark dP_dT3(x,mu,T,nodes) samples=100 seconds=10
#@show dP_dT4(x, mu, T, nodes)
#res = @benchmark dP_dT4(x,mu,T,nodes) samples=100 seconds=10
#display(res)

#fdm = central_fdm(5, 4)  # 五点中心差分，步长自动选择
function dP_dT4_direct(x, mu, T, nodes,fdm)
    f = T -> pressure_solve_core(x, mu, T, nodes)
    return fdm(f, T)
end
#@show dP_dT4_direct(x, mu, T, nodes,fdm)
#res = @benchmark dP_dT4_direct(x, mu, T, nodes,fdm) samples=100 seconds=10
#display(res)
Trho(100/hc,101/hc)