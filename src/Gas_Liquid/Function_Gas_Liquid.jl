include("Constants_Gas_Liquid.jl")
include("../init.jl")
using .Constants_Gas_Liquid: π, hc, m, mσ, mω, mρ, mδ, calculate_couplings
using .init: gauleg
using NLsolve
using FiniteDifferences

using BenchmarkTools

function get_nodes(n_p)
    p_nodes,p_weights = gauleg(0,20.0,n_p)
    coefficient = @. p_nodes^2 * p_weights / π^2
    return [p_nodes,coefficient]
end

@inline function fermion(E,μ,T)
    """费米子分布函数"""
    return 1 / (exp((E - μ) / T) + 1)
end

@inline function fermion_anti(E,μ,T)
    """反费米子分布函数"""
    return 1 / (exp((E + μ) / T) + 1)
end

@inline function calculate_log(E,μ,T)
    x = E - μ
    x_anti = E + μ
    term1 = 1 + exp(-x / T)
    term2 = 1 + exp(-x_anti / T)
    return log(term1) + log(term2)
end

@inline function calculate_mass(gσ,gδ)
    m_p = m-gσ-gδ
    m_n = m-gσ+gδ
    return m_p,m_n
end

@inline function calculate_energy(gσ,gδ,p_nodes)
    m_p,m_n = calculate_mass(gσ,gδ)
    p2 = @. p_nodes^2
    E_p = @. sqrt(p2 + m_p^2)
    E_n = @. sqrt(p2 + m_n^2)
    return E_p,E_n
end

@inline function calculate_ρ(E,μ,T,coef)
    # 避免创建中间数组，使用 map-reduce 模式
    return mapreduce(i -> (fermion(E[i],μ,T) - fermion_anti(E[i],μ,T))*coef[i], +, eachindex(E))
end

@inline function calculate_ρ_s(E,μ,T,coef,m)
    # 避免创建中间数组，使用 map-reduce 模式
    return mapreduce(i -> (fermion(E[i],μ,T) + fermion_anti(E[i],μ,T))*coef[i]*m/E[i], +, eachindex(E))
end

@inline function calculate_σ_term(gσ,ρ_ps,ρ_ns,couplings)
    fσ, _, _, _, b, c = couplings
    return - gσ + fσ * (ρ_ps + ρ_ns - b * m * gσ^2 - c * gσ^3)
end

@inline function calculate_δ_term(gδ,ρ_ps,ρ_ns,couplings)
    fδ = couplings[4]
    # 如果fδ = 0，则强制gδ = 0（δ介子耦合消失）
    if fδ == 0.0
        return -gδ  # 这将强制gδ = 0以满足约束条件
    else
        return -gδ + fδ * (ρ_ps - ρ_ns)
    end
end

@inline function calculate_ρ_term(gρ,ρ_p,ρ_n,couplings)
    fρ = couplings[3]
    return -gρ + fρ * (ρ_p - ρ_n)
end

@inline function calculate_ω_term(gω,ρ_p,ρ_n,couplings)
    fω = couplings[2]
    return -gω + fω * (ρ_p + ρ_n)
end

@inline function calculate_chemical_constraint(μ_B, μ_n,ρ_p, ρ_n,couplings)
    """计算重子化学势约束条件(同位旋不对称:质子中子数密度不等)"""
    gω = calculate_ω_term(0.0, ρ_p, ρ_n, couplings)
    gρ = calculate_ρ_term(0.0, ρ_p, ρ_n, couplings)
    return μ_B - μ_n - gω + gρ
end

@inline function calculate_asymmetry_constraint(ρ_n, ρ_p, target_asymmetry=0.198)
    """计算同位旋不对称度约束条件"""
    return target_asymmetry - (ρ_n - ρ_p)/(ρ_n + ρ_p)
end

function calculate_fun_constraint(x,nodes,couplings,params)
    """计算化学势约束条件下的残差方程"""
    gσ,gδ,μ_p,μ_n = x
    p_nodes, coefficient = nodes
    T = params[1]
    μ_B = params[2]
    
    m_p, m_n = calculate_mass(gσ, gδ)
    E_p, E_n = calculate_energy(gσ, gδ, p_nodes)
    ρ_p = calculate_ρ(E_p, μ_p, T, coefficient)
    ρ_n = calculate_ρ(E_n, μ_n, T, coefficient)
    ρ_ps = calculate_ρ_s(E_p, μ_p, T, coefficient, m_p)
    ρ_ns = calculate_ρ_s(E_n, μ_n, T, coefficient, m_n)

    σ_term = calculate_σ_term(gσ, ρ_ps, ρ_ns, couplings)
    δ_term = calculate_δ_term(gδ, ρ_ps, ρ_ns, couplings)
    # 计算重子化学势约束
    chem_constraint = calculate_chemical_constraint(μ_B, μ_n, ρ_p, ρ_n, couplings)
    
    # 计算同位旋不对称度约束
    asymmetry_constraint = calculate_asymmetry_constraint(ρ_n, ρ_p)

    return [σ_term, δ_term, chem_constraint, asymmetry_constraint]
end

function solve_fun_constraints(x0, nodes,couplings, params)
    """求解化学势约束条件"""
       
    # 使用NLsolve求解
    result = nlsolve(x -> calculate_fun_constraint(x, nodes, couplings, params), x0)
    
    return result.zero
end

@inline function calculate_init_term(E,μ,T,nodes)
    """计算积分项"""
    p, coef = nodes
    return mapreduce(i -> (fermion(E[i],μ,T) + fermion_anti(E[i],μ,T))*coef[i]*p[i]^2/E[i], +, eachindex(E))/3.0
end

@inline function calculate_pressure(gσ, gδ, gω, gρ, μ_p, μ_n, T, nodes, couplings)
    """计算压强"""
    fσ, fω, fρ, fδ, b, c = couplings
    p_nodes, _ = nodes
    
    # 计算质子和中子的有效质量和能量
    E_p, E_n = calculate_energy(gσ, gδ, p_nodes)

    # 计算p_p和p_n（动量积分项）
    p_p = calculate_init_term(E_p, μ_p, T, nodes)
    p_n = calculate_init_term(E_n, μ_n, T, nodes)
    
    # 计算压强
    pressure = -(1.0/3.0) * b * m * gσ^3 - 
               (1.0/4.0) * c * gσ^4 - 
               (1.0/(2.0*fσ)) * gσ^2 + 
               (1.0/(2.0*fω)) * gω^2 + 
               p_p + p_n + 
               (1.0/(2.0*fρ)) * gρ^2
    
    # 处理δ介子项：如果fδ = 0，则该项为0；否则正常计算
    if fδ != 0.0
        pressure -= (1.0/(2.0*fδ)) * gδ^2
    end
               
    return pressure
end

@inline function calculate_pressure_wrapper(x, nodes, couplings, params)
    """计算压强（自动计算场量gω和gρ）"""
    p_nodes, coef = nodes
    T, _ = params
    gσ, gδ, μ_p, μ_n = x

    # 计算密度用于确定场量
    E_p, E_n = calculate_energy(gσ, gδ, p_nodes)
    ρ_p = calculate_ρ(E_p, μ_p, T, coef)
    ρ_n = calculate_ρ(E_n, μ_n, T, coef)
    
    # 根据自洽条件计算场量
    gω = calculate_ω_term(0.0, ρ_p, ρ_n, couplings)
    gρ = calculate_ρ_term(0.0, ρ_p, ρ_n, couplings)

    return calculate_pressure(gσ, gδ, gω, gρ, μ_p, μ_n, T, nodes, couplings)
end

function calculate_pressure_solved(μ_B,T, x0, nodes, couplings)
    """计算压强，使用求解后的场量"""
    params = [T, μ_B]
    x = solve_fun_constraints(x0, nodes, couplings, params)
    return calculate_pressure_wrapper(x, nodes, couplings, params)
end
