include("Constants_PNJL.jl")
include("init.jl")
using .Constants_PNJL: hc, π, rho0, a0, a1, a2, b3, b4, T0, Nc, Lambda_f, G_f, K_f, m0, m0_q_f, m0_s_f
using .init: gauleg
using SpecialFunctions: log, exp

# 安全数学函数
@inline function safe_log(x; min_val=1e-16, handle_negative=:clamp)
    """
    安全对数函数，处理负值和接近零的值以确保数值稳定性
    
    参数：
    - x: 输入值
    - min_val: 最小值限制 (默认: 1e-16)
    - handle_negative: 处理负值的方式 (:clamp, :error, :nan)
    
    返回值：
    - 安全的对数值
    """
    if x <= 0
        if handle_negative == :error
            error("safe_log: Input must be positive, got $x")
        elseif handle_negative == :nan
            return NaN
        else  # :clamp
            return log(min_val)
        end
    elseif x < min_val
        return log(min_val)
    else
        return log(x)
    end
end
using ForwardDiff
using NLsolve

using BenchmarkTools
using StaticArrays
using FiniteDifferences

function get_nodes(p_num::Int, t_num::Int)
    # 动量节点和权重
    nodes1, weights1 = gauleg(0.0, Lambda_f, p_num)
    nodes2, weights2 = gauleg(0.0, 20.0, p_num)
    # 角度节点和权重（cosθ ∈ [0,1]）
    t_nodes, t_weights = gauleg(0.0, 1.0, t_num)

    # meshgrid (i,j) = (p, θ)
    p1_mesh = repeat(nodes1, 1, t_num)
    t1_mesh = repeat(t_nodes', p_num, 1)
    w1_mesh = weights1 * t_weights'  # 外积
    coefficient1 = w1_mesh .* p1_mesh .^ 2 ./ π^2  # 转换为球坐标后相同的系数p^2*4π

    p2_mesh = repeat(nodes2, 1, t_num)
    t2_mesh = t1_mesh
    w2_mesh = weights2 * t_weights'
    coefficient2 = w2_mesh .* p2_mesh .^ 2 ./ π^2  # 转换为球坐标后相同的系数p^2*4π

    nodes1 = [p1_mesh, t1_mesh, coefficient1]
    nodes2 = [p2_mesh, t2_mesh, coefficient2]
    return nodes1, nodes2
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
    # 使用安全对数函数避免负值或零值问题
    log_term = safe_log(value)
    U = T^4 * (-1 / 2 * Ta * Phi2 * Phi1 + Tb * log_term)  # 对数有效势
    return U
end

@inline function calculate_mass_vec(phi)
    """计算三种夸克的有效质量静态向量，兼容自动微分"""
    phiu, phid, phis = phi
    return SVector{3,eltype(phi)}(
        m0_q_f - 4 * G_f * phiu + 2 * K_f * phid * phis,
        m0_q_f - 4 * G_f * phid + 2 * K_f * phiu * phis,
        m0_s_f - 4 * G_f * phis + 2 * K_f * phiu * phid
    )
end

@inline function calculate_energy(mass_i, p, xi, t)
    p2 = p^2
    mass_i2 = mass_i^2
    term_xi = xi * (p * t)^2
    return sqrt(p2 + mass_i2 + term_xi)
end

@inline function calculate_log_term(E_i, mu_i, T, Phi1, Phi2)
    invT = 1.0 / T  # 预计算倒数
    x_i = (E_i - mu_i) * invT
    x_i_anti = (E_i + mu_i) * invT

    # 一次性计算所有指数项
    exp1 = exp(-x_i)
    exp2 = exp1 * exp1
    exp3 = exp1 * exp2
    exp1_anti = exp(-x_i_anti)
    exp2_anti = exp1_anti * exp1_anti
    exp3_anti = exp1_anti * exp2_anti

    f1_val = 1.0 + 3.0 * Phi1 * exp1 + 3.0 * Phi2 * exp2 + exp3
    f2_val = 1.0 + 3.0 * Phi2 * exp1_anti + 3.0 * Phi1 * exp2_anti + exp3_anti

    # 使用安全对数函数避免负值或零值问题
    return safe_log(f1_val) + safe_log(f2_val)
end

@inline function calculate_energy_sum(masses, p_nodes, coefficient, t_nodes, xi)
    """逐元素计算能量和"""
    total = 0.0
    # 完全展开嵌套循环，每次操作单个元素
    @inbounds for i in eachindex(masses)
        mass_i = masses[i]

        @inbounds @simd for j in eachindex(p_nodes)
            p = p_nodes[j]
            t = t_nodes[j]
            coefficient_j = coefficient[j]
            E = calculate_energy(mass_i, p, xi, t)
            total += E * coefficient_j
        end
    end
    return total * (-Nc)
end

@inline function calculate_log_sum(masses, p_nodes, Phi1, Phi2, mu, T, coefficient, t_nodes, xi)
    """逐元素计算对数项"""
    total = 0.0

    # 完全展开嵌套循环，每次操作单个元素
    @inbounds for i in eachindex(masses)
        mass_i = masses[i]
        mu_i = mu[i]  # 缓存数组元素
        @inbounds @simd for j in eachindex(p_nodes)
            p = p_nodes[j]
            t = t_nodes[j]
            coefficient_j = coefficient[j]
            E_i = calculate_energy(mass_i, p, xi, t)
            # 直接计算单个元素
            log_term = calculate_log_term(E_i, mu_i, T, Phi1, Phi2)
            total += log_term * coefficient_j
        end
    end
    return total * (-T)
end

function calculate_pressure(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi=0.0)
    """计算压力=-omega,T和mu传入前需归一化"""
    # 在函数开始时解包 nodes，并将数组部分转换为视图
    p_nodes1 = @view nodes_1[1][:]  # 假设 nodes[1] 是数组
    t_nodes1 = @view nodes_1[2][:]  # 假设 nodes[1] 是数组
    coef1 = @view nodes_1[3][:]  # 假设 nodes[1] 是数组
    p_nodes2 = @view nodes_2[1][:]  # 假设 nodes[2] 是数组
    t_nodes2 = @view nodes_2[2][:]  # 假设 nodes[2] 是数组
    coef2 = @view nodes_2[3][:]  # 假设 nodes[2] 是数组


    chi = calculate_chiral(phi)
    U = calculate_U(T, Phi1, Phi2)

    masses = calculate_mass_vec(phi)
    # 计算能量部分
    energy_sum = calculate_energy_sum(masses, p_nodes1, coef1, t_nodes1, xi)
    # 计算 log 部分
    log_sum = calculate_log_sum(masses, p_nodes2, Phi1, Phi2, mu, T, coef2, t_nodes2, xi)

    return -(chi + U + energy_sum + log_sum)
end

@inline function pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    phi = SVector{3}(x[1], x[2], x[3])
    Phi1, Phi2 = x[4], x[5]
    return calculate_pressure(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi)
end

function calculate_core(x, mu, T, nodes_1, nodes_2, xi)
    # 创建闭包函数，它捕获T、mu和nodes值
    f = x -> pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)

    return ForwardDiff.gradient(f, x)
end

@inline function calculate_rho(x, mu, T, nodes_1, nodes_2, xi)
    f_mu = mu -> pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    rho = ForwardDiff.gradient(f_mu, mu)
    return rho
end

@inline function calculate_thermo(x, mu, T, nodes_1, nodes_2, xi)
    rho = sum(calculate_rho(x, mu, T, nodes_1, nodes_2, xi)) / (3.0 * rho0)

    f_T = T -> pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    entropy = ForwardDiff.derivative(f_T, T)

    pressure = pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    energy = -pressure + sum(mu .* rho) + T * entropy  # 使用热力学关系计算能量

    return pressure, rho, entropy, energy
end

function calculate_t_rho(x, T, rho, nodes_1, nodes_2, xi, fvec=Vector{eltype(x)}(undef, 8))
    x_phi = SVector{5}(x[1:5])
    x_mu = SVector{3}(x[6:8])
    fvec[1:5] .= calculate_core(x_phi, x_mu, T, nodes_1, nodes_2, xi)
    fvec[6] = x_mu[1] - x_mu[2]  # μ_u - μ_d
    fvec[7] = x_mu[2] - x_mu[3]  # μ_d - μ_s
    fvec[8] = sum(calculate_rho(x_phi, x_mu, T, nodes_1, nodes_2, xi)) / (3.0 * rho0) - rho
    return fvec
end



function Trho(T_start, T_end; xi=0.0)
    # 获取节点（p_num=128, t_num=16）
    nodes_1, nodes_2 = get_nodes(1024, 32)

    # 输出目录和文件
    outdir = joinpath(@__DIR__, "..", "output")
    mkpath(outdir)
    outfile = joinpath(outdir, "trho_aniso.csv")

    # 初始x值
    #x_initial = [-1.8, -1.8, -2.1, 0.8, 0.8, 320 / hc, 320 / hc, 320 / hc]
    x_initial = [-0.07441, -0.07441, -1.86717, 0.02032, 0.02372, 2.15664, 2.15664, 2.15664]
    # 保存每个T下rho=3.00的解
    x_rho_3 = copy(x_initial)

    # 打开文件并写入表头，然后逐行追加扫描得到的解
    open(outfile, "w") do io
        println(io, "T,rho,phi_u,phi_d,phi_s,Phi1,Phi2,mu_u,mu_d,mu_s,pressure,entropy,energy,converged")

        # 主循环：按 T 扫描
        for T in T_start:-1/hc:T_end
            # 使用上一个T的rho=3.00的解作为初始值
            x = copy(x_rho_3)

            # 首先单独计算rho=3.00的情况
            rho = 3.00
            converged = false
            try
                res = nlsolve(x -> calculate_t_rho(x, T, rho, nodes_1, nodes_2, xi), x)
                converged = res.f_converged
                if converged
                    copyto!(x, res.zero)
                    # 保存当前T下rho=3.00的解，供下一个T循环使用
                    copyto!(x_rho_3, x)
                else
                    @warn "Root finding did not converge for T=$T and rho=$rho"
                end
            catch err
                @warn "Exception in root finding for T=$T and rho=$rho: $err"
                converged = false
            end

            # 计算热力学量（若收敛）并写入 csv 行
            if converged
                x_phi = SVector{5}(x[1:5])
                x_mu = SVector{3}(x[6:8])
                pressure, _, entropy, energy = calculate_thermo(x_phi, x_mu, T, nodes_1, nodes_2, 0.0)
            else
                pressure = NaN
                entropy = NaN
                energy = NaN
            end
            println(io, join([T * hc, rho, x..., pressure, entropy, energy, converged], ","))
            flush(io)

            # 然后计算剩余的rho值
            for rho in 2.99:-0.01:0.10
                try
                    res = nlsolve(x -> calculate_t_rho(x, T, rho, nodes_1, nodes_2, xi), x)
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
                    pressure, _, entropy, energy = calculate_thermo(x_phi, x_mu, T, nodes_1, nodes_2, 0.0)
                else
                    pressure = NaN
                    entropy = NaN
                    energy = NaN
                end
                println(io, join([T * hc, rho, x..., pressure, entropy, energy, converged], ","))
                flush(io)
            end
        end
    end

    return nothing
end


function pressure_solve_core(x, mu, T, nodes)
    X0_typed = convert.(promote_type(eltype(x), typeof(T)), x)
    res = nlsolve(x -> calculate_core(x, mu, T, nodes), X0_typed, autodiff=:forward)
    return pressure_wrapper(res.zero, mu, T, nodes)
end


function Tmu(; T_start, T_end, T_step, mu_start, mu_end, mu_step)
    # 节点
    nodes_1, nodes_2 = get_nodes(256, 16)

    # 输出文件
    outdir = joinpath(@__DIR__, "..", "output")
    mkpath(outdir)
    outfile = joinpath(outdir, "tmu_aniso.csv")

    # 初始 x（只有 5 个变量：phi_u, phi_d, phi_s, Phi1, Phi2）
    x_initial = [-1.8, -1.8, -2.1, 0.8, 0.8]
    x_prev = copy(x_initial)

    open(outfile, "w") do io
        println(io, "T,mu,phi_u,phi_d,phi_s,Phi1,Phi2,pressure,rho,entropy,energy,converged")

        for T in T_start:T_step:T_end
            # 对每个 T，在 mu 方向扫描
            for mu in mu_start:mu_step:mu_end
                x = copy(x_prev)
                converged = false
                pressure = NaN
                entropy = NaN
                energy = NaN
                rho = NaN
                try
                    # mu_vec 为三夸克相同的化学势
                    mu_vec = SVector{3}(mu, mu, mu)
                    # 使用 calculate_core 求解 5 个未知量
                    res = nlsolve(x -> calculate_core(x, mu_vec, T, nodes_1, nodes_2, 0.0), x; autodiff=:forward)
                    converged = res.f_converged
                    if converged
                        copyto!(x, res.zero)
                        x_prev .= x  # 用当前解作为下一点的初始值
                        # 计算热力学量（pressure, entropy, energy）
                        pressure, rho, entropy, energy = calculate_thermo(x, mu_vec, T, nodes_1, nodes_2, 0.0)
                    else
                        @warn "Root finding did not converge for T=$T and mu=$mu"
                    end
                catch err
                    @warn "Exception in root finding for T=$T and mu=$mu: $err"
                    converged = false
                end

                # 写入 csv：将 T, mu 转换回物理单位（乘 hc）以与 Trho 保持一致
                println(io, join([T * hc, mu * hc, x..., pressure, rho, entropy, energy, converged], ","))
                flush(io)
            end
        end
    end

    return nothing
end

#pressure_solve_core(x, mu, T, nodes)
#@show p = pressure_solve_core(x, mu, T, nodes)
#@code_warntype pressure_solve_core(x, mu, T, nodes)
#res = @benchmark pressure_solve_core(x, mu, T, nodes) samples=100 seconds=10
#display(res)
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
function dP_dT4_direct(x, mu, T, nodes, fdm)
    f = T -> pressure_solve_core(x, mu, T, nodes)
    return fdm(f, T)
end
#@show dP_dT4_direct(x, mu, T, nodes,fdm)
#res = @benchmark dP_dT4_direct(x, mu, T, nodes,fdm) samples=100 seconds=10
#display(res)

#Tmu(T_start=130/hc, T_end=131/hc, T_step=1/hc, mu_start=400/hc, mu_end=0.0, mu_step=-1/hc)
#Trho(T_start=100/hc, T_end=101/hc)

# ===================== find_cep 功能实现 =====================

function generate_single_temperature_data_aniso(T_target, xi=0.0)
    """为单个温度生成完整的rho-mu数据 (适用于aniso模型)"""
    nodes_1, nodes_2 = get_nodes(512, 32)

    # 初始猜测值 [phi_u, phi_d, phi_s, Phi1, Phi2, mu_u, mu_d, mu_s]
    #x_initial = [-0.4, -0.4, -1.9, 0.36, 0.38, 1.43, 1.43, 1.43]
    x_initial = [-0.07441, -0.07441, -1.86717, 0.02032, 0.02372, 2.15664, 2.15664, 2.15664]
    x = copy(x_initial)

    rho_mu_pairs = []

    # 从高密度到低密度扫描（与Trho函数一致）
    for rho in [3.00; collect(2.99:-0.01:0.10)]
        converged = false
        try
            res = nlsolve(x -> calculate_t_rho(x, T_target, rho, nodes_1, nodes_2, xi), x)
            converged = res.f_converged
            if converged
                copyto!(x, res.zero)
                # 对于aniso模型，取三个化学势的平均值作为代表性mu值
                mu_avg = (x[6] + x[7] + x[8]) / 3.0
                push!(rho_mu_pairs, (rho, mu_avg))
            else
                @debug "Root finding did not converge for T=$T_target and rho=$rho"
            end
        catch err
            @debug "Exception in root finding for T=$T_target and rho=$rho: $err"
            converged = false
        end
    end

    return rho_mu_pairs
end

function has_s_shape_aniso(T_mev, xi=0.0; min_negative_points=1, derivative_threshold=0.0)
    """
    检测给定温度下是否存在S形 (适用于aniso模型)
    
    参数:
    - T_mev: 温度 (MeV)
    - xi: 各向异性参数
    - min_negative_points: 至少需要多少个负导数点才认为存在S形
    - derivative_threshold: 负导数阈值
    
    返回: true表示存在S形，false表示不存在
    """
    T_julia = T_mev / hc  # 转换为Julia内部单位

    rho_mu_pairs = generate_single_temperature_data_aniso(T_julia, xi)

    if length(rho_mu_pairs) < 3
        @debug "温度 $T_mev MeV 处数据点太少，无法判断S形"
        return false
    end

    # 按rho排序（从小到大）
    sorted_pairs = sort(rho_mu_pairs, by=x -> x[1])

    # 计算数值导数 ∂mu/∂rho 并统计负值个数
    negative_derivative_count = 0
    max_negative_derivative = 0.0

    for i in 2:length(sorted_pairs)
        rho1, mu1 = sorted_pairs[i-1]
        rho2, mu2 = sorted_pairs[i]

        if rho2 != rho1  # 避免除零
            derivative = (mu2 - mu1) / (rho2 - rho1)
            if derivative < derivative_threshold
                negative_derivative_count += 1
                max_negative_derivative = min(max_negative_derivative, derivative)
            end
        end
    end

    has_s = negative_derivative_count >= min_negative_points

    if has_s
        println("  发现S形: $negative_derivative_count 个负导数点, 最小导数值 = $max_negative_derivative")
    end

    return has_s
end

function find_cep_aniso(xi, T_min=50.0, T_max=150.0, tolerance=0.5)
    """
    找到温度临界点T_cep，该温度下rho-mu关系正好不出现"S形" (适用于aniso模型)
    
    参数:
    - xi: 各向异性参数
    - T_min: 搜索的最低温度 (MeV)
    - T_max: 搜索的最高温度 (MeV) 
    - tolerance: 温度搜索精度 (MeV)
    
    返回:
    - T_cep: 临界温度 (MeV)，如果未找到相变则返回 NaN
    """

    scan_step = 5.0
    max_search_range = 200.0

    println("="^60)
    println("开始寻找临界温度 T_cep (PNJL aniso 模型)")
    println("各向异性参数 xi = $xi")
    println("搜索范围: [$T_min, $T_max] MeV")
    println("精度: $tolerance MeV")
    println("="^60)

    # 首先扫描整个范围，寻找S形存在的区域
    println("第一阶段: 扫描温度范围寻找S形现象...")

    scan_temps = collect(T_min:scan_step:T_max)
    s_shape_results = []

    for T_scan in scan_temps
        print("扫描 T=$T_scan MeV... ")
        has_s = has_s_shape_aniso(T_scan, xi)
        push!(s_shape_results, (T_scan, has_s))
        if !has_s
            println("无S形")
        end
    end

    # 寻找相变区间
    transition_found = false
    T_low = T_min
    T_high = T_max

    for i in 1:(length(s_shape_results)-1)
        T1, has_s1 = s_shape_results[i]
        T2, has_s2 = s_shape_results[i+1]

        if has_s1 != has_s2  # 找到相变
            T_low = has_s1 ? T1 : T2   # 有S形的温度作为下界
            T_high = has_s1 ? T2 : T1  # 无S形的温度作为上界
            transition_found = true
            println("✓ 发现相变区间: [$T_low, $T_high] MeV")
            break
        end
    end

    if !transition_found
        println("第二阶段: 扩展搜索范围...")

        # 分析初始扫描的结果
        all_have_s = all(result[2] for result in s_shape_results)
        all_no_s = all(!result[2] for result in s_shape_results)

        if all_have_s
            # 整个范围都有S形，说明临界点在更高温度
            println("初始范围内都有S形，向高温扩展搜索...")
            for T_test in (T_max+scan_step):scan_step:min(max_search_range, T_max + max_search_range)
                print("测试 T=$T_test MeV... ")
                current_has_s = has_s_shape_aniso(T_test, xi)

                if !current_has_s
                    println("无S形")
                    # 找到S形消失的点，临界点在T_max和T_test之间
                    T_low = T_max
                    T_high = T_test
                    transition_found = true
                    println("✓ 在高温区发现相变区间: [$T_low, $T_high] MeV")
                    break
                else
                    println("有S形")
                end
            end
        elseif all_no_s
            # 整个范围都无S形，说明临界点在更低温度
            println("初始范围内都无S形，向低温扩展搜索...")
            for T_test in (T_min-scan_step):-scan_step:max(5.0, T_min - max_search_range)
                print("测试 T=$T_test MeV... ")
                current_has_s = has_s_shape_aniso(T_test, xi)

                if current_has_s
                    println("有S形")
                    # 找到S形出现的点，临界点在T_test和T_min之间
                    T_low = T_test
                    T_high = T_min
                    transition_found = true
                    println("✓ 在低温区发现相变区间: [$T_low, $T_high] MeV")
                    break
                else
                    println("无S形")
                end
            end
        end
    end

    if !transition_found
        @warn "未找到S形相变，可能的原因："
        @warn "1. 当前物理参数下不存在S形相变"
        @warn "2. 相变温度超出搜索范围"
        @warn "3. 数值精度不足以检测S形"
        return NaN
    end

    # 二分搜索精确定位临界点
    println("第三阶段: 二分搜索精确定位临界点...")

    iteration = 0
    while (T_high - T_low) > tolerance && iteration < 15
        iteration += 1
        T_mid = (T_low + T_high) / 2

        print("第 $iteration 次迭代: 测试 T=$T_mid MeV... ")

        mid_has_s = has_s_shape_aniso(T_mid, xi)
        if !mid_has_s
            println("无S形")
        end

        if mid_has_s
            T_low = T_mid  # S形存在，临界点在更高温度
        else
            T_high = T_mid  # S形不存在，临界点在更低温度
        end

        println("    新区间: [$T_low, $T_high] MeV")
    end

    T_cep = (T_low + T_high) / 2
    println("="^60)
    println("✓ 找到临界温度: T_cep = $T_cep MeV")
    println("  最终搜索区间: [$T_low, $T_high] MeV")
    println("  区间宽度: $(T_high - T_low) MeV")
    println("="^60)

    # 记录结果到CSV文件
    if !isnan(T_cep)
        record_cep_result(xi, T_cep)
    end

    return T_cep
end

function record_cep_result(xi, T_cep)
    """记录CEP结果到CSV文件"""
    # 输出文件路径
    outdir = joinpath(@__DIR__, "..", "output")
    mkpath(outdir)
    outfile = joinpath(outdir, "aniso_cep.csv")

    # 检查文件是否存在，如果不存在则写入表头
    file_exists = isfile(outfile)

    open(outfile, "a") do io
        if !file_exists
            println(io, "xi,T_cep")
        end
        println(io, "$xi,$T_cep")
    end

    println("✓ 结果已记录到文件: $outfile")
    println("  xi = $xi, T_cep = $T_cep MeV")
end

# 便捷的测试函数
function test_find_cep_aniso()
    """测试find_cep_aniso功能的便捷函数"""
    println("测试find_cep_aniso功能...")
    println("注意: 如果当前参数设置下不存在S形相变，函数将返回NaN")

    # 使用较宽松的参数进行测试
    T_cep = find_cep_aniso(30.0, 100.0, 1.0)

    if isnan(T_cep)
        println("未找到临界温度，建议:")
        println("1. 调整物理参数 (xi, 初始条件等)")
        println("2. 扩大搜索范围")
        println("3. 调整S形检测敏感度")
    else
        println("找到临界温度: T_cep = $T_cep MeV")
    end

    return T_cep
end


#has_s_shape_aniso(140.0, xi)
#xi = 0.9
Trho(5.0/hc,5.0/hc,xi=0.9)
#find_cep_aniso(xi,10.0,20.0,0.01)
"""
for xi in 0.4:0.1:0.8
    find_cep_aniso(xi,90.0,115.0,0.01)
end
"""