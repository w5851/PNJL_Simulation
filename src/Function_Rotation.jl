include("Constants_Rotation.jl")
include("init.jl")
using .Constants_Rotation:hc, π, rho0, a0, a1, a2, a3, b3, b4, T0, Nc, Lambda_f, G_f, K_f, m0_q_f, m0_s_f, r0, coefficients
using .init:gauleg
using SpecialFunctions: log, exp, besselj
using ForwardDiff
using NLsolve

using BenchmarkTools
using StaticArrays
using FiniteDifferences

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

function Trho(T_start,T_end)
    nodes1, nodes2 = get_nodes(128,16)
    omega = 100/hc

    # 输出文件
    outdir = joinpath(@__DIR__, "..", "output")
    mkpath(outdir)
    outfile = joinpath(outdir, "trho_rotation.csv")

    x_initial = [-2.13,0.06,0.12, 310 / hc]
    x_rho_3 = copy(x_initial)

    open(outfile, "w") do io
        println(io, "T,rho,phi,Phi1,Phi2,mu,pressure,entropy,energy,converged")

        for T in T_start:1/hc:T_end
            x = copy(x_rho_3)
            rho = 6.00
            converged = false
            try
                res = nlsolve(x -> calculate_t_rho(x, T, rho, nodes1,omega), x)
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
                x_phi = SVector{3}(x[1:3])
                x_mu = x[4]
                pressure, _, entropy, energy = calculate_thermo(x_phi, x_mu, T, nodes1,omega)
            else
                pressure = NaN
                entropy = NaN
                energy = NaN
            end
            println(io, join([T*hc, rho, x..., pressure, entropy, energy, converged], ","))
            flush(io)

            for rho in 5.99:-0.01:0.10
                try
                    res = nlsolve(x -> calculate_t_rho(x, T, rho, nodes1,omega), x)
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
                    x_phi = SVector{3}(x[1:3])
                    x_mu = x[4]
                    pressure, _, entropy, energy = calculate_thermo(x_phi, x_mu, T, nodes1,omega)
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

function Tmu(;T_start, T_end, T_step, mu_start, mu_end, mu_step,omega=100/hc)
    # 节点
    nodes1, nodes2 = get_nodes(128, 16)
    #omega = 100/hc

    # 输出文件
    outdir = joinpath(@__DIR__, "..", "output")
    mkpath(outdir)
    outfile = joinpath(outdir, "tmu_rotation.csv")

    # 初始 x（只有 3 个变量：phi, Phi1, Phi2）
    x_initial = [-2.13, 0.06, 0.12]
    x_prev = copy(x_initial)

    open(outfile, "w") do io
        println(io, "T,mu,phi,Phi1,Phi2,mass,pressure,rho,entropy,energy,converged")

        for T in T_start:T_step:T_end
            # 对每个 T，在 mu 方向扫描
            for mu in mu_start:mu_step:mu_end
                x = copy(x_prev)
                converged = false
                pressure = NaN
                entropy = NaN
                energy = NaN
                rho = NaN
                mass = NaN
                try
                    # 使用 calculate_core 求解 3 个未知量
                    res = nlsolve(x -> calculate_core(x, mu, T, nodes1, omega), x; autodiff = :forward)
                    converged = res.f_converged
                    if converged
                        copyto!(x, res.zero)
                        x_prev .= x  # 用当前解作为下一点的初始值
                        # 计算热力学量（pressure, entropy, energy）
                        pressure, rho, entropy, energy = calculate_thermo(x, mu, T, nodes1, omega)
                        # 计算有效质量
                        mass = calculate_mass(x[1])  # x[1] 是 phi
                    else
                        @warn "Root finding did not converge for T=$T and mu=$mu"
                    end
                catch err
                    @warn "Exception in root finding for T=$T and mu=$mu: $err"
                    converged = false
                end

                # 写入 csv：将 T, mu 转换回物理单位（乘 hc）以与 Trho 保持一致
                println(io, join([T*hc, mu*hc, x..., mass, pressure, rho, entropy, energy, converged], ","))
                flush(io)
            end
        end
    end

    return nothing
end

# 示例调用（注释掉以避免自动运行）
Tmu(T_start=50/hc, T_end=350/hc, T_step=1/hc, mu_start=0/hc, mu_end=0.0, mu_step=-1/hc)
#Trho(10/hc,11/hc)

# ===================== find_cep 功能实现 =====================

function generate_single_temperature_data(T_target)
    """为单个温度生成完整的rho-mu数据"""
    nodes1, nodes2 = get_nodes(128, 16)
    omega = 100/hc
    
    # 初始猜测值
    x_initial = [-2.13, 0.06, 0.12, 310 / hc]
    x = copy(x_initial)
    
    rho_mu_pairs = []
    
    # 从高密度到低密度扫描（与Trho函数一致）
    for rho in [3.00; collect(2.99:-0.01:0.10)]
        converged = false
        try
            res = nlsolve(x -> calculate_t_rho(x, T_target, rho, nodes1, omega), x)
            converged = res.f_converged
            if converged
                copyto!(x, res.zero)
                mu_value = x[4]
                push!(rho_mu_pairs, (rho, mu_value))
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

function has_s_shape(T_mev; min_negative_points=3, derivative_threshold=-0.001)
    """
    检测给定温度下是否存在S形
    
    参数:
    - T_mev: 温度 (MeV)
    - min_negative_points: 至少需要多少个负导数点才认为存在S形
    - derivative_threshold: 负导数阈值
    
    返回: true表示存在S形，false表示不存在
    """
    T_julia = T_mev / hc  # 转换为Julia内部单位
    
    rho_mu_pairs = generate_single_temperature_data(T_julia)
    
    if length(rho_mu_pairs) < 3
        @debug "温度 $T_mev MeV 处数据点太少，无法判断S形"
        return false
    end
    
    # 按rho排序（从小到大）
    sorted_pairs = sort(rho_mu_pairs, by=x->x[1])
    
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

function find_cep(T_min=10.0, T_max=150.0, tolerance=0.5)
    """
    找到温度临界点T_cep，该温度下rho-mu关系正好不出现"S形"
    
    参数:
    - T_min: 搜索的最低温度 (MeV)
    - T_max: 搜索的最高温度 (MeV) 
    - tolerance: 温度搜索精度 (MeV)
    
    返回:
    - T_cep: 临界温度 (MeV)
    """
    
    println("开始寻找临界温度 T_cep...")
    println("搜索范围: [$T_min, $T_max] MeV")
    println("精度: $tolerance MeV")
    
    # 首先扫描整个范围，寻找S形存在的区域
    println("扫描温度范围寻找S形现象...")
    
    scan_temps = collect(T_min:5:T_max)
    s_shape_results = []
    
    for T_scan in scan_temps
        println("扫描 T=$T_scan MeV...")
        has_s = has_s_shape(T_scan)
        push!(s_shape_results, (T_scan, has_s))
        println("  $(has_s ? "有S形" : "无S形")")
    end
    
    # 寻找相变区间
    transition_found = false
    T_low = T_min
    T_high = T_max
    
    for i in 1:(length(s_shape_results)-1)
        T1, has_s1 = s_shape_results[i]
        T2, has_s2 = s_shape_results[i+1]
        
        if has_s1 != has_s2  # 找到相变
            T_low = has_s1 ? T2 : T1  # 有S形的温度作为下界
            T_high = has_s1 ? T1 : T2  # 无S形的温度作为上界
            transition_found = true
            println("发现相变区间: [$T_low, $T_high] MeV")
            break
        end
    end
    
    if !transition_found
        println("在扫描范围内未发现S形相变，尝试扩展搜索...")
        
        # 分析初始扫描的结果
        all_have_s = all(result[2] for result in s_shape_results)
        all_no_s = all(!result[2] for result in s_shape_results)
        
        if all_have_s
            # 整个范围都有S形，说明临界点在更高温度
            println("初始范围内都有S形，向高温扩展搜索...")
            for T_test in (T_max+5):10:200
                println("测试高温 T=$T_test MeV...")
                current_has_s = has_s_shape(T_test)
                
                if !current_has_s
                    println("  无S形")
                    # 找到S形消失的点，临界点在T_max和T_test之间
                    T_low = T_max
                    T_high = T_test
                    transition_found = true
                    println("在高温区发现相变区间: [$T_low, $T_high] MeV")
                    break
                else
                    println("  有S形")
                end
            end
        elseif all_no_s
            # 整个范围都无S形，说明临界点在更低温度
            println("初始范围内都无S形，向低温扩展搜索...")
            for T_test in (T_min-5):-10:10
                println("测试低温 T=$T_test MeV...")
                current_has_s = has_s_shape(T_test)
                
                if current_has_s
                    println("  有S形")
                    # 找到S形出现的点，临界点在T_test和T_min之间
                    T_low = T_test
                    T_high = T_min
                    transition_found = true
                    println("在低温区发现相变区间: [$T_low, $T_high] MeV")
                    break
                else
                    println("  无S形")
                end
            end
        end
    end
    
    if !transition_found
        @warn "未找到S形相变，返回搜索范围中点作为估计值"
        return (T_min + T_max) / 2
    end
    
    # 二分搜索精确定位临界点
    println("开始二分搜索...")
    
    iteration = 0
    while (T_high - T_low) > tolerance && iteration < 15
        iteration += 1
        T_mid = (T_low + T_high) / 2
        
        println("第 $iteration 次迭代: 测试 T=$T_mid MeV")
        
        mid_has_s = has_s_shape(T_mid)
        println("  结果: $(mid_has_s ? "有S形" : "无S形")")
        
        if mid_has_s
            T_low = T_mid  # S形存在，临界点在更高温度
        else
            T_high = T_mid  # S形不存在，临界点在更低温度
        end
        
        println("  新搜索区间: [$T_low, $T_high]")
    end
    
    T_cep = (T_low + T_high) / 2
    println("找到临界温度: T_cep = $T_cep MeV")
    
    return T_cep
end

# 便捷的测试函数
function test_find_cep()
    """测试find_cep功能的便捷函数"""
    println("测试find_cep功能...")
    println("注意: 如果当前参数设置下不存在S形相变，函数将返回NaN")
    
    # 使用较宽松的参数进行测试
    T_cep = find_cep(30.0, 100.0, 1.0)
    
    if isnan(T_cep)
        println("未找到临界温度，建议:")
        println("1. 调整物理参数 (omega, 初始条件等)")
        println("2. 扩大搜索范围")
        println("3. 调整S形检测敏感度")
    else
        println("找到临界温度: T_cep = $T_cep MeV")
    end
    
    return T_cep
end

"""
phi = -2.13
Phi1 = 0.06
Phi2 = 0.13
T=100/hc
mu = 308.3/hc
nodes1, nodes2 = get_nodes(128, 16)
omega = 100/hc
omega = 0.0
@show calculate_pressure(phi, Phi1, Phi2, mu, T, nodes1, omega)
x = SVector(phi, Phi1, Phi2)
@show calculate_core(x, mu, T, nodes1, omega)
@show calculate_thermo(x, mu, T, nodes1, omega)
x=SVector(phi, Phi1, Phi2, mu)
@show calculate_t_rho(x, T, 4.0, nodes1, omega)
"""