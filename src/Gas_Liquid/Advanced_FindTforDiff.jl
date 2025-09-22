# Advanced_FindTforDiff.jl
# 基于ForwardDiff自动微分的温度反向查找模块
# 给定固定的κ₃/κ₁和κ₄/κ₂值，寻找对应的温度
# 依赖于 Advanced_ForwardDiff.jl 中的扫描函数
# 使用简单的滑动数组检测和线性插值，无需外部依赖

include("Advanced_ForwardDiff.jl")
#include("Constants_Gas_Liquid.jl")
using .Constants_Gas_Liquid
using DataFrames, Statistics, Dates

function find_temperature_for_kappa_ratios(target_kappa3_kappa1, target_kappa4_kappa2, μ_B, 
                                          T_min, T_max, T_step_scan=1.0/hc;
                                          gsigma=1.25, gdelta=0.01,
                                          fs=15.0, fo=5.423, fr=0.95, fd=0.0,
                                          b=0.00692, c=-0.0048, n_nodes=256,
                                          verbose=true)
    """
    寻找给定κ₃/κ₁和κ₄/κ₂值对应的温度
    
    参数:
    - target_kappa3_kappa1: 目标κ₃/κ₁值
    - target_kappa4_kappa2: 目标κ₄/κ₂值
    - μ_B: 重子化学势
    - T_min: 温度搜索下限
    - T_max: 温度搜索上限
    - T_step_scan: 温度扫描步长 (用于初始扫描)
    - 其他参数: 模型参数
    - verbose: 是否打印详细信息
    
    返回:
    - (T_kappa3_kappa1, T_kappa4_kappa2): 对应的温度值 (取较大交点)
    """
    
    if verbose
        println("="^70)
        println("温度反向查找：基于ForwardDiff自动微分")
        println("="^70)
        println("目标值:")
        println("  κ₃/κ₁ = $target_kappa3_kappa1")
        println("  κ₄/κ₂ = $target_kappa4_kappa2")
        println("  μ_B = $(μ_B*hc) MeV")
        println("搜索范围:")
        println("  T: $(T_min*hc) - $(T_max*hc) MeV")
        println("  扫描步长: $(T_step_scan*hc) MeV")
    end
    
    # 设置模型参数
    nodes = get_nodes(n_nodes)
    couplings = [fs, fo, fr, fd, b, c]
    model_params = (nodes, couplings)
    
    # 第一步：粗扫描获取温度-κ比值关系
    if verbose
        println("\n第一步：执行温度扫描获取κ比值数据...")
    end
    
    T_array = T_min:T_step_scan:T_max
    kappa3_kappa1_array = Float64[]
    kappa4_kappa2_array = Float64[]
    T_valid_array = Float64[]
    
    for (i, T) in enumerate(T_array)
        if verbose && i % 20 == 1
            println("  扫描进度: $(i)/$(length(T_array)), T = $(round(T*hc, digits=1)) MeV")
        end
        
        try
            # 计算热力学涨落量
            κ1, κ2, κ3, κ4, κ3_κ1, κ4_κ2 = calculate_forwarddiff_thermodynamic_fluctuations(
                gsigma, gdelta, T, μ_B, model_params)
            
            # 只保存有效的数据点
            if isfinite(κ3_κ1) && isfinite(κ4_κ2)
                push!(T_valid_array, T)
                push!(kappa3_kappa1_array, κ3_κ1)
                push!(kappa4_kappa2_array, κ4_κ2)
            end
            
        catch e
            if verbose && i % 50 == 1
                println("    T = $(round(T*hc, digits=1)) MeV 计算失败: $e")
            end
        end
    end
    
    if length(T_valid_array) < 5
        error("有效数据点太少 ($(length(T_valid_array))个)，无法进行插值")
    end
    
    if verbose
        println("  扫描完成，获得 $(length(T_valid_array)) 个有效数据点")
        println("  κ₃/κ₁ 范围: $(round(minimum(kappa3_kappa1_array), digits=3)) - $(round(maximum(kappa3_kappa1_array), digits=3))")
        println("  κ₄/κ₂ 范围: $(round(minimum(kappa4_kappa2_array), digits=3)) - $(round(maximum(kappa4_kappa2_array), digits=3))")
    end
    
    # 第二步：寻找κ₃/κ₁对应的温度
    if verbose
        println("\n第二步：寻找κ₃/κ₁ = $target_kappa3_kappa1 对应的温度...")
    end
    
    T_kappa3_kappa1 = find_temperature_for_single_kappa_ratio(
        T_valid_array, kappa3_kappa1_array, target_kappa3_kappa1, "κ₃/κ₁", verbose)
    
    # 第三步：寻找κ₄/κ₂对应的温度
    if verbose
        println("\n第三步：寻找κ₄/κ₂ = $target_kappa4_kappa2 对应的温度...")
    end
    
    T_kappa4_kappa2 = find_temperature_for_single_kappa_ratio(
        T_valid_array, kappa4_kappa2_array, target_kappa4_kappa2, "κ₄/κ₂", verbose)
    
    # 结果总结
    if verbose
        println("\n" * "="^70)
        println("温度查找结果:")
        println("="^70)
        if !isnan(T_kappa3_kappa1)
            println("κ₃/κ₁ = $target_kappa3_kappa1 对应温度: $(round(T_kappa3_kappa1*hc, digits=2)) MeV")
        else
            println("κ₃/κ₁ = $target_kappa3_kappa1 未找到对应温度")
        end
        
        if !isnan(T_kappa4_kappa2)
            println("κ₄/κ₂ = $target_kappa4_kappa2 对应温度: $(round(T_kappa4_kappa2*hc, digits=2)) MeV")
        else
            println("κ₄/κ₂ = $target_kappa4_kappa2 未找到对应温度")
        end
    end
    
    return T_kappa3_kappa1, T_kappa4_kappa2
end

function find_temperature_for_kappa_ratios_with_optimization_params(
    target_kappa3_kappa1, target_kappa4_kappa2, μ_B, optimization_params,
    T_min, T_max, T_step_scan=1.0/hc;
    gsigma=1.25, gdelta=0.01, n_nodes=256, verbose=true)
    """
    通过优化参数计算给定κ₃/κ₁和κ₄/κ₂值对应的温度
    
    参数:
    - target_kappa3_kappa1: 目标κ₃/κ₁值
    - target_kappa4_kappa2: 目标κ₄/κ₂值
    - μ_B: 重子化学势
    - optimization_params: 优化参数元组 (ρ₀, B_A, K, m_ratio, E_sym)
        * ρ₀: 核饱和密度 (fm⁻³)
        * B_A: 结合能 (MeV)
        * K: 不可压缩模量 (MeV)
        * m_ratio: 有效质量比
        * E_sym: 对称能 (MeV)
    - T_min: 温度搜索下限
    - T_max: 温度搜索上限
    - T_step_scan: 温度扫描步长
    - gsigma, gdelta: 强子-夸克耦合参数
    - n_nodes: 积分节点数
    - verbose: 是否打印详细信息
    
    返回:
    - (T_kappa3_kappa1, T_kappa4_kappa2): 对应的温度值
    """
    
    # 解包优化参数
    ρ₀, B_A, K, m_ratio, E_sym = optimization_params
    
    if verbose
        println("="^80)
        println("基于优化参数的温度反向查找")
        println("="^80)
        println("优化参数:")
        println("  ρ₀ = $ρ₀ fm⁻³")
        println("  B_A = $B_A MeV")
        println("  K = $K MeV") 
        println("  m_ratio = $m_ratio")
        println("  E_sym = $E_sym MeV")
        println("目标κ比值:")
        println("  κ₃/κ₁ = $target_kappa3_kappa1")
        println("  κ₄/κ₂ = $target_kappa4_kappa2")
        println("  μ_B = $(μ_B*hc) MeV")
    end
    
    # 将优化参数转换为无量纲单位（除以hc）
    ρ₀_dimensionless = ρ₀
    B_A_dimensionless = B_A / hc
    K_dimensionless = K / hc  
    E_sym_dimensionless = E_sym / hc
    
    # 计算耦合常数
    if verbose
        println("\n计算耦合常数...")
    end
    
    try
        fσ, fω, fρ, fδ, b, c = calculate_couplings(
            ρ₀_dimensionless, B_A_dimensionless, K_dimensionless, m_ratio, E_sym_dimensionless)
        
        if verbose
            println("  耦合常数计算结果:")
            println("    fσ = $(round(fσ, digits=6))")
            println("    fω = $(round(fω, digits=6))")
            println("    fρ = $(round(fρ, digits=6))")
            println("    fδ = $(round(fδ, digits=6))")
            println("    b = $(round(b, digits=6))")
            println("    c = $(round(c, digits=6))")
        end
        
        # 调用原始的温度查找函数
        T_kappa3_kappa1, T_kappa4_kappa2 = find_temperature_for_kappa_ratios(
            target_kappa3_kappa1, target_kappa4_kappa2, μ_B, T_min, T_max, T_step_scan;
            gsigma=gsigma, gdelta=gdelta, fs=fσ, fo=fω, fr=fρ, fd=fδ,
            b=b, c=c, n_nodes=n_nodes, verbose=verbose)
        
        return T_kappa3_kappa1, T_kappa4_kappa2
        
    catch e
        if verbose
            println("✗ 耦合常数计算失败: $e")
        end
        return NaN, NaN
    end
end

function batch_find_temperatures_with_optimization_params(
    kappa_pairs, μ_B, optimization_params, T_min, T_max;
    T_step_scan=1.0/hc, gsigma=1.25, gdelta=0.01, n_nodes=256,
    output_file="output/Gas_Liquid/temperature_optimization_results.csv")
    """
    批量计算多组κ比值在给定优化参数下对应的温度
    
    参数:
    - kappa_pairs: κ比值对数组，格式 [(κ₃/κ₁, κ₄/κ₂), ...]
    - μ_B: 重子化学势
    - optimization_params: 优化参数元组 (ρ₀, B_A, K, m_ratio, E_sym)
    - T_min, T_max: 温度搜索范围
    - T_step_scan: 扫描步长
    - output_file: 输出文件路径
    - 其他参数: 模型参数
    
    返回:
    - DataFrame: 包含所有结果的数据框
    """
    
    # 解包优化参数
    ρ₀, B_A, K, m_ratio, E_sym = optimization_params
    
    println("="^80)
    println("基于优化参数的批量温度查找")
    println("="^80)
    println("优化参数: ρ₀=$ρ₀, B_A=$B_A MeV, K=$K MeV, m_ratio=$m_ratio, E_sym=$E_sym MeV")
    println("μ_B = $(μ_B*hc) MeV，共 $(length(kappa_pairs)) 组κ比值")
    println("温度搜索范围: $(T_min*hc) - $(T_max*hc) MeV")
    
    results = []
    
    for (i, (kappa3_kappa1, kappa4_kappa2)) in enumerate(kappa_pairs)
        println("\n处理第 $i/$(length(kappa_pairs)) 组: κ₃/κ₁ = $kappa3_kappa1, κ₄/κ₂ = $kappa4_kappa2")
        
        try
            T_k3k1, T_k4k2 = find_temperature_for_kappa_ratios_with_optimization_params(
                kappa3_kappa1, kappa4_kappa2, μ_B, optimization_params,
                T_min, T_max, T_step_scan;
                gsigma=gsigma, gdelta=gdelta, n_nodes=n_nodes, verbose=false)
            
            result_row = (
                kappa3_over_kappa1_target = kappa3_kappa1,
                kappa4_over_kappa2_target = kappa4_kappa2,
                T_for_kappa3_kappa1_MeV = isnan(T_k3k1) ? NaN : T_k3k1 * hc,
                T_for_kappa4_kappa2_MeV = isnan(T_k4k2) ? NaN : T_k4k2 * hc,
                temperature_difference_MeV = isnan(T_k3k1) || isnan(T_k4k2) ? NaN : abs(T_k3k1 - T_k4k2) * hc,
                mu_B_MeV = μ_B * hc,
                rho0 = ρ₀,
                B_A_MeV = B_A,
                K_MeV = K,
                m_ratio = m_ratio,
                E_sym_MeV = E_sym,
                status_kappa3_kappa1 = isnan(T_k3k1) ? "未找到" : "找到",
                status_kappa4_kappa2 = isnan(T_k4k2) ? "未找到" : "找到"
            )
            
            push!(results, result_row)
            
            temp_diff = isnan(T_k3k1) || isnan(T_k4k2) ? "N/A" : round(abs(T_k3k1 - T_k4k2) * hc, digits=2)
            println("  结果: T(κ₃/κ₁) = $(isnan(T_k3k1) ? "未找到" : round(T_k3k1*hc, digits=2))" * 
                   " MeV, T(κ₄/κ₂) = $(isnan(T_k4k2) ? "未找到" : round(T_k4k2*hc, digits=2)) MeV" *
                   ", |ΔT| = $temp_diff MeV")
            
        catch e
            println("  计算失败: $e")
            
            result_row = (
                kappa3_over_kappa1_target = kappa3_kappa1,
                kappa4_over_kappa2_target = kappa4_kappa2,
                T_for_kappa3_kappa1_MeV = NaN,
                T_for_kappa4_kappa2_MeV = NaN,
                temperature_difference_MeV = NaN,
                mu_B_MeV = μ_B * hc,
                rho0 = ρ₀,
                B_A_MeV = B_A,
                K_MeV = K,
                m_ratio = m_ratio,
                E_sym_MeV = E_sym,
                status_kappa3_kappa1 = "计算失败",
                status_kappa4_kappa2 = "计算失败"
            )
            
            push!(results, result_row)
        end
    end
    
    # 保存结果
    save_optimization_temperature_results(results, optimization_params, μ_B, T_min, T_max, 
                                        T_step_scan, output_file, gsigma, gdelta, n_nodes)
    
    return DataFrame(results)
end

function save_optimization_temperature_results(results, optimization_params, μ_B, T_min, T_max, 
                                             T_step_scan, output_file, gsigma, gdelta, n_nodes)
    """
    保存基于优化参数的温度查找结果到CSV文件
    """
    println("\n" * "="^60)
    println("保存优化参数温度查找结果")
    println("="^60)
    
    df = DataFrame(results)
    ρ₀, B_A, K, m_ratio, E_sym = optimization_params
    
    # 确保输出目录存在
    output_dir = dirname(output_file)
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    
    try
        open(output_file, "w") do io
            # 写入元数据头部
            println(io, "# Temperature Finder Results with Optimization Parameters")
            println(io, "# Generated on: $(Dates.now())")
            println(io, "# Optimization Parameters:")
            println(io, "# rho0 = $ρ₀ fm^-3")
            println(io, "# B_A = $B_A MeV")
            println(io, "# K = $K MeV")
            println(io, "# m_ratio = $m_ratio")
            println(io, "# E_sym = $E_sym MeV")
            println(io, "# Physical Parameters:")
            println(io, "# gsigma = $gsigma")
            println(io, "# gdelta = $gdelta")
            println(io, "# mu_B = $(μ_B*hc) MeV")
            println(io, "# T_search_range = $(T_min*hc) - $(T_max*hc) MeV")
            println(io, "# T_scan_step = $(T_step_scan*hc) MeV")
            println(io, "# nodes = $n_nodes")
            println(io, "#")
            
            # 写入CSV数据
            col_names = names(df)
            println(io, join(col_names, ","))
            
            for row in eachrow(df)
                values = [string(row[col]) for col in col_names]
                println(io, join(values, ","))
            end
        end
        
        println("✓ 结果已保存到: $output_file")
        
        # 显示统计信息
        successful_k3k1 = sum(df.status_kappa3_kappa1 .== "找到")
        successful_k4k2 = sum(df.status_kappa4_kappa2 .== "找到")
        both_successful = sum((df.status_kappa3_kappa1 .== "找到") .& (df.status_kappa4_kappa2 .== "找到"))
        total = nrow(df)
        
        println("\n统计信息:")
        println("总计算组数: $total")
        println("κ₃/κ₁ 成功找到温度: $successful_k3k1/$total")
        println("κ₄/κ₂ 成功找到温度: $successful_k4k2/$total")
        println("两个κ比值都找到温度: $both_successful/$total")
        
        if both_successful > 0
            valid_diffs = filter(!isnan, df.temperature_difference_MeV)
            if !isempty(valid_diffs)
                println("温度差统计 (|T_κ₃/κ₁ - T_κ₄/κ₂|):")
                println("  平均值: $(round(mean(valid_diffs), digits=2)) MeV")
                println("  中位数: $(round(median(valid_diffs), digits=2)) MeV")
                println("  最小值: $(round(minimum(valid_diffs), digits=2)) MeV")
                println("  最大值: $(round(maximum(valid_diffs), digits=2)) MeV")
            end
        end
        
    catch e
        println("✗ 保存文件失败: $e")
    end
end

function find_temperature_for_single_kappa_ratio(T_array, kappa_array, target_kappa, kappa_name, verbose=true)
    """
    为单个κ比值寻找对应的温度
    使用滑动数组检测峰值，并用线性插值找到交点
    
    参数:
    - T_array: 温度数组
    - kappa_array: κ比值数组
    - target_kappa: 目标κ比值
    - kappa_name: κ比值名称（用于打印）
    - verbose: 是否打印详细信息
    
    返回:
    - T_found: 找到的温度（取较大交点），如果未找到则返回NaN
    """
    
    # 检查目标值是否在数据范围内
    kappa_min, kappa_max = extrema(kappa_array)
    if target_kappa < kappa_min || target_kappa > kappa_max
        if verbose
            println("  ⚠️  目标值 $target_kappa 超出 $kappa_name 数据范围 [$kappa_min, $kappa_max]")
        end
        return NaN
    end
    
    # 寻找所有可能的交点
    crossing_points = Float64[]
    
    # 遍历相邻点对，寻找交点
    for i in 1:(length(kappa_array)-1)
        y1, y2 = kappa_array[i], kappa_array[i+1]
        T1, T2 = T_array[i], T_array[i+1]
        
        # 检查目标值是否在这两点之间
        if (y1 <= target_kappa <= y2) || (y2 <= target_kappa <= y1)
            # 线性插值找交点
            if abs(y2 - y1) > 1e-12  # 避免除零
                # 线性插值公式: T = T1 + (target_kappa - y1) * (T2 - T1) / (y2 - y1)
                T_interp = T1 + (target_kappa - y1) * (T2 - T1) / (y2 - y1)
                push!(crossing_points, T_interp)
                
                if verbose
                    println("    找到交点: T = $(round(T_interp*hc, digits=2)) MeV (在 $(round(T1*hc, digits=1)) - $(round(T2*hc, digits=1)) MeV 之间)")
                end
            end
        end
    end
    
    if isempty(crossing_points)
        if verbose
            println("  ❌ 未找到 $kappa_name = $target_kappa 的交点")
        end
        return NaN
    end
    
    # 如果有多个交点，寻找右侧峰值附近的交点（较大温度）
    if length(crossing_points) > 1
        if verbose
            println("  📍 找到 $(length(crossing_points)) 个交点，选择右侧交点")
        end
        
        # 寻找κ数组的峰值位置
        peak_idx = find_peak_index(kappa_array)
        peak_T = T_array[peak_idx]
        
        if verbose
            println("    峰值位置: T = $(round(peak_T*hc, digits=1)) MeV")
        end
        
        # 选择峰值右侧（较大温度）的交点
        right_crossings = filter(T -> T > peak_T, crossing_points)
        
        if !isempty(right_crossings)
            T_found = minimum(right_crossings)  # 峰值右侧最近的交点
            if verbose
                println("  ✅ 选择右侧交点: T = $(round(T_found*hc, digits=2)) MeV")
            end
        else
            # 如果峰值右侧没有交点，选择最大的交点
            T_found = maximum(crossing_points)
            if verbose
                println("  ⚠️  峰值右侧无交点，选择最大交点: T = $(round(T_found*hc, digits=2)) MeV")
            end
        end
    else
        T_found = crossing_points[1]
        if verbose
            println("  ✅ 找到唯一交点: T = $(round(T_found*hc, digits=2)) MeV")
        end
    end
    
    return T_found
end

function find_peak_index(data_array)
    """
    使用滑动窗口找到数组的峰值位置
    
    参数:
    - data_array: 数据数组
    
    返回:
    - peak_idx: 峰值位置索引
    """
    n = length(data_array)
    if n < 3
        return div(n, 2) + 1  # 数组太短，返回中点
    end
    
    # 寻找局部最大值
    max_val = -Inf
    peak_idx = 1
    
    for i in 2:(n-1)
        # 检查是否为局部最大值（比左右邻居都大）
        if data_array[i] > data_array[i-1] && data_array[i] > data_array[i+1]
            if data_array[i] > max_val
                max_val = data_array[i]
                peak_idx = i
            end
        end
    end
    
    # 如果没有找到局部最大值，找全局最大值
    if max_val == -Inf
        peak_idx = argmax(data_array)
    end
    
    return peak_idx
end

function batch_find_temperatures_for_kappa_ratios(kappa_pairs, μ_B, T_min, T_max;
                                                  T_step_scan=1.0/hc,
                                                  output_file="output/Gas_Liquid/temperature_finder_results.csv",
                                                  gsigma=1.25, gdelta=0.01,
                                                  fs=17.28476, fo=11.66174, fr=0.89363, fd=0.0,
                                                  b=0.00210, c=-0.00297, n_nodes=256)
    """
    批量寻找多组κ比值对应的温度
    
    参数:
    - kappa_pairs: κ比值对数组，格式 [(κ₃/κ₁, κ₄/κ₂), ...]
    - μ_B: 重子化学势
    - T_min, T_max: 温度搜索范围
    - T_step_scan: 扫描步长
    - output_file: 输出文件路径
    - 其他参数: 模型参数
    
    返回:
    - DataFrame: 包含所有结果的数据框
    """
    
    println("="^80)
    println("批量温度反向查找")
    println("="^80)
    println("μ_B = $(μ_B*hc) MeV，共 $(length(kappa_pairs)) 组κ比值")
    println("温度搜索范围: $(T_min*hc) - $(T_max*hc) MeV")
    
    results = []
    
    for (i, (kappa3_kappa1, kappa4_kappa2)) in enumerate(kappa_pairs)
        println("\n处理第 $i/$(length(kappa_pairs)) 组: κ₃/κ₁ = $kappa3_kappa1, κ₄/κ₂ = $kappa4_kappa2")
        
        try
            T_k3k1, T_k4k2 = find_temperature_for_kappa_ratios(
                kappa3_kappa1, kappa4_kappa2, μ_B, T_min, T_max, T_step_scan;
                gsigma=gsigma, gdelta=gdelta, fs=fs, fo=fo, fr=fr, fd=fd,
                b=b, c=c, n_nodes=n_nodes, verbose=false)
            
            result_row = (
                kappa3_over_kappa1_target = kappa3_kappa1,
                kappa4_over_kappa2_target = kappa4_kappa2,
                T_for_kappa3_kappa1_MeV = isnan(T_k3k1) ? NaN : T_k3k1 * hc,
                T_for_kappa4_kappa2_MeV = isnan(T_k4k2) ? NaN : T_k4k2 * hc,
                mu_B_MeV = μ_B * hc,
                status_kappa3_kappa1 = isnan(T_k3k1) ? "未找到" : "找到",
                status_kappa4_kappa2 = isnan(T_k4k2) ? "未找到" : "找到"
            )
            
            push!(results, result_row)
            
            println("  结果: T(κ₃/κ₁) = $(isnan(T_k3k1) ? "未找到" : round(T_k3k1*hc, digits=2))" * 
                   " MeV, T(κ₄/κ₂) = $(isnan(T_k4k2) ? "未找到" : round(T_k4k2*hc, digits=2)) MeV")
            
        catch e
            println("  计算失败: $e")
            
            result_row = (
                kappa3_over_kappa1_target = kappa3_kappa1,
                kappa4_over_kappa2_target = kappa4_kappa2,
                T_for_kappa3_kappa1_MeV = NaN,
                T_for_kappa4_kappa2_MeV = NaN,
                mu_B_MeV = μ_B * hc,
                status_kappa3_kappa1 = "计算失败",
                status_kappa4_kappa2 = "计算失败"
            )
            
            push!(results, result_row)
        end
    end
    
    # 保存结果
    save_temperature_finder_results(results, μ_B, T_min, T_max, T_step_scan, output_file,
                                   gsigma, gdelta, fs, fo, fr, fd, b, c, n_nodes)
    
    return DataFrame(results)
end

function save_temperature_finder_results(results, μ_B, T_min, T_max, T_step_scan, output_file,
                                        gsigma, gdelta, fs, fo, fr, fd, b, c, n_nodes)
    """
    保存温度查找结果到CSV文件
    """
    println("\n" * "="^60)
    println("保存温度查找结果")
    println("="^60)
    
    df = DataFrame(results)
    
    # 确保输出目录存在
    output_dir = dirname(output_file)
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    
    try
        open(output_file, "w") do io
            # 写入元数据头部
            println(io, "# Temperature Finder Results (ForwardDiff)")
            println(io, "# Generated on: $(Dates.now())")
            println(io, "# Model Parameters:")
            println(io, "# gsigma = $gsigma")
            println(io, "# gdelta = $gdelta")
            println(io, "# fs = $fs")
            println(io, "# fo = $fo")
            println(io, "# fr = $fr")
            println(io, "# fd = $fd")
            println(io, "# b = $b")
            println(io, "# c = $c")
            println(io, "# mu_B = $(μ_B*hc) MeV")
            println(io, "# T_search_range = $(T_min*hc) - $(T_max*hc) MeV")
            println(io, "# T_scan_step = $(T_step_scan*hc) MeV")
            println(io, "# nodes = $n_nodes")
            println(io, "#")
            
            # 写入CSV数据
            col_names = names(df)
            println(io, join(col_names, ","))
            
            for row in eachrow(df)
                values = [string(row[col]) for col in col_names]
                println(io, join(values, ","))
            end
        end
        
        println("✓ 结果已保存到: $output_file")
        
        # 显示统计信息
        successful_k3k1 = sum(df.status_kappa3_kappa1 .== "找到")
        successful_k4k2 = sum(df.status_kappa4_kappa2 .== "找到")
        total = nrow(df)
        
        println("\n统计信息:")
        println("总计算组数: $total")
        println("κ₃/κ₁ 成功找到温度: $successful_k3k1/$total")
        println("κ₄/κ₂ 成功找到温度: $successful_k4k2/$total")
        
    catch e
        println("✗ 保存文件失败: $e")
    end
end

function calculate_temperature_difference_sum_of_squares(
    kappa_pairs, μ_B_values, optimization_params, T_min, T_max;
    T_step_scan=1.0/hc, gsigma=1.25, gdelta=0.01, n_nodes=256,
    penalty_for_missing=1e6, verbose=false)
    """
    计算多组κ比值在给定优化参数下对应的温度差的平方和
    
    此函数基于 batch_find_temperatures_with_optimization_params 的逻辑，
    但不保存文件，直接返回优化目标值（所有组温度差的平方和）。
    适用于参数优化过程中的目标函数计算。
    
    参数:
    - kappa_pairs: κ比值对数组，格式 [(κ₃/κ₁, κ₄/κ₂), ...]
    - μ_B_values: 重子化学势数组，与 kappa_pairs 一一对应，格式 [μ_B1, μ_B2, ...]
    - optimization_params: 优化参数元组 (ρ₀, B_A, K, m_ratio, E_sym)
        * ρ₀: 核饱和密度 (fm⁻³)
        * B_A: 结合能 (MeV)
        * K: 不可压缩模量 (MeV)
        * m_ratio: 有效质量比
        * E_sym: 对称能 (MeV)
    - T_min, T_max: 温度搜索范围
    - T_step_scan: 扫描步长
    - gsigma, gdelta: 强子-夸克耦合参数
    - n_nodes: 积分节点数
    - penalty_for_missing: 当无法找到温度或计算失败时的惩罚值
    - verbose: 是否打印详细信息（建议在优化过程中设为false）
    
    返回:
    - sum_of_squares: 所有组温度差平方和 (单位: MeV²)
      对于每组 (κ₃/κ₁, κ₄/κ₂)，计算对应的温度 T₁, T₂，
      贡献 |T₁ - T₂|² 到总和中。
      如果某组计算失败或找不到温度，则贡献 penalty_for_missing。
    """
    
    # 验证输入参数
    if length(μ_B_values) != length(kappa_pairs)
        error("μ_B值数组长度 ($(length(μ_B_values))) 与κ比值对数组长度 ($(length(kappa_pairs))) 不匹配")
    end
    
    # 解包优化参数
    ρ₀, B_A, K, m_ratio, E_sym = optimization_params
    
    if verbose
        println("="^80)
        println("基于优化参数的温度差平方和计算")
        println("="^80)
        println("优化参数: ρ₀=$ρ₀, B_A=$B_A MeV, K=$K MeV, m_ratio=$m_ratio, E_sym=$E_sym MeV")
        println("共 $(length(kappa_pairs)) 组κ比值，每组对应不同的μ_B值")
        println("μ_B范围: $(round(minimum(μ_B_values)*hc, digits=1)) - $(round(maximum(μ_B_values)*hc, digits=1)) MeV")
        println("温度搜索范围: $(T_min*hc) - $(T_max*hc) MeV")
        println("无效值惩罚: $penalty_for_missing")
    end
    
    sum_of_squares = 0.0
    valid_pairs = 0
    failed_pairs = 0
    
    for (i, ((kappa3_kappa1, kappa4_kappa2), μ_B)) in enumerate(zip(kappa_pairs, μ_B_values))
        if verbose
            println("\n处理第 $i/$(length(kappa_pairs)) 组: κ₃/κ₁ = $kappa3_kappa1, κ₄/κ₂ = $kappa4_kappa2, μ_B = $(round(μ_B*hc, digits=1)) MeV")
        end
        
        try
            T_k3k1, T_k4k2 = find_temperature_for_kappa_ratios_with_optimization_params(
                kappa3_kappa1, kappa4_kappa2, μ_B, optimization_params,
                T_min, T_max, T_step_scan;
                gsigma=gsigma, gdelta=gdelta, n_nodes=n_nodes, verbose=false)
            
            # 检查是否成功找到两个温度
            if !isnan(T_k3k1) && !isnan(T_k4k2) && isfinite(T_k3k1) && isfinite(T_k4k2)
                # 计算温度差的平方 (转换为MeV单位)
                temp_diff_MeV = abs(T_k3k1 - T_k4k2) * hc
                contribution = temp_diff_MeV^2
                sum_of_squares += contribution
                valid_pairs += 1
                
                if verbose
                    println("  ✓ T(κ₃/κ₁) = $(round(T_k3k1*hc, digits=2)) MeV")
                    println("    T(κ₄/κ₂) = $(round(T_k4k2*hc, digits=2)) MeV")
                    println("    |ΔT| = $(round(temp_diff_MeV, digits=2)) MeV")
                    println("    贡献: $(round(contribution, digits=2)) MeV²")
                end
            else
                # 无法找到有效温度，应用惩罚
                sum_of_squares += penalty_for_missing
                failed_pairs += 1
                
                if verbose
                    t1_status = isnan(T_k3k1) ? "未找到" : "$(round(T_k3k1*hc, digits=2)) MeV"
                    t2_status = isnan(T_k4k2) ? "未找到" : "$(round(T_k4k2*hc, digits=2)) MeV"
                    println("  ✗ T(κ₃/κ₁) = $t1_status, T(κ₄/κ₂) = $t2_status")
                    println("    应用惩罚: $penalty_for_missing")
                end
            end
            
        catch e
            # 计算过程出错，应用惩罚
            sum_of_squares += penalty_for_missing
            failed_pairs += 1
            
            if verbose
                println("  ✗ 计算失败: $e")
                println("    应用惩罚: $penalty_for_missing")
            end
        end
    end
    
    if verbose
        println("\n" * "="^60)
        println("计算完成")
        println("="^60)
        println("总组数: $(length(kappa_pairs))")
        println("成功计算: $valid_pairs")
        println("失败/无效: $failed_pairs") 
        println("温度差平方和: $(round(sum_of_squares, digits=2)) MeV²")
        if valid_pairs > 0
            avg_contribution = (sum_of_squares - failed_pairs * penalty_for_missing) / valid_pairs
            println("平均有效贡献: $(round(avg_contribution, digits=2)) MeV²")
        end
    end
    
    return sum_of_squares
end

function calculate_temperature_difference_sum_of_squares_with_weights(
    kappa_pairs, weights, μ_B_values, optimization_params, T_min, T_max;
    T_step_scan=1.0/hc, gsigma=1.25, gdelta=0.01, n_nodes=256,
    penalty_for_missing=1e6, verbose=false)
    """
    计算多组κ比值在给定优化参数下对应的温度差的加权平方和
    
    参数:
    - kappa_pairs: κ比值对数组，格式 [(κ₃/κ₁, κ₄/κ₂), ...]
    - weights: 权重数组，与 kappa_pairs 长度相同
    - μ_B_values: 重子化学势数组，与 kappa_pairs 一一对应
    - 其他参数: 与 calculate_temperature_difference_sum_of_squares 相同
    
    返回:
    - weighted_sum_of_squares: 加权温度差平方和 (单位: MeV²)
    """
    
    if length(weights) != length(kappa_pairs)
        error("权重数组长度 ($(length(weights))) 与κ比值对数组长度 ($(length(kappa_pairs))) 不匹配")
    end
    
    if length(μ_B_values) != length(kappa_pairs)
        error("μ_B值数组长度 ($(length(μ_B_values))) 与κ比值对数组长度 ($(length(kappa_pairs))) 不匹配")
    end
    
    # 解包优化参数
    ρ₀, B_A, K, m_ratio, E_sym = optimization_params
    
    if verbose
        println("="^80)
        println("基于优化参数的加权温度差平方和计算")
        println("="^80)
        println("优化参数: ρ₀=$ρ₀, B_A=$B_A MeV, K=$K MeV, m_ratio=$m_ratio, E_sym=$E_sym MeV")
        println("共 $(length(kappa_pairs)) 组κ比值，每组对应不同的μ_B值")
        println("μ_B范围: $(round(minimum(μ_B_values)*hc, digits=1)) - $(round(maximum(μ_B_values)*hc, digits=1)) MeV")
        println("权重范围: $(round(minimum(weights), digits=3)) - $(round(maximum(weights), digits=3))")
    end
    
    weighted_sum_of_squares = 0.0
    valid_pairs = 0
    failed_pairs = 0
    
    for (i, ((kappa3_kappa1, kappa4_kappa2), weight, μ_B)) in enumerate(zip(kappa_pairs, weights, μ_B_values))
        if verbose
            println("\n处理第 $i/$(length(kappa_pairs)) 组: κ₃/κ₁ = $kappa3_kappa1, κ₄/κ₂ = $kappa4_kappa2, 权重 = $weight, μ_B = $(round(μ_B*hc, digits=1)) MeV")
        end
        
        try
            T_k3k1, T_k4k2 = find_temperature_for_kappa_ratios_with_optimization_params(
                kappa3_kappa1, kappa4_kappa2, μ_B, optimization_params,
                T_min, T_max, T_step_scan;
                gsigma=gsigma, gdelta=gdelta, n_nodes=n_nodes, verbose=false)
            
            # 检查是否成功找到两个温度
            if !isnan(T_k3k1) && !isnan(T_k4k2) && isfinite(T_k3k1) && isfinite(T_k4k2)
                # 计算加权温度差的平方 (转换为MeV单位)
                temp_diff_MeV = abs(T_k3k1 - T_k4k2) * hc
                contribution = weight * temp_diff_MeV^2
                weighted_sum_of_squares += contribution
                valid_pairs += 1
                
                if verbose
                    println("  ✓ |ΔT| = $(round(temp_diff_MeV, digits=2)) MeV")
                    println("    加权贡献: $(round(contribution, digits=2)) MeV²")
                end
            else
                # 无法找到有效温度，应用加权惩罚
                weighted_penalty = weight * penalty_for_missing
                weighted_sum_of_squares += weighted_penalty
                failed_pairs += 1
                
                if verbose
                    println("  ✗ 应用加权惩罚: $(round(weighted_penalty, digits=2))")
                end
            end
            
        catch e
            # 计算过程出错，应用加权惩罚
            weighted_penalty = weight * penalty_for_missing
            weighted_sum_of_squares += weighted_penalty
            failed_pairs += 1
            
            if verbose
                println("  ✗ 计算失败，应用加权惩罚: $(round(weighted_penalty, digits=2))")
            end
        end
    end
    
    if verbose
        println("\n" * "="^60)
        println("加权计算完成")
        println("="^60)
        println("总组数: $(length(kappa_pairs))")
        println("成功计算: $valid_pairs")
        println("失败/无效: $failed_pairs")
        println("加权温度差平方和: $(round(weighted_sum_of_squares, digits=2)) MeV²")
    end
    
    return weighted_sum_of_squares
end

function demo_temperature_difference_sum_of_squares()
    """
    演示函数：展示如何使用 calculate_temperature_difference_sum_of_squares 函数
    """
    
    println("="^80)
    println("演示：温度差平方和计算")
    println("="^80)
    
    # 定义测试参数
    kappa_pairs = [(1.2, 2.5), (1.5, 3.0), (1.8, 3.5)]  # 示例κ比值对
    μ_B_values = [300.0/hc, 320.0/hc, 340.0/hc]  # 每组对应不同的μ_B值
    optimization_params = (0.15, 16.0, 240.0, 0.7, 32.0)  # (ρ₀, B_A, K, m_ratio, E_sym)
    T_min = 80.0 / hc   # 80 MeV
    T_max = 200.0 / hc  # 200 MeV
    
    println("测试参数:")
    println("  κ比值对: $kappa_pairs")
    println("  μ_B值: $([round(μ*hc, digits=1) for μ in μ_B_values]) MeV")
    println("  优化参数: $optimization_params")
    println("  温度范围: $(T_min*hc) - $(T_max*hc) MeV")
    
    # 计算温度差平方和
    println("\n计算温度差平方和...")
    sum_of_squares = calculate_temperature_difference_sum_of_squares(
        kappa_pairs, μ_B_values, optimization_params, T_min, T_max;
        T_step_scan=2.0/hc, verbose=true, penalty_for_missing=1e4)
    
    println("\n最终结果:")
    println("温度差平方和 = $(round(sum_of_squares, digits=2)) MeV²")
    
    # 演示加权版本
    println("\n" * "="^60)
    println("演示：加权温度差平方和计算")
    println("="^60)
    
    weights = [1.0, 2.0, 0.5]  # 示例权重
    println("权重: $weights")
    
    weighted_sum = calculate_temperature_difference_sum_of_squares_with_weights(
        kappa_pairs, weights, μ_B_values, optimization_params, T_min, T_max;
        T_step_scan=2.0/hc, verbose=true, penalty_for_missing=1e4)
    
    println("\n最终结果:")
    println("加权温度差平方和 = $(round(weighted_sum, digits=2)) MeV²")
    
    return sum_of_squares, weighted_sum
end

function create_temperature_difference_objective(
    kappa_pairs, μ_B_values, T_min, T_max;
    T_step_scan=1.0/hc, gsigma=1.25, gdelta=0.01, n_nodes=256,
    penalty_for_missing=1e6, verbose=false)
    """
    创建温度差平方和目标函数的闭包
    
    此函数将实验确定的参数（kappa_pairs, μ_B_values, T_min, T_max）封装在闭包中，
    返回一个只需要优化参数作为输入的目标函数，适用于参数优化算法。
    
    参数:
    - kappa_pairs: κ比值对数组，格式 [(κ₃/κ₁, κ₄/κ₂), ...] (实验确定)
    - μ_B_values: 重子化学势数组，与 kappa_pairs 一一对应 (实验确定)
    - T_min, T_max: 温度搜索范围 (实验确定)
    - T_step_scan: 扫描步长
    - gsigma, gdelta: 强子-夸克耦合参数
    - n_nodes: 积分节点数
    - penalty_for_missing: 计算失败时的惩罚值
    - verbose: 是否显示详细信息（优化时建议设为false）
    
    返回:
    - objective_function: 闭包函数，签名为 f(optimization_params) -> Float64
      其中 optimization_params = (ρ₀, B_A, K, m_ratio, E_sym)
    
    使用示例:
    ```julia
    # 创建目标函数
    kappa_pairs = [(1.2, 2.5), (1.5, 3.0)]
    μ_B_values = [300.0/hc, 320.0/hc]  # 每组对应不同的μ_B
    T_min, T_max = 80.0/hc, 200.0/hc
    
    objective_func = create_temperature_difference_objective(
        kappa_pairs, μ_B, T_min, T_max; verbose=false)
    
    # 在优化中使用
    optimization_params = (0.15, 16.0, 240.0, 0.7, 32.0)
    result = objective_func(optimization_params)
    ```
    """
    
    function objective_function(optimization_params)
        """
        目标函数闭包：只需要优化参数输入
        
        参数:
        - optimization_params: 优化参数元组 (ρ₀, B_A, K, m_ratio, E_sym)
        
        返回:
        - sum_of_squares: 温度差平方和 (MeV²)
        """
        return calculate_temperature_difference_sum_of_squares(
            kappa_pairs, μ_B_values, optimization_params, T_min, T_max;
            T_step_scan=T_step_scan, gsigma=gsigma, gdelta=gdelta, 
            n_nodes=n_nodes, penalty_for_missing=penalty_for_missing, verbose=verbose)
    end
    
    return objective_function
end

function create_weighted_temperature_difference_objective(
    kappa_pairs, weights, μ_B_values, T_min, T_max;
    T_step_scan=1.0/hc, gsigma=1.25, gdelta=0.01, n_nodes=256,
    penalty_for_missing=1e6, verbose=false)
    """
    创建加权温度差平方和目标函数的闭包
    
    参数:
    - kappa_pairs: κ比值对数组
    - weights: 权重数组，与 kappa_pairs 长度相同
    - μ_B_values: 重子化学势数组，与 kappa_pairs 一一对应
    - 其他参数: 与 create_temperature_difference_objective 相同
    
    返回:
    - objective_function: 闭包函数，签名为 f(optimization_params) -> Float64
    """
    
    if length(weights) != length(kappa_pairs)
        error("权重数组长度 ($(length(weights))) 与κ比值对数组长度 ($(length(kappa_pairs))) 不匹配")
    end
    
    if length(μ_B_values) != length(kappa_pairs)
        error("μ_B值数组长度 ($(length(μ_B_values))) 与κ比值对数组长度 ($(length(kappa_pairs))) 不匹配")
    end
    
    function weighted_objective_function(optimization_params)
        """
        加权目标函数闭包：只需要优化参数输入
        
        参数:
        - optimization_params: 优化参数元组 (ρ₀, B_A, K, m_ratio, E_sym)
        
        返回:
        - weighted_sum_of_squares: 加权温度差平方和 (MeV²)
        """
        return calculate_temperature_difference_sum_of_squares_with_weights(
            kappa_pairs, weights, μ_B_values, optimization_params, T_min, T_max;
            T_step_scan=T_step_scan, gsigma=gsigma, gdelta=gdelta, 
            n_nodes=n_nodes, penalty_for_missing=penalty_for_missing, verbose=verbose)
    end
    
    return weighted_objective_function
end

function demo_objective_function_closure()
    """
    演示目标函数闭包的使用
    """
    
    println("="^80)
    println("演示：目标函数闭包")
    println("="^80)
    
    # 实验确定的参数
    kappa_pairs = [(1.2, 2.5), (1.5, 3.0), (1.8, 3.5)]
    μ_B_values = [300.0/hc, 320.0/hc, 340.0/hc]  # 每组对应不同的μ_B值
    T_min, T_max = 80.0/hc, 200.0/hc  # 80-200 MeV
    
    println("实验确定的参数:")
    println("  κ比值对: $kappa_pairs")
    println("  μ_B值: $([round(μ*hc, digits=1) for μ in μ_B_values]) MeV")
    println("  温度范围: $(T_min*hc) - $(T_max*hc) MeV")
    
    # 创建目标函数
    println("\n创建目标函数闭包...")
    objective_func = create_temperature_difference_objective(
        kappa_pairs, μ_B_values, T_min, T_max; 
        T_step_scan=2.0/hc, verbose=false, penalty_for_missing=1e4)
    
    # 测试不同的优化参数
    test_params = [
        (0.15, 16.0, 240.0, 0.7, 32.0),  # 参数组1
        (0.16, 15.5, 250.0, 0.75, 30.0), # 参数组2
        (0.14, 16.5, 230.0, 0.65, 34.0)  # 参数组3
    ]
    
    println("\n测试不同优化参数:")
    for (i, params) in enumerate(test_params)
        ρ₀, B_A, K, m_ratio, E_sym = params
        println("\n参数组 $i: ρ₀=$ρ₀, B_A=$B_A MeV, K=$K MeV, m_ratio=$m_ratio, E_sym=$E_sym MeV")
        
        # 使用闭包计算目标值
        result = objective_func(params)
        println("  目标函数值: $(round(result, digits=2)) MeV²")
    end
    
    # 演示加权版本
    println("\n" * "="^60)
    println("演示：加权目标函数闭包")
    println("="^60)
    
    weights = [1.0, 2.0, 0.5]
    println("权重: $weights")
    
    weighted_objective_func = create_weighted_temperature_difference_objective(
        kappa_pairs, weights, μ_B_values, T_min, T_max; 
        T_step_scan=2.0/hc, verbose=false, penalty_for_missing=1e4)
    
    println("\n测试加权目标函数:")
    for (i, params) in enumerate(test_params[1:2])  # 只测试前两组参数
        ρ₀, B_A, K, m_ratio, E_sym = params
        println("\n参数组 $i: ρ₀=$ρ₀, B_A=$B_A MeV, K=$K MeV, m_ratio=$m_ratio, E_sym=$E_sym MeV")
        
        result = weighted_objective_func(params)
        println("  加权目标函数值: $(round(result, digits=2)) MeV²")
    end
    
    println("\n目标函数闭包创建成功！可用于优化算法。")
    
    return objective_func, weighted_objective_func
end
