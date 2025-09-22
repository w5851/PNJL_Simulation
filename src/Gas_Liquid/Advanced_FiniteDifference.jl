# Advanced_Gas_Liquid.jl
# 高阶函数模块：包含压强导数计算、热力学涨落分析和批量处理功能
# 依赖于 Function_Gas_Liquid.jl 中的基础函数

include("Function_Gas_Liquid.jl")
include("Constants_Gas_Liquid.jl")
using .Constants_Gas_Liquid: hc
using FiniteDifferences

function calculate_pressure_derivatives(μ_B, T, x0, nodes, couplings; order=4, method=central_fdm(5, 1))
    """
    计算无量纲压强P/T⁴对无量纲化学势μ/T的一到四阶导数
    
    参数:
    - μ_B: 重子化学势
    - T: 温度
    - x0: 初始猜测值 [gσ, gδ, μ_p, μ_n]
    - nodes: 积分节点和权重
    - couplings: 耦合常数 [fσ, fω, fρ, fδ, b, c]
    - order: 计算导数的最高阶数 (默认4)
    - method: 有限差分方法 (默认5点中心差分)
    
    返回:
    - pressure_normalized: 无量纲压强 P/T⁴
    - derivatives: 包含一到四阶导数的数组
    """
    
    # 定义无量纲压强函数，只依赖于无量纲化学势
    pressure_normalized_func = μ_norm -> begin
        μ_physical = μ_norm * T  # 转换回物理化学势
        pressure_physical = calculate_pressure_solved(μ_physical, T, x0, nodes, couplings)
        return pressure_physical / (T^4)  # 返回无量纲压强
    end
    
    # 当前的无量纲化学势
    μ_B_normalized = μ_B / T
    
    # 计算无量纲压强值
    pressure_normalized = pressure_normalized_func(μ_B_normalized)
    
    # 计算各阶导数
    derivatives = zeros(order)
    
    for i in 1:order
        # 创建相应阶数的有限差分方法
        fdm_method = central_fdm(5, i)
        derivatives[i] = fdm_method(pressure_normalized_func, μ_B_normalized)
    end
    
    return pressure_normalized, derivatives
end

function calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings; h=1e-5)
    """
    高效计算无量纲压强P/T⁴对无量纲化学势μ/T的一到四阶导数
    这是热力学中的标准做法，可以得到正确的累积量
    
    参数:
    - μ_B: 重子化学势
    - T: 温度  
    - x0: 初始猜测值 [gσ, gδ, μ_p, μ_n]
    - nodes: 积分节点和权重
    - couplings: 耦合常数 [fσ, fω, fρ, fδ, b, c]
    - h: 步长 (默认1e-5)
    
    返回:
    - pressure_normalized: 无量纲压强 P/T⁴
    - dpre_dmu1: 一阶导数 ∂(P/T⁴)/∂(μ/T)
    - dpre_dmu2: 二阶导数 ∂²(P/T⁴)/∂(μ/T)²
    - dpre_dmu3: 三阶导数 ∂³(P/T⁴)/∂(μ/T)³
    - dpre_dmu4: 四阶导数 ∂⁴(P/T⁴)/∂(μ/T)⁴
    """
    
    # 定义无量纲压强函数 P/T⁴ 对无量纲化学势 μ/T 的函数
    # 注意：μ_normalized = μ/T，所以 μ = μ_normalized * T
    pressure_normalized_func = μ_norm -> begin
        μ_physical = μ_norm * T  # 转换回物理化学势
        pressure_physical = calculate_pressure_solved(μ_physical, T, x0, nodes, couplings)
        return pressure_physical / (T^4)  # 返回无量纲压强
    end
    
    # 当前的无量纲化学势
    μ_B_normalized = μ_B / T
    
    # 预定义有限差分方法
    fdm1 = central_fdm(5, 1)  # 一阶导数，5点中心差分
    fdm2 = central_fdm(5, 2)  # 二阶导数，5点中心差分
    fdm3 = central_fdm(7, 3)  # 三阶导数，7点中心差分
    fdm4 = central_fdm(7, 4)  # 四阶导数，7点中心差分
    
    # 计算无量纲压强和各阶导数
    pressure_normalized = pressure_normalized_func(μ_B_normalized)
    dpre_dmu1 = fdm1(pressure_normalized_func, μ_B_normalized)
    dpre_dmu2 = fdm2(pressure_normalized_func, μ_B_normalized)
    dpre_dmu3 = fdm3(pressure_normalized_func, μ_B_normalized)
    dpre_dmu4 = fdm4(pressure_normalized_func, μ_B_normalized)
    
    return pressure_normalized, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4
end

function calculate_thermodynamic_fluctuations(μ_B, T, x0, nodes, couplings)
    """
    计算热力学涨落相关量
    基于无量纲压强P/T⁴对无量纲化学势μ/T的导数计算累积量和涨落
    
    参数:
    - μ_B: 重子化学势
    - T: 温度
    - x0: 初始猜测值 [gσ, gδ, μ_p, μ_n] 
    - nodes: 积分节点和权重
    - couplings: 耦合常数 [fσ, fω, fρ, fδ, b, c]
    
    返回:
    - kappa1: 第一累积量 (重子数密度/T³)
    - kappa2: 第二累积量 (重子数涨落/T³)
    - kappa3: 第三累积量 (三阶累积量/T³)
    - kappa4: 第四累积量 (四阶累积量/T³)
    - fluctuation_ratios: 涨落比值 [κ₂/κ₁, κ₃/κ₂, κ₄/κ₂]
    """
    
    # 计算无量纲压强和导数
    pressure_normalized, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4 = 
        calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings)
    
    # 根据热力学关系计算累积量
    # 在相对论性情况下，累积量的定义为：
    # κₙ = ∂ⁿ(P/T⁴)/∂(μ/T)ⁿ
    kappa1 = dpre_dmu1  # 第一累积量：重子数密度/T³
    kappa2 = dpre_dmu2  # 第二累积量：重子数涨落/T³
    kappa3 = dpre_dmu3  # 第三累积量：三阶累积量/T³
    kappa4 = dpre_dmu4  # 第四累积量：四阶累积量/T³
    
    # 计算涨落比值
    fluctuation_ratios = [
        kappa2 / kappa1,   # κ₂/κ₁ (归一化方差)
        kappa3 / kappa2,   # κ₃/κ₂ (偏度相关)
        kappa4 / kappa2    # κ₄/κ₂ (峰度相关)
    ]
    
    return kappa1, kappa2, kappa3, kappa4, fluctuation_ratios
end

function calculate_derivatives_batch(μ_B_array, T, x0, nodes, couplings; save_results=false, output_file=joinpath(@__DIR__, "..", "..", "output", "Gas_Liquid", "derivatives_output.csv"))
    """
    批量计算多个化学势点的压强导数
    模拟Fortran代码中的循环计算过程
    
    参数:
    - μ_B_array: 化学势数组
    - T: 温度
    - x0: 初始猜测值 [gσ, gδ, μ_p, μ_n]
    - nodes: 积分节点和权重
    - couplings: 耦合常数
    - save_results: 是否保存结果到文件
    - output_file: 输出文件名
    
    返回:
    - results: 包含所有计算结果的NamedTuple
    """
    
    n_points = length(μ_B_array)
    
    # 预分配结果数组
    pressure_array = zeros(n_points)
    dpre_dmu1_array = zeros(n_points)
    dpre_dmu2_array = zeros(n_points)
    dpre_dmu3_array = zeros(n_points)
    dpre_dmu4_array = zeros(n_points)
    kappa1_array = zeros(n_points)
    kappa2_array = zeros(n_points)
    kappa3_array = zeros(n_points)
    kappa4_array = zeros(n_points)
    fluctuation_ratios_array = zeros(n_points, 3)
    
    println("开始批量计算 $(n_points) 个化学势点的压强导数...")
    
    # 循环计算每个化学势点
    for (i, μ_B) in enumerate(μ_B_array)
        try
            # 计算压强和导数
            pressure, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4 = 
                calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings)
            
            # 计算热力学涨落
            kappa1, kappa2, kappa3, kappa4, fluctuation_ratios = 
                calculate_thermodynamic_fluctuations(μ_B, T, x0, nodes, couplings)
            
            # 存储结果
            pressure_array[i] = pressure
            dpre_dmu1_array[i] = dpre_dmu1
            dpre_dmu2_array[i] = dpre_dmu2
            dpre_dmu3_array[i] = dpre_dmu3
            dpre_dmu4_array[i] = dpre_dmu4
            kappa1_array[i] = kappa1
            kappa2_array[i] = kappa2
            kappa3_array[i] = kappa3
            kappa4_array[i] = kappa4
            fluctuation_ratios_array[i, :] = fluctuation_ratios
            
            # 进度报告
            if i % 10 == 0 || i == n_points
                println("已完成: $(i)/$(n_points) ($(round(i/n_points*100, digits=1))%)")
            end
            
        catch e
            println("警告: 化学势 μ_B = $(μ_B*hc) MeV 处计算失败: $e")
            # 填入NaN值
            pressure_array[i] = NaN
            dpre_dmu1_array[i] = NaN
            dpre_dmu2_array[i] = NaN
            dpre_dmu3_array[i] = NaN
            dpre_dmu4_array[i] = NaN
            kappa1_array[i] = NaN
            kappa2_array[i] = NaN
            kappa3_array[i] = NaN
            kappa4_array[i] = NaN
            fluctuation_ratios_array[i, :] .= NaN
        end
    end
    
    # 创建结果NamedTuple
    results = (
        μ_B = μ_B_array,
        T = T,
        pressure = pressure_array,
        dpre_dmu1 = dpre_dmu1_array,
        dpre_dmu2 = dpre_dmu2_array,
        dpre_dmu3 = dpre_dmu3_array,
        dpre_dmu4 = dpre_dmu4_array,
        kappa1 = kappa1_array,
        kappa2 = kappa2_array,
        kappa3 = kappa3_array,
        kappa4 = kappa4_array,
        fluctuation_ratios = fluctuation_ratios_array
    )
    
    # 保存结果到文件
    if save_results
        save_derivatives_results(results, output_file)
        println("结果已保存到文件: $output_file")
    end
    
    println("批量计算完成!")
    return results
end

function save_derivatives_results(results, filename)
    """保存导数计算结果到CSV文件"""
    println("正在保存结果到CSV文件: $filename")
    
    # 确保输出目录存在
    output_dir = dirname(filename)
    if !isdir(output_dir)
        println("创建输出目录: $output_dir")
        mkpath(output_dir)
    end
    
    try
        open(filename, "w") do io
            # 写入CSV文件头
            write(io, "mu_B_MeV,T_MeV,pressure,dpre_dmu1,dpre_dmu2,dpre_dmu3,dpre_dmu4,kappa1,kappa2,kappa3,kappa4,kappa2_over_kappa1,kappa3_over_kappa2,kappa4_over_kappa2\n")
            
            # 写入数据
            for i in eachindex(results.μ_B)
                write(io, "$(results.μ_B[i]*hc),$(results.T*hc),$(results.pressure[i]),")
                write(io, "$(results.dpre_dmu1[i]),$(results.dpre_dmu2[i]),$(results.dpre_dmu3[i]),$(results.dpre_dmu4[i]),")
                write(io, "$(results.kappa1[i]),$(results.kappa2[i]),$(results.kappa3[i]),$(results.kappa4[i]),")
                write(io, "$(results.fluctuation_ratios[i,1]),$(results.fluctuation_ratios[i,2]),$(results.fluctuation_ratios[i,3])\n")
            end
        end
        println("CSV文件保存成功!")
    catch e
        println("保存CSV文件时出错: $e")
        rethrow(e)
    end
end

function calculate_fluctuation_ratios_vs_temperature(μ_B, T_min, T_max, x0, nodes, couplings; 
                                                   T_step=1.0/hc, save_results=false, 
                                                   output_file=joinpath(@__DIR__, "..", "..", "output", "Gas_Liquid", "fluctuation_ratios_vs_T.csv"))
    """
    固定重子化学势μ_B，改变温度T，计算不同温度下的κ₃/κ₁和κ₄/κ₂
    使用迭代初解策略：每个温度点使用上一个温度点的收敛解作为初始猜测值
    
    参数:
    - μ_B: 固定的重子化学势
    - T_min: 最小温度
    - T_max: 最大温度  
    - x0: 初始猜测值 [gσ, gδ, μ_p, μ_n]（仅用于第一个温度点）
    - nodes: 积分节点和权重
    - couplings: 耦合常数 [fσ, fω, fρ, fδ, b, c]
    - T_step: 温度步长 (默认1/hc)
    - save_results: 是否保存结果到文件
    - output_file: 输出文件名
    
    返回:
    - temperature_array: 温度数组
    - kappa3_over_kappa1: κ₃/κ₁数组
    - kappa4_over_kappa2: κ₄/κ₂数组
    - results_matrix: [温度, κ₃/κ₁, κ₄/κ₂] 格式的结果矩阵
    """
    
    # 生成温度数组
    T_array = T_min:T_step:T_max
    n_points = length(T_array)
    
    # 预分配结果数组
    temperature_array = zeros(n_points)
    kappa3_over_kappa1 = zeros(n_points)
    kappa4_over_kappa2 = zeros(n_points)
    
    # 存储每个温度点的收敛解，用于下一个点的初始猜测
    solution_history = Vector{Vector{Float64}}(undef, n_points)
    
    println("开始计算固定μ_B = $(μ_B*hc) MeV下，$(n_points) 个温度点的涨落比值...")
    println("温度范围: $(T_min*hc) - $(T_max*hc) MeV，步长: $(T_step*hc) MeV")
    println("使用迭代初解策略提高收敛性...")
    
    # 初始化第一个点的初始猜测值
    current_x0 = copy(x0)
    
    # 循环计算每个温度点
    for (i, T) in enumerate(T_array)
        try
            # 设置当前温度的参数
            params = [T, μ_B]
            
            # 求解自洽方程获得收敛解
            converged_solution = solve_fun_constraints(current_x0, nodes, couplings, params)
            
            # 存储收敛解
            solution_history[i] = copy(converged_solution)
            
            # 使用收敛解计算热力学涨落
            kappa1, kappa2, kappa3, kappa4, _ = 
                calculate_thermodynamic_fluctuations(μ_B, T, converged_solution, nodes, couplings)
            
            # 存储结果
            temperature_array[i] = T
            
            # 计算涨落比值，避免除零
            if abs(kappa1) > 1e-12
                kappa3_over_kappa1[i] = kappa3 / kappa1
            else
                kappa3_over_kappa1[i] = NaN
            end
            
            if abs(kappa2) > 1e-12
                kappa4_over_kappa2[i] = kappa4 / kappa2
            else
                kappa4_over_kappa2[i] = NaN
            end
            
            # 更新下一个温度点的初始猜测值
            if i < n_points
                current_x0 = copy(converged_solution)
                
                # 可选：对下一个温度点的初始猜测进行微调
                # 这里可以根据物理直觉对某些参数进行外推
                # 例如：假设场量随温度的变化是平滑的
                if i > 1
                    # 使用线性外推改善初始猜测
                    prev_solution = solution_history[i-1]
                    current_solution = solution_history[i]
                    extrapolated = current_solution .+ (current_solution .- prev_solution)
                    
                    # 混合使用外推值和当前解，避免过度外推
                    mix_factor = 0.3  # 外推权重
                    current_x0 = (1 - mix_factor) .* current_solution .+ mix_factor .* extrapolated
                end
            end
            
            # 进度报告
            if i % 10 == 0 || i == n_points
                println("已完成: $(i)/$(n_points) ($(round(i/n_points*100, digits=1))%) - T = $(T*hc) MeV")
                println("  当前解: [gσ=$(round(converged_solution[1], digits=4)), gδ=$(round(converged_solution[2], digits=4)), μ_p=$(round(converged_solution[3]*hc, digits=2)) MeV, μ_n=$(round(converged_solution[4]*hc, digits=2)) MeV]")
            end
            
        catch e
            println("警告: 温度 T = $(T*hc) MeV 处计算失败: $e")
            # 填入NaN值
            temperature_array[i] = T
            kappa3_over_kappa1[i] = NaN
            kappa4_over_kappa2[i] = NaN
            solution_history[i] = copy(current_x0)  # 保持当前初始猜测值不变
        end
    end
    
    # 创建结果矩阵 [温度, κ₃/κ₁, κ₄/κ₂]
    results_matrix = hcat(temperature_array .* hc, kappa3_over_kappa1, kappa4_over_kappa2)
    
    # 保存结果到文件
    if save_results
        save_fluctuation_ratios_results(results_matrix, μ_B, output_file)
        println("结果已保存到文件: $output_file")
    end
    
    println("温度扫描计算完成!")
    return temperature_array, kappa3_over_kappa1, kappa4_over_kappa2, results_matrix
end

function calculate_fluctuation_ratios_vs_temperature_advanced(μ_B, T_min, T_max, x0, nodes, couplings; 
                                                            T_step=1.0/hc, save_results=false, 
                                                            output_file=joinpath(@__DIR__, "..", "..", "output", "Gas_Liquid", "fluctuation_ratios_vs_T.csv"),
                                                            use_iterative_guess=true, extrapolation_weight=0.3,
                                                            return_solution_history=false)
    """
    固定重子化学势μ_B，改变温度T，计算不同温度下的κ₃/κ₁和κ₄/κ₂
    高级版本：提供更多控制选项
    
    参数:
    - μ_B: 固定的重子化学势
    - T_min: 最小温度
    - T_max: 最大温度  
    - x0: 初始猜测值 [gσ, gδ, μ_p, μ_n]
    - nodes: 积分节点和权重
    - couplings: 耦合常数 [fσ, fω, fρ, fδ, b, c]
    - T_step: 温度步长 (默认1/hc)
    - save_results: 是否保存结果到文件
    - output_file: 输出文件名
    - use_iterative_guess: 是否使用迭代初解策略 (默认true)
    - extrapolation_weight: 外推权重，控制线性外推的强度 (默认0.3)
    - return_solution_history: 是否返回每个温度点的收敛解历史 (默认false)
    
    返回:
    - temperature_array: 温度数组
    - kappa3_over_kappa1: κ₃/κ₁数组
    - kappa4_over_kappa2: κ₄/κ₂数组
    - results_matrix: [温度, κ₃/κ₁, κ₄/κ₂] 格式的结果矩阵
    - solution_history: 每个温度点的收敛解历史 (仅当return_solution_history=true时返回)
    """
    
    # 生成温度数组
    T_array = T_min:T_step:T_max
    n_points = length(T_array)
    
    # 预分配结果数组
    temperature_array = zeros(n_points)
    kappa3_over_kappa1 = zeros(n_points)
    kappa4_over_kappa2 = zeros(n_points)
    solution_history = Vector{Vector{Float64}}(undef, n_points)
    
    println("开始计算固定μ_B = $(μ_B*hc) MeV下，$(n_points) 个温度点的涨落比值...")
    println("温度范围: $(T_min*hc) - $(T_max*hc) MeV，步长: $(T_step*hc) MeV")
    
    if use_iterative_guess
        println("使用迭代初解策略，外推权重: $(extrapolation_weight)")
    else
        println("使用固定初解策略")
    end
    
    # 初始化第一个点的初始猜测值
    current_x0 = copy(x0)
    
    # 循环计算每个温度点
    for (i, T) in enumerate(T_array)
        try
            if use_iterative_guess
                # 使用迭代初解策略
                params = [T, μ_B]
                converged_solution = solve_fun_constraints(current_x0, nodes, couplings, params)
                solution_history[i] = copy(converged_solution)
                
                # 使用收敛解计算热力学涨落
                kappa1, kappa2, kappa3, kappa4, _ = 
                    calculate_thermodynamic_fluctuations(μ_B, T, converged_solution, nodes, couplings)
                
                # 更新下一个温度点的初始猜测值
                if i < n_points
                    current_x0 = copy(converged_solution)
                    
                    # 线性外推改善初始猜测
                    if i > 1 && extrapolation_weight > 0
                        prev_solution = solution_history[i-1]
                        current_solution = solution_history[i]
                        extrapolated = current_solution .+ (current_solution .- prev_solution)
                        current_x0 = (1 - extrapolation_weight) .* current_solution .+ extrapolation_weight .* extrapolated
                    end
                end
            else
                # 使用固定初解策略
                kappa1, kappa2, kappa3, kappa4, _ = 
                    calculate_thermodynamic_fluctuations(μ_B, T, x0, nodes, couplings)
                solution_history[i] = copy(x0)  # 保存初始猜测值作为"解"
            end
            
            # 存储结果
            temperature_array[i] = T
            
            # 计算涨落比值，避免除零
            if abs(kappa1) > 1e-12
                kappa3_over_kappa1[i] = kappa3 / kappa1
            else
                kappa3_over_kappa1[i] = NaN
            end
            
            if abs(kappa2) > 1e-12
                kappa4_over_kappa2[i] = kappa4 / kappa2
            else
                kappa4_over_kappa2[i] = NaN
            end
            
            # 进度报告
            if i % 10 == 0 || i == n_points
                println("已完成: $(i)/$(n_points) ($(round(i/n_points*100, digits=1))%) - T = $(T*hc) MeV")
                if use_iterative_guess
                    sol = solution_history[i]
                    println("  当前解: [gσ=$(round(sol[1], digits=4)), gδ=$(round(sol[2], digits=4)), μ_p=$(round(sol[3]*hc, digits=2)) MeV, μ_n=$(round(sol[4]*hc, digits=2)) MeV]")
                end
            end
            
        catch e
            println("警告: 温度 T = $(T*hc) MeV 处计算失败: $e")
            # 填入NaN值
            temperature_array[i] = T
            kappa3_over_kappa1[i] = NaN
            kappa4_over_kappa2[i] = NaN
            solution_history[i] = copy(current_x0)
        end
    end
    
    # 创建结果矩阵 [温度, κ₃/κ₁, κ₄/κ₂]
    results_matrix = hcat(temperature_array .* hc, kappa3_over_kappa1, kappa4_over_kappa2)
    
    # 保存结果到文件
    if save_results
        save_fluctuation_ratios_results(results_matrix, μ_B, output_file)
        println("结果已保存到文件: $output_file")
    end
    
    println("温度扫描计算完成!")
    
    if return_solution_history
        return temperature_array, kappa3_over_kappa1, kappa4_over_kappa2, results_matrix, solution_history
    else
        return temperature_array, kappa3_over_kappa1, kappa4_over_kappa2, results_matrix
    end
end

function save_fluctuation_ratios_results(results_matrix, μ_B, filename)
    """保存涨落比值计算结果到CSV文件"""
    println("正在保存涨落比值结果到CSV文件: $filename")
    
    # 确保输出目录存在
    output_dir = dirname(filename)
    if !isdir(output_dir)
        println("创建输出目录: $output_dir")
        mkpath(output_dir)
    end
    
    try
        open(filename, "w") do io
            # 写入CSV文件头
            write(io, "# 涨落比值随温度变化结果\n")
            write(io, "# 固定重子化学势: μ_B = $(μ_B*hc) MeV\n")
            write(io, "T_MeV,kappa3_over_kappa1,kappa4_over_kappa2\n")
            
            # 写入数据
            for i in 1:size(results_matrix, 1)
                write(io, "$(results_matrix[i,1]),$(results_matrix[i,2]),$(results_matrix[i,3])\n")
            end
        end
        println("CSV文件保存成功!")
    catch e
        println("保存CSV文件时出错: $e")
        rethrow(e)
    end
end