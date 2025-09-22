# test_find_temperature.jl
# 测试温度反向查找功能
# 测试Advanced_FindTforDiff.jl模块中的温度查找功能

# 激活项目环境
import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using Printf
using Dates

# 引入温度查找模块
include("../../src/Gas_Liquid/Advanced_FindTforDiff.jl")

println("="^80)
println("温度反向查找功能测试")
println("测试时间: $(Dates.now())")
println("="^80)

# 测试参数设置
μ_B_test = 697.0 / hc  # 697 MeV 重子化学势
T_min = 20.0 / hc      # 20 MeV 
T_max = 200.0 / hc     # 200 MeV
T_step = 2.0 / hc      # 2 MeV 步长（加快测试速度）

# 目标κ比值（示例值）
target_kappa3_kappa1 = 0.5  # 目标 κ₃/κ₁ 值
target_kappa4_kappa2 = 1.2  # 目标 κ₄/κ₂ 值

println("\n测试配置:")
println("μ_B = $(μ_B_test * hc) MeV")
println("T范围: $(T_min * hc) - $(T_max * hc) MeV")
println("T步长: $(T_step * hc) MeV")
println("目标 κ₃/κ₁ = $target_kappa3_kappa1")
println("目标 κ₄/κ₂ = $target_kappa4_kappa2")

# 执行温度查找测试
println("\n开始温度反向查找测试...")
try
    T_kappa3_kappa1, T_kappa4_kappa2 = find_temperature_for_kappa_ratios(
        target_kappa3_kappa1, target_kappa4_kappa2, μ_B_test, 
        T_min, T_max, T_step;
        gsigma=1.25, gdelta=0.01,
        fs=17.28476, fo=11.66174, fr=0.89363, fd=0.0,
        b=0.00210, c=-0.00297, n_nodes=128,  # 减少积分点以加快测试
        verbose=true
    )
    
    println("\n" * "="^80)
    println("温度反向查找测试完成!")
    println("="^80)
    
    # 输出结果
    println("\n最终结果:")
    if !isnan(T_kappa3_kappa1)
        @printf("κ₃/κ₁ = %.3f 对应温度: T = %.2f MeV\n", target_kappa3_kappa1, T_kappa3_kappa1 * hc)
    else
        println("κ₃/κ₁ = $target_kappa3_kappa1: 未找到对应温度")
    end
    
    if !isnan(T_kappa4_kappa2)
        @printf("κ₄/κ₂ = %.3f 对应温度: T = %.2f MeV\n", target_kappa4_kappa2, T_kappa4_kappa2 * hc)
    else
        println("κ₄/κ₂ = $target_kappa4_kappa2: 未找到对应温度")
    end
    
    # 计算温度差异
    if !isnan(T_kappa3_kappa1) && !isnan(T_kappa4_kappa2)
        T_diff = abs(T_kappa3_kappa1 - T_kappa4_kappa2) * hc
        @printf("\n两个κ比值对应温度差异: %.2f MeV\n", T_diff)
        
        if T_diff < 5.0  # 5 MeV容差
            println("✅ 两个κ比值在相近温度下同时满足")
        else
            println("⚠️  两个κ比值对应的温度差异较大")
        end
    end
    
catch e
    println("\n✗ 温度反向查找测试失败:")
    println("错误信息: $e")
    println("错误位置: $(stacktrace()[1])")
end

println("\n" * "="^80)
println("测试目标函数闭包")
println("="^80)

# 测试闭包函数
println("\n开始测试目标函数闭包...")
try
    # 实验确定的参数（这些参数将在实际使用时由用户填入）
    experimental_kappa_pairs = [(1.2, 2.5), (1.5, 3.0), (1.8, 3.5)]  # 实验κ比值对
    experimental_μ_B = 300.0 / hc  # 300 MeV，实验确定的重子化学势
    experimental_T_min = 80.0 / hc   # 80 MeV，实验温度下限
    experimental_T_max = 200.0 / hc  # 200 MeV，实验温度上限
    
    println("\n实验确定的参数:")
    println("  κ比值对: $experimental_kappa_pairs")
    println("  μ_B = $(experimental_μ_B*hc) MeV")
    println("  温度搜索范围: $(experimental_T_min*hc) - $(experimental_T_max*hc) MeV")
    
    # 创建目标函数闭包
    println("\n创建目标函数闭包...")
    objective_func = create_temperature_difference_objective(
        experimental_kappa_pairs, experimental_μ_B, experimental_T_min, experimental_T_max;
        T_step_scan=3.0/hc,  # 加快测试速度
        verbose=false,       # 减少输出
        penalty_for_missing=1e4,
        n_nodes=128         # 减少积分点以加快测试
    )
    
    println("✓ 目标函数闭包创建成功")
    
    # 测试不同的优化参数组合
    test_optimization_params = [
        (0.15, 16.0, 240.0, 0.7, 32.0),   # 标准参数组
        (0.16, 15.5, 250.0, 0.75, 30.0),  # 变化参数组1
        (0.14, 16.5, 230.0, 0.65, 34.0)   # 变化参数组2
    ]
    
    println("\n测试不同优化参数:")
    for (i, params) in enumerate(test_optimization_params)
        ρ₀, B_A, K, m_ratio, E_sym = params
        println("\n参数组 $i:")
        println("  ρ₀ = $ρ₀ fm⁻³")
        println("  B_A = $B_A MeV") 
        println("  K = $K MeV")
        println("  m_ratio = $m_ratio")
        println("  E_sym = $E_sym MeV")
        
        # 使用闭包计算目标函数值
        println("  计算目标函数值...")
        result = objective_func(params)
        
        if isfinite(result)
            @printf("  ✓ 目标函数值: %.2f MeV²\n", result)
        else
            println("  ✗ 目标函数值: $result (无效)")
        end
    end
    
    # 测试加权版本
    println("\n" * "-"^60)
    println("测试加权目标函数闭包")
    println("-"^60)
    
    weights = [1.0, 2.0, 0.5]  # 示例权重
    println("权重设置: $weights")
    
    weighted_objective_func = create_weighted_temperature_difference_objective(
        experimental_kappa_pairs, weights, experimental_μ_B, experimental_T_min, experimental_T_max;
        T_step_scan=3.0/hc,
        verbose=false,
        penalty_for_missing=1e4,
        n_nodes=128
    )
    
    println("✓ 加权目标函数闭包创建成功")
    
    # 测试加权版本（只测试前两个参数组以节省时间）
    println("\n测试加权目标函数:")
    for (i, params) in enumerate(test_optimization_params[1:2])
        ρ₀, B_A, K, m_ratio, E_sym = params
        println("\n参数组 $i:")
        @printf("  ρ₀=%.2f, B_A=%.1f MeV, K=%.1f MeV, m_ratio=%.2f, E_sym=%.1f MeV\n", 
                ρ₀, B_A, K, m_ratio, E_sym)
        
        weighted_result = weighted_objective_func(params)
        
        if isfinite(weighted_result)
            @printf("  ✓ 加权目标函数值: %.2f MeV²\n", weighted_result)
        else
            println("  ✗ 加权目标函数值: $weighted_result (无效)")
        end
    end
    
    println("\n✅ 目标函数闭包测试完成!")
    println("📝 注意: 请根据实际实验数据替换以下参数:")
    println("   - experimental_kappa_pairs: 实验观测的κ比值对")
    println("   - experimental_μ_B: 实验条件下的重子化学势")
    println("   - experimental_T_min, experimental_T_max: 实验温度范围")
    println("   - weights: 各组κ比值的实验权重")
    
catch e
    println("\n✗ 目标函数闭包测试失败:")
    println("错误信息: $e")
    println("错误位置: $(stacktrace()[1])")
end

println("\n测试结束时间: $(Dates.now())")
println("="^80)