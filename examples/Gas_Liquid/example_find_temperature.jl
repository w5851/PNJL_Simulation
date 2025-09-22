# example_find_temperature.jl
# 温度反向查找使用示例
# 演示如何使用Advanced_FindTforDiff.jl模块查找特定κ比值对应的温度

# 激活项目环境
import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using Printf

# 引入温度查找模块
include("../../src/Gas_Liquid/Advanced_FindTforDiff.jl")

println("="^70)
println("温度反向查找功能演示")
println("="^70)

# 示例1：查找单组κ比值对应的温度
println("\n【示例1】查找单组κ比值对应的温度")
println("-"^50)

# 设置物理参数
μ_B = 697.0 / hc        # 重子化学势 697 MeV
T_min = 30.0 / hc       # 最小温度 30 MeV  
T_max = 180.0 / hc      # 最大温度 180 MeV

# 目标κ比值
target_kappa3_kappa1 = 0.8
target_kappa4_kappa2 = 1.5

println("物理设置:")
println("  μ_B = $(μ_B * hc) MeV")
println("  搜索范围: $(T_min * hc) - $(T_max * hc) MeV")
println("  目标 κ₃/κ₁ = $target_kappa3_kappa1")
println("  目标 κ₄/κ₂ = $target_kappa4_kappa2")

# 执行查找
T_k31, T_k42 = find_temperature_for_kappa_ratios(
    target_kappa3_kappa1, target_kappa4_kappa2, μ_B, 
    T_min, T_max, 1.0/hc;  # 1 MeV 步长
    verbose=false  # 简化输出
)

println("\n查找结果:")
if !isnan(T_k31)
    @printf("  κ₃/κ₁ = %.1f → T = %.1f MeV\n", target_kappa3_kappa1, T_k31 * hc)
else
    println("  κ₃/κ₁ = $target_kappa3_kappa1 → 未找到对应温度")
end

if !isnan(T_k42)
    @printf("  κ₄/κ₂ = %.1f → T = %.1f MeV\n", target_kappa4_kappa2, T_k42 * hc)
else
    println("  κ₄/κ₂ = $target_kappa4_kappa2 → 未找到对应温度")
end

# 示例2：使用优化参数计算κ比值对应的温度
println("\n\n【示例2】使用优化参数计算κ比值对应的温度")
println("-"^50)

# 定义优化参数
optimization_params = (
    0.16,     # ρ₀: 核饱和密度 (fm⁻³)
    -16.0,    # B_A: 结合能 (MeV)
    240.0,    # K: 不可压缩模量 (MeV)
    0.75,     # m_ratio: 有效质量比
    31.3      # E_sym: 对称能 (MeV)
)

println("优化参数设置:")
println("  ρ₀ = $(optimization_params[1]) fm⁻³")
println("  B_A = $(optimization_params[2]) MeV")
println("  K = $(optimization_params[3]) MeV")
println("  m_ratio = $(optimization_params[4])")
println("  E_sym = $(optimization_params[5]) MeV")
println("  目标 κ₃/κ₁ = $target_kappa3_kappa1")
println("  目标 κ₄/κ₂ = $target_kappa4_kappa2")

# 执行查找
T_k31_opt, T_k42_opt = find_temperature_for_kappa_ratios_with_optimization_params(
    target_kappa3_kappa1, target_kappa4_kappa2, μ_B, optimization_params,
    T_min, T_max, 1.0/hc;
    verbose=false  # 简化输出
)

println("\n优化参数查找结果:")
if !isnan(T_k31_opt)
    @printf("  κ₃/κ₁ = %.1f → T = %.1f MeV\n", target_kappa3_kappa1, T_k31_opt * hc)
else
    println("  κ₃/κ₁ = $target_kappa3_kappa1 → 未找到对应温度")
end

if !isnan(T_k42_opt)
    @printf("  κ₄/κ₂ = %.1f → T = %.1f MeV\n", target_kappa4_kappa2, T_k42_opt * hc)
else
    println("  κ₄/κ₂ = $target_kappa4_kappa2 → 未找到对应温度")
end

# 计算温度差
if !isnan(T_k31_opt) && !isnan(T_k42_opt)
    temp_diff = abs(T_k31_opt - T_k42_opt) * hc
    @printf("  温度差 |ΔT| = %.2f MeV\n", temp_diff)
end

# 示例3：批量查找多组κ比值
println("\n\n【示例3】批量查找多组κ比值对应的温度")
println("-"^50)

# 定义多组目标κ比值
kappa_pairs = [
    (0.3, 0.8),   # (κ₃/κ₁, κ₄/κ₂)
    (0.5, 1.0),
    (0.7, 1.2),
    (1.0, 1.5)
]

println("批量查找设置:")
println("  κ比值对数: $(length(kappa_pairs)) 组")
println("  μ_B = $(μ_B * hc) MeV")

# 执行批量查找
try
    df_results = batch_find_temperatures_with_optimization_params(
        kappa_pairs, μ_B, optimization_params, T_min, T_max;
        T_step_scan=2.0/hc,  # 2 MeV 步长（加快速度）
        output_file=joinpath(@__DIR__, "../../output/Gas_Liquid/batch_temperature_finder_example.csv")
    )
    
    println("\n批量查找结果:")
    println("  κ₃/κ₁   κ₄/κ₂   T(κ₃/κ₁) [MeV]   T(κ₄/κ₂) [MeV]   ΔT [MeV]")
    println("  " * "-"^55)
    
    for row in eachrow(df_results)
        T_k31_str = isnan(row.T_for_kappa3_kappa1_MeV) ? "    N/A" : @sprintf("%7.1f", row.T_for_kappa3_kappa1_MeV)
        T_k42_str = isnan(row.T_for_kappa4_kappa2_MeV) ? "    N/A" : @sprintf("%7.1f", row.T_for_kappa4_kappa2_MeV)
        Delta_T_str = isnan(row.temperature_difference_MeV) ? "   N/A" : @sprintf("%6.1f", row.temperature_difference_MeV)
        
        @printf("  %5.1f   %5.1f   %s       %s       %s\n", 
                row.kappa3_over_kappa1_target, row.kappa4_over_kappa2_target, 
                T_k31_str, T_k42_str, Delta_T_str)
    end
    
    println("\n✅ 结果已保存到: output/Gas_Liquid/batch_temperature_finder_example.csv")
    
catch e
    println("\n❌ 批量查找失败: $e")
end

# 示例4：物理解释和应用建议
println("\n\n【示例4】物理意义和应用建议")
println("-"^50)

println("物理意义:")
println("  • κ₃/κ₁: 反映重子数涨落的偏度，是相变的重要信号")
println("  • κ₄/κ₂: 反映重子数涨落的峰度，临界点附近会有显著变化")
println("  • 温度查找: 帮助确定特定涨落特征对应的相变温度")

println("\n应用场景:")
println("  1. 相变研究: 寻找临界点对应的温度")
println("  2. 实验对比: 将理论计算与实验数据进行匹配")
println("  3. 参数扫描: 在参数空间中寻找特定物理条件")
println("  4. 相图绘制: 确定相边界位置")

println("\n使用建议:")
println("  • 温度搜索范围应覆盖感兴趣的相变区域")
println("  • 扫描步长影响精度：1-2 MeV适合大多数应用")
println("  • 如果存在多个交点，物理上通常选择右侧（高温）交点")
println("  • 建议先做粗扫描确定大致范围，再精细搜索")
println("  • 优化参数方法可用于贝叶斯优化等参数调节任务")

println("\n" * "="^70)
println("演示完成！")
println("="^70)
