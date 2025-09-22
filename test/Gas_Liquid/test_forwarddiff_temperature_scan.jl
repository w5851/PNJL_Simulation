# test_forwarddiff_temperature_scan.jl
# ForwardDiff温度扫描功能测试
# 测试Advanced_ForwardDiff.jl模块中的ForwardDiff自动微分功能

# 激活项目环境
import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using Printf
using Dates

# 引入Advanced_ForwardDiff模块和相关依赖
include("../../src/Gas_Liquid/Advanced_ForwardDiff.jl")

# 主程序：ForwardDiff温度扫描测试
println("="^80)
println("ForwardDiff自动微分温度扫描测试")
println("测试时间: $(Dates.now())")
println("="^80)


# 固定μ_B值，使用自然单位
μ_B_fixed = 697.0 / hc  # 转换为自然单位：697 MeV

# 温度扫描参数
T_min = 20.0 / hc        # 20 MeV
T_max = 200.0 / hc      # 200 MeV
T_step = 1.0 / hc      # 1 MeV步长

# 输出文件 - 使用绝对路径避免相对路径混淆
project_root = dirname(dirname(@__DIR__))  # 获取项目根目录
output_file = joinpath(project_root, "output", "Gas_Liquid", "forwarddiff_temperature_scan.csv")

println("\n扫描配置:")
println("μ_B = $(μ_B_fixed * hc) MeV")
println("T范围: $(T_min * hc) - $(T_max * hc) MeV")
println("T步长: $(T_step * hc) MeV")
println("输出文件: $output_file")

# 执行ForwardDiff温度扫描
println("\n开始ForwardDiff温度扫描...")
try
    df_results = forwarddiff_temperature_scan(
        μ_B_fixed, T_min, T_max, T_step, output_file;
        gsigma=1.25, gdelta=0.01,      # 场初值
        fs=15.0,                    # σ耦合
        fo=5.423,                    # ω耦合
        fr=0.95,                     # ρ耦合
        fd=0.0,                        # δ耦合
        b=0.00692,                     # 三阶耦合
        c=-0.0048,                    # 四阶耦合
        n_nodes=256                    # 积分节点数
    )
    
    println("\n" * "="^80)
    println("ForwardDiff温度扫描测试完成!")
    println("="^80)
    
    # 输出最终统计信息
    total_points = nrow(df_results)
    successful_points = sum(.!isnan.(df_results.P_T4))
    success_rate = round(successful_points / total_points * 100, digits=1)
    
    println("\n最终统计:")
    println("总计算点数: $total_points")
    println("成功计算点数: $successful_points")
    println("成功率: $success_rate%")
    
    if successful_points > 0
        # 显示一些关键结果
        println("\n关键结果预览 (前5个成功点):")
        valid_rows = filter(row -> !isnan(row.P_T4), df_results)
        for (i, row) in enumerate(eachrow(valid_rows))
            if i <= 5
                @printf("T=%.1f MeV: P/T⁴=%.6f, κ₃/κ₁=%.6f, κ₄/κ₂=%.6f\n", 
                        row.T_MeV, row.P_T4, 
                        isfinite(row.kappa3_over_kappa1) ? row.kappa3_over_kappa1 : NaN,
                        isfinite(row.kappa4_over_kappa2) ? row.kappa4_over_kappa2 : NaN)
            end
        end
    end
    
catch e
    println("\n✗ ForwardDiff温度扫描测试失败:")
    println("错误信息: $e")
    println("错误位置: $(stacktrace()[1])")
end

println("\n测试结束时间: $(Dates.now())")
println("="^80)