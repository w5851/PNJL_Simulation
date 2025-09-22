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

println("\n\n【示例2】使用优化参数计算κ比值对应的温度")
println("-"^50)
# NOTICE: This Julia example has been moved to examples/Gas_Liquid/example_find_temperature.jl
println("This example was moved to examples/Gas_Liquid/example_find_temperature.jl")