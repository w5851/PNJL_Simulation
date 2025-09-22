# test_bayesian_optimization.jl
# 测试Advanced_BayesianOptimization.jl模块的贝叶斯优化功能

# 激活项目环境
import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

# 引入贝叶斯优化模块
include("../../src/Gas_Liquid/Advanced_BayesianOptimization.jl")

demo_bayesian_optimization_with_warmup()