# 激活项目环境
import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))


# 引入温度查找模块
include("../../src/Gas_Liquid/Advanced_FindTforDiff.jl")

res = find_temperature_for_kappa_ratios(1.11111023684003, 0.224522832511389, 697.0/hc, 20.0/hc, 200.0/hc, 2.0/hc;
    gsigma=1.25, gdelta=0.01,
    fs=15.0,                    # σ耦合
    fo=5.423,                    # ω耦合
    fr=0.95,                     # ρ耦合
    fd=0.0,                        # δ耦合
    b=0.00692,                     # 三阶耦合
    c=-0.0048,
    n_nodes=256,
    verbose=true)
#666-1.06152332992368-0.164279260625683
res2 = find_temperature_for_kappa_ratios(1.06152332992368, 0.164279260625683, 666.0/hc, 20.0/hc, 200.0/hc, 2.0/hc;
    gsigma=1.25, gdelta=0.01,
    fs=15.0,                    # σ耦合
    fo=5.423,                    # ω耦合
    fr=0.95,                     # ρ耦合
    fd=0.0,                        # δ耦合
    b=0.00692,                     # 三阶耦合
    c=-0.0048,
    n_nodes=256,
    verbose=true)
#632-1.09031788496341--0.28904867673079
res3 = find_temperature_for_kappa_ratios(1.09031788496341, -0.28904867673079, 632.0/hc, 20.0/hc, 200.0/hc, 2.0/hc;
    gsigma=1.25, gdelta=0.01,
    fs=15.0,                    # σ耦合
    fo=5.423,                    # ω耦合
    fr=0.95,                     # ρ耦合
    fd=0.0,                        # δ耦合
    b=0.00692,                     # 三阶耦合
    c=-0.0048,
    n_nodes=256,
    verbose=true)
println("结果: $res")
println("结果: $res2")
println("结果: $res3")