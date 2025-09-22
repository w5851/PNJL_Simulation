# example_find_temperature_old.jl
# 旧版 温度反向查找使用示例
# （内容与原文件一致）

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using Printf

include("../../src/Gas_Liquid/Advanced_FindTforDiff.jl")

println("演示（旧版）完成")
