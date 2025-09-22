module Parameter_Rotation

include("Constants_Rotation.jl")
using .Constants_Rotation: hc
# 初始参数（以自然单位 fm^{-1} 为单位，与代码中使用的一致）
# 对应 `Tmu_scan.jl` 中的 x_initial = [-2.13, 0.06, 0.12]
const x_initial_Tmu = [-2.13, 0.06, 0.12]

# 对应 `Trho_scan.jl` 中的 x_initial = [-2.13,0.06,0.12, 310 / hc]
# 注意最后一项是 mu（以 MeV 为单位除以 hc 转换为 fm^{-1}）
const x_initial_Trho = [-2.13, 0.06, 0.12, 310 / hc]

export x_initial_Tmu, x_initial_Trho

end # module Parameter_Rotation
