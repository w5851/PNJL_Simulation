module Constants_Rotation

# PNJL 常数
const π = 3.141592653
const hc = 197.33
const rho0 = 0.16
const T0 = 270 / hc

# PNJL 模型常数
const Nc = 3.0  # 颜色数 Nc

const Lambda = 650.0    # !(MeV)
const G_Lam2 = 4.93 / 1e6  # (MeV), 自耦合
const K_Lam5 = 12.36      # (MeV), 三味耦合

const m0_q = 5.0       # (MeV)
const m0_s = 5.0
const m0_q_f = m0_q / hc   # (fm⁻¹)
const m0_s_f = m0_s / hc   # (fm⁻¹)

const Lambda_f = Lambda / hc      # !(fm**(-1))
const G_f = G_Lam2 * hc^2   # (fm²)
const K_f = K_Lam5 / Lambda_f^5  # (fm⁵)

# Polyakov loop 参数
const a0 = 6.75
const a1 = -1.95
const a2 = 2.625
const a3 = -7.44
const b3 = 0.75
const b4 = 7.5

# 强旋转常数
const r0 = 0.1 / 1000 * hc
const C = 4.0

# f(T, ω) 的系数表
const coefficients = Dict(
    "a" => [0.0454431, -1.27942e-5, -5.43339e-9],
    "b" => [46.8263, -0.0210165, -2.15394e-5],
    "c" => [1.00298, 1.55157e-4, -5.99032e-8],
    "d" => [0.0600157, -5.74388e-6, -8.24192]
)
export hc, π, rho0, a0, a1, a2, b3, b4, T0, Nc, Lambda_f, G_f, K_f, m0_q_f, m0_s_f, r0, coefficients
end # module Constants_Rotation
