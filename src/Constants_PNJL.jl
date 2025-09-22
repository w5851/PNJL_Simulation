module Constants_PNJL

# 热力学常数
const π = 3.141592653
const hc = 197.33
const rho0 = 0.16


# PNJL 模型常数
const T0 = 210 / hc
const Nc = 3.0

const Lambda = 602.3    # !(MeV)
const G_Lam2 = 1.835
const K_Lam5 = 12.36

const m0_q = 5.5       # !(MeV)
const m0_s = 140.7

const Lambda_f = Lambda / hc      # !(fm**(-1))
const G_f = G_Lam2 / Lambda_f^2   # !(fm**(2))
const K_f = K_Lam5 / Lambda_f^5   # !(fm**(5))

const m0_q_f = m0_q / hc       # !(fm**(-1))
const m0_s_f = m0_s / hc

const m0 = [m0_q_f, m0_q_f, m0_s_f]

# emm常数
const a0 = 3.51
const a1 = -2.47
const a2 = 15.2
const b3 = -1.75
const b4 = 7.555
export hc, π, rho0, a0, a1, a2, b3, b4, T0, Nc, Lambda_f, G_f, K_f, m0_q_f, m0_s_f
end # module Constants_PNJL