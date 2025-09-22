module Constants_Gas_Liquid
using SpecialFunctions: log, exp

const π = 3.141592653
const hc = 197.33
const m = 939.0 / hc  # MeV to fm⁻¹
const mσ = 550.0 / hc  # MeV to fm⁻¹
const mω = 783.0 / hc  # MeV to fm⁻¹
const mρ = 775.0 / hc  # MeV to fm⁻¹
const mδ = 980.0 / hc  # MeV to fm⁻¹
"""
const rho0 = 0.16
const B_A = -16.0 / hc  # MeV to fm⁻¹
const K = 240.0 / hc  # MeV to fm⁻¹
const E_sym = 31.3 / hc  # MeV to fm⁻¹
const ratio = 0.75  # 有效质量比
"""
function calculate_couplings(ρ0,B_A,K,m_ratio,E_sym)
    """输入参数的单位已经转换(/hc)"""     
    # 初始化变量
    meff = m_ratio * m  # 有效质量 (MeV)
    fδ = 0.0  # δ 耦合常数 fδ=(gδ/mδ)^2(固定为0.0)
        
    # 计算 Fermi 动量 kF 和 Fermi 能量 EF
    kF = (1.5π^2 * ρ0)^(1/3)
    EF = sqrt(kF^2 + meff^2)
    
    # 计算 sigma 场强度 gsigma
    gσ = m - meff
    
    # 计算 omega 耦合常数 fω
    fω = (m + B_A - EF) / ρ0
    
    # 计算中间变量 x, t, I1, I2, I3
    x = kF / meff
    t =  sqrt(1 + x^2)

    term1 =  0.5x * t + x / t - 1.5log(x + t)
    I1 =  (2 / π^2) * meff^2 * term1
    
    term2 =  0.25(x * t^3 - 0.5x * t - 0.5log(x + t))
    I2 =  (2 / π^2) * meff^4 * term2
    
    term3 =  0.5(x * t - log(x + t))
    I3 =  (2 / π^2) * meff^3 * term3
    
    # 计算 α, β, γ, δ 系数
    α1 =  K - fω * (6kF^3 / π^2) - 3kF^2 / EF
    β1 =  2m * gσ * α1
    γ1 =  3gσ^2 * α1
    δ1 =  -(6kF^3 / π^2) * (meff / EF)^2 - α1 * I1
    
    α2 =  0.5gσ^2
    β2 =  (1/3)m * gσ^3
    γ2 =  0.25gσ^4
    δ2 =  ρ0 * (m + B_A) - 0.5fω * ρ0^2 - I2
    
    α3 = gσ
    β3 =  m * gσ^2
    γ3 =  gσ^3
    δ3 = I3
    
    # 计算参数 c 和 b
    denom1 =  (γ1 * α2 - γ2 * α1) * (β2 * α3 - β3 * α2) -
             (γ2 * α3 - γ3 * α2) * (β1 * α2 - β2 * α1)
             
    num1 =  (δ1 * α2 - δ2 * α1) * (β2 * α3 - β3 * α2) -
           (δ2 * α3 - δ3 * α2) * (β1 * α2 - β2 * α1)
           
    c =  num1 / denom1
    
    denom2 =  β1 * α2 - β2 * α1
    num2 =  (δ1 * α2 - δ2 * α1) - (γ1 * α2 - γ2 * α1) * c
    b =  num2 / denom2
    
    # 计算 sigma 耦合常数 fσ
    fσ =  α1 / (δ1 - β1 * b - γ1 * c)
    
    # 计算 rho 耦合常数 fρ
    term_fρ1 =  (meff / EF)^2
    term_fρ2 =  1 + 0.25fδ * I1
    fρ =  (2 / ρ0) * (E_sym - (1/6)*(kF^2 / EF) + (ρ0 / 8) * fδ * term_fρ1 / term_fρ2)
    
    return (fσ, fω, fρ, fδ, b, c)
end

export π, hc, m, mσ, mω, mρ, mδ, calculate_couplings
end # module Constants_Gas_Liquid