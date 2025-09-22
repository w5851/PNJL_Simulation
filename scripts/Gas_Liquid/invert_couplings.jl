# NOTICE: invert_couplings.jl moved to examples/Gas_Liquid/invert_couplings.jl
println("This script was moved to examples/Gas_Liquid/invert_couplings.jl")

# 函数 wrapper：给定 x 返回 (fs,fo,fr,fd,b,c)
function F_from_x(x)
    ρ0 = x[1]
    B_A = x[2]
    K = x[3]
    m_ratio = x[4]
    E_sym = x[5]
    return calculate_couplings(ρ0, B_A, K, m_ratio, E_sym)
end

# 返回不包含 fδ 的向量 [fs, fo, fr, b, c]
function F_from_x5(x)
    tup = calculate_couplings(x[1], x[2], x[3], x[4], x[5])
    return Float64[tup[1], tup[2], tup[3], tup[5], tup[6]]
end

# 残差函数 r(x) of length 5
function resid!(res, x)
    # domain checks and finite checks
    if !(all(isfinite, x))
        fill!(res, 1e6); return
    end
    ρ0 = x[1]
    m_ratio = x[4]
    if ρ0 <= 0 || m_ratio <= 0
        fill!(res, 1e6); return
    end
    # compute outputs safely (ignore fδ)
    try
        yvec = F_from_x5(x)
    catch e
        fill!(res, 1e6); return
    end
    if any(z -> !isfinite(z), yvec)
        fill!(res, 1e6); return
    end
    # yvec = [fs, fo, fr, b, c]
    res[1] = yvec[1] - Ytarget[1]
    res[2] = yvec[2] - Ytarget[2]
    res[3] = yvec[3] - Ytarget[3]
    res[4] = yvec[4] - Ytarget[4]  # b
    res[5] = yvec[5] - Ytarget[5]  # c
end

# wrapper 返回残差向量（供 nlsolve 在 autodiff=false 时使用）
function resid_vec(x)
    # self-contained residual vector to avoid closure/scope issues
    if !(all(isfinite, x))
        return fill(1e6, 5)
    end
    ρ0 = x[1]
    m_ratio = x[4]
    if ρ0 <= 0 || m_ratio <= 0
        return fill(1e6, 5)
    end
    y = try
        F_from_x5(x)
    catch
        return fill(1e6, 5)
    end
    if any(z -> !isfinite(z), y)
        return fill(1e6, 5)
    end
    return y .- Ytarget
end

# 初值（物理常见值，单位注意与 module 内一致）
# use constants from module
x0 = [0.16, -16.0/hc, 240.0/hc, 0.75, 31.3/hc]
println("初始猜测 x0 = ", x0)
# 诊断：打印初值的输出与残差（安全评估）
y0 = try
    F_from_x5(x0)
catch e
    println("评估 F(x0) 失败: ", e)
    nothing
end

if y0 !== nothing
    println("F(x0) = ", y0)
    println("resid(x0) = ", resid_vec(x0))
else
    println("无法评估初值的输出/残差")
end

# 目标函数（最小化残差平方和）
obj(x) = sum(resid_vec(x) .^ 2)

# 设定变量上下界（使用与 x0 同阶的合理范围）
hconst = Constants_Gas_Liquid.hc
lower = [1e-6, -500.0/hconst, 1.0/hconst, 0.01, 0.1/hconst]
upper = [1.0, 500.0/hconst, 1000.0/hconst, 10.0, 500.0/hconst]

function run_inversion(; nstarts=200, rng_seed=1234)
    Random.seed!(rng_seed)
    best_rnorm = Inf
    best_x = copy(x0)
    for i in 1:nstarts
        init = lower .+ rand(5) .* (upper .- lower)
        # 先快速评估起点残差，跳过明显不可行的点
        rtest = resid_vec(init)
        if any(isnan, rtest) || any(!isfinite, rtest) || sqrt(sum(abs2, rtest)) > 1e5
            continue
        end
        try
            # 使用有界优化（Fminbox + BFGS），让实现自身做有限差分梯度
            resopt = optimize(obj, lower, upper, init, Fminbox(BFGS()))
            xopt = Optim.minimizer(resopt)
            r = resid_vec(xopt)
            rnorm = sqrt(sum(abs2, r))
            println("尝试起点 #", i, " -> rnorm=", rnorm)
            if rnorm < best_rnorm
                best_rnorm = rnorm
                best_x = xopt
            end
        catch e
            println("优化起点出错(#", i, "): ", e)
        end
    end

    println("最佳残差范数 = ", best_rnorm)
    println("最佳 x = ", best_x)
    return best_x, best_rnorm
end

best_x, best_rnorm = run_inversion(nstarts=200, rng_seed=20250918)

# 对最佳解做诊断输出
try
    yopt = F_from_x5(best_x)
    println("对应输出 (fs,fo,fr,b,c) = ", yopt)
    println("目标 = ", Ytarget)
    println("残差向量 = ", yopt - Ytarget)
    println("残差范数 = ", sqrt(sum(abs2, yopt - Ytarget)))

    ρ0,B_A,K,m_ratio,E_sym = best_x
    C = Constants_Gas_Liquid
    meff = m_ratio * C.m
    kF = (1.5 * C.π^2 * ρ0)^(1/3)
    EF = sqrt(kF^2 + meff^2)
    gσ = C.m - meff
    fω = (C.m + B_A - EF) / ρ0
    x = kF / meff
    t = sqrt(1 + x^2)
    term1 = 0.5x * t + x / t - 1.5 * log(x + t)
    I1 = (2 / C.π^2) * meff^2 * term1
    term2 = 0.25 * (x * t^3 - 0.5x * t - 0.5 * log(x + t))
    I2 = (2 / C.π^2) * meff^4 * term2
    term3 = 0.5 * (x * t - log(x + t))
    I3 = (2 / C.π^2) * meff^3 * term3
    α1 = K - fω * (6kF^3 / C.π^2) - 3kF^2 / EF
    β1 = 2C.m * gσ * α1
    γ1 = 3gσ^2 * α1
    δ1 = -(6kF^3 / C.π^2) * (meff / EF)^2 - α1 * I1
    α2 = 0.5gσ^2
    β2 = (1/3)C.m * gσ^3
    γ2 = 0.25gσ^4
    δ2 = ρ0 * (C.m + B_A) - 0.5fω * ρ0^2 - I2
    α3 = gσ
    β3 = C.m * gσ^2
    γ3 = gσ^3
    δ3 = I3
    denom1 = (γ1 * α2 - γ2 * α1) * (β2 * α3 - β3 * α2) - (γ2 * α3 - γ3 * α2) * (β1 * α2 - β2 * α1)
    denom2 = β1 * α2 - β2 * α1
    println("诊断: ρ0=", ρ0, ", m_ratio=", m_ratio, ", EF=", EF)
    println("denom1=", denom1)
    println("denom2=", denom2)
catch e
    println("无法对最佳解计算诊断: ", e)
end
