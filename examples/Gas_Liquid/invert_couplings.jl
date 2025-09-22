# invert_couplings.jl
# 尝试用 NLsolve 求解输入 (ρ0,B_A,K,m_ratio,E_sym) 使得
# calculate_couplings(...) 的输出 (fσ,fω,fρ,fδ,b,c) 在 fδ=0 条件下，匹配目标 (fs,fo,fr,fd=0,b,c)

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using NLsolve
using LinearAlgebra
using Random
try
    using Optim
catch
    println("Optim 未安装，尝试安装...")
    Pkg.add("Optim")
    using Optim
end

include(joinpath(@__DIR__, "../..", "src", "Gas_Liquid", "Constants_Gas_Liquid.jl"))
using .Constants_Gas_Liquid

# 目标输出（用户给定）
fs_target = 15.0
fo_target = 5.423
fr_target = 0.95
fd_target = 0.0
b_target = 0.00692
c_target = -0.0048

Ytarget = [fs_target, fo_target, fr_target, b_target, c_target]

function F_from_x(x)
    ρ0 = x[1]
    B_A = x[2]
    K = x[3]
    m_ratio = x[4]
    E_sym = x[5]
    return calculate_couplings(ρ0, B_A, K, m_ratio, E_sym)
end

function F_from_x5(x)
    tup = calculate_couplings(x[1], x[2], x[3], x[4], x[5])
    return Float64[tup[1], tup[2], tup[3], tup[5], tup[6]]
end

function resid!(res, x)
    if !(all(isfinite, x))
        fill!(res, 1e6); return
    end
    ρ0 = x[1]
    m_ratio = x[4]
    if ρ0 <= 0 || m_ratio <= 0
        fill!(res, 1e6); return
    end
    try
        yvec = F_from_x5(x)
    catch e
        fill!(res, 1e6); return
    end
    if any(z -> !isfinite(z), yvec)
        fill!(res, 1e6); return
    end
    res[1] = yvec[1] - Ytarget[1]
    res[2] = yvec[2] - Ytarget[2]
    res[3] = yvec[3] - Ytarget[3]
    res[4] = yvec[4] - Ytarget[4]
    res[5] = yvec[5] - Ytarget[5]
end

function resid_vec(x)
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

x0 = [0.16, -16.0/hc, 240.0/hc, 0.75, 31.3/hc]
println("初始猜测 x0 = ", x0)

# 诊断...

obj(x) = sum(resid_vec(x) .^ 2)

hconst = Constants_Gas_Liquid.hc
lower = [1e-6, -500.0/hconst, 1.0/hconst, 0.01, 0.1/hconst]
upper = [1.0, 500.0/hconst, 1000.0/hconst, 10.0, 500.0/hconst]

function run_inversion(; nstarts=200, rng_seed=1234)
    Random.seed!(rng_seed)
    best_rnorm = Inf
    best_x = copy(x0)
    for i in 1:nstarts
        init = lower .+ rand(5) .* (upper .- lower)
        rtest = resid_vec(init)
        if any(isnan, rtest) || any(!isfinite, rtest) || sqrt(sum(abs2, rtest)) > 1e5
            continue
        end
        try
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

# 诊断输出
try
    yopt = F_from_x5(best_x)
    println("对应输出 (fs,fo,fr,b,c) = ", yopt)
catch e
    println("无法对最佳解计算诊断: ", e)
end
