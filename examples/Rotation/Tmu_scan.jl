#!/usr/bin/env julia
# Moved Tmu functionality from src/Rotation/Function_Rotation.jl
# Writes output to results/output/Rotation

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using Printf
using NLsolve
using StaticArrays

include("../../src/Rotation/Function_Rotation.jl")
using .Function_Rotation: get_nodes, calculate_core, calculate_thermo, calculate_mass, hc
include("../../src/Rotation/Parameter_Rotation.jl")
using .Parameter_Rotation: x_initial_Tmu

function Tmu(T_start; T_end=T_start, T_step=1/hc, mu_start=0/hc, mu_end=0.0, mu_step=-1/hc, omega=100/hc)
    # 节点
    nodes1, nodes2 = get_nodes(128, 16)

    # 输出文件（写入 results/output/Rotation）
    outdir = joinpath(@__DIR__, "..", "..", "results", "output", "Rotation")
    mkpath(outdir)
    outfile = joinpath(outdir, "tmu_rotation.csv")

    # 初始 x（只有 3 个变量：phi, Phi1, Phi2），从模块中读取
    x_initial = x_initial_Tmu
    x_prev = copy(x_initial)

    open(outfile, "w") do io
        println(io, "T,mu,phi,Phi1,Phi2,mass,pressure,rho,entropy,energy,converged")

        for T in T_start:T_step:T_end
            # 对每个 T，在 mu 方向扫描
            for mu in mu_start:mu_step:mu_end
                x = copy(x_prev)
                converged = false
                pressure = NaN
                entropy = NaN
                energy = NaN
                rho = NaN
                mass = NaN
                try
                    # 使用 calculate_core 求解 3 个未知量
                    res = nlsolve(x -> calculate_core(x, mu, T, nodes1, omega), x; autodiff = :forward)
                    converged = res.f_converged
                    if converged
                        copyto!(x, res.zero)
                        x_prev .= x  # 用当前解作为下一点的初始值
                        # 计算热力学量（pressure, entropy, energy）
                        pressure, rho, entropy, energy = calculate_thermo(x, mu, T, nodes1, omega)
                        # 计算有效质量
                        mass = calculate_mass(x[1])  # x[1] 是 phi
                    else
                        @warn "Root finding did not converge for T=$T and mu=$mu"
                    end
                catch err
                    @warn "Exception in root finding for T=$T and mu=$mu: $err"
                    converged = false
                end

                # 写入 csv：将 T, mu 转换回物理单位（乘 hc）以与 Trho 保持一致
                println(io, join([T*hc, mu*hc, x..., mass, pressure, rho, entropy, energy, converged], ","))
                flush(io)
            end
        end
    end

    return nothing
end

# 示例直接运行（注释掉以避免长时间运行）
Tmu(50/hc; T_end=60/hc, T_step=1/hc, mu_start=0/hc, mu_end=0.0, mu_step=-1/hc)
