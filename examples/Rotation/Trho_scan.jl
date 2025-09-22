#!/usr/bin/env julia
# Moved Trho functionality from src/Rotation/Function_Rotation.jl
# Writes output to results/output/Rotation

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using Printf
using NLsolve
using StaticArrays

include("../../src/Rotation/Function_Rotation.jl")
using Function_Rotation: get_nodes, calculate_t_rho, calculate_thermo, hc, SVector

function Trho(T_start,T_end)
    nodes1, nodes2 = get_nodes(128,16)
    omega = 100/hc

    # 输出文件（写入 results/output/Rotation）
    outdir = joinpath(@__DIR__, "..", "..", "results", "output", "Rotation")
    mkpath(outdir)
    outfile = joinpath(outdir, "trho_rotation.csv")

    x_initial = [-2.13,0.06,0.12, 310 / hc]
    x_rho_3 = copy(x_initial)

    open(outfile, "w") do io
        println(io, "T,rho,phi,Phi1,Phi2,mu,pressure,entropy,energy,converged")

        for T in T_start:1/hc:T_end
            x = copy(x_rho_3)
            rho = 6.00
            converged = false
            try
                res = nlsolve(x -> calculate_t_rho(x, T, rho, nodes1,omega), x)
                converged = res.f_converged
                if converged
                    copyto!(x, res.zero)
                    copyto!(x_rho_3, x)
                else
                    @warn "Root finding did not converge for T=$T and rho=$rho"
                end
            catch err
                @warn "Exception in root finding for T=$T and rho=$rho: $err"
                converged = false
            end

            if converged
                x_phi = SVector{3}(x[1:3])
                x_mu = x[4]
                pressure, _, entropy, energy = calculate_thermo(x_phi, x_mu, T, nodes1,omega)
            else
                pressure = NaN
                entropy = NaN
                energy = NaN
            end
            println(io, join([T*hc, rho, x..., pressure, entropy, energy, converged], ","))
            flush(io)

            for rho in 5.99:-0.01:0.10
                try
                    res = nlsolve(x -> calculate_t_rho(x, T, rho, nodes1,omega), x)
                    converged = res.f_converged
                    if converged
                        copyto!(x, res.zero)
                    else
                        @warn "Root finding did not converge for T=$T and rho=$rho"
                    end
                catch err
                    @warn "Exception in root finding for T=$T and rho=$rho: $err"
                    converged = false
                end

                if converged
                    x_phi = SVector{3}(x[1:3])
                    x_mu = x[4]
                    pressure, _, entropy, energy = calculate_thermo(x_phi, x_mu, T, nodes1,omega)
                else
                    pressure = NaN
                    entropy = NaN
                    energy = NaN
                end
                println(io, join([T*hc, rho, x..., pressure, entropy, energy, converged], ","))
                flush(io)
            end
        end
    end

    return nothing
end

# If run directly, provide a small example (commented to avoid accidental long runs)
Trho(20/hc,21/hc)
