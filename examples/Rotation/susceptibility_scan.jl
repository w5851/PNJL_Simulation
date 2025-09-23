#!/usr/bin/env julia
# Scan n-th generalized susceptibility over mu or omega (or both) for the Rotation model
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
include("../../src/Rotation/Advanced_ForwardDiff.jl")

"""
susceptibility_scan(T; kwargs...)

Scans generalized susceptibilities (w.r.t mu_B or omega) over a parameter range.

Keyword args:
- mode = :mu or :omega or :both (what to scan)
- n = integer order of susceptibility (default 1)
- mu_range = (mu_start, mu_end, mu_step)
- omega_range = (omega_start, omega_end, omega_step)
- T_range = (T_start, T_end, T_step)
- nodes = (p_num, t_num)
- outname = filename suffix (default "susceptibility_scan.csv")

Behavior:
- x is iterated similarly to Tmu: for adjacent parameter steps, the previous solution is used as initial guess.
- Writes CSV with columns: T, mu, omega, n, susceptibility, converged, phi, Phi1, Phi2
"""
function susceptibility_scan(T_start; T_end=T_start, T_step=1/hc,
                              mode::Symbol = :mu,
                              n::Integer = 1,
                              mu_range = (0/hc, 0.0, -1/hc),
                              omega_range = (100/hc, 100/hc, 0.0),
                              nodes=(128,16),
                              outname = "susceptibility_scan.csv")

    # nodes
    nodes1, nodes2 = get_nodes(nodes[1], nodes[2])

    # output
    outdir = joinpath(@__DIR__, "..", "..", "results", "output", "Rotation")
    mkpath(outdir)
    outfile = joinpath(outdir, outname)

    # initial x
    x_initial = x_initial_Tmu
    x_prev = copy(x_initial)

    mu_start, mu_end, mu_step = mu_range
    omega_start, omega_end, omega_step = omega_range

    open(outfile, "w") do io
        println(io, "T,mu,omega,n,susceptibility,converged,phi,Phi1,Phi2")

        for T in T_start:T_step:T_end
            # depending on mode, choose loops
            mu_vals = mu_start:mu_step:mu_end
            omega_vals = omega_start:omega_step:omega_end

            # if mode is :mu only, omega_vals is single value; vice versa for :omega
            for mu in mu_vals
                for omega in omega_vals
                    x = copy(x_prev)
                    converged = false
                    susc = NaN
                    try
                        # solve core (3 variables) using current params
                        res = nlsolve(x -> calculate_core(x, mu, T, nodes1, omega), x; autodiff=:forward)
                        converged = res.f_converged
                        if converged
                            copyto!(x, res.zero)
                            x_prev .= x

                            # compute susceptibility depending on mode
                            if mode == :mu || mode == :both
                                susc = generalized_susceptibility_mu(n, x, mu, T, nodes1, omega)
                            elseif mode == :omega
                                susc = generalized_susceptibility_omega(n, x, mu, T, nodes1, omega)
                            else
                                @warn "Unknown mode: $mode"
                            end
                        else
                            @warn "Root finding did not converge for T=$T, mu=$mu, omega=$omega"
                        end
                    catch err
                        @warn "Exception in root finding or susceptibility eval for T=$T, mu=$mu, omega=$omega: $err"
                        converged = false
                    end

                    println(io, join([T*hc, mu*hc, omega*hc, n, susc, converged, x...,], ","))
                    flush(io)
                end
            end
        end
    end

    return nothing
end

# example invocation (commented out to avoid long runs)
susceptibility_scan(80/hc; T_end=200/hc, mode=:omega, n=1, mu_range=(200/hc, 200/hc, -1/hc), omega_range=(0.02,0.02,0.1))
