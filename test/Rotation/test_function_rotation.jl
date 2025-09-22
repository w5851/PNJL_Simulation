import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using StaticArrays
using Test

# 以非模块方式包含实现（保持与原项目相同的加载风格）
include(joinpath(@__DIR__, "..", "..", "src", "Rotation", "Function_Rotation.jl"))

# 由于源文件是 module，我们需要使用模块名来访问导出的函数
using .Function_Rotation: get_nodes, calculate_pressure, calculate_core, calculate_thermo, calculate_t_rho, hc

@testset "Function_Rotation low-level" begin
    # 使用较小节点以加快测试速度
    nodes1, nodes2 = get_nodes(128, 16)
    omega = 0.0

    phi = -2.13
    Phi1 = 0.06
    Phi2 = 0.13
    T = 100 / hc
    mu = 308.3 / hc

    @testset "calculate_pressure" begin
        p = calculate_pressure(phi, Phi1, Phi2, mu, T, nodes1, omega)
        @show p
        @test isa(p, Number)
        @test !isnan(p)
    end

    @testset "calculate_core and thermo" begin
        x = SVector(phi, Phi1, Phi2)
        core = calculate_core(x, mu, T, nodes1, omega)
        @show core
        @test length(core) == 3
        pressure, rho, entropy, energy = calculate_thermo(x, mu, T, nodes1, omega)
        @show pressure
        @show rho
        @show entropy
        @show energy
        @test isa(pressure, Number)
        @test isa(rho, Number)
        @test isa(entropy, Number)
        @test isa(energy, Number)
    end

    @testset "calculate_t_rho" begin
        x_full = SVector(phi, Phi1, Phi2, mu)
        fvec = calculate_t_rho(x_full, T, 4.0, nodes1, omega)
        @show fvec
        @test length(fvec) == 4
    end
end
