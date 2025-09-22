import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using Test

# 包含要测试的实现
include(joinpath(@__DIR__, "..", "..", "src", "Rotation", "Advanced_FindCep.jl"))

# 把原来在源文件中的测试函数实现移动到这里（源文件不是模块）
function test_find_cep()
    """便捷的测试函数"""
    println("测试find_cep功能...")
    println("注意: 如果当前参数设置下不存在S形相变，函数将返回NaN")

    # 使用较宽松的参数进行测试（可调整以减小运行量）
    T_cep = find_cep(30.0, 100.0, 1.0)

    if isnan(T_cep)
        println("未找到临界温度，建议:")
        println("1. 调整物理参数 (omega, 初始条件等)")
        println("2. 扩大搜索范围")
        println("3. 调整S形检测敏感度")
    else
        println("找到临界温度: T_cep = $T_cep MeV")
    end

    return T_cep
end


# 使用 Test 框架包装调用，主要验证函数可执行且返回一个数值
@testset "Advanced_FindCep basic" begin
    val = test_find_cep()
    @test isfinite(val) || isnan(val)
end
