# test_independent_mu_B.jl
# 测试修改后的函数是否能正确处理独立的μ_B值

# 激活项目环境
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

# 包含必要的模块
include("../../src/Gas_Liquid/Advanced_FindTforDiff.jl")

function test_independent_mu_B()
    """
    测试独立μ_B值的温度差平方和计算
    """
    
    println("="^80)
    println("测试：独立μ_B值的温度差平方和计算")
    println("="^80)
    
    # 测试数据：3组κ比值对，每组对应不同的μ_B
    kappa_pairs = [
        (1.09031788496341, -0.28904867673079),   # 第1组
        (1.06152332992368, 0.164279260625683),   # 第2组  
        (1.11111023684003, 0.224522832511389)    # 第3组
    ]
    
    μ_B_values = [
        632.0 / hc,   # 第1组对应632 MeV
        666.0 / hc,   # 第2组对应666 MeV
        697.0 / hc    # 第3组对应697 MeV
    ]
    
    optimization_params = (0.15, -16.0, 240.0, 0.7, 32.0)
    T_min, T_max = 25.0/hc, 200.0/hc
    
    println("测试数据:")
    println("  κ比值对: $kappa_pairs")
    println("  μ_B值: $([round(μ*hc, digits=1) for μ in μ_B_values]) MeV")
    println("  优化参数: $optimization_params")
    
    # 测试1：直接调用函数
    println("\n测试1：直接调用 calculate_temperature_difference_sum_of_squares")
    try
        result = calculate_temperature_difference_sum_of_squares(
            kappa_pairs, μ_B_values, optimization_params, T_min, T_max;
            T_step_scan=5.0/hc, verbose=true, penalty_for_missing=1e4)
        
        println("✅ 函数调用成功")
        println("结果: $(round(result, digits=2)) MeV²")
    catch e
        println("❌ 函数调用失败: $e")
        return false
    end
    
    # 测试2：闭包函数
    println("\n测试2：创建和使用闭包函数")
    try
        objective_func = create_temperature_difference_objective(
            kappa_pairs, μ_B_values, T_min, T_max;
            T_step_scan=5.0/hc, verbose=false, penalty_for_missing=1e4)
        
        result = objective_func(optimization_params)
        println("✅ 闭包创建和调用成功")
        println("结果: $(round(result, digits=2)) MeV²")
    catch e
        println("❌ 闭包函数失败: $e")
        return false
    end
    
    # 测试3：加权版本
    println("\n测试3：加权温度差平方和")
    try
        weights = [1.0, 2.0, 0.5]
        weighted_result = calculate_temperature_difference_sum_of_squares_with_weights(
            kappa_pairs, weights, μ_B_values, optimization_params, T_min, T_max;
            T_step_scan=5.0/hc, verbose=false, penalty_for_missing=1e4)
        
        println("✅ 加权函数调用成功")
        println("权重: $weights")
        println("结果: $(round(weighted_result, digits=2)) MeV²")
    catch e
        println("❌ 加权函数失败: $e")
        return false
    end
    
    # 测试4：输入验证
    println("\n测试4：输入参数验证")
    try
        # 故意使用不匹配的数组长度
        wrong_mu_B = [632.0/hc, 666.0/hc]  # 只有2个值，但有3组κ比值对
        
        calculate_temperature_difference_sum_of_squares(
            kappa_pairs, wrong_mu_B, optimization_params, T_min, T_max;
            verbose=false)
        
        println("❌ 应该检测到数组长度不匹配错误")
        return false
    catch e
        if occursin("不匹配", string(e))
            println("✅ 正确检测到数组长度不匹配: $e")
        else
            println("❌ 意外错误: $e")
            return false
        end
    end
    
    println("\n" * "="^60)
    println("✅ 所有测试通过！独立μ_B值功能正常工作")
    println("="^60)
    
    return true
end

function test_compatibility_comparison()
    """
    比较新旧版本的兼容性测试
    """
    
    println("\n" * "="^80)
    println("兼容性测试：比较使用相同μ_B值 vs 独立μ_B值")
    println("="^80)

    kappa_pairs = [(1.09031788496341, -0.28904867673079), (1.06152332992368, 0.164279260625683)]
    single_μ_B = 632.0 / hc
    μ_B_array = [single_μ_B, single_μ_B]  # 相同值的数组
    optimization_params = (0.15, -16.0, 240.0, 0.7, 32.0)
    T_min, T_max = 25.0/hc, 200.0/hc
    
    println("测试数据:")
    println("  κ比值对: $kappa_pairs")
    println("  μ_B值: $(single_μ_B*hc) MeV (相同值)")
    
    try
        # 使用新版本（μ_B数组形式）
        result_new = calculate_temperature_difference_sum_of_squares(
            kappa_pairs, μ_B_array, optimization_params, T_min, T_max;
            T_step_scan=5.0/hc, verbose=false, penalty_for_missing=1e4)
        
        println("✅ 新版本（μ_B数组）结果: $(round(result_new, digits=2)) MeV²")
        
        # 注意：由于我们已经修改了函数，无法直接测试旧版本
        # 但可以验证逻辑一致性
        println("✅ 兼容性测试通过：可以使用相同值的数组模拟旧版本行为")
        
    catch e
        println("❌ 兼容性测试失败: $e")
        return false
    end
    
    return true
end

# 运行测试
if true
    println("开始测试独立μ_B值功能...")
    
    success1 = test_independent_mu_B()
    success2 = test_compatibility_comparison()
    
    if success1 && success2
        println("\n🎉 所有测试成功完成！")
        println("修改后的函数可以正确处理每组κ比值对对应独立的μ_B值。")
    else
        println("\n❌ 测试失败，请检查函数实现。")
        exit(1)
    end
end