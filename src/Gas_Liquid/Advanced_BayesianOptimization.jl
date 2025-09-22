# Advanced_BayesianOptimization.jl
# 基于贝叶斯优化的PNJL模型参数优化模块
# 使用 BayesianOptimization.jl 优化温度差平方和目标函数
# 依赖于 Advanced_FindTforDiff.jl 中的闭包函数

include("Advanced_FindTforDiff.jl")


using BayesianOptimization
using GaussianProcesses
using Printf
using Dates
using CSV
using DataFrames

# 导入必要的函数用于高斯过程更新
import GaussianProcesses: update!, initialise_model!

"""
优化参数结构体，定义五个PNJL模型参数
"""
struct PNJLOptimizationParams
    ρ₀::Float64      # 核饱和密度 (fm⁻³)
    B_A::Float64     # 结合能 (MeV)
    K::Float64       # 不可压缩模量 (MeV)
    m_ratio::Float64 # 有效质量比
    E_sym::Float64   # 对称能 (MeV)
end

# 参数转换工具函数
params_to_vector(p::PNJLOptimizationParams) = [p.ρ₀, p.B_A, p.K, p.m_ratio, p.E_sym]
vector_to_params(v::Vector{Float64}) = PNJLOptimizationParams(v[1], v[2], v[3], v[4], v[5])

function create_bayesian_objective_function(
    kappa_pairs, μ_B_values, T_min, T_max;
    T_step_scan=1.0/hc, gsigma=1.25, gdelta=0.01, n_nodes=256,
    penalty_for_missing=1e6, verbose=false)
    """
    创建用于贝叶斯优化的目标函数闭包
    
    参数:
    - kappa_pairs: κ比值对数组，格式 [(κ₃/κ₁, κ₄/κ₂), ...]
    - μ_B_values: 重子化学势数组，与 kappa_pairs 一一对应
    - T_min, T_max: 温度搜索范围
    - T_step_scan: 扫描步长
    - 其他参数: 模型计算参数
    
    返回:
    - objective_function: 接受参数向量 [ρ₀, B_A, K, m_ratio, E_sym] 的目标函数
    """
    
    println("="^80)
    println("创建贝叶斯优化目标函数")
    println("="^80)
    println("实验数据配置:")
    println("  κ比值对数量: $(length(kappa_pairs))")
    println("  μ_B值: $([round(μ*hc, digits=1) for μ in μ_B_values]) MeV")
    println("  温度范围: $(T_min*hc) - $(T_max*hc) MeV")
    println("  扫描步长: $(T_step_scan*hc) MeV")
    println("  惩罚值: $penalty_for_missing")
    
    # 创建基础闭包函数
    base_objective = create_temperature_difference_objective(
        kappa_pairs, μ_B_values, T_min, T_max;
        T_step_scan=T_step_scan, gsigma=gsigma, gdelta=gdelta, 
        n_nodes=n_nodes, penalty_for_missing=penalty_for_missing, verbose=verbose)
    
    # 包装为接受向量参数的函数（贝叶斯优化要求）
    function bayesian_objective(param_vector::Vector{Float64})
        try
            # 将向量转换为参数元组
            ρ₀, B_A, K, m_ratio, E_sym = param_vector
            optimization_params = (ρ₀, B_A, K, m_ratio, E_sym)
            
            # 调用基础目标函数
            result = base_objective(optimization_params)
            
            # 确保返回有限值
            if !isfinite(result)
                return penalty_for_missing
            end
            
            return result
            
        catch e
            # 如果计算失败，返回惩罚值
            if verbose
                println("目标函数评估失败: $e")
            end
            return penalty_for_missing
        end
    end
    
    return bayesian_objective
end

function warmup_objective_function(
    objective_function, lower_bounds, upper_bounds;
    n_warmup_samples=3, verbose=true)
    """
    预热目标函数：通过预计算加速后续贝叶斯优化
    
    参数:
    - objective_function: 目标函数
    - lower_bounds: 参数下界
    - upper_bounds: 参数上界
    - n_warmup_samples: 预热样本数
    - verbose: 是否显示详细信息
    
    返回:
    - warmup_results: 预热计算结果
    - estimated_time_per_eval: 估算的单次评估时间
    """
    
    if verbose
        println("\n" * "="^60)
        println("目标函数预热计算")
        println("="^60)
        println("预热样本数: $n_warmup_samples")
    end
    
    warmup_results = []
    total_time = 0.0
    
    for i in 1:n_warmup_samples
        # 生成随机参数（在边界内）
        params = lower_bounds .+ (upper_bounds .- lower_bounds) .* rand(length(lower_bounds))
        
        if verbose
            println("\n预热 $i/$n_warmup_samples:")
            println("  参数: [$(join([round(p, digits=4) for p in params], ", "))]")
        end
        
        # 计时评估
        start_time = time()
        try
            result = objective_function(params)
            eval_time = time() - start_time
            total_time += eval_time
            
            push!(warmup_results, (params=params, result=result, time=eval_time))
            
            if verbose
                println("  结果: $(round(result, digits=2))")
                println("  用时: $(round(eval_time, digits=2)) 秒")
            end
            
        catch e
            eval_time = time() - start_time
            total_time += eval_time
            
            push!(warmup_results, (params=params, result=Inf, time=eval_time, error=e))
            
            if verbose
                println("  结果: 计算失败 - $e")
                println("  用时: $(round(eval_time, digits=2)) 秒")
            end
        end
    end
    
    estimated_time_per_eval = total_time / n_warmup_samples
    
    if verbose
        println("\n预热完成!")
        println("平均评估时间: $(round(estimated_time_per_eval, digits=2)) 秒")
        
        # 统计成功率
        successful_evals = sum([r.result != Inf for r in warmup_results])
        success_rate = successful_evals / n_warmup_samples * 100
        println("成功率: $(round(success_rate, digits=1))% ($successful_evals/$n_warmup_samples)")
        
        if successful_evals > 0
            successful_results = [r.result for r in warmup_results if r.result != Inf]
            println("成功结果范围: $(round(minimum(successful_results), digits=2)) - $(round(maximum(successful_results), digits=2))")
        end
        
        # 时间估算
        println("\n时间估算 (基于预热结果):")
        for (desc, iters) in [("10次迭代", 10), ("25次迭代", 25), ("50次迭代", 50)]
            estimated_total = estimated_time_per_eval * iters
            if estimated_total < 60
                println("  $desc: ~$(round(estimated_total, digits=1)) 秒")
            elseif estimated_total < 3600
                println("  $desc: ~$(round(estimated_total/60, digits=1)) 分钟")
            else
                println("  $desc: ~$(round(estimated_total/3600, digits=1)) 小时")
            end
        end
    end
    
    return warmup_results, estimated_time_per_eval
end

function optimize_pnjl_parameters_with_warmup(
    kappa_pairs, μ_B_values, T_min, T_max, lower_bounds, upper_bounds;
    maxiterations=50, initial_samples=10, T_step_scan=2.0/hc,
    acquisition_function=:expected_improvement, model_optimization_frequency=5,
    penalty_for_missing=1e6, verbose_objective=false, verbosity_level=:progress,
    warmup_samples=3, skip_warmup=false,
    output_file="Rotation_PNJL/output/Gas_Liquid/bayesian_optimization_results.csv")
    """
    带预热的PNJL模型参数贝叶斯优化
    
    参数:
    - 所有原有参数...
    - warmup_samples: 预热样本数 (默认3)
    - skip_warmup: 是否跳过预热 (默认false)
    
    返回:
    - optimization_result: 优化结果结构体
    """
    
    println("="^100)
    println("PNJL模型参数贝叶斯优化 (带预热)")
    println("="^100)
    println("开始时间: $(Dates.now())")
    
    # 验证输入
    if length(μ_B_values) != length(kappa_pairs)
        error("μ_B值数组长度与κ比值对数组长度不匹配")
    end
    
    if length(lower_bounds) != 5 || length(upper_bounds) != 5
        error("参数边界必须是5维向量 [ρ₀, B_A, K, m_ratio, E_sym]")
    end
    
    # 显示优化配置
    param_names = ["ρ₀ (fm⁻³)", "B_A (MeV)", "K (MeV)", "m_ratio", "E_sym (MeV)"]
    println("\n优化配置:")
    println("  实验数据点: $(length(kappa_pairs)) 组")
    println("  参数边界:")
    for i in 1:5
        println("    $(param_names[i]): $(lower_bounds[i]) - $(upper_bounds[i])")
    end
    println("  最大迭代: $maxiterations")
    println("  初始采样: $initial_samples")
    println("  采集函数: $acquisition_function")
    println("  预热样本: $warmup_samples")
    
    # 创建目标函数
    println("\n创建目标函数...")
    objective_function = create_bayesian_objective_function(
        kappa_pairs, μ_B_values, T_min, T_max;
        T_step_scan=T_step_scan, penalty_for_missing=penalty_for_missing,
        verbose=verbose_objective)
    
    # 预热阶段
    warmup_time = 0.0
    estimated_time_per_eval = 0.0
    
    if !skip_warmup
        println("\n" * "="^60)
        println("第一阶段：目标函数预热")
        println("="^60)
        
        warmup_start = time()
        warmup_results, estimated_time_per_eval = warmup_objective_function(
            objective_function, lower_bounds, upper_bounds;
            n_warmup_samples=warmup_samples, verbose=true)
        warmup_time = time() - warmup_start
        
        println("\n预热阶段完成，用时: $(round(warmup_time, digits=2)) 秒")
    else
        println("\n跳过预热阶段")
    end
    
    # 设置优化器
    println("\n" * "="^60)
    println("第二阶段：设置贝叶斯优化器")
    println("="^60)
    optimizer = setup_bayesian_optimizer(
        objective_function, lower_bounds, upper_bounds;
        maxiterations=maxiterations, initial_samples=initial_samples,
        acquisition_function=acquisition_function,
        model_optimization_frequency=model_optimization_frequency,
        verbosity_level=verbosity_level)
    
    # 执行优化
    println("\n" * "="^60)
    println("第三阶段：贝叶斯优化")
    println("="^60)
    
    if estimated_time_per_eval > 0
        estimated_total_time = estimated_time_per_eval * maxiterations
        if estimated_total_time < 60
            println("预计优化时间: ~$(round(estimated_total_time, digits=1)) 秒")
        elseif estimated_total_time < 3600
            println("预计优化时间: ~$(round(estimated_total_time/60, digits=1)) 分钟")
        else
            println("预计优化时间: ~$(round(estimated_total_time/3600, digits=1)) 小时")
        end
    end
    
    optimization_start = time()
    
    try
        result = boptimize!(optimizer)
        
        optimization_time = time() - optimization_start
        total_time = warmup_time + optimization_time
        
        println("\n" * "="^80)
        println("贝叶斯优化完成!")
        println("="^80)
        println("预热时间: $(round(warmup_time, digits=2)) 秒")
        println("优化时间: $(round(optimization_time, digits=2)) 秒")
        println("总用时: $(round(total_time, digits=2)) 秒")
        
        # 提取最优结果
        best_params = result.observed_optimizer
        best_value = result.observed_optimum
        
        println("\n最优参数:")
        for i in 1:5
            println("  $(param_names[i]) = $(round(best_params[i], digits=6))")
        end
        println("\n最优目标函数值: $(round(best_value, digits=4)) MeV²")
        
        # 保存结果（包含预热信息）
        save_optimization_results_with_warmup(optimizer, result, kappa_pairs, μ_B_values, 
                                             T_min, T_max, lower_bounds, upper_bounds,
                                             maxiterations, initial_samples, T_step_scan,
                                             acquisition_function, output_file, 
                                             total_time, warmup_time, optimization_time,
                                             warmup_samples, estimated_time_per_eval)
        
        return result
        
    catch e
        println("\n❌ 贝叶斯优化失败: $e")
        println("错误位置: $(stacktrace()[1])")
        return nothing
    end
end

function save_optimization_results_with_warmup(
    optimizer, result, kappa_pairs, μ_B_values, T_min, T_max,
    lower_bounds, upper_bounds, maxiterations, initial_samples, T_step_scan,
    acquisition_function, output_file, total_time, warmup_time, optimization_time,
    warmup_samples, estimated_time_per_eval)
    """
    保存包含预热信息的贝叶斯优化结果
    """
    
    println("\n" * "="^60)
    println("保存优化结果")
    println("="^60)
    
    # 确保输出目录存在
    output_dir = dirname(output_file)
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    
    try
        # 提取优化历史
        n_evals = length(optimizer.model.y)
        param_names = ["rho0", "B_A", "K", "m_ratio", "E_sym"]
        
        # 创建结果DataFrame
        results_data = []
        
        for i in 1:n_evals
            x_i = optimizer.model.x[:, i]  # 第i个评估点的参数
            y_i = optimizer.model.y[i]     # 第i个评估点的目标函数值
            
            row = Dict(
                "iteration" => i,
                "objective_value" => y_i
            )
            
            # 添加参数值
            for (j, param_name) in enumerate(param_names)
                row[param_name] = x_i[j]
            end
            
            push!(results_data, row)
        end
        
        df = DataFrame(results_data)
        
        # 写入CSV文件（带预热元数据）
        open(output_file, "w") do io
            # 写入元数据
            println(io, "# Bayesian Optimization Results for PNJL Model Parameters (With Warmup)")
            println(io, "# Generated on: $(Dates.now())")
            println(io, "# Optimization Configuration:")
            println(io, "# maxiterations = $maxiterations")
            println(io, "# initial_samples = $initial_samples")
            println(io, "# acquisition_function = $acquisition_function")
            println(io, "# T_step_scan = $(T_step_scan*hc) MeV")
            println(io, "#")
            println(io, "# Timing Information:")
            println(io, "# warmup_samples = $warmup_samples")
            println(io, "# warmup_time = $(round(warmup_time, digits=2)) seconds")
            println(io, "# optimization_time = $(round(optimization_time, digits=2)) seconds")
            println(io, "# total_time = $(round(total_time, digits=2)) seconds")
            println(io, "# estimated_time_per_eval = $(round(estimated_time_per_eval, digits=2)) seconds")
            println(io, "#")
            println(io, "# Parameter Bounds:")
            param_names_full = ["ρ₀ (fm⁻³)", "B_A (MeV)", "K (MeV)", "m_ratio", "E_sym (MeV)"]
            for i in 1:5
                println(io, "# $(param_names_full[i]): $(lower_bounds[i]) - $(upper_bounds[i])")
            end
            println(io, "#")
            println(io, "# Experimental Data:")
            println(io, "# kappa_pairs = $kappa_pairs")
            println(io, "# mu_B_values = $([round(μ*hc, digits=1) for μ in μ_B_values]) MeV")
            println(io, "# T_range = $(T_min*hc) - $(T_max*hc) MeV")
            println(io, "#")
            println(io, "# Best Result:")
            println(io, "# best_objective = $(result.observed_optimum)")
            println(io, "# best_params = $(result.observed_optimizer)")
            println(io, "#")
            
            # 写入CSV数据
            CSV.write(io, df)
        end
        
        println("✅ 优化结果已保存到: $output_file")
        
        # 显示统计信息
        best_iteration = argmin(df.objective_value)
        worst_value = maximum(df.objective_value)
        
        println("\n优化统计:")
        println("  总评估次数: $n_evals")
        println("  最佳迭代: $best_iteration")
        println("  最优值: $(round(result.observed_optimum, digits=4))")
        println("  最差值: $(round(worst_value, digits=4))")
        println("  改进倍数: $(round(worst_value / result.observed_optimum, digits=2))x")
        println("  平均评估时间: $(round(estimated_time_per_eval, digits=2)) 秒")
        
    catch e
        println("❌ 保存结果失败: $e")
    end
end

function setup_bayesian_optimizer(
    objective_function, lower_bounds, upper_bounds;
    maxiterations=50, initial_samples=10, 
    acquisition_function=:expected_improvement,
    model_optimization_frequency=5,
    verbosity_level=:progress)
    """
    设置贝叶斯优化器
    
    参数:
    - objective_function: 目标函数（接受向量参数）
    - lower_bounds: 参数下界向量 [ρ₀_min, B_A_min, K_min, m_ratio_min, E_sym_min]
    - upper_bounds: 参数上界向量 [ρ₀_max, B_A_max, K_max, m_ratio_max, E_sym_max]
    - maxiterations: 最大迭代次数
    - initial_samples: 初始采样点数
    - acquisition_function: 采集函数类型
    - model_optimization_frequency: 模型超参数优化频率
    - verbosity_level: 输出详细程度
    
    返回:
    - BOpt: 配置好的贝叶斯优化器
    """
    
    # 参数维度
    n_dims = length(lower_bounds)
    
    println("\n贝叶斯优化器配置:")
    println("  参数维度: $n_dims")
    println("  参数下界: $lower_bounds")
    println("  参数上界: $upper_bounds")
    println("  最大迭代: $maxiterations")
    println("  初始采样: $initial_samples")
    
    # 验证边界
    for i in 1:n_dims
        if lower_bounds[i] >= upper_bounds[i]
            error("参数 $i 的下界 ($(lower_bounds[i])) 必须小于上界 ($(upper_bounds[i]))")
        end
    end
    
    # 创建高斯过程模型
    println("\n创建高斯过程模型...")
    model = ElasticGPE(
        n_dims,                           # 输入维度
        mean = MeanConst(0.0),           # 常数均值函数
        kernel = SEArd(ones(n_dims), 0.0), # 自动相关确定的平方指数核
        logNoise = 0.,                 # 对数噪声
        capacity = 3000                  # 容量
    )
    
    # 选择采集函数
    if acquisition_function == :expected_improvement
        acquisition = ExpectedImprovement()
    elseif acquisition_function == :upper_confidence_bound
        acquisition = UpperConfidenceBound()
    elseif acquisition_function == :probability_of_improvement
        acquisition = ProbabilityOfImprovement()
    else
        acquisition = ExpectedImprovement()  # 默认
    end
    
    # 设置模型优化器
    println("设置模型超参数优化器...")
    
    # 为kernbounds设置合理的边界
    # 对于SEArd核：需要 n_dims + 1 个参数（length scales + signal variance）
    kern_lower = [-2 * ones(n_dims); -2]  # [长度尺度下界..., 信号方差下界] - 缩小范围
    kern_upper = [2 * ones(n_dims); 1]    # [长度尺度上界..., 信号方差上界] - 缩小范围
    
    # 使用更保守的模型优化设置，减少优化失败
    modeloptimizer = MAPGPOptimizer(
        every = model_optimization_frequency * 2, # 减少超参数优化频率
        noisebounds = [-3, 2],                   # 缩小噪声边界范围
        kernbounds = [kern_lower, kern_upper],    # 使用更保守的核参数边界
        maxeval = 50                             # 减少超参数优化的评估次数
    )
    
    # 设置输出详细程度
    if verbosity_level == :silent
        verbosity = Silent
    elseif verbosity_level == :timings
        verbosity = Timings
    else
        verbosity = Progress  # 默认
    end
    
    # 设置采集函数优化选项，避免FORCED_STOP错误
    acquisition_options = (
        method = :GD_STOGO_RAND,     # 使用L-BFGS优化方法
        restarts = 20,           # 减少重启次数
        maxeval = 10000,          # 减少最大评估次数
        maxtime = 30.0         # 设置时间限制（30秒）
        #ftol_rel = 1e-3,        # 设置相对容差
        #ftol_abs = 1e-4         # 设置绝对容差
    )
    
    # 创建贝叶斯优化器
    println("创建贝叶斯优化器...")
    optimizer = BOpt(
        objective_function,              # 目标函数
        model,                          # 高斯过程模型
        acquisition,                    # 采集函数
        modeloptimizer,                 # 模型优化器
        lower_bounds,                   # 下界
        upper_bounds,                   # 上界
        sense = Min,                    # 最小化目标函数
        maxiterations = maxiterations,   # 最大迭代次数
        verbosity = verbosity,          # 输出详细程度
        initializer_iterations = initial_samples,  # 初始采样点数
        acquisitionoptions = acquisition_options   # 采集函数优化选项
    )
    
    println("✅ 贝叶斯优化器创建完成")
    
    return optimizer
end

function optimize_pnjl_parameters(
    kappa_pairs, μ_B_values, T_min, T_max, lower_bounds, upper_bounds;
    maxiterations=50, initial_samples=10, T_step_scan=2.0/hc,
    acquisition_function=:expected_improvement, model_optimization_frequency=5,
    penalty_for_missing=1e6, verbose_objective=false, verbosity_level=:progress,
    output_file="Rotation_PNJL/output/Gas_Liquid/bayesian_optimization_results.csv")
    """
    使用贝叶斯优化来优化PNJL模型参数
    
    参数:
    - kappa_pairs: κ比值对数组
    - μ_B_values: 重子化学势数组
    - T_min, T_max: 温度搜索范围
    - lower_bounds: 参数下界 [ρ₀_min, B_A_min, K_min, m_ratio_min, E_sym_min]
    - upper_bounds: 参数上界 [ρ₀_max, B_A_max, K_max, m_ratio_max, E_sym_max]
    - maxiterations: 最大迭代次数
    - initial_samples: 初始采样点数
    - T_step_scan: 温度扫描步长
    - acquisition_function: 采集函数类型
    - model_optimization_frequency: 模型超参数优化频率
    - penalty_for_missing: 计算失败惩罚值
    - verbose_objective: 目标函数是否显示详细信息
    - verbosity_level: 优化过程输出详细程度
    - output_file: 结果保存文件
    
    返回:
    - optimization_result: 优化结果结构体
    """
    
    println("="^100)
    println("PNJL模型参数贝叶斯优化")
    println("="^100)
    println("开始时间: $(Dates.now())")
    
    # 验证输入
    if length(μ_B_values) != length(kappa_pairs)
        error("μ_B值数组长度与κ比值对数组长度不匹配")
    end
    
    if length(lower_bounds) != 5 || length(upper_bounds) != 5
        error("参数边界必须是5维向量 [ρ₀, B_A, K, m_ratio, E_sym]")
    end
    
    # 显示优化配置
    param_names = ["ρ₀ (fm⁻³)", "B_A (MeV)", "K (MeV)", "m_ratio", "E_sym (MeV)"]
    println("\n优化配置:")
    println("  实验数据点: $(length(kappa_pairs)) 组")
    println("  参数边界:")
    for i in 1:5
        println("    $(param_names[i]): $(lower_bounds[i]) - $(upper_bounds[i])")
    end
    println("  最大迭代: $maxiterations")
    println("  初始采样: $initial_samples")
    println("  采集函数: $acquisition_function")
    
    # 创建目标函数
    println("\n创建目标函数...")
    objective_function = create_bayesian_objective_function(
        kappa_pairs, μ_B_values, T_min, T_max;
        T_step_scan=T_step_scan, penalty_for_missing=penalty_for_missing,
        verbose=verbose_objective)
    
    # 设置优化器
    println("\n设置贝叶斯优化器...")
    optimizer = setup_bayesian_optimizer(
        objective_function, lower_bounds, upper_bounds;
        maxiterations=maxiterations, initial_samples=initial_samples,
        acquisition_function=acquisition_function,
        model_optimization_frequency=model_optimization_frequency,
        verbosity_level=verbosity_level)
    
    # 执行优化
    println("\n" * "="^80)
    println("开始贝叶斯优化...")
    println("="^80)
    
    start_time = time()
    
    try
        result = boptimize!(optimizer)
        
        end_time = time()
        elapsed_time = end_time - start_time
        
        println("\n" * "="^80)
        println("贝叶斯优化完成!")
        println("="^80)
        println("优化用时: $(round(elapsed_time, digits=2)) 秒")
        
        # 提取最优结果
        best_params = result.observed_optimizer
        best_value = result.observed_optimum
        
        println("\n最优参数:")
        for i in 1:5
            println("  $(param_names[i]) = $(round(best_params[i], digits=6))")
        end
        println("\n最优目标函数值: $(round(best_value, digits=4)) MeV²")
        
        # 保存结果
        save_optimization_results(optimizer, result, kappa_pairs, μ_B_values, 
                                 T_min, T_max, lower_bounds, upper_bounds,
                                 maxiterations, initial_samples, T_step_scan,
                                 acquisition_function, output_file, elapsed_time)
        
        return result
        
    catch e
        println("\n❌ 贝叶斯优化失败: $e")
        println("错误位置: $(stacktrace()[1])")
        return nothing
    end
end

function save_optimization_results(
    optimizer, result, kappa_pairs, μ_B_values, T_min, T_max,
    lower_bounds, upper_bounds, maxiterations, initial_samples, T_step_scan,
    acquisition_function, output_file, elapsed_time)
    """
    保存贝叶斯优化结果到CSV文件
    """
    
    println("\n" * "="^60)
    println("保存优化结果")
    println("="^60)
    
    # 确保输出目录存在
    output_dir = dirname(output_file)
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    
    try
        # 提取优化历史
        n_evals = length(optimizer.model.y)
        param_names = ["rho0", "B_A", "K", "m_ratio", "E_sym"]
        
        # 创建结果DataFrame
        results_data = []
        
        for i in 1:n_evals
            x_i = optimizer.model.x[:, i]  # 第i个评估点的参数
            y_i = optimizer.model.y[i]     # 第i个评估点的目标函数值
            
            row = Dict(
                "iteration" => i,
                "objective_value" => y_i
            )
            
            # 添加参数值
            for (j, param_name) in enumerate(param_names)
                row[param_name] = x_i[j]
            end
            
            push!(results_data, row)
        end
        
        df = DataFrame(results_data)
        
        # 写入CSV文件（带元数据头部）
        open(output_file, "w") do io
            # 写入元数据
            println(io, "# Bayesian Optimization Results for PNJL Model Parameters")
            println(io, "# Generated on: $(Dates.now())")
            println(io, "# Optimization Configuration:")
            println(io, "# maxiterations = $maxiterations")
            println(io, "# initial_samples = $initial_samples")
            println(io, "# acquisition_function = $acquisition_function")
            println(io, "# T_step_scan = $(T_step_scan*hc) MeV")
            println(io, "# elapsed_time = $(round(elapsed_time, digits=2)) seconds")
            println(io, "#")
            println(io, "# Parameter Bounds:")
            param_names_full = ["ρ₀ (fm⁻³)", "B_A (MeV)", "K (MeV)", "m_ratio", "E_sym (MeV)"]
            for i in 1:5
                println(io, "# $(param_names_full[i]): $(lower_bounds[i]) - $(upper_bounds[i])")
            end
            println(io, "#")
            println(io, "# Experimental Data:")
            println(io, "# kappa_pairs = $kappa_pairs")
            println(io, "# mu_B_values = $([round(μ*hc, digits=1) for μ in μ_B_values]) MeV")
            println(io, "# T_range = $(T_min*hc) - $(T_max*hc) MeV")
            println(io, "#")
            println(io, "# Best Result:")
            println(io, "# best_objective = $(result.observed_optimum)")
            println(io, "# best_params = $(result.observed_optimizer)")
            println(io, "#")
            
            # 写入CSV数据
            CSV.write(io, df)
        end
        
        println("✅ 优化结果已保存到: $output_file")
        
        # 显示统计信息
        best_iteration = argmin(df.objective_value)
        worst_value = maximum(df.objective_value)
        
        println("\n优化统计:")
        println("  总评估次数: $n_evals")
        println("  最佳迭代: $best_iteration")
        println("  最优值: $(round(result.observed_optimum, digits=4))")
        println("  最差值: $(round(worst_value, digits=4))")
        println("  改进倍数: $(round(worst_value / result.observed_optimum, digits=2))x")
        
    catch e
        println("❌ 保存结果失败: $e")
    end
end

function demo_bayesian_optimization()
    """
    演示贝叶斯优化功能
    使用测试文件中的实验数据
    """
    
    println("="^100)
    println("演示：PNJL模型参数贝叶斯优化")
    println("="^100)
    
    # 使用测试文件中的实验数据
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
    
    T_min, T_max = 25.0/hc, 200.0/hc
    
    # 设置参数边界（请根据物理约束调整）
    lower_bounds = [0.145, -17.0, 212.0, 0.55, 26.1]   # [ρ₀, B_A, K, m_ratio, E_sym]
    upper_bounds = [0.170, -15.6, 401.0, 0.75, 44.0]   # [ρ₀, B_A, K, m_ratio, E_sym]
    
    println("演示配置:")
    println("  实验数据: $(length(kappa_pairs)) 组κ比值对")
    println("  μ_B值: $([round(μ*hc, digits=1) for μ in μ_B_values]) MeV")
    println("  温度范围: $(T_min*hc) - $(T_max*hc) MeV")
    println("  参数边界: $lower_bounds - $upper_bounds")
    
    # 执行贝叶斯优化（演示版本使用较少迭代）
    result = optimize_pnjl_parameters(
        kappa_pairs, μ_B_values, T_min, T_max, lower_bounds, upper_bounds;
        maxiterations=25,              # 演示用较少迭代
        initial_samples=8,             # 演示用较少初始样本
        T_step_scan=3.0/hc,           # 演示用粗扫描
        acquisition_function=:expected_improvement,
        model_optimization_frequency=5,
        penalty_for_missing=1e5,
        verbose_objective=false,
        verbosity_level=:timings,      # 使用timings减少输出
        output_file="Rotation_PNJL/output/Gas_Liquid/demo_bayesian_optimization.csv"
    )
    
    if result !== nothing
        println("\n" * "="^80)
        println("✅ 演示成功完成!")
        println("="^80)
        println("请检查输出文件了解详细结果。")
        println("您可以调整参数边界和迭代次数来获得更好的优化结果。")
    else
        println("\n❌ 演示失败。")
    end
    
    return result
end

function demo_bayesian_optimization_with_warmup()
    """
    演示带预热的贝叶斯优化功能
    使用测试文件中的实验数据
    """
    
    println("="^100)
    println("演示：PNJL模型参数贝叶斯优化 (带预热)")
    println("="^100)
    
    # 使用测试文件中的实验数据
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
    
    T_min, T_max = 70.0/hc, 120.0/hc
    
    # 设置参数边界（请根据物理约束调整）
    lower_bounds = [0.145, -17.0, 212.0, 0.55, 26.1]   # [ρ₀, B_A, K, m_ratio, E_sym]
    upper_bounds = [0.170, -15.6, 401.0, 0.75, 44.0]   # [ρ₀, B_A, K, m_ratio, E_sym]
    
    println("演示配置:")
    println("  实验数据: $(length(kappa_pairs)) 组κ比值对")
    println("  μ_B值: $([round(μ*hc, digits=1) for μ in μ_B_values]) MeV")
    println("  温度范围: $(T_min*hc) - $(T_max*hc) MeV")
    println("  参数边界: $lower_bounds - $upper_bounds")
    
    # 执行带预热的贝叶斯优化
    result = optimize_pnjl_parameters_with_warmup(
        kappa_pairs, μ_B_values, T_min, T_max, lower_bounds, upper_bounds;
        maxiterations=20,              # 减少迭代次数
        initial_samples=10,             # 减少初始样本
        T_step_scan=2.0/hc,           # 使用更粗的扫描步长
        acquisition_function=:expected_improvement,
        model_optimization_frequency=3, # 减少模型优化频率
        penalty_for_missing=1e3,       # 降低惩罚值
        verbose_objective=false,
        verbosity_level=:Progress,      # 使用Progress显示进度
        warmup_samples=1,              # 预热样本数
        skip_warmup=true,             # 跳过预热
        output_file="Rotation_PNJL/output/Gas_Liquid/demo_bayesian_optimization_warmup.csv"
    )
    
    if result !== nothing
        println("\n" * "="^80)
        println("✅ 带预热的演示成功完成!")
        println("="^80)
        println("预热功能的优势:")
        println("- Julia函数编译优化（JIT预热）")
        println("- 目标函数计算时间估算")
        println("- 优化时间预测")
        println("- 潜在问题早期发现")
        println("\n请检查输出文件了解详细结果。")
    else
        println("\n❌ 带预热的演示失败。")
    end
    
    return result
end

function demo_simple_bayesian_optimization_with_warmup()
    """
    简化版演示：带预热的快速贝叶斯优化测试
    """
    
    println("="^80)
    println("简化演示：PNJL模型参数贝叶斯优化（带预热的快速测试版）")
    println("="^80)
    
    # 使用更少的实验数据进行快速测试
    kappa_pairs = [
        (1.09031788496341, -0.28904867673079),   # 第1组
        (1.06152332992368, 0.164279260625683),   # 第2组  
        (1.11111023684003, 0.224522832511389)    # 第3组
    ]
    
    μ_B_values = [
        632.0 / hc,       # 第1组对应632 MeV
        666.0 / hc        # 第2组对应666 MeV
    ]
    
    T_min, T_max = 50.0/hc, 120.0/hc  # 缩小温度范围
    
    # 设置较宽松的参数边界用于测试
    lower_bounds = [0.13, -17.5, 220.0, 0.65, 30.0]   # [ρ₀, B_A, K, m_ratio, E_sym]
    upper_bounds = [0.17, -15.5, 260.0, 0.75, 34.0]   # [ρ₀, B_A, K, m_ratio, E_sym]
    
    println("简化演示配置:")
    println("  实验数据: $(length(kappa_pairs)) 组κ比值对")
    println("  μ_B值: $([round(μ*hc, digits=1) for μ in μ_B_values]) MeV")
    println("  温度范围: $(T_min*hc) - $(T_max*hc) MeV（缩小范围）")
    println("  参数边界: $lower_bounds - $upper_bounds")
    
    # 执行简化的带预热贝叶斯优化
    result = optimize_pnjl_parameters_with_warmup(
        kappa_pairs, μ_B_values, T_min, T_max, lower_bounds, upper_bounds;
        maxiterations=6,               # 非常少的迭代
        initial_samples=3,             # 非常少的初始样本
        T_step_scan=8.0/hc,           # 粗扫描
        acquisition_function=:expected_improvement,
        model_optimization_frequency=10, # 很少优化模型
        penalty_for_missing=5e3,       # 较低惩罚值
        verbose_objective=false,
        verbosity_level=:timings,      # 简化输出
        warmup_samples=2,              # 很少的预热样本
        skip_warmup=false,
        output_file="Rotation_PNJL/output/Gas_Liquid/demo_simple_bayesian_warmup.csv"
    )
    
    if result !== nothing
        println("\n" * "="^60)
        println("✅ 简化的带预热演示成功完成!")
        println("="^60)
        println("这是一个超快速测试版本，使用了:")
        println("- 较少的实验数据点 ($(length(kappa_pairs))组)")
        println("- 较小的温度范围 ($(T_min*hc)-$(T_max*hc) MeV)")
        println("- 较粗的扫描步长 (8 MeV)")
        println("- 较少的优化迭代 (6次)")
        println("- 预热功能 (2个样本)")
        println("\n如需完整优化，请使用 demo_bayesian_optimization_with_warmup() 函数")
    else
        println("\n❌ 简化的带预热演示失败。")
    end
    
    return result
end

function quick_test_objective()
    """
    快速测试目标函数是否正常工作
    """
    
    println("="^60)
    println("快速测试目标函数")
    println("="^60)
    
    # 测试数据
    kappa_pairs = [
        (1.09031788496341, -0.28904867673079),   # 第1组
        (1.06152332992368, 0.164279260625683),   # 第2组  
        (1.11111023684003, 0.224522832511389)    # 第3组
    ]
    μ_B_values = [632.0/hc, 666.0/hc]
    T_min, T_max = 50.0/hc, 150.0/hc
    
    # 创建目标函数
    objective_func = create_bayesian_objective_function(
        kappa_pairs, μ_B_values, T_min, T_max;
        T_step_scan=5.0/hc, penalty_for_missing=1e4, verbose=false)
    
    # 测试几个参数点
    test_params = [
        [0.15, -16.0, 240.0, 0.7, 32.0],   # 测试参数1
        [0.12, -18.0, 220.0, 0.6, 30.0],   # 测试参数2
        [0.18, -14.0, 260.0, 0.8, 35.0]    # 测试参数3
    ]
    
    println("测试目标函数...")
    for (i, params) in enumerate(test_params)
        try
            result = objective_func(params)
            println("  测试 $i: 参数 $params → 目标值 $(round(result, digits=2))")
        catch e
            println("  测试 $i: 失败 - $e")
        end
    end
    
    println("✅ 目标函数测试完成")
end

function test_warmup_only()
    """
    仅测试预热功能，不进行完整优化
    """
    
    println("="^80)
    println("预热功能单独测试")
    println("="^80)
    
    # 测试数据
    kappa_pairs = [
        (1.09031788496341, -0.28904867673079),   # 第1组
        (1.06152332992368, 0.164279260625683),   # 第2组  
        (1.11111023684003, 0.224522832511389)    # 第3组
    ]  # 只用一组数据
    μ_B_values = [632.0/hc]
    T_min, T_max = 60.0/hc, 100.0/hc  # 很小的温度范围
    
    # 创建目标函数
    println("创建目标函数...")
    objective_func = create_bayesian_objective_function(
        kappa_pairs, μ_B_values, T_min, T_max;
        T_step_scan=10.0/hc, penalty_for_missing=1e4, verbose=false)
    
    # 设置测试边界
    lower_bounds = [0.14, -17.0, 230.0, 0.65, 30.0]
    upper_bounds = [0.16, -15.0, 250.0, 0.75, 34.0]
    
    # 执行预热测试
    println("\n开始预热测试...")
    warmup_results, estimated_time = warmup_objective_function(
        objective_func, lower_bounds, upper_bounds;
        n_warmup_samples=3, verbose=true)
    
    println("\n预热测试完成!")
    println("现在您可以选择是否继续进行完整的贝叶斯优化。")
    
    return warmup_results, estimated_time
end

function quick_warmup_optimization()
    """
    快速预热+最小化优化测试
    """
    
    println("="^80)
    println("快速预热+最小化优化测试")
    println("="^80)
    
    # 使用最简配置
    kappa_pairs = [
        (1.09031788496341, -0.28904867673079),   # 第1组
        (1.06152332992368, 0.164279260625683),   # 第2组  
        (1.11111023684003, 0.224522832511389)    # 第3组
    ]
    μ_B_values = [632.0/hc]
    T_min, T_max = 80.0/hc, 120.0/hc  # 很窄的温度范围
    lower_bounds = [0.14, -17.0, 235.0, 0.68, 31.0]
    upper_bounds = [0.16, -15.0, 245.0, 0.72, 33.0]  # 很窄的参数范围
    
    println("超快速测试配置:")
    println("  数据: 1组κ比值")
    println("  温度: $(T_min*hc)-$(T_max*hc) MeV")
    println("  边界: 很窄的参数空间")
    
    # 执行超快速优化
    result = optimize_pnjl_parameters_with_warmup(
        kappa_pairs, μ_B_values, T_min, T_max, lower_bounds, upper_bounds;
        maxiterations=3,               # 超少迭代
        initial_samples=2,             # 超少样本
        T_step_scan=10.0/hc,          # 超粗扫描
        acquisition_function=:expected_improvement,
        model_optimization_frequency=20, # 基本不优化模型
        penalty_for_missing=1e3,
        verbose_objective=false,
        verbosity_level=:timings,
        warmup_samples=2,              # 超少预热
        skip_warmup=false,
        output_file="Rotation_PNJL/output/Gas_Liquid/quick_warmup_test.csv"
    )
    
    if result !== nothing
        println("\n✅ 快速预热+优化测试成功!")
        println("这验证了整个流程可以正常工作。")
    else
        println("\n❌ 快速测试失败。")
    end
    
    return result
end

function load_previous_optimization_results(csv_file::String)
    """
    从CSV文件加载前一次优化的结果
    
    参数:
    - csv_file: CSV结果文件路径
    
    返回:
    - (X_observed, y_observed): 观测点和对应的目标函数值
    - best_params: 最优参数
    - best_value: 最优目标函数值
    """
    
    if !isfile(csv_file)
        println("⚠️  未找到文件 $(csv_file)，将从头开始优化")
        return nothing, nothing, nothing, nothing
    end
    
    try
        df = CSV.read(csv_file, DataFrame)
        
        # 提取参数列（按顺序：ρ₀, B_A, K, m_ratio, E_sym）
        X_observed = Matrix{Float64}(df[:, [:rho0, :B_A, :K, :m_ratio, :E_sym]])
        y_observed = Vector{Float64}(abs.(df.objective_value))  # 取绝对值确保正值
        
        # 找到最优值
        best_idx = argmin(y_observed)
        best_params = X_observed[best_idx, :]
        best_value = y_observed[best_idx]
        
        println("✅ 成功加载 $(size(X_observed, 1)) 个历史优化点")
        println("   最优参数: [$(join([round(p, digits=4) for p in best_params], ", "))]")
        println("   最优值: $(round(best_value, digits=2))")
        
        return X_observed, y_observed, best_params, best_value
        
    catch e
        println("❌ 加载CSV文件失败: $e")
        return nothing, nothing, nothing, nothing
    end
end

function setup_bayesian_optimizer_with_prior_data(
    objective_function, lower_bounds, upper_bounds, X_prior=nothing, y_prior=nothing;
    maxiterations=50, initial_samples=10, 
    acquisition_function=:expected_improvement,
    model_optimization_frequency=5,
    verbosity_level=:progress)
    """
    设置包含先验数据的贝叶斯优化器
    
    参数:
    - X_prior: 先验观测点 (n_points × n_dims)
    - y_prior: 先验目标函数值 (n_points,)
    - 其他参数与原函数相同
    """
    
    n_dims = length(lower_bounds)
    
    println("\n设置带先验数据的贝叶斯优化器:")
    println("  参数维度: $n_dims")
    println("  先验数据点: $(X_prior === nothing ? 0 : size(X_prior, 1))")
    println("  最大迭代: $maxiterations")
    println("  初始采样: $initial_samples")
    
    # 验证边界
    for i in 1:n_dims
        if upper_bounds[i] <= lower_bounds[i]
            error("参数 $i 的上界必须大于下界: $(lower_bounds[i]) < $(upper_bounds[i])")
        end
    end
    
    # 创建高斯过程模型
    println("\n创建高斯过程模型...")
    model = ElasticGPE(
        n_dims,
        mean = MeanConst(0.),
        kernel = SEArd(ones(n_dims), 0.0),  # 与原始函数保持一致
        logNoise = log(0.01),
        capacity = 3000
    )
    
    # 选择采集函数
    if acquisition_function == :expected_improvement
        acquisition = ExpectedImprovement()
    elseif acquisition_function == :upper_confidence_bound
        acquisition = UpperConfidenceBound()
    elseif acquisition_function == :probability_of_improvement
        acquisition = ProbabilityOfImprovement()
    else
        acquisition = ExpectedImprovement()
    end
    
    # 设置模型优化器
    kern_lower = [-2 * ones(n_dims); -2]
    kern_upper = [2 * ones(n_dims); 1]
    
    modeloptimizer = MAPGPOptimizer(
        every = model_optimization_frequency,
        noisebounds = [-4, 1],
        kernbounds = (kern_lower, kern_upper),
        maxeval = 30
    )
    
    # 设置输出详细程度
    if verbosity_level == :silent
        verbosity = Silent
    elseif verbosity_level == :timings
        verbosity = Timings
    else
        verbosity = Progress
    end
    
    # 设置采集函数优化选项
    acquisition_options = (
        method = :GN_DIRECT_L,
        restarts = 10,
        maxtime = 30.0,
        maxeval = 10000
        #xtol_abs = 1e-7,
        #ftol_abs = 1e-9
    )
    
    # 创建贝叶斯优化器
    println("创建贝叶斯优化器...")
    optimizer = BOpt(
        objective_function,
        model,
        acquisition,
        modeloptimizer,
        lower_bounds, upper_bounds,
        repetitions = 1,
        maxiterations = maxiterations,
        sense = Min,
        verbosity = verbosity,
        acquisitionoptions = acquisition_options,
        initializer_iterations = initial_samples
    )
    
    # 如果有先验数据，添加到优化器中
    if X_prior !== nothing && y_prior !== nothing
        println("添加先验数据到优化器...")
        
        # 验证数据维度
        if size(X_prior, 2) != n_dims
            error("先验数据维度 $(size(X_prior, 2)) 与参数维度 $n_dims 不匹配")
        end
        
        # 过滤并准备有效的先验数据
        valid_X = Float64[]
        valid_y = Float64[]
        valid_count = 0
        
        for i in 1:size(X_prior, 1)
            x_point = X_prior[i, :]
            y_point = y_prior[i]
            
            # 验证点在边界内且有效
            if all(lower_bounds .<= x_point .<= upper_bounds) && isfinite(y_point) && y_point > 0
                if valid_count == 0
                    valid_X = reshape(x_point, 1, :)
                    valid_y = [y_point]
                else
                    valid_X = vcat(valid_X, reshape(x_point, 1, :))
                    push!(valid_y, y_point)
                end
                valid_count += 1
            end
        end
        
        if valid_count > 0
            println("✅ 筛选出 $valid_count 个有效的先验数据点")
            
            # 通过直接初始化模型来使用先验数据
            try
                # 准备数据格式用于高斯过程初始化
                GP_X = valid_X'  # 转置为 n_dims × n_points 格式
                GP_y = valid_y
                
                # 初始化高斯过程模型（这会替换空模型）
                initialise_model!(optimizer.model; GP_X=GP_X, GP_y=GP_y)
                
                println("✅ 成功使用 $valid_count 个先验观测点初始化模型")
            catch e
                println("⚠️  添加先验数据失败: $e")
                println("将不使用先验数据，从头开始优化")
            end
        else
            println("⚠️  没有找到有效的先验数据点")
        end
    end
    
    println("✅ 带先验数据的贝叶斯优化器创建完成")
    
    return optimizer
end

function continue_bayesian_optimization_from_csv(
    csv_file::String, kappa_pairs, μ_B_values, T_min, T_max, lower_bounds, upper_bounds;
    additional_iterations=25, T_step_scan=1.0/hc,
    acquisition_function=:expected_improvement, model_optimization_frequency=5,
    penalty_for_missing=1e6, verbose_objective=false, verbosity_level=:progress,
    output_file="Rotation_PNJL/output/Gas_Liquid/continued_bayesian_optimization.csv")
    """
    从CSV文件继续贝叶斯优化
    
    参数:
    - csv_file: 前一次优化结果的CSV文件
    - additional_iterations: 额外的迭代次数
    - 其他参数与原优化函数相同
    """
    
    println("="^100)
    println("从CSV文件继续贝叶斯优化")
    println("="^100)
    println("开始时间: $(Dates.now())")
    println("前次结果文件: $csv_file")
    println("额外迭代: $additional_iterations")
    
    # 加载先验数据
    X_prior, y_prior, best_params, best_value = load_previous_optimization_results(csv_file)
    
    if X_prior === nothing
        println("❌ 无法加载先验数据，改为从头开始优化")
        return optimize_pnjl_parameters_with_warmup(
            kappa_pairs, μ_B_values, T_min, T_max, lower_bounds, upper_bounds;
            maxiterations=additional_iterations, output_file=output_file)
    end
    
    # 创建目标函数
    println("\n创建目标函数...")
    objective_function = create_bayesian_objective_function(
        kappa_pairs, μ_B_values, T_min, T_max;
        T_step_scan=T_step_scan, penalty_for_missing=penalty_for_missing, verbose=verbose_objective)
    
    # 设置带先验数据的优化器
    println("\n设置带先验数据的优化器...")
    optimizer = setup_bayesian_optimizer_with_prior_data(
        objective_function, lower_bounds, upper_bounds, X_prior, y_prior;
        maxiterations=additional_iterations, initial_samples=5,  # 减少初始采样，因为有先验数据
        acquisition_function=acquisition_function,
        model_optimization_frequency=model_optimization_frequency,
        verbosity_level=verbosity_level)
    
    # 执行优化
    println("\n" * "="^80)
    println("开始继续优化...")
    println("="^80)
    println("基于 $(size(X_prior, 1)) 个历史数据点")
    println("当前最优值: $(round(best_value, digits=2))")
    
    start_time = time()
    
    try
        boptimize!(optimizer)
        elapsed_time = time() - start_time
        
        # 获取最终结果
        final_best_params = optimizer.observed_optimizer
        final_best_value = optimizer.observed_optimum
        
        println("\n" * "="^80)
        println("继续优化完成!")
        println("="^80)
        println("用时: $(round(elapsed_time, digits=2)) 秒")
        println("总评估次数: $(length(optimizer.observed_optimum_trace))")
        
        # 比较改进
        improvement = best_value - final_best_value
        improvement_percent = (improvement / best_value) * 100
        
        println("\n结果比较:")
        println("  初始最优值: $(round(best_value, digits=4))")
        println("  最终最优值: $(round(final_best_value, digits=4))")
        if improvement > 0
            println("  改进: $(round(improvement, digits=4)) ($(round(improvement_percent, digits=2))%)")
        else
            println("  未发现更好的解")
        end
        
        println("\n最优参数:")
        param_names = ["ρ₀ (fm⁻³)", "B_A (MeV)", "K (MeV)", "m_ratio", "E_sym (MeV)"]
        for (i, name) in enumerate(param_names)
            println("  $name = $(round(final_best_params[i], digits=6))")
        end
        
        # 保存结果
        result = (optimizer=optimizer, best_params=final_best_params, 
                 best_value=final_best_value, elapsed_time=elapsed_time)
        
        save_optimization_results(
            optimizer, result, kappa_pairs, μ_B_values, T_min, T_max,
            lower_bounds, upper_bounds, additional_iterations, 5, T_step_scan,
            acquisition_function, output_file, elapsed_time)
        
        println("\n✅ 继续优化成功完成!")
        println("请检查输出文件: $output_file")
        
        return result
        
    catch e
        println("❌ 继续优化失败: $e")
        return nothing
    end
end

function demo_continue_optimization()
    """
    演示如何从CSV文件继续优化
    """
    
    println("="^100)
    println("演示：从CSV文件继续贝叶斯优化")
    println("="^100)
    
    # 检查是否存在之前的结果
    previous_csv = "Rotation_PNJL/output/Gas_Liquid/demo_bayesian_optimization_warmup.csv"
    
    # 实验数据配置
    kappa_pairs = [
        (1.09031788496341, -0.28904867673079),   # 第1组
        (1.06152332992368, 0.164279260625683),   # 第2组  
        (1.11111023684003, 0.224522832511389)    # 第3组
    ]
    
    μ_B_values = [
        632.0/hc,
        666.0/hc,
        697.0/hc
    ]
    
    T_min, T_max = 25.0/hc, 200.0/hc
    
    # 参数边界 - 可以根据第一次的结果调整边界
    lower_bounds = [0.145, -17.0, 212.0, 0.55, 26.1]
    upper_bounds = [0.170, -15.6, 401.0, 0.75, 44.0]
    
    println("继续优化配置:")
    println("  前次结果文件: $(previous_csv)")
    println("  额外迭代: 25 次")
    println("  采集函数: upper_confidence_bound (更好的探索)")
    
    # 从CSV继续优化
    result = continue_bayesian_optimization_from_csv(
        previous_csv, kappa_pairs, μ_B_values, T_min, T_max, lower_bounds, upper_bounds;
        additional_iterations=25,
        acquisition_function=:upper_confidence_bound,  # 使用UCB进行更好的探索
        model_optimization_frequency=3,
        verbosity_level=:timings,
        output_file="Rotation_PNJL/output/Gas_Liquid/continued_optimization_demo.csv"
    )
    
    if result !== nothing
        println("\n✅ 继续优化演示成功!")
        println("建议:")
        println("  1. 比较新旧结果文件，观察改进情况")
        println("  2. 可以多次运行以继续改进")
        println("  3. 调整参数边界以探索更广的空间")
    else
        println("\n❌ 继续优化演示失败")
    end
    
    return result
end

# =============================================================================
# 简化版贝叶斯优化接口 - 仿照用户提供的简洁风格
# =============================================================================

"""
简化版贝叶斯优化 - 带先验数据支持

仿照 BayesianOptimization.jl 的简洁风格，专为 PNJL 模型参数优化设计。
支持从先验数据开始优化，减少不必要的复杂性。
"""

function create_pnjl_objective(kappa_pairs, μ_B_values, T_min, T_max; 
                              T_step=1.0/hc, penalty=1e6)
    """
    创建PNJL模型参数优化目标函数
    
    参数:
    - kappa_pairs: κ比值对 [(κ₃/κ₁, κ₄/κ₂), ...]
    - μ_B_values: 重子化学势数组 [μ_B1, μ_B2, ...]  
    - T_min, T_max: 温度搜索范围
    - T_step: 温度扫描步长
    - penalty: 计算失败惩罚值
    
    返回:
    - objective: 目标函数 f(x) where x = [ρ₀, B_A, K, m_ratio, E_sym]
    """
    
    # 创建基础闭包函数
    base_func = create_temperature_difference_objective(
        kappa_pairs, μ_B_values, T_min, T_max;
        T_step_scan=T_step, penalty_for_missing=penalty, verbose=false)
    
    # 包装为向量输入格式
    function objective(x::Vector{Float64})
        try
            params = (x[1], x[2], x[3], x[4], x[5])  # ρ₀, B_A, K, m_ratio, E_sym
            result = base_func(params)
            return isfinite(result) ? result : penalty
        catch
            return penalty
        end
    end
    
    return objective
end

function optimize_with_prior!(gp, bounds, method; max_iterations=15, verbosity=2)
    """
    使用已有高斯过程模型进行贝叶斯优化
    
    参数:
    - gp: 已拟合的高斯过程模型
    - bounds: 参数边界 [(min1, max1), (min2, max2), ...]
    - method: 采集函数方法 (如 EI())
    - max_iterations: 最大迭代次数
    - verbosity: 输出详细程度
    
    返回:
    - result: 优化结果
    """
    
    # 将 GP 包装为 BayesianOptimization 模型
    model = GaussianProcessModel(gp)
    objective = BayesianOptimization.Objective(model=model, sense=MinSense)
    
    # 运行优化
    result = optimize(objective, bounds, method, max_iterations=max_iterations, verbosity=verbosity)
    
    return result
end

function bayesian_optimize_pnjl_simple(kappa_pairs, μ_B_values, T_min, T_max, bounds;
                                      prior_X=nothing, prior_y=nothing,
                                      max_iterations=15, T_step=1.0/hc,
                                      kernel_params=(0.0, 0.0), verbosity=2)
    """
    简化的PNJL模型贝叶斯优化接口
    
    参数:
    - kappa_pairs: κ比值对数组
    - μ_B_values: 重子化学势数组
    - T_min, T_max: 温度范围
    - bounds: 参数边界 [(min1, max1), (min2, max2), ...]
    - prior_X: 先验观测点 [[x1, x2, x3, x4, x5], ...]
    - prior_y: 先验目标值 [y1, y2, ...]
    - max_iterations: 优化迭代次数
    - T_step: 温度扫描步长
    - kernel_params: 核函数初始参数
    - verbosity: 输出详细程度
    
    返回:
    - result: 优化结果
    """
    
    println("="^60)
    println("PNJL模型贝叶斯优化 (简化版)")
    println("="^60)
    
    # 1. 创建目标函数
    objective = create_pnjl_objective(kappa_pairs, μ_B_values, T_min, T_max; 
                                    T_step=T_step, penalty=1e6)
    
    # 2. 设置高斯过程
    n_dims = length(bounds)
    mean = MeanZero()
    kernel = SE(kernel_params[1], kernel_params[2])  # 平方指数核
    
    # 3. 创建和拟合高斯过程模型
    if prior_X !== nothing && prior_y !== nothing
        println("使用 $(length(prior_y)) 个先验数据点")
        
        # 转换数据格式
        X_matrix = hcat(prior_X...)  # 转为矩阵格式，每列一个样本
        y_vector = map(y -> objective([y...]), eachcol(X_matrix))  # 重新计算以确保一致性
        
        # 创建并拟合高斯过程
        gp = GP(X_matrix, y_vector, mean, kernel)
        optimize!(gp)  # 优化超参数
        
    else
        println("从头开始优化（无先验数据）")
        
        # 生成初始样本点
        n_initial = 5
        X_initial = []
        for i in 1:n_initial
            x = [bounds[j][1] + (bounds[j][2] - bounds[j][1]) * rand() for j in 1:n_dims]
            push!(X_initial, x)
        end
        
        X_matrix = hcat(X_initial...)
        y_vector = map(x -> objective(x), X_initial)
        
        # 创建并拟合高斯过程
        gp = GP(X_matrix, y_vector, mean, kernel)
        optimize!(gp)
    end
    
    # 4. 设置优化方法
    model = GaussianProcessModel(gp)
    objective_wrapped = BayesianOptimization.Objective(model=model, sense=MinSense)
    method = BayesianOptimization.EI()  # 期望提升
    
    # 5. 运行优化
    println("开始贝叶斯优化，最大迭代: $max_iterations")
    result = optimize(objective_wrapped, bounds, method, 
                     max_iterations=max_iterations, verbosity=verbosity)
    
    # 6. 输出结果
    println("\n" * "="^50)
    println("优化完成!")
    println("="^50)
    best_params = result.model_args.X[:, end]
    best_value = result.model_args.y[end]
    
    param_names = ["ρ₀", "B_A", "K", "m_ratio", "E_sym"]
    println("最佳参数:")
    for i in 1:length(param_names)
        println("  $(param_names[i]) = $(round(best_params[i], digits=6))")
    end
    println("最佳目标值: $(round(best_value, digits=4))")
    
    return result
end

function demo_simple_bayesian_optimization()
    """
    演示简化版贝叶斯优化的使用
    """
    
    println("="^80)
    println("演示：简化版PNJL贝叶斯优化")
    println("="^80)
    
    # 实验数据
    kappa_pairs = [
        (1.09031788496341, -0.28904867673079),
        (1.06152332992368, 0.164279260625683),
        (1.11111023684003, 0.224522832511389)
    ]
    
    μ_B_values = [632.0/hc, 666.0/hc, 697.0/hc]
    T_min, T_max = 50.0/hc, 150.0/hc
    
    # 参数边界
    bounds = [
        (0.145, 0.170),   # ρ₀
        (-17.0, -15.6),   # B_A  
        (212.0, 401.0),   # K
        (0.55, 0.75),     # m_ratio
        (26.1, 44.0)      # E_sym
    ]
    
    # 示例：创建一些先验数据
    prior_X = [
        [0.155, -16.2, 250.0, 0.65, 32.0],
        [0.160, -16.5, 280.0, 0.68, 35.0],
        [0.150, -16.0, 230.0, 0.62, 30.0]
    ]
    
    println("使用先验数据的优化:")
    # 先计算先验点的目标函数值
    objective = create_pnjl_objective(kappa_pairs, μ_B_values, T_min, T_max; T_step=5.0/hc)
    prior_y = [objective(x) for x in prior_X]
    
    println("先验数据点:")
    for (i, (x, y)) in enumerate(zip(prior_X, prior_y))
        println("  点 $i: $(round.(x, digits=3)) → $(round(y, digits=2))")
    end
    
    # 执行优化
    result = bayesian_optimize_pnjl_simple(
        kappa_pairs, μ_B_values, T_min, T_max, bounds;
        prior_X=prior_X, prior_y=prior_y,
        max_iterations=10, T_step=5.0/hc, verbosity=1
    )
    
    println("\n✅ 简化版贝叶斯优化演示完成!")
    println("对比复杂版本，这个实现:")
    println("  - 代码更简洁，易于理解")
    println("  - 支持先验数据")
    println("  - 减少了不必要的功能")
    println("  - 保持了核心优化能力")
    
    return result
end

function load_prior_data_from_csv(csv_file::String)
    """
    从CSV文件加载先验数据
    
    参数:
    - csv_file: CSV文件路径
    
    返回:
    - (prior_X, prior_y): 先验观测点和目标值，如果失败返回(nothing, nothing)
    """
    
    if !isfile(csv_file)
        println("⚠️  未找到文件: $csv_file")
        return nothing, nothing
    end
    
    try
        df = CSV.read(csv_file, DataFrame)
        
        # 提取参数列
        param_cols = [:rho0, :B_A, :K, :m_ratio, :E_sym]
        if !all(col in names(df) for col in param_cols)
            println("❌ CSV文件缺少必要的参数列")
            return nothing, nothing
        end
        
        prior_X = []
        prior_y = Float64[]
        
        for row in eachrow(df)
            x = [row[col] for col in param_cols]
            y = abs(row.objective_value)  # 确保为正值
            
            if isfinite(y) && y > 0
                push!(prior_X, x)
                push!(prior_y, y)
            end
        end
        
        if length(prior_X) > 0
            println("✅ 从 $csv_file 加载了 $(length(prior_X)) 个有效数据点")
            return prior_X, prior_y
        else
            println("⚠️  CSV文件中没有有效的数据点")
            return nothing, nothing
        end
        
    catch e
        println("❌ 读取CSV文件失败: $e")
        return nothing, nothing
    end
end

function bayesian_optimize_from_csv(csv_file::String, kappa_pairs, μ_B_values, T_min, T_max, bounds;
                                   max_iterations=15, T_step=1.0/hc, verbosity=2)
    """
    从CSV文件加载先验数据并继续贝叶斯优化
    
    参数:
    - csv_file: 包含先验数据的CSV文件
    - 其他参数与 bayesian_optimize_pnjl_simple 相同
    
    返回:
    - result: 优化结果
    """
    
    println("="^80)
    println("从CSV文件继续贝叶斯优化 (简化版)")
    println("="^80)
    println("CSV文件: $csv_file")
    
    # 加载先验数据
    prior_X, prior_y = load_prior_data_from_csv(csv_file)
    
    # 执行优化
    result = bayesian_optimize_pnjl_simple(
        kappa_pairs, μ_B_values, T_min, T_max, bounds;
        prior_X=prior_X, prior_y=prior_y,
        max_iterations=max_iterations, T_step=T_step, verbosity=verbosity
    )
    
    return result
end

function demo_bayesian_with_csv()
    """
    演示从CSV文件加载先验数据的贝叶斯优化
    """
    
    println("="^80)
    println("演示：从CSV文件加载先验数据的贝叶斯优化")
    println("="^80)
    
    # 实验数据
    kappa_pairs = [
        (1.09031788496341, -0.28904867673079),
        (1.06152332992368, 0.164279260625683)
    ]
    
    μ_B_values = [632.0/hc, 666.0/hc]
    T_min, T_max = 50.0/hc, 150.0/hc
    
    # 参数边界
    bounds = [
        (0.145, 0.170),   # ρ₀
        (-17.0, -15.6),   # B_A  
        (212.0, 401.0),   # K
        (0.55, 0.75),     # m_ratio
        (26.1, 44.0)      # E_sym
    ]
    
    # 尝试从现有CSV文件加载数据
    csv_file = "Rotation_PNJL/output/Gas_Liquid/demo_bayesian_optimization_warmup.csv"
    
    result = bayesian_optimize_from_csv(
        csv_file, kappa_pairs, μ_B_values, T_min, T_max, bounds;
        max_iterations=10, T_step=5.0/hc, verbosity=1
    )
    
    println("\n✅ 从CSV加载的贝叶斯优化演示完成!")
    
    return result
end

# =============================================================================
# 简化版使用指南和总结
# =============================================================================

"""
# 简化版贝叶斯优化使用指南

## 主要改进

相比原来臃肿的实现，新的简化版本提供了：

1. **简洁的API** - 仿照 BayesianOptimization.jl 的风格
2. **先验数据支持** - 可以从已有数据继续优化  
3. **减少复杂性** - 去除了预热、详细日志等非核心功能
4. **保持向后兼容** - 原有复杂函数仍然可用

## 核心函数

### 1. 基础优化
```julia
result = bayesian_optimize_pnjl_simple(
    kappa_pairs, μ_B_values, T_min, T_max, bounds;
    max_iterations=15, T_step=1.0/hc, verbosity=2
)
```

### 2. 使用先验数据
```julia
prior_X = [[0.155, -16.2, 250.0, 0.65, 32.0], ...]
prior_y = [125.3, 143.7, ...]

result = bayesian_optimize_pnjl_simple(
    kappa_pairs, μ_B_values, T_min, T_max, bounds;
    prior_X=prior_X, prior_y=prior_y,
    max_iterations=15
)
```

### 3. 从CSV文件继续优化
```julia
result = bayesian_optimize_from_csv(
    "previous_results.csv", 
    kappa_pairs, μ_B_values, T_min, T_max, bounds;
    max_iterations=15
)
```

## 参数说明

- `kappa_pairs`: κ比值对 [(κ₃/κ₁, κ₄/κ₂), ...]
- `μ_B_values`: 重子化学势数组（以 1/hc 为单位）
- `T_min, T_max`: 温度搜索范围（以 1/hc 为单位）
- `bounds`: 参数边界 [(min1,max1), (min2,max2), ...]，顺序为 [ρ₀, B_A, K, m_ratio, E_sym]
- `prior_X`: 先验观测点列表 [[x1,x2,x3,x4,x5], ...]
- `prior_y`: 先验目标值列表 [y1, y2, ...]
- `max_iterations`: 优化迭代次数
- `T_step`: 温度扫描步长
- `verbosity`: 输出详细程度 (0=静默, 1=基本, 2=详细)

## 演示函数

- `demo_simple_bayesian_optimization()` - 基础演示
- `demo_bayesian_with_csv()` - CSV文件演示

## 与原有代码的关系

简化版本位于文件末尾，不影响原有复杂功能。你可以：
- 使用简化版进行快速开发和测试
- 需要高级功能时继续使用原有函数
- 逐步迁移到简化版本

## 使用建议

1. 新项目优先使用简化版本
2. 先用少量迭代测试，确认目标函数正常工作
3. 有先验数据时一定要使用，可以显著提高效率
4. 调整 T_step 平衡计算速度和精度
"""

function usage_example_simplified()
    """
    简化版贝叶斯优化的完整使用示例
    """
    
    println("="^80)
    println("简化版贝叶斯优化完整使用示例")
    println("="^80)
    
    # 1. 定义实验数据
    kappa_pairs = [
        (1.09031788496341, -0.28904867673079),
        (1.06152332992368, 0.164279260625683)
    ]
    μ_B_values = [632.0/hc, 666.0/hc]
    T_min, T_max = 60.0/hc, 140.0/hc
    
    # 2. 定义参数边界
    bounds = [
        (0.150, 0.165),   # ρ₀ (fm⁻³)
        (-16.8, -16.0),   # B_A (MeV)  
        (230.0, 280.0),   # K (MeV)
        (0.60, 0.70),     # m_ratio
        (30.0, 36.0)      # E_sym (MeV)
    ]
    
    println("步骤 1: 从头开始优化（无先验数据）")
    result1 = bayesian_optimize_pnjl_simple(
        kappa_pairs, μ_B_values, T_min, T_max, bounds;
        max_iterations=8, T_step=8.0/hc, verbosity=1
    )
    
    # 3. 模拟从第一次优化中提取先验数据
    println("\n步骤 2: 使用先验数据继续优化")
    
    # 提取前几个最好的点作为先验数据
    all_X = [result1.model_args.X[:, i] for i in 1:size(result1.model_args.X, 2)]
    all_y = result1.model_args.y
    
    # 选择目标值最小的3个点作为先验
    sorted_indices = sortperm(all_y)
    best_indices = sorted_indices[1:min(3, length(sorted_indices))]
    
    prior_X = [all_X[i] for i in best_indices]
    prior_y = [all_y[i] for i in best_indices]
    
    println("使用的先验数据:")
    for (i, (x, y)) in enumerate(zip(prior_X, prior_y))
        println("  先验点 $i: [$(join(round.(x, digits=3), ", "))] → $(round(y, digits=2))")
    end
    
    # 用先验数据再次优化
    result2 = bayesian_optimize_pnjl_simple(
        kappa_pairs, μ_B_values, T_min, T_max, bounds;
        prior_X=prior_X, prior_y=prior_y,
        max_iterations=5, T_step=8.0/hc, verbosity=1
    )
    
    # 4. 比较结果
    println("\n结果比较:")
    println("  第一次优化最优值: $(round(minimum(result1.model_args.y), digits=3))")
    println("  第二次优化最优值: $(round(minimum(result2.model_args.y), digits=3))")
    
    improvement = minimum(result1.model_args.y) - minimum(result2.model_args.y)
    if improvement > 0
        println("  改进: $(round(improvement, digits=3)) ($(round(improvement/minimum(result1.model_args.y)*100, digits=1))%)")
    else
        println("  第二次优化未发现更好的解")
    end
    
    println("\n✅ 完整使用示例完成!")
    println("\n总结:")
    println("  - 简化版本大大减少了代码复杂性")
    println("  - 支持先验数据可以加速收敛")
    println("  - API设计更加直观易用")
    println("  - 适合快速原型开发和测试")
    
    return result1, result2
end

# =============================================================================
# 快速测试函数
# =============================================================================

function quick_test_simplified_api()
    """
    快速测试简化版API是否工作正常
    """
    
    println("="^60)
    println("快速测试简化版贝叶斯优化API")
    println("="^60)
    
    # 使用最小配置进行测试
    kappa_pairs = [(1.09, -0.29)]
    μ_B_values = [632.0/hc]
    T_min, T_max = 80.0/hc, 120.0/hc
    
    bounds = [
        (0.150, 0.160),   # ρ₀
        (-16.5, -16.0),   # B_A  
        (240.0, 260.0),   # K
        (0.65, 0.70),     # m_ratio
        (31.0, 33.0)      # E_sym
    ]
    
    println("测试配置:")
    println("  数据点: 1 组")
    println("  温度范围: $(T_min*hc)-$(T_max*hc) MeV")
    println("  迭代次数: 3")
    
    try
        result = bayesian_optimize_pnjl_simple(
            kappa_pairs, μ_B_values, T_min, T_max, bounds;
            max_iterations=3, T_step=10.0/hc, verbosity=0
        )
        
        println("✅ 简化版API测试成功!")
        println("   总评估: $(length(result.model_args.y))")
        println("   最优值: $(round(minimum(result.model_args.y), digits=2))")
        
        return true
        
    catch e
        println("❌ 简化版API测试失败: $e")
        return false
    end
end
