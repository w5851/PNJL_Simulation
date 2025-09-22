include("Function_Rotation.jl")
using .Function_Rotation: get_nodes, calculate_t_rho, nlsolve, hc

function generate_single_temperature_data(T_target)
    """为单个温度生成完整的rho-mu数据
    T_target: Julia 内部单位的温度
    返回: rho_mu_pairs (数组 of (rho, mu))
    """
    nodes1, nodes2 = get_nodes(128, 16)
    omega = 100/hc
    
    # 初始猜测值
    x_initial = [-2.13, 0.06, 0.12, 310 / hc]
    x = copy(x_initial)
    
    rho_mu_pairs = []
    
    # 从高密度到低密度扫描（与Trho函数一致）
    for rho in [3.00; collect(2.99:-0.01:0.10)]
        converged = false
        try
            res = nlsolve(x -> calculate_t_rho(x, T_target, rho, nodes1, omega), x)
            converged = res.f_converged
            if converged
                copyto!(x, res.zero)
                mu_value = x[4]
                push!(rho_mu_pairs, (rho, mu_value))
            else
                @debug "Root finding did not converge for T=$T_target and rho=$rho"
            end
        catch err
            @debug "Exception in root finding for T=$T_target and rho=$rho: $err"
            converged = false
        end
    end
    
    return rho_mu_pairs
end

function has_s_shape(T_mev; min_negative_points=3, derivative_threshold=-0.001)
    """
    检测给定温度下是否存在S形
    参数:
    - T_mev: 温度 (MeV)
    - min_negative_points: 至少需要多少个负导数点才认为存在S形
    - derivative_threshold: 负导数阈值
    返回: true表示存在S形，false表示不存在
    """
    T_julia = T_mev / hc  # 转换为Julia内部单位
    
    rho_mu_pairs = generate_single_temperature_data(T_julia)
    
    if length(rho_mu_pairs) < 3
        @debug "温度 $T_mev MeV 处数据点太少，无法判断S形"
        return false
    end
    
    # 按rho排序（从小到大）
    sorted_pairs = sort(rho_mu_pairs, by=x->x[1])
    
    # 计算数值导数 ∂mu/∂rho 并统计负值个数
    negative_derivative_count = 0
    max_negative_derivative = 0.0
    
    for i in 2:length(sorted_pairs)
        rho1, mu1 = sorted_pairs[i-1]
        rho2, mu2 = sorted_pairs[i]
        
        if rho2 != rho1  # 避免除零
            derivative = (mu2 - mu1) / (rho2 - rho1)
            if derivative < derivative_threshold
                negative_derivative_count += 1
                max_negative_derivative = min(max_negative_derivative, derivative)
            end
        end
    end
    
    has_s = negative_derivative_count >= min_negative_points
    
    if has_s
        println("  发现S形: $negative_derivative_count 个负导数点, 最小导数值 = $max_negative_derivative")
    end
    
    return has_s
end

function find_cep(T_min=10.0, T_max=150.0, tolerance=0.5)
    """
    找到温度临界点T_cep，该温度下rho-mu关系正好不出现"S形"
    返回: T_cep (MeV)
    """
    println("开始寻找临界温度 T_cep...")
    println("搜索范围: [$T_min, $T_max] MeV")
    println("精度: $tolerance MeV")
    
    # 首先扫描整个范围，寻找S形存在的区域
    println("扫描温度范围寻找S形现象...")
    
    scan_temps = collect(T_min:5:T_max)
    s_shape_results = []
    
    for T_scan in scan_temps
        println("扫描 T=$T_scan MeV...")
        has_s = has_s_shape(T_scan)
        push!(s_shape_results, (T_scan, has_s))
        println("  $(has_s ? "有S形" : "无S形")")
    end
    
    # 寻找相变区间
    transition_found = false
    T_low = T_min
    T_high = T_max
    
    for i in 1:(length(s_shape_results)-1)
        T1, has_s1 = s_shape_results[i]
        T2, has_s2 = s_shape_results[i+1]
        
        if has_s1 != has_s2  # 找到相变
            T_low = has_s1 ? T2 : T1  # 有S形的温度作为下界
            T_high = has_s1 ? T1 : T2  # 无S形的温度作为上界
            transition_found = true
            println("发现相变区间: [$T_low, $T_high] MeV")
            break
        end
    end
    
    if !transition_found
        println("在扫描范围内未发现S形相变，尝试扩展搜索...")
        
        # 分析初始扫描的结果
        all_have_s = all(result[2] for result in s_shape_results)
        all_no_s = all(!result[2] for result in s_shape_results)
        
        if all_have_s
            # 整个范围都有S形，说明临界点在更高温度
            println("初始范围内都有S形，向高温扩展搜索...")
            for T_test in (T_max+5):10:200
                println("测试高温 T=$T_test MeV...")
                current_has_s = has_s_shape(T_test)
                
                if !current_has_s
                    println("  无S形")
                    # 找到S形消失的点，临界点在T_max和T_test之间
                    T_low = T_max
                    T_high = T_test
                    transition_found = true
                    println("在高温区发现相变区间: [$T_low, $T_high] MeV")
                    break
                else
                    println("  有S形")
                end
            end
        elseif all_no_s
            # 整个范围都无S形，说明临界点在更低温度
            println("初始范围内都无S形，向低温扩展搜索...")
            for T_test in (T_min-5):-10:10
                println("测试低温 T=$T_test MeV...")
                current_has_s = has_s_shape(T_test)
                
                if current_has_s
                    println("  有S形")
                    # 找到S形出现的点，临界点在T_test和T_min之间
                    T_low = T_test
                    T_high = T_min
                    transition_found = true
                    println("在低温区发现相变区间: [$T_low, $T_high] MeV")
                    break
                else
                    println("  无S形")
                end
            end
        end
    end
    
    if !transition_found
        @warn "未找到S形相变，返回搜索范围中点作为估计值"
        return (T_min + T_max) / 2
    end
    
    # 二分搜索精确定位临界点
    println("开始二分搜索...")
    
    iteration = 0
    while (T_high - T_low) > tolerance && iteration < 15
        iteration += 1
        T_mid = (T_low + T_high) / 2
        
        println("第 $iteration 次迭代: 测试 T=$T_mid MeV")
        
        mid_has_s = has_s_shape(T_mid)
        println("  结果: $(mid_has_s ? "有S形" : "无S形")")
        
        if mid_has_s
            T_low = T_mid  # S形存在，临界点在更高温度
        else
            T_high = T_mid  # S形不存在，临界点在更低温度
        end
        
        println("  新搜索区间: [$T_low, $T_high]")
    end
    
    T_cep = (T_low + T_high) / 2
    println("找到临界温度: T_cep = $T_cep MeV")
    
    return T_cep
end


