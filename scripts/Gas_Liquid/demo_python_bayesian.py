def setup_chinese_font():
#!/usr/bin/env python3
"""
NOTICE: This demo script has been moved.

The original file has been relocated to:
  examples/Gas_Liquid/demo_python_bayesian.py

Please edit and run the copy in the `examples/Gas_Liquid` directory.
"""

import sys

print("This demo has been moved to examples/Gas_Liquid/demo_python_bayesian.py")
sys.exit(0)

def check_dependencies():
    """检查所需依赖"""
    print("检查Python依赖...")
    
    missing_deps = []
    
    try:
        import numpy as np
        print("✅ numpy")
    except ImportError:
        missing_deps.append("numpy")
        print("❌ numpy")
    
    try:
        import pandas as pd
        print("✅ pandas")
    except ImportError:
        missing_deps.append("pandas")
        print("❌ pandas")
    
    try:
        print("✅ matplotlib")
    except ImportError:
        missing_deps.append("matplotlib")
        print("❌ matplotlib")
    
    try:
        from skopt import gp_minimize
        print("✅ scikit-optimize")
    except ImportError:
        missing_deps.append("scikit-optimize")
        print("❌ scikit-optimize")
    
    try:
        import julia
        print("✅ julia")
    except ImportError:
        missing_deps.append("julia")
        print("❌ julia")
    
    if missing_deps:
        print(f"\n❌ 缺少依赖: {missing_deps}")
        print("请运行以下命令安装:")
        print(f"pip install {' '.join(missing_deps)}")
        
        if 'julia' in missing_deps:
            print("\n对于julia，还需要运行:")
            print("python -c \"import julia; julia.install()\"")
        
        return False
    
    print("\n✅ 所有依赖都已安装")
    return True

def demo_simple_optimization():
    """简单的优化演示"""
    
    if not check_dependencies():
        return None
    
    print("\n" + "="*80)
    print("PNJL模型Python贝叶斯优化演示")
    print("="*80)
    
    import numpy as np
    import pandas as pd
    from skopt import gp_minimize
    from skopt.space import Real
    import time
    from datetime import datetime
    
    # 设置中文字体
    font_name = setup_chinese_font()
    use_chinese = font_name is not None
    
    # Julia接口
    try:
        import julia
        from julia import Main as jl
        import os
        
        # 设置Julia环境
        print("初始化Julia环境...")
        
        # 确保在正确的项目目录下
        script_path = os.path.dirname(os.path.abspath(__file__))
        julia_project_path = os.path.dirname(os.path.dirname(script_path))  # 从scripts/Gas_Liquid回到项目根目录
        print(f"Julia项目路径: {julia_project_path}")
        
        # 切换到Julia项目目录
        original_cwd = os.getcwd()
        os.chdir(julia_project_path)
        print(f"当前工作目录: {os.getcwd()}")
        
        # 激活Julia项目环境
        jl.eval('using Pkg; Pkg.activate(".")')
        
        # 检查文件是否存在
        julia_file_path = "src/Gas_Liquid/Advanced_FindTforDiff.jl"
        full_path = os.path.join(julia_project_path, julia_file_path)
        print(f"尝试加载Julia文件: {full_path}")
        
        if not os.path.exists(full_path):
            print(f"❌ 文件不存在: {full_path}")
            print("可用的文件:")
            gas_liquid_dir = "src/Gas_Liquid"
            if os.path.exists(gas_liquid_dir):
                for file in os.listdir(gas_liquid_dir):
                    if file.endswith('.jl'):
                        print(f"  - {file}")
            os.chdir(original_cwd)  # 恢复原始目录
            return None
        
        # 加载Julia文件和必要的包
        jl.eval('using Pkg')
        jl.eval('Pkg.instantiate()')  # 确保所有依赖都已安装
        jl.eval(f'include("{julia_file_path}")')
        
        # 检查函数是否可用
        try:
            jl.eval('typeof(create_temperature_difference_objective)')
            print("✅ 目标函数 create_temperature_difference_objective 已加载")
        except Exception as e:
            print(f"❌ 无法找到目标函数: {e}")
            os.chdir(original_cwd)
            return None
            
        print("✅ Julia环境初始化成功")
        
        # 恢复原始工作目录
        os.chdir(original_cwd)
        
    except Exception as e:
        print(f"❌ Julia初始化失败: {e}")
        return None
    
    # 常数
    hc = 197.327  # MeV·fm
    
    # 实验数据（与Julia版本相同）
    kappa_pairs = [
        (1.09031788496341, -0.28904867673079),   # 第1组
        (1.06152332992368, 0.164279260625683),   # 第2组
        (1.11111023684003, 0.224522832511389)    # 第3组
    ]
    
    mu_B_values = [632.0, 666.0, 697.0]  # MeV
    T_min, T_max = 70.0, 120.0           # MeV
    T_step_scan = 2.0                    # MeV
    
    # 转换到Julia单位
    mu_B_julia = [mu / hc for mu in mu_B_values]
    T_min_julia = T_min / hc
    T_max_julia = T_max / hc
    T_step_julia = T_step_scan / hc
    
    print(f"实验数据: {len(kappa_pairs)} 组κ比值对")
    print(f"μ_B值: {mu_B_values} MeV")
    print(f"温度范围: {T_min} - {T_max} MeV")
    
    # 在Julia端创建并重用目标函数闭包，避免在每次评估中重复 include/eval（会导致性能下降）
    jl_base_objective = None
    try:
        # 注意：此处使用Python端的列表，pyjulia会进行类型转换
        jl_base_objective = jl.create_temperature_difference_objective(
            kappa_pairs, [mu / hc for mu in mu_B_values], T_min_julia, T_max_julia,
            T_step_scan=T_step_julia, penalty_for_missing=1e6, verbose=False)
        print("✅ 已在Julia端创建目标函数闭包 base_objective（将被重用）")
    except Exception as e:
        print(f"⚠️ 无法在Julia端创建目标函数闭包（将回退到eval）：{e}")
    
    def objective_function(params):
        """目标函数：优先直接调用在Julia端创建并重用的闭包，回退到eval以保证兼容性。"""
        try:
            rho0, B_A, K, m_ratio, E_sym = params

            if jl_base_objective is not None:
                # 直接调用Julia闭包（pyjulia会做类型转换）
                params_tuple = (rho0, B_A, K, m_ratio, E_sym)
                jl_result = jl_base_objective(params_tuple)
            else:
                # 回退到原先的eval方式（较慢）以保证功能完整
                kappa_pairs_str = str(kappa_pairs).replace('(', '[').replace(')', ']')
                mu_B_julia_str = str(mu_B_julia)
                jl_code = f"""
base_objective = create_temperature_difference_objective({kappa_pairs_str}, {mu_B_julia_str}, {T_min_julia}, {T_max_julia}; T_step_scan={T_step_julia}, penalty_for_missing=1e6, verbose=false)
params_tuple = ({rho0}, {B_A}, {K}, {m_ratio}, {E_sym})
base_objective(params_tuple)
"""
                jl_result = jl.eval(jl_code)

            return float(jl_result) if np.isfinite(jl_result) else 1e6

        except Exception as e:
            print(f"目标函数评估失败: {e}")
            return 1e6
    
    # 参数边界（与Julia版本相同）
    param_bounds = [
        (0.145, 0.170),    # ρ₀ (fm⁻³)
        (-17.0, -15.6),    # B_A (MeV)
        (212.0, 401.0),    # K (MeV)
        (0.55, 0.75),      # m_ratio
        (26.1, 44.0)       # E_sym (MeV)
    ]
    
    # 参数名称（根据字体支持选择中文或英文）
    if use_chinese:
        param_names = ["ρ₀ (fm⁻³)", "B_A (MeV)", "K (MeV)", "m_ratio", "E_sym (MeV)"]
        title1 = "优化收敛历史"
        title2 = "累积最优值"
        xlabel = "迭代次数"
        ylabel1 = "目标函数值"
        ylabel2 = "最优目标函数值"
        legend1 = "最优值"
        legend2 = "累积最优"
    else:
        param_names = ["rho0 (fm^-3)", "B_A (MeV)", "K (MeV)", "m_ratio", "E_sym (MeV)"]
        title1 = "Optimization Convergence History"
        title2 = "Cumulative Best Value"
        xlabel = "Iteration"
        ylabel1 = "Objective Function Value"
        ylabel2 = "Best Objective Function Value"
        legend1 = "Best Value"
        legend2 = "Cumulative Best"
    
    print("\n参数边界:")
    for name, (low, high) in zip(param_names, param_bounds):
        print(f"  {name}: {low} - {high}")
    
    # 设置优化空间
    dimensions = [Real(low, high) for low, high in param_bounds]
    
    # 预热测试
    print("\n预热目标函数...")
    test_params = [0.155, -16.2, 250.0, 0.65, 32.0]
    start_time = time.time()
    test_result = objective_function(test_params)
    eval_time = time.time() - start_time
    
    print(f"预热结果: {test_result:.4f} (用时: {eval_time:.2f}秒)")
    
    # 执行贝叶斯优化
    print(f"\n开始贝叶斯优化...")
    print(f"开始时间: {datetime.now()}")
    
    max_iterations = 150  # 减少迭代次数以便快速演示
    
    optimization_start = time.time()
    
    try:
        result = gp_minimize(
            func=objective_function,
            dimensions=dimensions,
            n_calls=max_iterations,
            n_initial_points=8,
            acq_func='EI',  # Expected Improvement
            random_state=42,
            verbose=True
        )
        
        optimization_time = time.time() - optimization_start
        
        print("\n" + "="*80)
        print("✅ 贝叶斯优化完成!")
        print("="*80)
        print(f"优化时间: {optimization_time:.2f} 秒")
        print(f"总评估次数: {len(result.func_vals)}")
        
        # 最优结果
        best_params = result.x
        best_value = result.fun
        
        print(f"\n最优参数:")
        for name, value in zip(param_names, best_params):
            print(f"  {name} = {value:.6f}")
        print(f"\n最优目标函数值: {best_value:.4f} MeV²")
        
        # 简单的收敛图
        plt.figure(figsize=(12, 6))
        
        plt.subplot(1, 2, 1)
        plt.plot(result.func_vals, 'b-o', markersize=4)
        plt.axhline(y=best_value, color='r', linestyle='--', alpha=0.7, label=legend1)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel1)
        plt.title(title1)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 累积最优值
        cumulative_best = []
        current_best = float('inf')
        for val in result.func_vals:
            if val < current_best:
                current_best = val
            cumulative_best.append(current_best)
        
        plt.subplot(1, 2, 2)
        plt.plot(cumulative_best, 'g-', linewidth=2, label=legend2)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel2)
        plt.title(title2)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # 确保输出目录存在
        output_dir = os.path.join(julia_project_path, "output", "Gas_Liquid")
        os.makedirs(output_dir, exist_ok=True)
        
        plot_path = os.path.join(output_dir, 'demo_optimization_result.png')
        csv_path = os.path.join(output_dir, 'demo_optimization_result.csv')
        
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close()  # 关闭图形以节省内存
        
        # 保存结果到CSV
        optimization_data = []
        for i, (params, value) in enumerate(zip(result.x_iters, result.func_vals)):
            row = {
                'iteration': i + 1,
                'rho0': params[0],
                'B_A': params[1], 
                'K': params[2],
                'm_ratio': params[3],
                'E_sym': params[4],
                'objective_value': value,
                'is_best': (i == np.argmin(result.func_vals))
            }
            optimization_data.append(row)
        
        df = pd.DataFrame(optimization_data)
        df.to_csv(csv_path, index=False)
        print(f"\n✅ 结果已保存到: {csv_path}")
        print(f"✅ 收敛图已保存到: {plot_path}")
        
        return result
        
    except Exception as e:
        print(f"\n❌ 优化失败: {e}")
        return None

def main():
    """主函数"""
    print("PNJL模型Python贝叶斯优化演示")
    print("="*50)
    
    result = demo_simple_optimization()
    
    if result is not None:
        print("\n🎉 演示成功完成!")
        print("\n接下来您可以:")
        print("1. 查看 output/Gas_Liquid/demo_optimization_result.csv 了解详细结果")
        print("2. 查看 output/Gas_Liquid/demo_optimization_result.png 了解收敛情况")
        print("3. 使用完整版本 src/Gas_Liquid/Advanced_Bayesian.py 进行更复杂的优化")
        print("4. 修改参数进行更多实验")
    else:
        print("\n❌ 演示失败，请检查依赖和环境配置")

if __name__ == "__main__":
    main()