"""
Advanced Bayesian Optimization for PNJL Model Parameters
使用Python的贝叶斯优化库调用Julia函数实现PNJL模型参数优化

基于Julia中的demo_bayesian_optimization_with_warmup()功能
依赖：
- scikit-optimize: Python贝叶斯优化库
- julia: Python调用Julia的接口
- numpy, pandas: 数据处理
- matplotlib: 可视化

作者：AI Assistant
日期：2025年9月19日
"""

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import time
from datetime import datetime
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# 贝叶斯优化相关库
try:
    from skopt import gp_minimize
    from skopt.space import Real
    from skopt.acquisition import gaussian_ei, gaussian_lcb, gaussian_pi
    from skopt.utils import use_named_args
    from skopt import dump, load
    print("✅ scikit-optimize库导入成功")
except ImportError as e:
    print(f"❌ 请安装scikit-optimize: pip install scikit-optimize")
    print(f"导入错误: {e}")

# Julia接口
try:
    import julia
    from julia import Main as jl
    print("✅ PyJulia库导入成功")
except ImportError as e:
    print(f"❌ 请安装PyJulia: pip install julia")
    print(f"导入错误: {e}")
    jl = None

class PNJLBayesianOptimizer:
    """
    PNJL模型参数贝叶斯优化器
    
    使用Python的scikit-optimize库进行贝叶斯优化，
    调用Julia函数进行PNJL模型的物理计算
    """
    
    def __init__(self, julia_project_path="d:/Desktop/Julia/Rotation_PNJL", verbose=True):
        """
        初始化优化器
        
        参数:
        - julia_project_path: Julia项目路径
        - verbose: 是否显示详细信息
        """
        self.julia_project_path = Path(julia_project_path)
        self.verbose = verbose
        self.hc = 197.327  # MeV·fm conversion factor
        
        # 参数名称和边界
        self.param_names = ["ρ₀ (fm⁻³)", "B_A (MeV)", "K (MeV)", "m_ratio", "E_sym (MeV)"]
        self.param_symbols = ["ρ₀", "B_A", "K", "m_ratio", "E_sym"]
        
        # 优化历史
        self.optimization_history = []
        self.warmup_results = []
        
        # 初始化Julia环境
        self._setup_julia_environment()
    
    def _setup_julia_environment(self):
        """设置Julia环境和加载必要的模块"""
        if jl is None:
            raise ImportError("PyJulia未安装，无法调用Julia函数")
        
        if self.verbose:
            print("="*80)
            print("初始化Julia环境")
            print("="*80)
            print(f"Julia项目路径: {self.julia_project_path}")
        
        try:
            # 激活Julia项目环境
            jl.eval(f'using Pkg; Pkg.activate("{str(self.julia_project_path).replace(chr(92), "/")}")')            # 加载必要的Julia模块
            julia_file_path = self.julia_project_path / "src" / "Gas_Liquid" / "Advanced_FindTforDiff.jl"
            
            if self.verbose:
                print(f"加载Julia文件: {julia_file_path}")
            
            # 加载Julia文件
            jl.eval(f'include("{str(julia_file_path).replace(chr(92), "/")}")')
            
            # 测试Julia函数是否可用
            test_result = jl.eval("typeof(find_temperature_for_kappa_ratios)")
            
            if self.verbose:
                print(f"✅ Julia环境初始化成功")
                print(f"find_temperature_for_kappa_ratios函数类型: {test_result}")
            
        except Exception as e:
            print(f"❌ Julia环境初始化失败: {e}")
            raise e
    
    def create_objective_function(self, kappa_pairs, mu_B_values, T_min, T_max,
                                T_step_scan=None, penalty_for_missing=1e6):
        """
        创建PNJL模型参数优化的目标函数
        
        参数:
        - kappa_pairs: κ比值对列表 [(κ₃/κ₁, κ₄/κ₂), ...]
        - mu_B_values: 重子化学势数组 [μ_B1, μ_B2, ...] (MeV单位)
        - T_min, T_max: 温度搜索范围 (MeV单位)
        - T_step_scan: 温度扫描步长 (MeV单位)
        - penalty_for_missing: 计算失败时的惩罚值
        
        返回:
        - objective_function: 目标函数 f(params) -> float
        """
        
        if T_step_scan is None:
            T_step_scan = 1.0  # 默认1 MeV
        
        # 转换单位到Julia的自然单位 (1/hc)
        mu_B_julia = [mu / self.hc for mu in mu_B_values]
        T_min_julia = T_min / self.hc
        T_max_julia = T_max / self.hc
        T_step_julia = T_step_scan / self.hc
        
        if self.verbose:
            print("="*80)
            print("创建贝叶斯优化目标函数")
            print("="*80)
            print("实验数据配置:")
            print(f"  κ比值对数量: {len(kappa_pairs)}")
            print(f"  μ_B值: {mu_B_values} MeV")
            print(f"  温度范围: {T_min} - {T_max} MeV")
            print(f"  扫描步长: {T_step_scan} MeV")
            print(f"  惩罚值: {penalty_for_missing}")
        
        def objective_function(params):
            """
            目标函数：计算给定PNJL参数下的温度差平方和
            
            参数:
            - params: [ρ₀, B_A, K, m_ratio, E_sym]
            
            返回:
            - 温度差平方和 (较小值表示更好的拟合)
            """
            try:
                rho0, B_A, K, m_ratio, E_sym = params
                
                # 调用Julia函数计算温度差
                jl_result = jl.eval(f"""
                    # 设置PNJL模型参数
                    rho0 = {rho0}
                    B_A = {B_A}
                    K = {K}
                    m_ratio = {m_ratio}
                    E_sym = {E_sym}
                    
                    # 实验数据
                    kappa_pairs = {kappa_pairs}
                    mu_B_values = {mu_B_julia}
                    T_min = {T_min_julia}
                    T_max = {T_max_julia}
                    T_step = {T_step_julia}
                    
                    # 创建温度差计算函数
                    base_objective = create_temperature_difference_objective(
                        kappa_pairs, mu_B_values, T_min, T_max;
                        T_step_scan=T_step, penalty_for_missing={penalty_for_missing}, verbose=false)
                    
                    # 计算目标函数值
                    params_tuple = (rho0, B_A, K, m_ratio, E_sym)
                    result = base_objective(params_tuple)
                    result
                """)
                
                # 确保返回有限值
                if not np.isfinite(jl_result):
                    return penalty_for_missing
                
                return float(jl_result)
                
            except Exception as e:
                if self.verbose:
                    print(f"目标函数评估失败: {e}")
                return penalty_for_missing
        
        return objective_function
    
    def warmup_objective_function(self, objective_function, param_bounds, n_warmup_samples=3):
        """
        预热目标函数，估算计算时间
        
        参数:
        - objective_function: 目标函数
        - param_bounds: 参数边界 [(min1, max1), (min2, max2), ...]
        - n_warmup_samples: 预热样本数
        
        返回:
        - warmup_results: 预热计算结果列表
        - estimated_time_per_eval: 估算的单次评估时间 (秒)
        """
        
        if self.verbose:
            print("="*60)
            print("目标函数预热")
            print("="*60)
            print(f"预热样本数: {n_warmup_samples}")
        
        warmup_results = []
        total_time = 0.0
        
        for i in range(n_warmup_samples):
            # 生成随机参数
            test_params = []
            for (low, high) in param_bounds:
                test_params.append(np.random.uniform(low, high))
            
            if self.verbose:
                print(f"\n预热样本 {i+1}/{n_warmup_samples}")
                print(f"  测试参数: {[round(p, 4) for p in test_params]}")
            
            # 计算目标函数
            start_time = time.time()
            try:
                result = objective_function(test_params)
                eval_time = time.time() - start_time
                total_time += eval_time
                
                warmup_results.append({
                    'params': test_params.copy(),
                    'value': result,
                    'time': eval_time,
                    'success': True
                })
                
                if self.verbose:
                    print(f"  结果: {result:.4f}")
                    print(f"  用时: {eval_time:.2f} 秒")
                    
            except Exception as e:
                eval_time = time.time() - start_time
                total_time += eval_time
                
                warmup_results.append({
                    'params': test_params.copy(),
                    'value': None,
                    'time': eval_time,
                    'success': False,
                    'error': str(e)
                })
                
                if self.verbose:
                    print(f"  失败: {e}")
                    print(f"  用时: {eval_time:.2f} 秒")
        
        estimated_time_per_eval = total_time / n_warmup_samples
        
        if self.verbose:
            print(f"\n预热阶段完成")
            print(f"总用时: {total_time:.2f} 秒")
            print(f"平均单次评估时间: {estimated_time_per_eval:.2f} 秒")
            
            success_count = sum(1 for r in warmup_results if r['success'])
            print(f"成功率: {success_count}/{n_warmup_samples} ({100*success_count/n_warmup_samples:.1f}%)")
            
            if success_count > 0:
                successful_values = [r['value'] for r in warmup_results if r['success']]
                print(f"目标函数值范围: {min(successful_values):.2f} - {max(successful_values):.2f}")
        
        self.warmup_results = warmup_results
        return warmup_results, estimated_time_per_eval
    
    def optimize_with_warmup(self, kappa_pairs, mu_B_values, T_min, T_max, param_bounds,
                           max_iterations=20, initial_samples=10, T_step_scan=2.0,
                           acquisition_function='EI', warmup_samples=3, skip_warmup=False,
                           output_file=None, random_state=42):
        """
        带预热的PNJL模型参数贝叶斯优化
        
        参数:
        - kappa_pairs: κ比值对列表
        - mu_B_values: 重子化学势数组 (MeV)
        - T_min, T_max: 温度范围 (MeV)
        - param_bounds: 参数边界 [(min1, max1), ...]
        - max_iterations: 最大迭代次数
        - initial_samples: 初始随机采样数
        - T_step_scan: 温度扫描步长 (MeV)
        - acquisition_function: 采集函数 ('EI', 'LCB', 'PI')
        - warmup_samples: 预热样本数
        - skip_warmup: 是否跳过预热
        - output_file: 结果保存文件
        - random_state: 随机种子
        
        返回:
        - optimization_result: 优化结果字典
        """
        
        print("="*100)
        print("PNJL模型参数贝叶斯优化 (带预热)")
        print("="*100)
        print(f"开始时间: {datetime.now()}")
        
        # 验证输入
        if len(mu_B_values) != len(kappa_pairs):
            raise ValueError("μ_B值数组长度与κ比值对数组长度不匹配")
        
        if len(param_bounds) != 5:
            raise ValueError("参数边界必须是5维向量 [ρ₀, B_A, K, m_ratio, E_sym]")
        
        # 显示优化配置
        print("\n优化配置:")
        print(f"  实验数据点: {len(kappa_pairs)} 组")
        print("  参数边界:")
        for i, (name, (low, high)) in enumerate(zip(self.param_names, param_bounds)):
            print(f"    {name}: {low} - {high}")
        print(f"  最大迭代: {max_iterations}")
        print(f"  初始采样: {initial_samples}")
        print(f"  采集函数: {acquisition_function}")
        print(f"  预热样本: {warmup_samples}")
        
        # 创建目标函数
        print("\n创建目标函数...")
        objective_function = self.create_objective_function(
            kappa_pairs, mu_B_values, T_min, T_max, T_step_scan)
        
        # 预热阶段
        warmup_time = 0.0
        estimated_time_per_eval = 0.0
        
        if not skip_warmup:
            print("\n" + "="*60)
            print("第一阶段：目标函数预热")
            print("="*60)
            
            warmup_start = time.time()
            warmup_results, estimated_time_per_eval = self.warmup_objective_function(
                objective_function, param_bounds, n_warmup_samples=warmup_samples)
            warmup_time = time.time() - warmup_start
            
            print(f"\n预热阶段完成，用时: {warmup_time:.2f} 秒")
        else:
            print("\n跳过预热阶段")
        
        # 执行贝叶斯优化
        print("\n" + "="*60)
        print("第二阶段：贝叶斯优化")
        print("="*60)
        
        if estimated_time_per_eval > 0:
            estimated_total_time = estimated_time_per_eval * max_iterations
            if estimated_total_time < 60:
                print(f"预计优化时间: ~{estimated_total_time:.1f} 秒")
            elif estimated_total_time < 3600:
                print(f"预计优化时间: ~{estimated_total_time/60:.1f} 分钟")
            else:
                print(f"预计优化时间: ~{estimated_total_time/3600:.1f} 小时")
        
        # 设置优化空间
        dimensions = [Real(low, high, name=name) for (low, high), name in 
                     zip(param_bounds, self.param_symbols)]
        
        # 选择采集函数
        acq_func_map = {
            'EI': 'EI',        # Expected Improvement
            'LCB': 'LCB',      # Lower Confidence Bound  
            'PI': 'PI'         # Probability of Improvement
        }
        acq_func = acq_func_map.get(acquisition_function, 'EI')
        
        optimization_start = time.time()
        
        try:
            # 执行贝叶斯优化
            result = gp_minimize(
                func=objective_function,
                dimensions=dimensions,
                n_calls=max_iterations,
                n_initial_points=initial_samples,
                acq_func=acq_func,
                random_state=random_state,
                verbose=self.verbose
            )
            
            optimization_time = time.time() - optimization_start
            total_time = warmup_time + optimization_time
            
            print("\n" + "="*80)
            print("贝叶斯优化完成!")
            print("="*80)
            print(f"预热时间: {warmup_time:.2f} 秒")
            print(f"优化时间: {optimization_time:.2f} 秒")
            print(f"总用时: {total_time:.2f} 秒")
            
            # 提取最优结果
            best_params = result.x
            best_value = result.fun
            
            print("\n最优参数:")
            for i, (name, value) in enumerate(zip(self.param_names, best_params)):
                print(f"  {name} = {value:.6f}")
            print(f"\n最优目标函数值: {best_value:.4f} MeV²")
            
            # 构建结果字典
            optimization_result = {
                'best_params': best_params,
                'best_value': best_value,
                'optimization_result': result,
                'warmup_results': getattr(self, 'warmup_results', []),
                'warmup_time': warmup_time,
                'optimization_time': optimization_time,
                'total_time': total_time,
                'estimated_time_per_eval': estimated_time_per_eval,
                'config': {
                    'kappa_pairs': kappa_pairs,
                    'mu_B_values': mu_B_values,
                    'T_min': T_min,
                    'T_max': T_max,
                    'param_bounds': param_bounds,
                    'max_iterations': max_iterations,
                    'initial_samples': initial_samples,
                    'T_step_scan': T_step_scan,
                    'acquisition_function': acquisition_function,
                    'warmup_samples': warmup_samples
                }
            }
            
            # 保存结果
            if output_file:
                self.save_optimization_results(optimization_result, output_file)
            
            return optimization_result
            
        except Exception as e:
            print(f"\n❌ 贝叶斯优化失败: {e}")
            return None

    def save_optimization_results(self, result, output_file):
        """
        保存优化结果到CSV文件
        
        参数:
        - result: 优化结果字典
        - output_file: 输出文件路径
        """
        
        # 如果没有指定绝对路径，默认保存到output/Gas_Liquid目录
        if not os.path.isabs(output_file):
            project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
            output_dir = os.path.join(project_root, "output", "Gas_Liquid")
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, output_file)
        
        print(f"\n保存优化结果到: {output_file}")
        
        # 确保输出目录存在
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            # 创建结果DataFrame
            optimization_data = []
            
            skopt_result = result['optimization_result']
            for i, (params, value) in enumerate(zip(skopt_result.x_iters, skopt_result.func_vals)):
                row = {
                    'iteration': i + 1,
                    'rho0': params[0],
                    'B_A': params[1], 
                    'K': params[2],
                    'm_ratio': params[3],
                    'E_sym': params[4],
                    'objective_value': value,
                    'is_best': (i == np.argmin(skopt_result.func_vals))
                }
                optimization_data.append(row)
            
            df = pd.DataFrame(optimization_data)
            
            # 保存到CSV
            df.to_csv(output_file, index=False)
            print(f"✅ 优化历史已保存到: {output_file}")
            
            # 保存配置信息到单独文件
            config_file = output_path.with_suffix('.config.json')
            import json
            with open(config_file, 'w', encoding='utf-8') as f:
                # 处理不可序列化的对象
                config_copy = result['config'].copy()
                config_copy['optimization_summary'] = {
                    'best_params': result['best_params'],
                    'best_value': result['best_value'],
                    'total_time': result['total_time'],
                    'warmup_time': result['warmup_time'],
                    'optimization_time': result['optimization_time']
                }
                json.dump(config_copy, f, indent=2, ensure_ascii=False)
            print(f"✅ 配置信息已保存到: {config_file}")
            
        except Exception as e:
            print(f"❌ 保存结果失败: {e}")

    def plot_optimization_history(self, result, save_path=None):
        """
        绘制优化历史图
        
        参数:
        - result: 优化结果字典
        - save_path: 图片保存路径，如果为None则保存到默认路径
        """
        
        # 设置中文字体
        try:
            plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
            plt.rcParams['axes.unicode_minus'] = False
        except:
            pass  # 如果字体设置失败，使用默认字体
        
        skopt_result = result['optimization_result']
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle('PNJL模型参数贝叶斯优化结果', fontsize=16)
        
        # 优化历史
        axes[0, 0].plot(skopt_result.func_vals, 'b-o', markersize=4)
        axes[0, 0].axhline(y=result['best_value'], color='r', linestyle='--', alpha=0.7)
        axes[0, 0].set_xlabel('迭代次数')
        axes[0, 0].set_ylabel('目标函数值')
        axes[0, 0].set_title('优化收敛历史')
        axes[0, 0].grid(True, alpha=0.3)
        
        # 参数演化
        param_history = np.array(skopt_result.x_iters)
        
        for i, name in enumerate(self.param_symbols[:5]):
            row, col = divmod(i+1, 3)
            if row < 2 and col < 3:
                axes[row, col].plot(param_history[:, i], 'g-o', markersize=3)
                axes[row, col].axhline(y=result['best_params'][i], color='r', linestyle='--', alpha=0.7)
                axes[row, col].set_xlabel('迭代次数')
                axes[row, col].set_ylabel(name)
                axes[row, col].set_title(f'{self.param_names[i]} 演化')
                axes[row, col].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # 如果没有指定保存路径，使用默认路径
        if save_path is None:
            project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
            output_dir = os.path.join(project_root, "output", "Gas_Liquid")
            os.makedirs(output_dir, exist_ok=True)
            save_path = os.path.join(output_dir, "optimization_history.png")
        
        # 确保输出目录存在
        save_dir = os.path.dirname(save_path)
        os.makedirs(save_dir, exist_ok=True)
        
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"✅ 优化历史图已保存到: {save_path}")
        
        plt.close()  # 关闭图形以节省内存

def demo_bayesian_optimization_with_warmup(optimizer='auto'):
    """
    演示：PNJL模型参数贝叶斯优化 (带预热)
    与Julia版本的demo_bayesian_optimization_with_warmup()功能相同
    
    参数:
    - optimizer: 优化库选择 ('auto', 'skopt', 'optuna', 'hyperopt', 'bayes_opt')
                'auto' 自动选择最佳可用库
    """
    
    print("="*100)
    print("演示：PNJL模型参数贝叶斯优化 (带预热) - Python版本")
    print("="*100)
    
    # 自动选择最佳优化库
    if optimizer == 'auto':
        # 优先级：Optuna > skopt > hyperopt > bayes_opt
        if optuna_available:
            optimizer = 'optuna'
            print("🚀 自动选择: Optuna (高性能，现代化)")
        elif skopt_available:
            optimizer = 'skopt'
            print("🔬 自动选择: scikit-optimize (科学计算友好)")
        elif hyperopt_available:
            optimizer = 'hyperopt'
            print("🏃 自动选择: Hyperopt (TPE算法)")
        elif bayesopt_available:
            optimizer = 'bayes_opt'
            print("📖 自动选择: bayes_opt (简单易用)")
        else:
            print("❌ 没有可用的优化库，请安装至少一个:")
            print("  pip install optuna          # 推荐")
            print("  pip install scikit-optimize")
            print("  pip install hyperopt")
            print("  pip install bayesian-optimization")
            return None
    
    print(f"使用优化库: {optimizer}")
    
    # 创建优化器
    try:
        pnjl_optimizer = PNJLBayesianOptimizer(verbose=True)
    except Exception as e:
        print(f"❌ 优化器初始化失败: {e}")
        return None
    
    # 使用与Julia版本相同的实验数据
    kappa_pairs = [
        (1.09031788496341, -0.28904867673079),   # 第1组
        (1.06152332992368, 0.164279260625683),   # 第2组  
        (1.11111023684003, 0.224522832511389)    # 第3组
    ]
    
    mu_B_values = [
        632.0,   # 第1组对应632 MeV
        666.0,   # 第2组对应666 MeV
        697.0    # 第3组对应697 MeV
    ]
    
    T_min, T_max = 70.0, 120.0  # MeV
    
    # 设置参数边界（与Julia版本相同）
    param_bounds = [
        (0.145, 0.170),    # ρ₀ (fm⁻³)
        (-17.0, -15.6),    # B_A (MeV)
        (212.0, 401.0),    # K (MeV)
        (0.55, 0.75),      # m_ratio
        (26.1, 44.0)       # E_sym (MeV)
    ]
    
    print("演示配置:")
    print(f"  实验数据: {len(kappa_pairs)} 组κ比值对")
    print(f"  μ_B值: {mu_B_values} MeV")
    print(f"  温度范围: {T_min} - {T_max} MeV")
    print(f"  参数边界: {param_bounds}")
    
    # 执行带预热的贝叶斯优化（使用与Julia版本相似的参数）
    result = optimizer.optimize_with_warmup(
        kappa_pairs=kappa_pairs,
        mu_B_values=mu_B_values,
        T_min=T_min,
        T_max=T_max,
        param_bounds=param_bounds,
        max_iterations=20,              # 与Julia版本相同
        initial_samples=10,             # 与Julia版本相同
        T_step_scan=2.0,               # 与Julia版本相同，使用更粗的扫描步长
        acquisition_function='EI',      # Expected Improvement
        warmup_samples=1,              # 与Julia版本相同
        skip_warmup=True,              # 与Julia版本相同，跳过预热
        output_file="Rotation_PNJL/output/Gas_Liquid/demo_bayesian_optimization_warmup_python.csv",
        random_state=42
    )
    
    if result is not None:
        print("\n" + "="*80)
        print("✅ 带预热的演示成功完成!")
        print("="*80)
        print("Python版本的优势:")
        print("- 丰富的贝叶斯优化库生态")
        print("- 更好的可视化支持")
        print("- 更灵活的数据处理")
        print("- 与Julia函数的无缝集成")
        print("\n请检查输出文件了解详细结果。")
        
        # 绘制优化历史
        optimizer.plot_optimization_history(
            result, 
            save_path="Rotation_PNJL/output/Gas_Liquid/optimization_history_python.png"
        )
        
    else:
        print("\n❌ 带预热的演示失败。")
    
    return result

if __name__ == "__main__":
    # 运行演示
    result = demo_bayesian_optimization_with_warmup()

# =============================================================================
# 额外功能：从CSV继续优化、高级可视化等
# =============================================================================

def load_previous_optimization_results(csv_file):
    """
    从CSV文件加载前一次优化的结果
    
    参数:
    - csv_file: CSV结果文件路径
    
    返回:
    - (X_observed, y_observed): 观测点和对应的目标函数值
    - best_params: 最优参数
    - best_value: 最优目标函数值
    """
    
    if not Path(csv_file).exists():
        print(f"❌ 文件不存在: {csv_file}")
        return None, None, None, None
    
    try:
        df = pd.read_csv(csv_file)
        
        # 提取观测点和目标函数值
        param_cols = ['rho0', 'B_A', 'K', 'm_ratio', 'E_sym']
        X_observed = df[param_cols].values.tolist()
        y_observed = df['objective_value'].values.tolist()
        
        # 找到最优点
        best_idx = df['objective_value'].idxmin()
        best_params = df.loc[best_idx, param_cols].values.tolist()
        best_value = df.loc[best_idx, 'objective_value']
        
        print(f"✅ 成功加载前次优化结果:")
        print(f"  观测点数: {len(X_observed)}")
        print(f"  最优值: {best_value:.4f}")
        print(f"  最优参数: {[round(p, 4) for p in best_params]}")
        
        return X_observed, y_observed, best_params, best_value
        
    except Exception as e:
        print(f"❌ 加载文件失败: {e}")
        return None, None, None, None

def continue_optimization_from_csv(csv_file, kappa_pairs, mu_B_values, T_min, T_max, param_bounds,
                                 additional_iterations=25, T_step_scan=1.0,
                                 acquisition_function='LCB', output_file=None):
    """
    从CSV文件继续贝叶斯优化
    
    参数:
    - csv_file: 前一次优化结果的CSV文件
    - additional_iterations: 额外的迭代次数
    - 其他参数与原优化函数相同
    """
    
    print("="*100)
    print("从CSV文件继续贝叶斯优化")
    print("="*100)
    print(f"开始时间: {datetime.now()}")
    print(f"前次结果文件: {csv_file}")
    print(f"额外迭代: {additional_iterations}")
    
    # 加载先验数据
    X_prior, y_prior, best_params, best_value = load_previous_optimization_results(csv_file)
    
    if X_prior is None:
        print("❌ 无法加载先验数据，启动全新优化...")
        optimizer = PNJLBayesianOptimizer(verbose=True)
        return optimizer.optimize_with_warmup(
            kappa_pairs, mu_B_values, T_min, T_max, param_bounds,
            max_iterations=additional_iterations, output_file=output_file)
    
    # 创建优化器
    optimizer = PNJLBayesianOptimizer(verbose=True)
    
    # 创建目标函数
    objective_function = optimizer.create_objective_function(
        kappa_pairs, mu_B_values, T_min, T_max, T_step_scan)
    
    # 设置优化空间
    dimensions = [Real(low, high, name=name) for (low, high), name in 
                 zip(param_bounds, optimizer.param_symbols)]
    
    print(f"\n基于 {len(X_prior)} 个历史数据点继续优化")
    print(f"当前最优值: {best_value:.4f}")
    
    start_time = time.time()
    
    try:
        # 使用先验数据进行贝叶斯优化
        result = gp_minimize(
            func=objective_function,
            dimensions=dimensions,
            n_calls=additional_iterations,
            x0=X_prior,  # 先验观测点
            y0=y_prior,  # 先验目标函数值
            acq_func=acquisition_function,
            random_state=42,
            verbose=True
        )
        
        optimization_time = time.time() - start_time
        
        print("\n" + "="*80)
        print("继续优化完成!")
        print("="*80)
        print(f"优化时间: {optimization_time:.2f} 秒")
        
        # 提取最优结果
        new_best_params = result.x
        new_best_value = result.fun
        
        print(f"\n新的最优参数:")
        for i, (name, value) in enumerate(zip(optimizer.param_names, new_best_params)):
            print(f"  {name} = {value:.6f}")
        print(f"\n新的最优目标函数值: {new_best_value:.4f} MeV²")
        
        if new_best_value < best_value:
            improvement = best_value - new_best_value
            print(f"🎉 改进了 {improvement:.4f} MeV²!")
        else:
            print("⚠️ 未发现更好的解")
        
        # 构建结果字典
        optimization_result = {
            'best_params': new_best_params,
            'best_value': new_best_value,
            'optimization_result': result,
            'previous_best_value': best_value,
            'optimization_time': optimization_time,
            'prior_data_points': len(X_prior)
        }
        
        # 保存结果
        if output_file:
            optimizer.save_optimization_results(optimization_result, output_file)
        
        return optimization_result
        
    except Exception as e:
        print(f"\n❌ 继续优化失败: {e}")
        return None

def compare_optimization_methods():
    """
    比较不同采集函数的优化效果
    """
    
    print("="*100)
    print("比较不同贝叶斯优化采集函数")
    print("="*100)
    
    # 实验数据（使用较少数据进行快速测试）
    kappa_pairs = [
        (1.09031788496341, -0.28904867673079),   # 第1组
        (1.06152332992368, 0.164279260625683)    # 第2组
    ]
    
    mu_B_values = [632.0, 666.0]  # MeV
    T_min, T_max = 80.0, 120.0     # 缩小温度范围
    
    # 参数边界
    param_bounds = [
        (0.14, 0.16),      # ρ₀ (缩小范围)
        (-17.0, -15.5),    # B_A
        (230.0, 250.0),    # K (缩小范围)
        (0.65, 0.75),      # m_ratio
        (30.0, 34.0)       # E_sym (缩小范围)
    ]
    
    # 测试不同的采集函数
    acquisition_functions = ['EI', 'LCB', 'PI']
    results = {}
    
    for acq_func in acquisition_functions:
        print(f"\n{'='*60}")
        print(f"测试采集函数: {acq_func}")
        print(f"{'='*60}")
        
        optimizer = PNJLBayesianOptimizer(verbose=False)
        
        result = optimizer.optimize_with_warmup(
            kappa_pairs=kappa_pairs,
            mu_B_values=mu_B_values,
            T_min=T_min,
            T_max=T_max,
            param_bounds=param_bounds,
            max_iterations=15,
            initial_samples=5,
            T_step_scan=2.0,
            acquisition_function=acq_func,
            skip_warmup=True,
            output_file=f"Rotation_PNJL/output/Gas_Liquid/compare_{acq_func.lower()}.csv"
        )
        
        if result:
            results[acq_func] = {
                'best_value': result['best_value'],
                'best_params': result['best_params'],
                'total_time': result['total_time']
            }
            print(f"✅ {acq_func}: {result['best_value']:.4f} (用时: {result['total_time']:.1f}s)")
        else:
            print(f"❌ {acq_func}: 优化失败")
    
    # 比较结果
    if results:
        print(f"\n{'='*80}")
        print("采集函数比较结果")
        print(f"{'='*80}")
        
        best_method = min(results.keys(), key=lambda k: results[k]['best_value'])
        
        for method, result in results.items():
            marker = "🏆" if method == best_method else "  "
            print(f"{marker} {method}: {result['best_value']:.4f} MeV² (用时: {result['total_time']:.1f}s)")
        
        print(f"\n🏆 最佳采集函数: {best_method}")
    
    return results

def advanced_visualization(optimization_result, save_dir="Rotation_PNJL/output/Gas_Liquid/"):
    """
    高级可视化功能
    
    参数:
    - optimization_result: 优化结果字典
    - save_dir: 图片保存目录
    """
    
    save_path = Path(save_dir)
    save_path.mkdir(parents=True, exist_ok=True)
    
    skopt_result = optimization_result['optimization_result']
    
    # 1. 优化收敛分析
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('PNJL贝叶斯优化详细分析', fontsize=16)
    
    # 收敛历史
    axes[0, 0].plot(skopt_result.func_vals, 'b-o', markersize=4, alpha=0.7)
    
    # 计算累积最优值
    cumulative_best = []
    current_best = float('inf')
    for val in skopt_result.func_vals:
        if val < current_best:
            current_best = val
        cumulative_best.append(current_best)
    
    axes[0, 0].plot(cumulative_best, 'r-', linewidth=2, label='累积最优')
    axes[0, 0].set_xlabel('迭代次数')
    axes[0, 0].set_ylabel('目标函数值')
    axes[0, 0].set_title('优化收敛历史')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 改进分析
    improvements = []
    for i in range(1, len(cumulative_best)):
        if cumulative_best[i] < cumulative_best[i-1]:
            improvements.append(i)
    
    axes[0, 1].bar(range(len(skopt_result.func_vals)), 
                   [1 if i in improvements else 0 for i in range(len(skopt_result.func_vals))],
                   alpha=0.7, color='green')
    axes[0, 1].set_xlabel('迭代次数')
    axes[0, 1].set_ylabel('是否改进')
    axes[0, 1].set_title(f'改进次数: {len(improvements)}')
    axes[0, 1].grid(True, alpha=0.3)
    
    # 参数相关性分析
    param_history = np.array(skopt_result.x_iters)
    
    # 选择最重要的两个参数进行2D可视化
    param_names = ["ρ₀", "B_A", "K", "m_ratio", "E_sym"]
    
    # 散点图：参数1 vs 参数2，颜色表示目标函数值
    scatter = axes[1, 0].scatter(param_history[:, 0], param_history[:, 1], 
                               c=skopt_result.func_vals, cmap='viridis_r', alpha=0.7)
    axes[1, 0].set_xlabel(param_names[0])
    axes[1, 0].set_ylabel(param_names[1])
    axes[1, 0].set_title('参数空间探索')
    plt.colorbar(scatter, ax=axes[1, 0], label='目标函数值')
    
    # 参数方差分析
    param_stds = np.std(param_history, axis=0)
    axes[1, 1].bar(param_names, param_stds, alpha=0.7, color='orange')
    axes[1, 1].set_ylabel('标准差')
    axes[1, 1].set_title('参数探索广度')
    axes[1, 1].tick_params(axis='x', rotation=45)
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path / 'advanced_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 2. 参数演化详细图
    fig, axes = plt.subplots(3, 2, figsize=(12, 15))
    fig.suptitle('PNJL模型参数演化分析', fontsize=16)
    
    for i, name in enumerate(param_names):
        row, col = divmod(i, 2)
        if row < 3:
            axes[row, col].plot(param_history[:, i], 'b-o', markersize=3, alpha=0.7)
            axes[row, col].axhline(y=optimization_result['best_params'][i], 
                                 color='r', linestyle='--', alpha=0.8, label='最优值')
            axes[row, col].set_xlabel('迭代次数')
            axes[row, col].set_ylabel(name)
            axes[row, col].set_title(f'{name} 演化历史')
            axes[row, col].legend()
            axes[row, col].grid(True, alpha=0.3)
    
    # 删除多余的子图
    if len(param_names) % 2 == 1:
        fig.delaxes(axes[2, 1])
    
    plt.tight_layout()
    plt.savefig(save_path / 'parameter_evolution.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"✅ 高级可视化图已保存到: {save_path}")

def quick_optimization_test():
    """
    快速优化测试 - 用于验证功能
    """
    
    print("="*80)
    print("快速优化测试")
    print("="*80)
    
    # 最小化配置
    kappa_pairs = [(1.09031788496341, -0.28904867673079)]
    mu_B_values = [632.0]
    T_min, T_max = 80.0, 120.0
    
    param_bounds = [
        (0.14, 0.16),     # ρ₀ 
        (-17.0, -15.5),   # B_A
        (235.0, 245.0),   # K
        (0.68, 0.72),     # m_ratio
        (31.0, 33.0)      # E_sym
    ]
    
    optimizer = PNJLBayesianOptimizer(verbose=True)
    
    result = optimizer.optimize_with_warmup(
        kappa_pairs=kappa_pairs,
        mu_B_values=mu_B_values,
        T_min=T_min,
        T_max=T_max,
        param_bounds=param_bounds,
        max_iterations=10,
        initial_samples=5,
        T_step_scan=3.0,
        skip_warmup=True,
        output_file="Rotation_PNJL/output/Gas_Liquid/quick_test.csv"
    )
    
    if result:
        print("✅ 快速测试成功!")
        print("这验证了Python-Julia接口和贝叶斯优化流程都工作正常。")
        
        # 简单可视化
        advanced_visualization(result)
    else:
        print("❌ 快速测试失败")
    
    return result

# =============================================================================
# 使用说明和示例
# =============================================================================

def print_usage_guide():
    """
    打印使用说明
    """
    
    print("="*100)
    print("PNJL模型Python贝叶斯优化使用指南")
    print("="*100)
    
    print("""
🎯 主要功能

1. 基础优化：
   optimizer = PNJLBayesianOptimizer()
   result = optimizer.optimize_with_warmup(kappa_pairs, mu_B_values, T_min, T_max, bounds)

2. 从CSV继续优化：
   result = continue_optimization_from_csv(csv_file, kappa_pairs, mu_B_values, T_min, T_max, bounds)

3. 方法比较：
   results = compare_optimization_methods()

4. 高级可视化：
   advanced_visualization(result)

📊 参数说明

- kappa_pairs: κ比值对 [(κ₃/κ₁, κ₄/κ₂), ...]
- mu_B_values: 重子化学势 [μ_B1, μ_B2, ...] (MeV单位)
- T_min, T_max: 温度范围 (MeV单位)
- param_bounds: 参数边界 [(min1,max1), ...] 对应 [ρ₀, B_A, K, m_ratio, E_sym]

🚀 快速开始

1. 运行演示: demo_bayesian_optimization_with_warmup()
2. 快速测试: quick_optimization_test()
3. 查看使用指南: print_usage_guide()

📈 与Julia版本的对比优势

✅ Python优势:
- 更丰富的优化库生态 (scikit-optimize, optuna等)
- 更好的数据可视化 (matplotlib, seaborn等)
- 更强的数据处理能力 (pandas, numpy等)
- 更容易的结果分析和后处理

✅ Julia优势:
- 高性能的数值计算
- 专业的物理建模能力
- JIT编译优化

🔧 依赖安装

pip install scikit-optimize julia numpy pandas matplotlib

然后设置PyJulia:
python -c "import julia; julia.install()"

""")

# 主程序入口扩展
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='PNJL模型贝叶斯优化')
    parser.add_argument('--mode', choices=['demo', 'test', 'compare', 'help'], 
                       default='demo', help='运行模式')
    parser.add_argument('--csv', type=str, help='从CSV文件继续优化')
    
    args = parser.parse_args()
    
    if args.mode == 'demo':
        print("🚀 运行演示...")
        result = demo_bayesian_optimization_with_warmup()
    elif args.mode == 'test':
        print("🧪 运行快速测试...")
        result = quick_optimization_test()
    elif args.mode == 'compare':
        print("📊 比较不同方法...")
        result = compare_optimization_methods()
    elif args.mode == 'help':
        print_usage_guide()
    
    if args.csv:
        print(f"📁 从CSV文件继续优化: {args.csv}")
        # 这里需要用户提供其他参数，暂时使用默认值
        kappa_pairs = [(1.09031788496341, -0.28904867673079), (1.06152332992368, 0.164279260625683)]
        mu_B_values = [632.0, 666.0]
        T_min, T_max = 70.0, 120.0
        param_bounds = [(0.145, 0.170), (-17.0, -15.6), (212.0, 401.0), (0.55, 0.75), (26.1, 44.0)]
        
        result = continue_optimization_from_csv(
            args.csv, kappa_pairs, mu_B_values, T_min, T_max, param_bounds,
            output_file="Rotation_PNJL/output/Gas_Liquid/continued_optimization.csv"
        )