#!/usr/bin/env python3
"""
测试Python贝叶斯优化功能
确保Python调用Julia函数和贝叶斯优化都能正常工作

运行方式:
cd d:/Desktop/Julia/Rotation_PNJL
python test/Gas_Liquid/test_python_bayesian.py
"""

import sys
import os
from pathlib import Path

def get_project_paths():
    """获取项目路径"""
    script_dir = Path(__file__).parent
    project_dir = script_dir.parent.parent  # 从test/Gas_Liquid回到项目根目录
    src_dir = project_dir / "src"
    gas_liquid_dir = src_dir / "Gas_Liquid" 
    output_dir = project_dir / "output" / "Gas_Liquid"
    
    # 添加src路径到Python路径
    sys.path.insert(0, str(src_dir))
    sys.path.insert(0, str(gas_liquid_dir))
    
    return project_dir, gas_liquid_dir, output_dir

# 在文件开头就设置路径
project_dir, src_dir, output_dir = get_project_paths()

# 尝试不同的导入方式
try:
    # 方式1：直接导入模块
    import Advanced_Bayesian
    PNJLBayesianOptimizer = Advanced_Bayesian.PNJLBayesianOptimizer
    ADVANCED_BAYESIAN_AVAILABLE = True
except ImportError:
    try:
        # 方式2：从包导入
        from Gas_Liquid.Advanced_Bayesian import PNJLBayesianOptimizer
        ADVANCED_BAYESIAN_AVAILABLE = True
    except ImportError:
        try:
            # 方式3：直接from导入
            from Advanced_Bayesian import PNJLBayesianOptimizer
            ADVANCED_BAYESIAN_AVAILABLE = True
        except ImportError:
            # 最后的备选方案
            PNJLBayesianOptimizer = None
            ADVANCED_BAYESIAN_AVAILABLE = False

# 在模块加载时立即确保将 src/Gas_Liquid 添加到 sys.path，以便编辑器和运行时都能解析本地模块 Advanced_Bayesian
project_dir, src_dir, output_dir = get_project_paths()

def test_imports():
    """测试所有必要的库导入"""
    print("="*60)
    print("测试库导入")
    print("="*60)
    
    try:
        # 测试基础库
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        print("✅ 基础库 (numpy, pandas, matplotlib) 导入成功")
    except ImportError as e:
        print(f"❌ 基础库导入失败: {e}")
        return False
    
    try:
        # 测试scikit-optimize
        from skopt import gp_minimize
        from skopt.space import Real
        print("✅ scikit-optimize 导入成功")
    except ImportError as e:
        print(f"❌ scikit-optimize 导入失败: {e}")
        print("请运行: pip install scikit-optimize")
        return False
    
    try:
        # 测试PyJulia
        import julia
        print("✅ PyJulia 导入成功")
    except ImportError as e:
        print(f"❌ PyJulia 导入失败: {e}")
        print("请运行: pip install julia")
        print("然后运行: python -c \"import julia; julia.install()\"")
        return False
    
    try:
        # 测试我们的模块
        if not ADVANCED_BAYESIAN_AVAILABLE:
            print("❌ PNJLBayesianOptimizer 在启动时导入失败")
            return False
        print("✅ PNJLBayesianOptimizer 导入成功")
    except Exception as e:
        print(f"❌ PNJLBayesianOptimizer 测试失败: {e}")
        return False
    
    return True

def test_julia_connection():
    """测试Julia连接和函数调用"""
    print("\n" + "="*60)
    print("测试Julia连接")
    print("="*60)
    
    try:
        project_dir, _, _ = get_project_paths()
        
        # 使用全局导入的类
        if not ADVANCED_BAYESIAN_AVAILABLE:
            print("❌ PNJLBayesianOptimizer 不可用")
            return False
        
        # 创建优化器（这会初始化Julia环境）
        optimizer = PNJLBayesianOptimizer(
            julia_project_path=str(project_dir),
            verbose=True
        )
        
        print("✅ Julia环境初始化成功")
        return True
        
    except Exception as e:
        print(f"❌ Julia连接失败: {e}")
        return False

def test_objective_function():
    """测试目标函数创建和调用"""
    print("\n" + "="*60)
    print("测试目标函数")
    print("="*60)
    
    try:
        project_dir, _, _ = get_project_paths()
        
        # 使用全局导入的类
        if not ADVANCED_BAYESIAN_AVAILABLE:
            print("❌ PNJLBayesianOptimizer 不可用")
            return False
        
        optimizer = PNJLBayesianOptimizer(
            julia_project_path=str(project_dir), 
            verbose=False
        )
        
        # 简单的测试数据
        kappa_pairs = [(1.09, -0.29)]
        mu_B_values = [632.0]
        T_min, T_max = 80.0, 120.0
        
        # 创建目标函数
        objective_func = optimizer.create_objective_function(
            kappa_pairs, mu_B_values, T_min, T_max, T_step_scan=3.0
        )
        
        # 测试目标函数调用
        test_params = [0.155, -16.2, 240.0, 0.70, 32.0]
        result = objective_func(test_params)
        
        print(f"✅ 目标函数调用成功")
        print(f"  测试参数: {test_params}")
        print(f"  函数值: {result}")
        
        return True
        
    except Exception as e:
        print(f"❌ 目标函数测试失败: {e}")
        return False

def test_quick_optimization():
    """测试快速贝叶斯优化"""
    print("\n" + "="*60)
    print("测试快速贝叶斯优化")
    print("="*60)
    
    try:
        project_dir, _, output_dir = get_project_paths()
        
        # 使用全局导入的类
        if not ADVANCED_BAYESIAN_AVAILABLE:
            print("❌ PNJLBayesianOptimizer 不可用")
            return False
        
        optimizer = PNJLBayesianOptimizer(
            julia_project_path=str(project_dir),
            verbose=False
        )
        
        # 最小化测试配置
        kappa_pairs = [(1.09, -0.29)]
        mu_B_values = [632.0]
        T_min, T_max = 80.0, 120.0
        
        param_bounds = [
            (0.14, 0.16),     # ρ₀ 
            (-17.0, -15.5),   # B_A
            (235.0, 245.0),   # K
            (0.68, 0.72),     # m_ratio
            (31.0, 33.0)      # E_sym
        ]
        
        print("运行微型优化测试 (5次迭代)...")
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        test_output_file = output_dir / "test_optimization_result.csv"
        
        result = optimizer.optimize_with_warmup(
            kappa_pairs=kappa_pairs,
            mu_B_values=mu_B_values,
            T_min=T_min,
            T_max=T_max,
            param_bounds=param_bounds,
            max_iterations=5,
            initial_samples=3,
            T_step_scan=3.0,
            skip_warmup=True,
            output_file=str(test_output_file)
        )
        
        if result:
            print("✅ 快速优化测试成功!")
            print(f"  最优值: {result['best_value']:.4f}")
            print(f"  最优参数: {[round(p, 4) for p in result['best_params']]}")
            print(f"  结果文件: {test_output_file}")
            return True
        else:
            print("❌ 快速优化返回None")
            return False
            
    except Exception as e:
        print(f"❌ 快速优化测试失败: {e}")
        return False

def test_output_directory():
    """测试输出目录创建和文件写入"""
    print("\n" + "="*60)
    print("测试输出目录和文件操作")
    print("="*60)
    
    try:
        project_dir, _, output_dir = get_project_paths()
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        print(f"✅ 输出目录创建成功: {output_dir}")
        
        # 测试文件写入
        test_file = output_dir / "test_file.txt"
        with open(test_file, 'w', encoding='utf-8') as f:
            f.write("测试文件写入\nTest file writing\n")
        
        # 测试文件读取
        with open(test_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        print(f"✅ 文件读写测试成功")
        print(f"  测试文件: {test_file}")
        
        # 清理测试文件
        test_file.unlink()
        print("✅ 测试文件已清理")
        
        return True
        
    except Exception as e:
        print(f"❌ 输出目录测试失败: {e}")
        return False

def main():
    """主测试流程"""
    print("="*80)
    print("PNJL Python贝叶斯优化功能测试")
    print("="*80)
    
    # 检查工作目录和路径
    project_dir, src_dir, output_dir = get_project_paths()
    print(f"当前工作目录: {os.getcwd()}")
    print(f"脚本位置: {__file__}")
    print(f"项目根目录: {project_dir}")
    print(f"源代码目录: {src_dir}")
    print(f"输出目录: {output_dir}")
    
    tests = [
        ("库导入测试", test_imports),
        ("输出目录测试", test_output_directory),
        ("Julia连接测试", test_julia_connection),
        ("目标函数测试", test_objective_function),
        ("快速优化测试", test_quick_optimization)
    ]
    
    results = []
    
    for test_name, test_func in tests:
        try:
            success = test_func()
            results.append((test_name, success))
            
            if not success:
                print(f"\n⚠️ {test_name}失败，停止后续测试")
                break
                
        except Exception as e:
            print(f"\n💥 {test_name}出现异常: {e}")
            results.append((test_name, False))
            break
    
    # 汇总结果
    print("\n" + "="*80)
    print("测试结果汇总")
    print("="*80)
    
    for test_name, success in results:
        status = "✅ 通过" if success else "❌ 失败"
        print(f"{test_name}: {status}")
    
    success_count = sum(success for _, success in results)
    total_count = len(results)
    
    print(f"\n总体结果: {success_count}/{total_count} 通过")
    
    if success_count == total_count:
        print("\n🎉 所有测试通过！Python贝叶斯优化功能可以正常使用。")
        print("\n接下来可以运行:")
        print("  python scripts/Gas_Liquid/demo_python_bayesian.py")
        print("  python src/Gas_Liquid/Advanced_Bayesian.py")
    else:
        print("\n⚠️ 部分测试失败，请检查依赖和环境配置。")

if __name__ == "__main__":
    main()