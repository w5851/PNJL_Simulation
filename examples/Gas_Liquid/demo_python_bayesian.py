#!/usr/bin/env python3
"""
PNJL模型Python贝叶斯优化演示脚本
不依赖复杂的模块导入，直接运行演示

使用方法:
cd d:/Desktop/Julia/Rotation_PNJL
python scripts/Gas_Liquid/demo_python_bayesian.py
"""

import matplotlib
# 设置中文字体和后端
matplotlib.use('Agg')  # 使用非交互式后端
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# 尝试设置中文字体
def setup_chinese_font():
    """设置中文字体，解决中文显示问题"""
    try:
        # 尝试常见的中文字体
        chinese_fonts = [
            'SimHei',           # 黑体 (Windows)
            'Microsoft YaHei',  # 微软雅黑 (Windows)
            'Arial Unicode MS', # Arial Unicode (macOS/Windows)
            'DejaVu Sans',      # DejaVu (Linux)
            'WenQuanYi Micro Hei' # 文泉驿微米黑 (Linux)
        ]
        
        available_fonts = [f.name for f in fm.fontManager.ttflist]
        
        for font in chinese_fonts:
            if font in available_fonts:
                plt.rcParams['font.sans-serif'] = [font]
                plt.rcParams['axes.unicode_minus'] = False
                print(f"✅ 使用字体: {font}")
                return font
        
        # 如果没有找到中文字体，使用英文标签
        print("⚠️ 未找到中文字体，将使用英文标签")
        return None
        
    except Exception as e:
        print(f"⚠️ 字体设置失败: {e}，将使用英文标签")
        return None

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
    jl.eval('using Pkg; Pkg.activate("")')
    
    # 检查文件是否存在
    julia_file_path = "src/Gas_Liquid/Advanced_FindTforDiff.jl"
    full_path = os.path.join(julia_project_path, julia_file_path)
    print(f"尝试加载Julia文件: {full_path}")
    
    if not os.path.exists(full_path):
        print(f"❌ 文件不存在: {full_path}")
        os.chdir(original_cwd)  # 恢复原始目录
        raise FileNotFoundError(full_path)
    
    jl.eval('using Pkg')
    jl.eval('Pkg.instantiate()')  # 确保所有依赖都已安装
    jl.eval(f'include("{julia_file_path}")')
    
    try:
        jl.eval('typeof(create_temperature_difference_objective)')
        print("✅ 目标函数 create_temperature_difference_objective 已加载")
    except Exception as e:
        print(f"❌ 无法找到目标函数: {e}")
        os.chdir(original_cwd)
        raise
    
    print("✅ Julia环境初始化成功")
    os.chdir(original_cwd)
    
except Exception as e:
    print(f"❌ Julia初始化失败: {e}")
    # 继续执行（脚本可回退或仅做演示）

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
    jl_base_objective = jl.create_temperature_difference_objective(
        kappa_pairs, [mu / hc for mu in mu_B_values], T_min_julia, T_max_julia,
        T_step_scan=T_step_julia, penalty_for_missing=1e6, verbose=False)
    print("✅ 已在Julia端创建目标函数闭包 base_objective（将被重用）")
except Exception as e:
    print(f"⚠️ 无法在Julia端创建目标函数闭包（将回退到eval）：{e}")

def objective_function(params):
    try:
        rho0, B_A, K, m_ratio, E_sym = params

        if jl_base_objective is not None:
            params_tuple = (rho0, B_A, K, m_ratio, E_sym)
            jl_result = jl_base_objective(params_tuple)
        else:
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

# 以下省略：优化与绘图逻辑与原始脚本一致

def main():
    print('Demo main placeholder - original script continues with optimization and plotting logic')

if __name__ == '__main__':
    main()
