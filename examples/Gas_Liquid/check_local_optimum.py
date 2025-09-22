#!/usr/bin/env python3
"""
在给定最优点附近均匀采样并评估目标函数，判断是否存在更好的点。
"""
import os
import sys
import time
import numpy as np
import pandas as pd

# 参数：最优点（从你提供的结果）
best_params = [0.17, -16.613855, 401.0, 0.75, 26.1]
best_value_reported = 270.2564

# 采样设置：每个参数在中心点附近取 n_steps 点（包含中心），使用相对或绝对步长
n_steps = 3  # 每个参数沿轴的点数（奇数以保证中心包含），3^5=243 点，外加随机点，实用运行时间
rel_frac = [0.01, 0.01, 0.005, 0.01, 0.01]  # 每个参数的相对变化幅度（fraction of param）

# Julia 闭包创建设置（与 demo 脚本一致）
from julia import Main as jl

# 初始化Julia环境（与 demo 脚本相同），切换到项目根并 include 目标文件
import os
original_cwd = os.getcwd()
script_path = os.path.dirname(os.path.abspath(__file__))
julia_project_path = os.path.dirname(os.path.dirname(script_path))
os.chdir(julia_project_path)
try:
    jl.eval('using Pkg; Pkg.activate(".")')
    jl.eval('using Pkg')
    jl.eval('Pkg.instantiate()')
    jl.eval('include("src/Gas_Liquid/Advanced_FindTforDiff.jl")')
    print('✅ 已初始化Julia项目并include文件')
except Exception as e:
    print('⚠️ 初始化Julia项目或include文件失败:', e)
    os.chdir(original_cwd)
    raise
finally:
    # 恢复原始cwd（create_temperature_difference_objective 已被包含到 Main）
    os.chdir(original_cwd)

# 实验数据（与 demo 脚本匹配）
hc = 197.327
kappa_pairs = [
    (1.09031788496341, -0.28904867673079),
    (1.06152332992368, 0.164279260625683),
    (1.11111023684003, 0.224522832511389)
]
mu_B_values = [632.0, 666.0, 697.0]
T_min, T_max = 70.0, 120.0
T_step_scan = 2.0

mu_B_julia = [mu / hc for mu in mu_B_values]
T_min_julia = T_min / hc
T_max_julia = T_max / hc
T_step_julia = T_step_scan / hc

# 创建Julia目标闭包（一次）
try:
    jl_base_objective = jl.create_temperature_difference_objective(
        kappa_pairs, mu_B_julia, T_min_julia, T_max_julia,
        T_step_scan=T_step_julia, penalty_for_missing=1e6, verbose=False)
except Exception as e:
    print('无法在Julia端创建闭包:', e)
    sys.exit(1)

# 生成网格采样点（每个参数在中心附近线性取值）
param_grid = []
for i, p in enumerate(best_params):
    frac = rel_frac[i]
    if p != 0:
        step = abs(p) * frac
    else:
        step = frac
    vals = np.linspace(p - step*(n_steps//2), p + step*(n_steps//2), n_steps)
    param_grid.append(vals)

# 构建笛卡尔积（如果全网格太大，则限制）
from itertools import product
points = list(product(*param_grid))

# 额外加些随机扰动点
rng = np.random.default_rng(42)
n_random = 50
for _ in range(n_random := n_random if (n_random:=n_random) else 0):
    pass
# 直接用固定数量
for _ in range(n_random):
    perturb = [best_params[i] + rng.normal(0, abs(best_params[i])*rel_frac[i]/2 if best_params[i]!=0 else rel_frac[i]/2) for i in range(len(best_params))]
    points.append(tuple(perturb))

print(f"将评估 {len(points)} 个点（包括网格与随机点）...")

results = []
start = time.time()
for idx, pt in enumerate(points, 1):
    t0 = time.time()
    val = jl_base_objective(tuple(pt))
    t1 = time.time()
    results.append({
        'idx': idx,
        'rho0': pt[0],
        'B_A': pt[1],
        'K': pt[2],
        'm_ratio': pt[3],
        'E_sym': pt[4],
        'objective_value': float(val),
        'eval_time_s': t1 - t0
    })
    if idx % 10 == 0 or idx == len(points):
        print(f"  进度 {idx}/{len(points)} — 当前值 {float(val):.4f} — 用时 {t1-t0:.2f}s")

total_time = time.time() - start
print(f"全部评估完成（{len(points)} 点），总用时 {total_time:.2f}s")

# 保存为CSV
out_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'output', 'Gas_Liquid')
os.makedirs(out_dir, exist_ok=True)
out_file = os.path.join(out_dir, 'local_optimum_scan.csv')

pd.DataFrame(results).to_csv(out_file, index=False)
print('结果已保存到', out_file)

# 分析结果：是否存在比 reported 更小的点
vals = np.array([r['objective_value'] for r in results])
min_val = vals.min()
min_idx = vals.argmin()
print(f"采样最小值: {min_val:.4f}（reported {best_value_reported:.4f}）")
if min_val < best_value_reported - 1e-6:
    print('存在比报告的最优值更小的点：可能说明报告点不是全局最优或附近有更好点')
else:
    print('在采样范围内未发现更小的值：报告点在此邻域看起来为局部最优')

# 给出最小点细节
print('最小点信息:')
print(pd.DataFrame(results).iloc[min_idx])

# 任务完成
print('完成')
