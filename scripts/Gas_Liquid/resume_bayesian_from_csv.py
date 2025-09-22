#!/usr/bin/env python3
"""
NOTICE: This script has been moved to examples/Gas_Liquid/resume_bayesian_from_csv.py

Please use the copy in the `examples/Gas_Liquid` directory.
"""

import sys

print('This script has been moved to examples/Gas_Liquid/resume_bayesian_from_csv.py')
sys.exit(0)
#!/usr/bin/env python3
"""
从已有的 CSV （local_optimum_scan.csv）加载评估点并继续执行基于 gp_minimize 的贝叶斯优化。

用法：调整参数后直接运行脚本。脚本会：
- 加载 CSV 的 (x, y) 数据作为已有评估点。
- 使用 skopt.gp_minimize 并传入 x0, y0 继续优化 n_extra_calls 次。
- 合并旧/新结果并保存 CSV 和收敛图。
"""
import os
import sys
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 配置
# project_root should be the repository root (Rotation_PNJL). __file__ is scripts/Gas_Liquid/...
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
csv_path = os.path.join(project_root, 'scripts', 'output', 'Gas_Liquid', 'local_optimum_scan.csv')
output_dir = os.path.join(project_root, 'output', 'Gas_Liquid')
os.makedirs(output_dir, exist_ok=True)

# 优化设置
n_extra_calls = 100  # 额外的评估次数
n_initial_random = 0  # gp_minimize 的 n_initial_points，当提供 x0 时可以为0

# 参数边界（与 demo 保持一致）
param_bounds = [
    (0.145, 0.170),    # rho0
    (-17.0, -15.6),    # B_A
    (212.0, 401.0),    # K
    (0.55, 0.75),      # m_ratio
    (26.1, 44.0)       # E_sym
]

# 读取 CSV
if not os.path.exists(csv_path):
    print('CSV 文件未找到:', csv_path)
    sys.exit(1)

df = pd.read_csv(csv_path)
# CSV 列在 check_local_optimum.py 中是：rho0,B_A,K,m_ratio,E_sym,objective_value
expected_cols = {'rho0','B_A','K','m_ratio','E_sym','objective_value'}
if not expected_cols.issubset(set(df.columns)):
    print('CSV 列不包含预期列，看到的列:', df.columns.tolist())
    sys.exit(1)

# 准备 x0, y0
full_x = df[['rho0','B_A','K','m_ratio','E_sym']].values
full_y = df['objective_value'].values

# 过滤不在 bounds 内的点
def in_bounds(x):
    for xi, (low, high) in zip(x, param_bounds):
        if not (low <= xi <= high):
            return False
    return True

filtered_x = []
filtered_y = []
for xi, yi in zip(full_x, full_y):
    if in_bounds(xi):
        filtered_x.append(list(xi))
        filtered_y.append(float(yi))

print(f'加载 {len(full_x)} 条历史评估点，{len(filtered_x)} 条点位于参数边界内并将被用于继续优化')

# 如果没有可用的历史点，则正常开始
if len(filtered_x) == 0:
    x0 = None
    y0 = None
else:
    x0 = filtered_x
    y0 = filtered_y

# Julia 初始化 & 闭包创建（复用 demo 脚本逻辑）
from julia import Main as jl
hc = 197.327
# Try to load experimental data from data/raw/Gas_Liquid; fall back to hard-coded constants
data_dir = os.path.join(project_root, 'data', 'raw', 'Gas_Liquid')

# Default (hard-coded) values
default_kappa_pairs = [
    (1.09031788496341, -0.28904867673079),
    (1.06152332992368, 0.164279260625683),
    (1.11111023684003, 0.224522832511389)
]
default_mu_B_values = [632.0, 666.0, 697.0]

def try_load_combined_muB_kappa(path):
    csv_file = os.path.join(path, 'muB_kappa_pairs.csv')
    if not os.path.exists(csv_file):
        return None
    try:
        df = pd.read_csv(csv_file)
        # Expect columns: mu_B,kappa3_1,kappa4_2
        if {'mu_B', 'kappa3_1', 'kappa4_2'}.issubset(set(df.columns)):
            mu_vals = list(df['mu_B'].astype(float).values)
            kappa_pairs = list(df[['kappa3_1', 'kappa4_2']].itertuples(index=False, name=None))
            return mu_vals, kappa_pairs
        # try to infer: first numeric is mu_B, next two numeric are kappas
        numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        if len(numeric_cols) >= 3:
            mu_vals = list(df[numeric_cols[0]].astype(float).values)
            kappa_pairs = list(df[numeric_cols[1:3]].itertuples(index=False, name=None))
            return mu_vals, kappa_pairs
    except Exception as e:
        print('Warning: 无法解析 muB_kappa_pairs.csv:', e)
    return None

def try_load_kappa_pairs(path):
    # keep previous behavior for backward compatibility
    csv_file = os.path.join(path, 'kappa_pairs.csv')
    if not os.path.exists(csv_file):
        return None
    try:
        df_k = pd.read_csv(csv_file)
        if {'k1', 'k2'}.issubset(set(df_k.columns)):
            return list(df_k[['k1', 'k2']].itertuples(index=False, name=None))
        numeric_cols = [c for c in df_k.columns if pd.api.types.is_numeric_dtype(df_k[c])]
        if len(numeric_cols) >= 2:
            return list(df_k[numeric_cols[:2]].itertuples(index=False, name=None))
    except Exception as e:
        print('Warning: 无法解析 kappa_pairs.csv:', e)
    return None

def try_load_mu_B_values(path):
    csv_file = os.path.join(path, 'mu_B_values.csv')
    if not os.path.exists(csv_file):
        return None
    try:
        df_m = pd.read_csv(csv_file)
        if 'mu_B' in df_m.columns:
            return list(df_m['mu_B'].astype(float).values)
        numeric_cols = [c for c in df_m.columns if pd.api.types.is_numeric_dtype(df_m[c])]
        if len(numeric_cols) >= 1:
            return list(df_m[numeric_cols[0]].astype(float).values)
    except Exception as e:
        print('Warning: 无法解析 mu_B_values.csv:', e)
    return None

# Load combined first; fall back to separate files; otherwise use defaults
combined = try_load_combined_muB_kappa(data_dir)
if combined is not None:
    mu_B_values, kappa_pairs = combined
else:
    # try separate files
    kappa_pairs = try_load_kappa_pairs(data_dir) or default_kappa_pairs
    mu_B_values = try_load_mu_B_values(data_dir) or default_mu_B_values

print(f'Using kappa_pairs (count={len(kappa_pairs)}):', kappa_pairs)
print(f'Using mu_B_values (count={len(mu_B_values)}):', mu_B_values)
T_min, T_max = 70.0, 120.0
T_step_scan = 2.0
mu_B_julia = [mu / hc for mu in mu_B_values]
T_min_julia = T_min / hc
T_max_julia = T_max / hc
T_step_julia = T_step_scan / hc

# ensure Julia file included
import os
original_cwd = os.getcwd()
script_path = os.path.dirname(os.path.abspath(__file__))
julia_project_path = os.path.dirname(os.path.dirname(script_path))
os.chdir(julia_project_path)
try:
    # Allow skipping Julia initialization for testing data reading only
    if os.environ.get('SKIP_JULIA', '0') != '1':
        jl.eval('using Pkg; Pkg.activate(".")')
        jl.eval('using Pkg')
        jl.eval('Pkg.instantiate()')
        jl.eval('include("src/Gas_Liquid/Advanced_FindTforDiff.jl")')
        print('✅ 已初始化Julia项目并include文件')
    else:
        print('ℹ️ SKIP_JULIA=1: 已跳过 Julia 初始化（仅测试数据读取）')
finally:
    os.chdir(original_cwd)

# create objective closure once
if os.environ.get('SKIP_JULIA', '0') != '1':
    jl_base_objective = jl.create_temperature_difference_objective(
        kappa_pairs, mu_B_julia, T_min_julia, T_max_julia,
        T_step_scan=T_step_julia, penalty_for_missing=1e6, verbose=False)
else:
    jl_base_objective = None

# wrapper for skopt
def objective_wrapper(x):
    if jl_base_objective is None:
        raise RuntimeError('Julia objective is not initialized (SKIP_JULIA=1).')
    return float(jl_base_objective(tuple(x)))

# run gp_minimize with initial data
from skopt import gp_minimize
from skopt.space import Real

space = [Real(low, high) for (low, high) in param_bounds]

print('开始继续贝叶斯优化: 额外调用 =', n_extra_calls)
start_time = time.time()
res = gp_minimize(func=objective_wrapper,
                  dimensions=space,
                  n_calls=n_extra_calls,
                  x0=x0,
                  y0=y0,
                  n_initial_points=n_initial_random,
                  random_state=42,
                  verbose=True)
end_time = time.time()
print('继续优化完成，用时:', end_time-start_time)

# 合并历史与新结果
# res.x_iters 与 res.func_vals 包含从 x0/y0 开始的完整列表（如果x0/y0提供，会在前面）
all_x = res.x_iters
all_y = res.func_vals

# 保存合并结果
out_csv = os.path.join(output_dir, 'resume_optimization_result.csv')
rows = []
for i, (x, y) in enumerate(zip(all_x, all_y), 1):
    rows.append({'iteration': i,
                 'rho0': x[0], 'B_A': x[1], 'K': x[2], 'm_ratio': x[3], 'E_sym': x[4],
                 'objective_value': y})
pd.DataFrame(rows).to_csv(out_csv, index=False)
print('合并结果已保存到', out_csv)

# 画收敛图
plt.figure(figsize=(10,5))
plt.subplot(1,2,1)
plt.plot(all_y, 'b-o', markersize=3)
plt.xlabel('Iteration')
plt.ylabel('Objective')
plt.grid(True)

cumulative = []
best = float('inf')
for v in all_y:
    best = min(best, v)
    cumulative.append(best)

plt.subplot(1,2,2)
plt.plot(cumulative, 'g-')
plt.xlabel('Iteration')
plt.ylabel('Cumulative Best')
plt.grid(True)
plt.tight_layout()
plot_path = os.path.join(output_dir, 'resume_optimization_convergence.png')
plt.savefig(plot_path, dpi=150)
print('收敛图已保存到', plot_path)

print('完成')
