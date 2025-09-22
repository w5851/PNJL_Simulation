#!/usr/bin/env python3
"""
从 CSV 恢复贝叶斯优化的示例脚本（完整版本）。
支持从 data/raw/Gas_Liquid/muB_kappa_pairs.csv 自动加载 mu_B 与 kappa 对。
在环境变量 SKIP_JULIA=1 时，会跳过 Julia 初始化，仅进行数据加载与打印（便于测试）。
"""

import os
import numpy as np
import pandas as pd
from skopt import gp_minimize
from skopt.space import Real
import json

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_DIR = os.path.join(PROJECT_ROOT, 'data', 'raw', 'Gas_Liquid')

def try_load_combined_muB_kappa(path):
	if not os.path.exists(path):
		return None
	try:
		df = pd.read_csv(path)
		expected_cols = {'mu_B', 'kappa3_1', 'kappa4_2'}
		if not expected_cols.issubset(df.columns):
			print(f"CSV 缺少列，期望列: {expected_cols}, 实际列: {df.columns}")
			return None
		mu_B = df['mu_B'].to_numpy()
		kappa_pairs = list(zip(df['kappa3_1'].to_numpy(), df['kappa4_2'].to_numpy()))
		return mu_B, kappa_pairs
	except Exception as e:
		print(f"读取合并 CSV 失败: {e}")
		return None

def try_load_kappa_pairs(path):
	if not os.path.exists(path):
		return None
	try:
		df = pd.read_csv(path)
		if not {'k1', 'k2'}.issubset(df.columns):
			return None
		return list(zip(df['k1'].to_numpy(), df['k2'].to_numpy()))
	except Exception as e:
		print(f"读取 kappa_pairs 失败: {e}")
		return None

def try_load_mu_B_values(path):
	if not os.path.exists(path):
		return None
	try:
		df = pd.read_csv(path)
		if 'mu_B' not in df.columns:
			return None
		return df['mu_B'].to_numpy()
	except Exception as e:
		print(f"读取 mu_B_values 失败: {e}")
		return None

combined_csv = os.path.join(DATA_DIR, 'muB_kappa_pairs.csv')
kappa_csv = os.path.join(DATA_DIR, 'kappa_pairs.csv')
muB_csv = os.path.join(DATA_DIR, 'mu_B_values.csv')

loaded = try_load_combined_muB_kappa(combined_csv)
if loaded is not None:
	mu_B_values, kappa_pairs = loaded
	print(f"Loaded combined CSV: {combined_csv}")
else:
	kappa_pairs = try_load_kappa_pairs(kappa_csv)
	mu_B_values = try_load_mu_B_values(muB_csv)

if kappa_pairs is None or mu_B_values is None:
	print("未找到合适的 CSV 数据，使用默认示例数据")
	kappa_pairs = [(1.09031788496341, -0.28904867673079), (1.06152332992368, 0.164279260625683), (1.11111023684003, 0.224522832511389)]
	mu_B_values = np.array([632.0, 666.0, 697.0])

print(f"使用的 kappa_pairs (count={len(kappa_pairs)}): {kappa_pairs}")
print(f"使用的 mu_B_values (count={len(mu_B_values)}): {mu_B_values}")

SKIP_JULIA = os.environ.get('SKIP_JULIA', '0') == '1'
if SKIP_JULIA:
	print('ℹ️ SKIP_JULIA=1: 已跳过 Julia 初始化（仅测试数据读取）')
	jl_base_objective = None
else:
	try:
		from julia import Main as jl
		# 初始化 Julia 并加载项目（略）
		jl.eval('using Pkg; Pkg.activate("."); Pkg.instantiate()')
		jl.eval('include("src/Gas_Liquid/Advanced_FindTforDiff.jl")')
		# 调用创建目标函数的工厂函数
		jl_base_objective = jl.create_temperature_difference_objective(kappa_pairs, [mu/197.327 for mu in mu_B_values], 70.0/197.327, 120.0/197.327, T_step_scan=2.0/197.327, penalty_for_missing=1e6, verbose=False)
		print('✅ Julia 目标函数已初始化')
	except Exception as e:
		print(f"⚠️ 初始化 Julia 失败: {e}")
		jl_base_objective = None

def objective_wrapper(x):
	if SKIP_JULIA or jl_base_objective is None:
		raise RuntimeError('Julia objective is not initialized (SKIP_JULIA=1).')
	return float(jl_base_objective(tuple(x)))

space = [Real(0.145, 0.170, name='rho0'), Real(-17.0, -15.6, name='B_A'), Real(212.0, 401.0, name='K'), Real(0.55, 0.75, name='m_ratio'), Real(26.1, 44.0, name='E_sym')]

def main():
	try:
		print('开始贝叶斯优化（演示）...')
		res = gp_minimize(objective_wrapper, space, n_calls=10, random_state=42)
		print('优化完成，结果: ', res.x, res.fun)
	except RuntimeError as e:
		print('运行时错误:', e)

if __name__ == '__main__':
	main()
