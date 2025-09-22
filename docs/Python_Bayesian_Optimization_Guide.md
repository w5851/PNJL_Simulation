# PNJL模型Python贝叶斯优化安装和使用指南

## 概述

本项目实现了使用Python的贝叶斯优化库调用Julia函数进行PNJL模型参数优化的功能，与Julia版本的`demo_bayesian_optimization_with_warmup()`提供相同的功能，但具有以下优势：

- 🚀 更丰富的优化库生态（scikit-optimize, optuna等）
- 📊 更好的数据可视化能力
- 🔧 更灵活的数据处理和分析
- 🧪 更方便的实验和调试

## 安装依赖

### 1. Python依赖

```bash
# 基础科学计算库
pip install numpy pandas matplotlib

# 贝叶斯优化库
pip install scikit-optimize

# Python调用Julia接口
pip install julia
```

### 2. 设置PyJulia

```bash
# 初始化PyJulia
python -c "import julia; julia.install()"
```

### 3. 验证Julia环境

确保Julia项目环境正确配置：

```bash
cd d:/Desktop/Julia/Rotation_PNJL
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

## 快速开始

### 1. 运行功能测试

```bash
cd d:/Desktop/Julia/Rotation_PNJL
python scripts/test_python_bayesian.py
```

### 2. 运行演示

```bash
# 基础演示
python src/Gas_Liquid/Advanced_Bayesian.py --mode demo

# 快速测试
python src/Gas_Liquid/Advanced_Bayesian.py --mode test

# 方法比较
python src/Gas_Liquid/Advanced_Bayesian.py --mode compare

# 查看帮助
python src/Gas_Liquid/Advanced_Bayesian.py --mode help
```

### 3. 从CSV文件继续优化

```bash
python src/Gas_Liquid/Advanced_Bayesian.py --csv path/to/previous_results.csv
```

## 主要功能

### 1. 基础优化功能

```python
from src.Gas_Liquid.Advanced_Bayesian import PNJLBayesianOptimizer

# 创建优化器
optimizer = PNJLBayesianOptimizer()

# 设置实验数据
kappa_pairs = [
    (1.09031788496341, -0.28904867673079),   # 第1组
    (1.06152332992368, 0.164279260625683),   # 第2组
    (1.11111023684003, 0.224522832511389)    # 第3组
]

mu_B_values = [632.0, 666.0, 697.0]  # MeV
T_min, T_max = 70.0, 120.0           # MeV

# 参数边界 [ρ₀, B_A, K, m_ratio, E_sym]
param_bounds = [
    (0.145, 0.170),    # ρ₀ (fm⁻³)
    (-17.0, -15.6),    # B_A (MeV)
    (212.0, 401.0),    # K (MeV)
    (0.55, 0.75),      # m_ratio
    (26.1, 44.0)       # E_sym (MeV)
]

# 执行优化
result = optimizer.optimize_with_warmup(
    kappa_pairs=kappa_pairs,
    mu_B_values=mu_B_values,
    T_min=T_min,
    T_max=T_max,
    param_bounds=param_bounds,
    max_iterations=20,
    initial_samples=10,
    T_step_scan=2.0,
    acquisition_function='EI',
    output_file="output/optimization_result.csv"
)
```

### 2. 从先前结果继续优化

```python
from src.Gas_Liquid.Advanced_Bayesian import continue_optimization_from_csv

result = continue_optimization_from_csv(
    csv_file="previous_results.csv",
    kappa_pairs=kappa_pairs,
    mu_B_values=mu_B_values,
    T_min=T_min,
    T_max=T_max,
    param_bounds=param_bounds,
    additional_iterations=25,
    output_file="continued_optimization.csv"
)
```

### 3. 高级可视化

```python
from src.Gas_Liquid.Advanced_Bayesian import advanced_visualization

# 生成详细的优化分析图
advanced_visualization(result, save_dir="output/plots/")
```

## 与Julia版本对比

| 特性 | Julia版本 | Python版本 |
|------|-----------|-------------|
| 优化算法 | BayesianOptimization.jl | scikit-optimize |
| 性能 | 高（JIT编译） | 中等（但足够用） |
| 可视化 | 基础 | 丰富（matplotlib） |
| 数据处理 | 基础 | 强大（pandas） |
| 生态系统 | 科学计算专用 | 通用机器学习 |
| 学习曲线 | 陡峭 | 平缓 |
| 调试便利性 | 一般 | 优秀 |

## 核心优势

### 1. 更丰富的优化选项
- 支持多种采集函数（EI, LCB, PI等）
- 可以轻松切换到其他优化库（如optuna）
- 更好的超参数控制

### 2. 更好的数据分析
- 自动保存优化历史到CSV
- 支持从之前的结果继续优化
- 丰富的可视化选项

### 3. 更好的用户体验
- 清晰的进度显示
- 详细的错误处理
- 灵活的配置选项

## 技术实现

### 核心架构

```
Python (贝叶斯优化) ←→ PyJulia接口 ←→ Julia (PNJL物理计算)
     ↓                                      ↑
scikit-optimize                    Advanced_FindTforDiff.jl
     ↓                                      ↑  
 优化决策                              目标函数评估
```

### 关键组件

1. **PNJLBayesianOptimizer**: 主优化器类
2. **create_objective_function**: 目标函数构造器
3. **warmup_objective_function**: 预热功能
4. **advanced_visualization**: 高级可视化

## 故障排除

### 常见问题

1. **PyJulia导入失败**
   ```bash
   pip install julia
   python -c "import julia; julia.install()"
   ```

2. **Julia函数调用失败**
   - 检查Julia项目路径
   - 确认Julia环境已激活
   - 验证必要的Julia包已安装

3. **优化收敛慢**
   - 增加初始采样点数
   - 调整参数边界范围
   - 尝试不同的采集函数

4. **内存不足**
   - 减少最大迭代次数
   - 增大温度扫描步长
   - 减少实验数据点数

## 性能优化建议

1. **调整扫描精度**: 增大`T_step_scan`可以显著提高速度
2. **合理设置边界**: 窄的参数边界有助于更快收敛
3. **预热功能**: 用于估算计算时间，可以跳过以节省时间
4. **并行化**: 未来可以扩展为并行评估多个参数点

## 扩展性

本实现设计为可扩展的：

- 可以轻松添加新的采集函数
- 支持集成其他优化库（optuna, hyperopt等）
- 可以扩展到多目标优化
- 支持约束优化

## 参与贡献

欢迎改进和扩展此实现：

1. 添加更多的优化算法
2. 改进可视化功能
3. 优化性能
4. 添加更多的物理约束

## 联系和支持

如有问题或建议，请参考：
- Julia版本实现：`Advanced_BayesianOptimization.jl`
- 技术文档：项目docs目录
- 测试用例：`test_python_bayesian.py`