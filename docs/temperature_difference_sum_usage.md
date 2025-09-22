# 温度差平方和计算函数使用说明

## 概述
基于 `Advanced_FindTforDiff.jl` 中的批量计算逻辑，新增了两个专门用于优化的函数，用于计算多组κ比值对应的温度差的平方和。

## 主要函数

### 1. `calculate_temperature_difference_sum_of_squares`

计算多组κ比值在给定优化参数下对应的温度差的平方和，返回标量值（MeV²单位）。

#### 函数签名
```julia
calculate_temperature_difference_sum_of_squares(
    kappa_pairs, μ_B, optimization_params, T_min, T_max;
    T_step_scan=1.0/hc, gsigma=1.25, gdelta=0.01, n_nodes=256,
    penalty_for_missing=1e6, verbose=false)
```

### 2. `calculate_temperature_difference_sum_of_squares_with_weights`

计算加权温度差平方和，允许为不同的κ比值对分配不同权重。

### 3. `create_temperature_difference_objective` ⭐ **推荐用于优化**

创建温度差平方和目标函数的闭包，将实验确定的参数封装起来，返回只需要优化参数作为输入的函数。

#### 函数签名
```julia
create_temperature_difference_objective(
    kappa_pairs, μ_B, T_min, T_max;
    T_step_scan=1.0/hc, gsigma=1.25, gdelta=0.01, n_nodes=256,
    penalty_for_missing=1e6, verbose=false)
```

#### 返回值
- `objective_function`: 闭包函数，签名为 `f(optimization_params) -> Float64`

### 4. `create_weighted_temperature_difference_objective`

创建加权版本的目标函数闭包。

## 使用示例

### 基本使用
```julia
# 包含模块
# include("src/Gas_Liquid/Advanced_FindTforDiff.jl")

# 定义输入参数
kappa_pairs = [(1.2, 2.5), (1.5, 3.0), (1.8, 3.5)]
μ_B = 300.0 / hc  # 300 MeV
optimization_params = (0.15, 16.0, 240.0, 0.7, 32.0)  # (ρ₀, B_A, K, m_ratio, E_sym)
T_min, T_max = 80.0/hc, 200.0/hc  # 80-200 MeV

# 直接计算温度差平方和
sum_sq = calculate_temperature_difference_sum_of_squares(
    kappa_pairs, μ_B, optimization_params, T_min, T_max,
    verbose=false)

println("温度差平方和: $sum_sq MeV²")
```

### 推荐用法：使用闭包函数 ⭐

```julia
# 实验确定的参数（固定）
kappa_pairs = [(1.2, 2.5), (1.5, 3.0), (1.8, 3.5)]
μ_B_values = [300.0/hc, 320.0/hc, 340.0/hc]  # 每组对应不同的μ_B值 ⚠️ 重要变更
T_min, T_max = 80.0/hc, 200.0/hc  # 80-200 MeV

# 创建目标函数闭包
objective_func = create_temperature_difference_objective(
    kappa_pairs, μ_B_values, T_min, T_max; 
    verbose=false, penalty_for_missing=1e6)

# 现在目标函数只需要优化参数输入
optimization_params = (0.15, 16.0, 240.0, 0.7, 32.0)
result = objective_func(optimization_params)

println("目标函数值: $result MeV²")

# 可以直接用于优化算法
using Optim
initial_params = [0.15, 16.0, 240.0, 0.7, 32.0]
result = optimize(params -> objective_func(tuple(params...)), initial_params)
```

### 加权版本使用
```julia
# 使用权重版本
weights = [1.0, 2.0, 0.5]
weighted_objective_func = create_weighted_temperature_difference_objective(
    kappa_pairs, weights, μ_B_values, T_min, T_max; verbose=false)

weighted_result = weighted_objective_func(optimization_params)
println("加权目标函数值: $weighted_result MeV²")
```

### ⚠️ 重要变更说明

**从单一μ_B到μ_B数组**：
- **旧版本**: 所有κ比值对使用相同的μ_B值
- **新版本**: 每组κ比值对对应独立的μ_B值，更符合实验实际

**函数签名变更**：
```julia
# 旧版本
calculate_temperature_difference_sum_of_squares(kappa_pairs, μ_B, ...)

# 新版本  
calculate_temperature_difference_sum_of_squares(kappa_pairs, μ_B_values, ...)
```

**输入数据格式**：
```julia
# 示例：3组实验数据
kappa_pairs = [(1.2, 2.5), (1.5, 3.0), (1.8, 3.5)]
μ_B_values = [300.0/hc, 320.0/hc, 340.0/hc]  # 一一对应

# 数组长度必须相等
assert(length(kappa_pairs) == length(μ_B_values))
```

## 运行演示

```julia
# 运行基本演示函数
demo_temperature_difference_sum_of_squares()

# 运行闭包演示函数 ⭐ 推荐
demo_objective_function_closure()
```

## 闭包函数的优势

### 1. **简化优化接口**
```julia
# 传统方式：需要传递多个固定参数
function objective(params)
    return calculate_temperature_difference_sum_of_squares(
        kappa_pairs, μ_B, params, T_min, T_max; verbose=false)
end

# 闭包方式：一次设置，重复使用
objective_func = create_temperature_difference_objective(
    kappa_pairs, μ_B, T_min, T_max; verbose=false)
# 现在 objective_func 只需要 optimization_params
```

### 2. **与优化库完美集成**
```julia
using Optim, NLopt

# 直接用于Optim.jl
result = optimize(objective_func, initial_guess)

# 直接用于NLopt.jl
opt = Opt(:LN_NELDERMEAD, 5)  # 5个优化参数
opt.min_objective = (params, grad) -> objective_func(params)
```

### 3. **参数验证和错误处理**
闭包在创建时验证参数有效性，避免运行时错误。

## 错误处理

函数具有完善的错误处理机制：
1. **计算失败**: 当某组κ比值的温度计算失败时，自动应用惩罚值
2. **无效温度**: 当找不到有效温度时，同样应用惩罚值
3. **优化友好**: 即使部分计算失败，函数仍能返回有效的目标值用于优化

## 优化应用

这些函数特别适用于：
1. **参数优化**: 作为优化算法的目标函数
2. **敏感性分析**: 评估参数变化对温度一致性的影响
3. **模型验证**: 检验理论模型的内在一致性

## 性能考虑

- 建议在优化过程中设置 `verbose=false` 以减少输出
- 可以调整 `penalty_for_missing` 值来控制优化行为
- `T_step_scan` 参数影响计算精度和速度的平衡