# Rotation模型的Tmu函数说明

## 新增函数

### 主要函数

1. **`Tmu(;T_start, T_end, T_step, mu_start, mu_end, mu_step)`**
   - 功能：在温度-化学势(T-μ)平面内进行二维扫描，输出状态方程
   - 输出文件：`output/tmu_rotation.csv`
   - 输出列：`T,mu,phi,Phi1,Phi2,mass,pressure,rho,entropy,energy,converged`

2. **`pressure_solve_core(x, mu, T, nodes1, omega)`**
   - 功能：给定μ和T，求解平衡态并返回压力
   - 用途：辅助函数，用于压力导数计算

### 辅助分析函数

3. **`dP_dT_rotation(x, mu, T, nodes1, omega)`** - 压力对温度的一阶导数
4. **`dP_dT2_rotation(x, mu, T, nodes1, omega)`** - 压力对温度的二阶导数
5. **`dP_dT3_rotation(x, mu, T, nodes1, omega)`** - 压力对温度的三阶导数
6. **`dP_dT4_rotation(x, mu, T, nodes1, omega)`** - 压力对温度的四阶导数

## 与PNJL_aniso模型的主要差异

### 变量数量
- **Rotation模型**：3个变量 `[phi, Phi1, Phi2]`
- **PNJL_aniso模型**：5个变量 `[phi_u, phi_d, phi_s, Phi1, Phi2]`

### 函数参数
- **Rotation模型**：使用 `nodes1, omega`
- **PNJL_aniso模型**：使用 `nodes_1, nodes_2, xi`

### 化学势处理
- **Rotation模型**：单一化学势 `mu`
- **PNJL_aniso模型**：三夸克化学势向量，但在Tmu中设为相同值

### 初始值设置
- **Rotation模型**：`[-2.13, 0.06, 0.12]`
- **PNJL_aniso模型**：`[-1.8, -1.8, -2.1, 0.8, 0.8]`

## 使用示例

```julia
# 扫描T从130到131 MeV，μ从400到0 MeV
Tmu(T_start=130/hc, T_end=131/hc, T_step=1/hc, 
    mu_start=400/hc, mu_end=0.0, mu_step=-1/hc)
```

## 输出文件格式

生成的CSV文件包含以下列：
- `T`: 温度 (MeV)
- `mu`: 化学势 (MeV)
- `phi`: 手征凝聚
- `Phi1`, `Phi2`: Polyakov环期望值
- `mass`: 有效质量 (通过calculate_mass(phi)计算)
- `pressure`: 压力
- `rho`: 密度
- `entropy`: 熵
- `energy`: 能量密度
- `converged`: 收敛状态 (true/false)

## 注意事项

1. 所有内部计算都使用自然单位，输出时乘以`hc`转换为物理单位
2. 函数会自动创建`output`目录
3. 使用与现有`Trho`函数相同的节点设置（128×16）
4. 包含异常处理和收敛检查
