# Gas_Liquid 脚本文件夹

这个文件夹包含了与气液相变程序相关的应用层脚本和可视化工具。

> **🏗️ 架构说明**: 核心ForwardDiff自动微分计算功能已迁移到 `src/Gas_Liquid/Advanced_ForwardDiff.jl` 模块中。本文件夹现在专注于应用层脚本和数据可视化。

## � 文件说明

### 1. `plot_temperature_scan.jl` - 数据可视化脚本
**功能**: 读取ForwardDiff温度扫描结果并生成可视化图表

**主要特性**:
- 📈 **数据可视化**: 绘制κ₃/κ₁和κ₄/κ₂随温度的变化曲线
- 📋 **元数据解析**: 自动读取CSV文件中的运行参数元数据
- 🖼️ **图像保存**: 生成高分辨率PNG图像文件
- 📊 **多图模式**: 支持比值图和各个κ值的独立展示

### 2. ~~`forwarddiff_temperature_scan.jl`~~ (已迁移)
**🔄 已迁移到**: `src/Gas_Liquid/Advanced_ForwardDiff.jl` 和 `test/Gas_Liquid/test_forwarddiff_temperature_scan.jl`
- ✅ **函数定义**: 移动到 `src/Gas_Liquid/Advanced_ForwardDiff.jl` 模块
- ✅ **测试程序**: 移动到 `test/Gas_Liquid/test_forwarddiff_temperature_scan.jl`

## 🚀 快速开始

### 使用新的模块化架构运行ForwardDiff计算

#### 第一步：运行ForwardDiff温度扫描测试
```bash
# 在项目根目录下运行新的测试文件
cd D:\Desktop\Julia\Rotation_PNJL
julia --project=. test/Gas_Liquid/test_forwarddiff_temperature_scan.jl
```

**计算配置**:
- 重子化学势: μ_B = 697 MeV (固定)
- 温度范围: 20 - 200 MeV
- 温度步长: 1 MeV
- 输出文件: `output/Gas_Liquid/forwarddiff_temperature_scan.csv`

#### 第二步：生成可视化图表
```bash
# 绘制结果图表
julia --project=. scripts/Gas_Liquid/plot_temperature_scan.jl
```

**输出图像**:
- `output/Gas_Liquid/kappa_ratios_temperature_scan.png` - κ比值图
- `output/Gas_Liquid/individual_kappas_temperature_scan.png` - 各κ值图

### 使用模块化API进行自定义计算

#### 方法1：直接调用Advanced_ForwardDiff模块
```julia
# 引入Advanced_ForwardDiff模块
include("src/Gas_Liquid/Advanced_ForwardDiff.jl")

# 设置模型参数
nodes = get_nodes(256)
couplings = [17.28476, 11.66174, 0.89363, 0.0, 0.00210, -0.00297]
model_params = (nodes, couplings)

# 计算单点热力学涨落
T = 100.0 / hc  # 100 MeV
μ_B = 697.0 / hc  # 697 MeV
κ1, κ2, κ3, κ4, κ3_κ1, κ4_κ2 = calculate_forwarddiff_thermodynamic_fluctuations(
    1.25, 0.01, T, μ_B, model_params)

println("κ₃/κ₁ = $κ3_κ1")
println("κ₄/κ₂ = $κ4_κ2")
```

#### 方法2：使用温度扫描函数
```julia
include("src/Gas_Liquid/Advanced_ForwardDiff.jl")

# 自定义温度扫描
df_results = forwarddiff_temperature_scan(
    697.0/hc,           # μ_B
    50.0/hc,            # T_min  
    150.0/hc,           # T_max
    2.0/hc,             # T_step
    "output/Gas_Liquid/custom_scan.csv";  # 输出文件
    gsigma=1.25,        # 场初值
    gdelta=0.01,
    n_nodes=128         # 减少积分点以提高速度
)
```

**输出图像**:
- `output/Gas_Liquid/kappa_ratios_temperature_scan.png` - κ比值图
- `output/Gas_Liquid/individual_kappas_temperature_scan.png` - 各κ值图

## 📋 详细使用说明

### `plot_temperature_scan.jl` 使用方法

#### 基本运行（自动模式）
```julia
julia --project=. scripts/Gas_Liquid/plot_temperature_scan.jl
```
自动查找 `output/Gas_Liquid/forwarddiff_temperature_scan.csv` 文件并生成图表。

#### 指定输入文件
```julia
julia --project=. scripts/Gas_Liquid/plot_temperature_scan.jl custom_data.csv
```

#### 函数调用方式
```julia
# 在Julia REPL中使用
include("scripts/Gas_Liquid/plot_temperature_scan.jl")

# 绘制比值图
p1 = plot_temperature_scan_results("output/Gas_Liquid/forwarddiff_temperature_scan.csv", 
                                   "my_ratios.png")

# 绘制各个κ值
p2 = plot_individual_kappas("output/Gas_Liquid/forwarddiff_temperature_scan.csv", 
                            "my_kappas.png")
```

### ForwardDiff计算模块使用方法

> **📍 核心计算功能位置**: `src/Gas_Liquid/Advanced_ForwardDiff.jl`

#### 基本ForwardDiff计算
```julia
include("src/Gas_Liquid/Advanced_ForwardDiff.jl")

# 设置模型参数
nodes = get_nodes(256)
couplings = [17.28476, 11.66174, 0.89363, 0.0, 0.00210, -0.00297]
model_params = (nodes, couplings)

# 计算单点导数
T = 100.0 / hc
μ_B = 697.0 / hc
d1, d2, d3, d4 = calculate_forwarddiff_derivatives(1.25, 0.01, T, μ_B, model_params)
```

#### 模型参数配置
```julia
# PNJL模型参数可以在调用时指定：
df_results = forwarddiff_temperature_scan(
    μ_B, T_min, T_max, T_step, output_file;
    gsigma=1.25,        # sigma场初值
    gdelta=0.01,        # delta场初值  
    fs=17.28476,        # sigma耦合常数
    fo=11.66174,        # omega耦合常数
    fr=0.89363,         # rho耦合常数
    fd=0.0,             # delta耦合常数
    b=0.00210,          # 三次项系数
    c=-0.00297,         # 四次项系数
    n_nodes=256         # 积分节点数
)
```

## 📊 架构变更说明

### 🔄 迁移前后对比

| 迁移前 | 迁移后 |
|--------|--------|
| `scripts/Gas_Liquid/forwarddiff_temperature_scan.jl` | `src/Gas_Liquid/Advanced_ForwardDiff.jl` (函数定义) |
| 单文件包含函数+主程序 | `test/Gas_Liquid/test_forwarddiff_temperature_scan.jl` (测试程序) |
| 脚本层面调用 | 模块化API调用 |

### ✅ 新架构优势

1. **模块化设计**: 
   - 函数定义与测试程序分离
   - 可重用的API接口
   - 更清晰的代码结构

2. **更好的可维护性**:
   - 函数集中在src目录便于管理
   - 测试文件独立，便于验证
   - 文档结构更加清晰

3. **灵活的调用方式**:
   - 可以直接调用单个函数
   - 支持批量计算和自定义参数
   - 便于集成到其他脚本中

## 📊 输出文件格式

### CSV数据文件
**文件位置**: `output/Gas_Liquid/forwarddiff_temperature_scan.csv`

**元数据头部** (以 # 开头的注释行):
```
# ForwardDiff Temperature Scan Results
# Generated on: 2024-01-15T10:30:00.000
# Model Parameters:
# gsigma = 1.25
# gdelta = 0.01
# fs = 17.28476
# fo = 11.66174
# fr = 0.89363
# fd = 0.0
# b = 0.00210
# c = -0.00297
# mu_B = 697.0 MeV
# T_range = 20.0 - 200.0 MeV
# T_step = 1.0 MeV
# nodes = 256
```

**数据列说明**:
| 列名 | 单位 | 说明 |
|------|------|------|
| `T_MeV` | MeV | 温度 |
| `P_T4` | 无量纲 | 归一化压强 P/T⁴ |
| `kappa1` | 无量纲 | 一阶累积量 κ₁ |
| `kappa2` | 无量纲 | 二阶累积量 κ₂ |
| `kappa3` | 无量纲 | 三阶累积量 κ₃ |
| `kappa4` | 无量纲 | 四阶累积量 κ₄ |
| `kappa3_over_kappa1` | 无量纲 | 热力学涨落比值 κ₃/κ₁ |
| `kappa4_over_kappa2` | 无量纲 | 热力学涨落比值 κ₄/κ₂ |
| `mu_over_T` | 无量纲 | 化学势温度比 μ_B/T |

### 图像文件
**输出位置**: `output/Gas_Liquid/`

1. **`kappa_ratios_temperature_scan.png`**
   - 主要结果图：κ₃/κ₁ 和 κ₄/κ₂ 随温度变化
   - 包含水平参考线 (比值 = 1)
   - 高分辨率 (300 DPI)

2. **`individual_kappas_temperature_scan.png`**
   - 各个κ值的独立展示
   - 对数坐标显示
   - 用于诊断和详细分析

## ⚙️ 技术细节

### ForwardDiff自动微分
脚本使用Julia的ForwardDiff.jl包进行自动微分：

1. **一阶导数**: `ForwardDiff.derivative()`
2. **高阶导数**: 递归应用自动微分
3. **类型安全**: 使用类型提升确保数值稳定性
4. **增强四阶导数**: 使用5点中心差分提高精度

### 数值方法特点
- **约束方程求解**: 使用NLsolve.jl的Newton方法
- **类型保护**: 动态类型提升避免类型不匹配
- **错误处理**: 计算失败时返回NaN而不是中断
- **内存优化**: 预分配结果数组减少内存分配

### 物理意义
- **κ₃/κ₁**: 反映重子数涨落的偏度，相变信号
- **κ₄/κ₂**: 反映重子数涨落的峰度，临界点标识
- **μ_B/T**: 化学势温度比，控制相图位置

## 🔧 故障排除

### 常见问题及解决方案

**1. 包依赖问题**
```julia
# 确保安装了必要的包
using Pkg
Pkg.add(["NLsolve", "ForwardDiff", "CSV", "DataFrames", "Plots"])
```

**2. 输出目录不存在**
```bash
# 手动创建输出目录
mkdir -p output/Gas_Liquid
```

**3. ForwardDiff计算收敛失败**
- 调整模型参数 (gsigma, gdelta)  
- 减小温度步长
- 检查化学势范围的合理性
- 参考 `src/Gas_Liquid/Advanced_ForwardDiff.jl` 中的错误处理

**4. 模块导入失败**
```julia
# 确保使用正确的相对路径
include("src/Gas_Liquid/Advanced_ForwardDiff.jl")  # 从项目根目录运行
include("../../src/Gas_Liquid/Advanced_ForwardDiff.jl")  # 从test目录运行
```

**5. 绘图失败 (SSH环境)**
```julia
# 确保使用GR后端且保存文件而不显示
using Plots
gr()  # 设置后端
# 脚本会自动保存PNG文件
```

**5. 内存不足**
- 减少温度扫描范围
- 增大温度步长
- 监控系统内存使用

### 性能优化建议

**1. 计算性能**
- 使用较少的积分节点 (如128而不是256)
- 增大温度步长减少计算点数
- 在多核系统上可考虑并行化

**2. 数据管理**
- 定期清理output目录中的旧文件
- 使用压缩格式保存大数据集
- 备份重要的计算结果

## 📚 相关文档

- **主项目文档**: `../../docs/`
- **ForwardDiff模块**: `../../src/Gas_Liquid/Advanced_ForwardDiff.jl`
- **模块分离说明**: `../../src/Gas_Liquid/MODULE_SEPARATION_README.md`
- **核心函数文档**: `../../src/Gas_Liquid/Function_Gas_Liquid.jl`
- **测试文件**: `../../test/Gas_Liquid/test_forwarddiff_temperature_scan.jl`
- **使用指南**: `../../USAGE_GUIDE.md`

## 🤝 贡献指南

如需修改或扩展这些脚本：

1. **备份原始文件**
2. **更新模块化代码**: 修改 `src/Gas_Liquid/Advanced_ForwardDiff.jl` 中的函数
3. **更新测试文件**: 修改 `test/Gas_Liquid/test_forwarddiff_temperature_scan.jl`
4. **测试修改后的功能**  
5. **更新文档和注释**
6. **保持API兼容性**

---
**最后更新**: 2025年9月14日  
**架构版本**: v2.0 (模块化)  
**维护者**: Rotation_PNJL项目组
4. **添加适当的错误处理**
5. **保持代码风格一致性**

---
**最后更新**: 2025年9月14日  
**版本**: v1.0  
**维护者**: Rotation_PNJL项目组
