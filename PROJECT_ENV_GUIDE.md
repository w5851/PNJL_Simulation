# 🚀 项目环境使用指南

## 📦 一键安装和配置

在项目根目录运行：

```bash
julia install.jl
```

这将自动：
- ✅ 激活项目环境
- ✅ 安装所有依赖包
- ✅ 验证包的正确性
- ✅ 创建测试快捷脚本

## 🎯 快速开始

### 方式1: 使用项目运行器（推荐）
```bash
julia run_project.jl
```
提供交互式菜单，包含所有项目功能。

### 方式2: 直接运行功能
```bash
# 贝叶斯优化测试
julia test/test_bayesian_fixed.jl

# 可视化测试
julia test/run_visualization_test.jl

# 气液相变测试
julia test/Gas_Liquid/test_find_temperature.jl
```

### 方式3: Julia REPL中使用
```julia
julia> import Pkg
julia> Pkg.activate(".")
julia> using BayesianOptimization, Plots
```

## 📋 项目结构

```
📁 Rotation_PNJL/
├── 🔧 install.jl           # 一键安装脚本
├── 🎯 run_project.jl       # 统一运行器
├── 📦 Project.toml         # 项目依赖配置
├── 📚 src/                 # 核心源代码
├── 🧪 test/               # 测试文件
├── 📊 scripts/            # 演示脚本
└── 📈 output/             # 输出结果
```

## 🔍 环境激活模式

项目中的所有脚本已统一使用项目环境：

### 自动激活（推荐）
所有脚本已内置环境激活：
```julia
import Pkg
Pkg.activate(".")  # 或相对路径
```

### 手动激活
如需在REPL中手动激活：
```julia
# 在项目根目录
import Pkg; Pkg.activate(".")

# 在子目录（如test/）
import Pkg; Pkg.activate("..")
```

## 📦 包依赖列表

项目环境包含以下包：
- **BayesianOptimization** - 贝叶斯优化
- **GaussianProcesses** - 高斯过程
- **ForwardDiff** - 自动微分
- **NLsolve** - 非线性求解器
- **SpecialFunctions** - 特殊函数
- **StaticArrays** - 静态数组
- **CSV & DataFrames** - 数据处理
- **Plots** - 数据可视化
- **FastGaussQuadrature** - 快速高斯积分
- **FiniteDifferences** - 有限差分
- **BenchmarkTools** - 性能测试

## 🛠️ 故障排除

### 包缺失问题
```bash
# 重新安装环境
julia install.jl

# 手动安装单个包
julia -e "import Pkg; Pkg.activate('.'); Pkg.add('PackageName')"
```

### 环境冲突问题
```bash
# 清理并重建环境
julia -e "import Pkg; Pkg.activate('.'); Pkg.resolve(); Pkg.instantiate()"
```

### 路径问题
确保在正确的目录运行脚本：
```bash
cd "path/to/Rotation_PNJL"
julia script_name.jl
```

## 💡 最佳实践

1. **总是使用项目环境**：确保所有脚本激活项目环境
2. **统一入口点**：推荐使用`run_project.jl`统一运行
3. **环境隔离**：不要在全局环境安装项目专用包
4. **版本固定**：依赖特定版本的包通过Project.toml管理

## 🎉 现在开始使用

```bash
# 第一次使用
julia install.jl

# 日常使用
julia run_project.jl
```

享受Rotation_PNJL项目的强大功能！🚀