Usage Guide — Rotation_PNJL
===========================

前提
-----
- 需要安装 Julia（建议 1.8+ 或 1.11.x），并在 REPL 中能载入本目录的依赖包（若使用 CSV.jl、ForwardDiff、NLsolve 等）。

快速运行示例（PowerShell）
----------------------------
在项目根路径下运行：

```powershell
# 交互式 include 后运行函数
julia -e 'include("d:/Desktop/Julia/Rotation_PNJL/src/Function_PNJL_aniso.jl"); Trho(100/hc, 110/hc)'

# 或运行 Tmu 扫描（使用位置参数或关键字参数，视函数签名而定）
julia -e 'include("d:/Desktop/Julia/Rotation_PNJL/src/Function_PNJL_aniso.jl"); Tmu(130/hc,140/hc,1/hc,0.0,400/hc,10/hc)'
```

输出
-----
- 生成的 CSV 文件位于 `Rotation_PNJL/output/`，文件名例如 `trho_scan.csv`, `tmu_scan.csv` 或 `trho_scan_rotation.csv`。

注意事项
---------
- 文件中可能含有示例调用；在 include 时会执行这些调用。建议在开发期间将示例调用注释掉，或在交互式会话中手动运行。
- 对于大型扫描，建议在脚本中设置较小的节点数或更小的网格范围进行测试后再跑全域扫描。

安装依赖（推荐）
-----------------
本项目包含 `Project.toml`，推荐使用 Julia 的包环境来安装和管理依赖。以下为在 Windows PowerShell 下的示例：

1. 激活并安装 Project.toml 中列出的依赖（会在该目录创建/使用本地环境）：

```powershell
julia -e "using Pkg; Pkg.activate('d:/Desktop/Julia/Rotation_PNJL'); Pkg.instantiate()"
```

2. 如果需要单独添加某个包（例如手工添加 CSV 或特定版本）：

```powershell
julia -e "using Pkg; Pkg.activate('d:/Desktop/Julia/Rotation_PNJL'); Pkg.add('CSV')"
```

3. 进入交互式 REPL 进行调试时，可在 Julia REPL 内执行：

```julia
using Pkg
Pkg.activate("d:/Desktop/Julia/Rotation_PNJL")
Pkg.instantiate()
```

注意：Pkg.instantiate() 会根据 `Project.toml`（及可能存在的 `Manifest.toml`）安装依赖，执行前请确保网络可用。对于首次大型依赖安装请耐心等待。
