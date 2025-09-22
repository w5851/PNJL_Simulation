# Function_Rotation.jl 说明文档

文件路径: `src/Function_Rotation.jl`

## 一句话概述
为旋转（rotation）模型计算热力学量、密度扫描与生成 `output/trho_rotation.csv` 的主计算脚本，包含节点生成、压力/能量/熵的计算以及对温度-密度路径的根求解 (nlsolve)。

## 主要用途
- 生成 T–rho 扫描并将结果写入 `output/trho_rotation.csv`。
- 提供压力、能量、熵及相空间积分相关的工具函数，供数值扫描与最小化/求根使用。

## 主要依赖
- 本地模块/文件：`Constants_Rotation.jl`、`init3.jl`（提供常量与 `gauleg` 节点）
- 外部包：`SpecialFunctions`, `ForwardDiff`, `NLsolve`, `BenchmarkTools`, `StaticArrays`, `FiniteDifferences`

请确保在运行前在对应项目环境中安装并激活上述包（Project.toml 已在仓内）。

## 输出
- 路径: `output/trho_rotation.csv`（相对 `Rotation_PNJL` 根目录）
- 列名（文件首行）: T,rho,phi,Phi1,Phi2,mu,pressure,entropy,energy,converged
  - 注意：代码里写入的 T 值为 `T*hc`（见文件实现），mu 在内部以能量单位（除以 hc）处理。

## 文件中重要函数与变量（摘要）
- get_nodes(p_num::Int, t_num::Int)
  - 作用：构造 (p, θ, n) 三维离散网格与权重，并返回两个节点集合 `nodes1, nodes2`（不同动量区间）。
  - 输出：`nodes1 = [p_mesh, n_mesh, coefficient1]`，`nodes2 = [p2_mesh, n2_mesh, coefficient2]`。

- init_bessel(p,theta,n,w)
  - 作用：计算与贝塞尔函数相关的系数项，返回整合后权重因子。

- calculate_chiral(phi)
  - 作用：返回手征相关项：G_f * mean(phi.^2)

- calculate_U(T, Phi1, Phi2)
  - 作用：计算并返回 U(Phi)（多项式形式，乘以 T^4）。

- calculate_mass(phi)
  - 作用：基于 phi 返回粒子有效质量 m = m0_q_f - 2 G_f phi。

- calculate_energy(mass, p, n, omega)
  - 作用：返回单粒子能量项 sqrt(p^2 + mass^2) - (0.5 + n) * omega

- AA, AAbar, calculate_log_term
  - 作用：封装分布函数相关表达式并返回对数项用于压力积分。

- calculate_log_sum(masses, p_nodes, Phi1, Phi2, mu, T, coefficient, n_nodes, omega)
  - 作用：对网格逐元素累加对数项，返回 -T * sum(log_terms * weights)。
  - 性能说明：内部使用双层循环 + @simd，针对元素计算进行了局部缓存优化。

- calculate_pressure(phi,Phi1,Phi2,mu,T,nodes1,omega)
  - 作用：将 chi、U 与对数积分项合成并返回压力（注意函数定义中返回值已取负号以实现特定最小化/根求策略）。

- calculate_core(x, mu, T, nodes1, omega)
  - 作用：对 x = (phi1, Phi1, Phi2) 使用 ForwardDiff 计算梯度（用于 nlsolve 的方程组）。

- calculate_rho / calculate_thermo
  - 作用：计算密度、熵、能量等热力学量（利用 ForwardDiff 求导数）。

- calculate_t_rho(x,T,rho,nodes1,omega)
  - 作用：构造 nlsolve 所需的 4 元方程，前三为梯度（gap 方程），第四为密度差（rho - 目标 rho）。

- Trho(T_start,T_end)
  - 作用：主控制流程，生成节点（get_nodes(128,16)）、循环温度 T，并对每个 T 做 rho 从 3.00 递减到 0.10 的扫面，调用 nlsolve 记录结果到 CSV。
  - 注意：文件末尾当前会调用 `Trho(100/hc,101/hc)` 以示例运行。

## 关键常量与假设（来自 `Constants_Rotation.jl` 的引用）
- hc, π, rho0, a0..a3, b3, b4, T0, Nc, Lambda_f, G_f, K_f, m0_q_f, m0_s_f, r0, coefficients
- 假设：`rho0` 为归一化密度基准，`Lambda_f` 为动量上限，`r0` 用于贝塞尔项的空间尺度。

## 已识别的关键信息项（需要在开始前确认）
1. `Constants_Rotation.jl` 的精确常量值与 `coefficients` 的结构（代码引用 `Constants.coefficients`，但文件中使用 `Constants_Rotation`，需确认命名空间一致性）。
2. `init3.jl` 中的 `gauleg` 实现细节（区间、精度），以确认节点/权重的数值意义与单位。
3. `hc` 的单位与 T、mu 的物理单位换算（代码在写入 CSV 时使用 `T*hc`）。
4. 是否期望 `Trho` 在模块加载时自动运行（文件当前尾部调用 `Trho(...)`），或仅作为 API 被调用。

缺失项标注为问题（一次只问一个）：
- 问：您确认希望文档保存为 `src/Function_Rotation.md` 并以中文撰写吗？（是/否）

## 最小可行测试（smoke test）
目的：验证能在本机用最小配置成功运行并生成 CSV 文件头与部分数据。

步骤（Windows PowerShell，可复制执行）：

```powershell
cd D:\Desktop\Julia\Rotation_PNJL
# 激活项目环境（如需要）并运行文件
julia --project=. src\Function_Rotation.jl
# 查看生成的文件头与前5行数据
Get-Content .\output\trho_rotation.csv -TotalCount 6
```

预期输出：CSV 第一行为列头 `T,rho,phi,Phi1,Phi2,mu,pressure,entropy,energy,converged`，后续行包含数值（T 单位与代码实现相关）。若 `Trho(...)` 被注释或删除，则不会生成文件。

## 已知风险与改进建议
- 收敛失败：使用 `nlsolve` 时可能在某些 (T,rho) 不收敛，当前策略是写入 NaN 并记录 warning。建议增加更稳定的初值策略（例如使用相邻 mu 的解或前一步解作为猜测，这在仓内其他示例可能已有实现）。
- 性能：节点数 (128,16) + bessel 计算开销大，可通过：
  - 预计算并缓存常用项，避免在循环内重复计算相同 bessel 值；
  - 减小 p, t 网格进行快速 smoke 测试；
  - 使用多线程或并行化（需注意 ForwardDiff 与线程的兼容）。
- 数值稳定性：对数项中 f1,f2 等可能接近 0，建议在临界区域加入小量截断以避免 log(0)。

## 建议的后续小任务
- 将 `Trho(...)` 的自动执行改为可选（例如包外部调用），以避免在导入模块时立即执行长时间任务。
- 编写单元测试：针对 `get_nodes`（尺寸与权重和）、`calculate_log_term`（已知输入输出）等提供快速单元测试。

---

## 变更记录
- 文档作者：agent 自动生成（中文）
- 说明：此文档仅为代码说明与运行指南，不修改源代码。



