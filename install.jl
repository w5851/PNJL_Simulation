# install.jl
# 一键安装和配置Rotation_PNJL项目环境

println("🚀 Rotation_PNJL 项目环境安装器")
println("=" ^ 60)

# 检查当前目录是否为项目根目录
if !isfile("Project.toml")
    println("❌ 错误: 未找到Project.toml文件")
    println("   请确保在项目根目录运行此脚本")
    exit(1)
end

println("📁 当前工作目录: $(pwd())")
println("✅ 找到项目配置文件: Project.toml")

# 激活项目环境
println("\n🔧 激活项目环境...")
import Pkg
Pkg.activate(".")
println("✅ 项目环境已激活")

# 读取并显示项目信息
println("\n📋 项目信息:")
try
    project_toml = Pkg.project()
    if haskey(project_toml, "name")
        println("   项目名称: $(project_toml.name)")
    end
    if haskey(project_toml, "version")
        println("   项目版本: $(project_toml.version)")
    end
    println("   项目路径: $(project_toml.path)")
catch e
    println("   无法读取项目信息: $e")
end

# 检查Julia版本兼容性
println("\n🔍 检查Julia版本兼容性...")
julia_version = VERSION
println("   当前Julia版本: $julia_version")

# 读取Project.toml中的compat信息
try
    using TOML
    project_data = TOML.parsefile("Project.toml")
    if haskey(project_data, "compat") && haskey(project_data["compat"], "julia")
        required_version = project_data["compat"]["julia"]
        println("   项目要求Julia版本: $required_version")
        if julia_version >= v"1.6"
            println("   ✅ Julia版本兼容")
        else
            println("   ⚠️  Julia版本可能不兼容，建议升级到1.6+")
        end
    end
catch e
    println("   ⚠️  无法检查版本兼容性: $e")
end

# 更新注册表
println("\n📦 更新包注册表...")
try
    Pkg.Registry.update()
    println("✅ 注册表更新完成")
catch e
    println("⚠️  注册表更新失败: $e")
    println("   继续安装包...")
end

# 安装依赖包
println("\n📚 安装项目依赖包...")
println("   这可能需要几分钟时间，请耐心等待...")

try
    # 首先解析依赖关系
    Pkg.resolve()
    println("✅ 依赖关系解析完成")
    
    # 安装包
    Pkg.instantiate()
    println("✅ 依赖包安装完成")
    
    # 预编译包
    println("\n🔨 预编译包...")
    Pkg.precompile()
    println("✅ 包预编译完成")
    
catch e
    println("❌ 包安装失败: $e")
    println("\n🔧 尝试手动安装关键包...")
    
    # 手动安装关键包
    key_packages = [
        "BayesianOptimization",
        "GaussianProcesses", 
        "ForwardDiff",
        "NLsolve",
        "SpecialFunctions",
        "StaticArrays",
        "CSV",
        "DataFrames",
        "Plots",
        "FastGaussQuadrature",
        "FiniteDifferences",
        "BenchmarkTools"
    ]
    
    for pkg in key_packages
        try
            println("   安装 $pkg...")
            Pkg.add(pkg)
            println("   ✅ $pkg 安装成功")
        catch pkg_error
            println("   ❌ $pkg 安装失败: $pkg_error")
        end
    end
end

# 验证关键包
println("\n🧪 验证关键包安装...")
test_packages = [
    ("BayesianOptimization", "贝叶斯优化"),
    ("GaussianProcesses", "高斯过程"),
    ("ForwardDiff", "自动微分"),
    ("NLsolve", "非线性求解器"),
    ("SpecialFunctions", "特殊函数"),
    ("StaticArrays", "静态数组"),
    ("CSV", "CSV文件处理"),
    ("DataFrames", "数据框处理"),
    ("Plots", "数据可视化"),
    ("FastGaussQuadrature", "快速高斯积分"),
    ("FiniteDifferences", "有限差分"),
    ("BenchmarkTools", "性能测试")
]

success_count = 0
for (pkg_name, desc) in test_packages
    try
        eval(:(using $(Symbol(pkg_name))))
        println("   ✅ $desc ($pkg_name)")
        success_count += 1
    catch e
        println("   ❌ $desc ($pkg_name) - 无法加载: $e")
    end
end

# 显示安装状态
println("\n📊 安装状态:")
Pkg.status()

# 创建测试脚本快捷方式
println("\n🎯 创建测试快捷方式...")
try
    # 检查test目录下是否有测试文件
    test_files = []
    if isdir("test")
        for file in readdir("test")
            if endswith(file, ".jl") && occursin("test", lowercase(file))
                push!(test_files, file)
            end
        end
    end
    
    if !isempty(test_files)
        println("   发现测试文件:")
        for file in test_files
            println("     - test/$file")
        end
        
        # 创建运行测试的快捷脚本
        test_script = """
# run_tests.jl
# 快速运行项目测试

import Pkg
Pkg.activate(".")

println("🧪 运行项目测试...")
include("test/runtests.jl")
"""
        
        open("run_tests.jl", "w") do f
            write(f, test_script)
        end
        println("   ✅ 创建测试运行脚本: run_tests.jl")
        
        # 如果有BayesianOptimization测试，创建专门的运行脚本
        if any(occursin("bayesian", lowercase(file)) for file in test_files)
            bayesian_script = """
# run_bayesian_test.jl
# 快速运行贝叶斯优化测试

import Pkg
Pkg.activate(".")

println("🧠 运行贝叶斯优化测试...")
for file in readdir("test")
    if occursin("bayesian", lowercase(file)) && endswith(file, ".jl")
        println("执行: test/\$file")
        include("test/\$file")
        break
    end
end
"""
            open("run_bayesian_test.jl", "w") do f
                write(f, bayesian_script)
            end
            println("   ✅ 创建贝叶斯测试运行脚本: run_bayesian_test.jl")
        end
    end
catch e
    println("   ⚠️  创建测试快捷方式失败: $e")
end

# 最终结果
println("\n" * "=" ^ 60)
println("🎉 安装完成!")
println("=" ^ 60)

if success_count >= 6
    println("✅ 环境安装成功! ($success_count/$(length(test_packages))个关键包已安装)")
    println("\n📖 使用说明:")
    println("1. 运行项目测试: julia run_tests.jl")
    if isfile("run_bayesian_test.jl")
        println("2. 运行贝叶斯测试: julia run_bayesian_test.jl")
    end
    println("3. 进入Julia REPL:")
    println("   julia --project=.")
    println("4. 或在脚本中激活环境:")
    println("   import Pkg; Pkg.activate(\".\")")
    
    println("\n🚀 项目已准备就绪!")
else
    println("⚠️  环境安装部分成功 ($success_count/$(length(test_packages))个包)")
    println("   某些包可能需要手动安装或检查网络连接")
    println("   运行 'julia install.jl' 重试安装")
end

println("\n💡 提示: 如遇到问题，请检查:")
println("   - 网络连接是否正常")
println("   - Julia版本是否为1.6+")
println("   - 是否有足够的磁盘空间")
println("=" ^ 60)