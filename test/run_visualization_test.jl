# run_visualization_test.jl
# 测试可视化功能

import Pkg
# 获取当前文件的目录，然后向上一级到项目根目录
project_root = joinpath(@__DIR__, "..")
Pkg.activate(project_root)
println("✅ 激活项目环境: ", project_root)

println("🎨 测试数据可视化功能")
println("=" ^ 60)

# 测试基础包加载
println("\n📦 加载基础包...")
try
    using Plots, CSV, DataFrames
    println("✅ Plots.jl 加载成功")
    println("✅ CSV.jl 加载成功") 
    println("✅ DataFrames.jl 加载成功")
catch e
    println("❌ 包加载失败: $e")
    println("请运行: julia install.jl")
    exit(1)
end

# 设置绘图后端
println("\n🎨 设置绘图后端...")
try
    gr()
    println("✅ GR后端设置成功")
catch e
    println("⚠️  GR后端设置失败，尝试其他后端...")
    try
        plotlyjs()
        println("✅ PlotlyJS后端设置成功")
    catch e2
        println("❌ 所有后端都失败: $e2")
    end
end

# 创建测试数据
println("\n📊 创建测试数据...")
x = 0:0.1:4π
y1 = sin.(x)
y2 = cos.(x)

# 测试基本绘图
println("\n🎯 测试基本绘图功能...")
try
    p = plot(x, y1, 
             label="sin(x)", 
             linewidth=2,
             title="测试绘图",
             xlabel="x", 
             ylabel="y")
    
    plot!(p, x, y2, 
          label="cos(x)", 
          linewidth=2, 
          linestyle=:dash)
    
    println("✅ 基本绘图功能正常")
    
    # 保存图片
    savefig(p, "test_plot.png")
    println("✅ 图片保存为: test_plot.png")
    
catch e
    println("❌ 绘图测试失败: $e")
end

# 测试数据处理
println("\n📋 测试数据处理功能...")
try
    # 创建测试DataFrame
    df = DataFrame(
        x = x[1:10:end],
        sin_x = sin.(x[1:10:end]),
        cos_x = cos.(x[1:10:end])
    )
    
    println("✅ DataFrame创建成功")
    println("   数据维度: $(size(df))")
    
    # 保存CSV
    CSV.write("test_data.csv", df)
    println("✅ CSV保存成功: test_data.csv")
    
    # 读取CSV
    df_read = CSV.read("test_data.csv", DataFrame)
    println("✅ CSV读取成功")
    
catch e
    println("❌ 数据处理测试失败: $e")
end

# 测试项目中的绘图脚本
println("\n🔍 测试项目绘图脚本...")
plot_script = "scripts/Gas_Liquid/plot_temperature_scan.jl"
if isfile(plot_script)
    println("   找到绘图脚本: $plot_script")
    try
        include(plot_script)
        println("✅ 绘图脚本加载成功")
    catch e
        println("⚠️  绘图脚本加载失败: $e")
        println("   这可能是因为缺少特定的数据文件")
    end
else
    println("   ⚠️  未找到项目绘图脚本")
end

# 清理测试文件
println("\n🧹 清理测试文件...")
for file in ["test_plot.png", "test_data.csv"]
    if isfile(file)
        rm(file)
        println("   删除: $file")
    end
end

println("\n" ^ 60)
println("🎉 可视化功能测试完成!")
println("=" ^ 60)

println("\n💡 可视化功能使用说明:")
println("1. 基本绘图: using Plots; plot(x, y)")
println("2. 数据处理: using DataFrames, CSV")
println("3. 项目绘图: include(\"scripts/Gas_Liquid/plot_temperature_scan.jl\")")
println("4. 后端选择: gr(), plotlyjs(), pyplot()")

println("\n📁 相关文件:")
println("   - scripts/Gas_Liquid/plot_temperature_scan.jl")
println("   - output/*.csv (数据文件)")
println("=" ^ 60)