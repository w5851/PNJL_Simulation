# run_project.jl
# Rotation_PNJL项目统一运行器
# 提供项目中各种功能的统一入口点

# 激活项目环境
println("🚀 激活Rotation_PNJL项目环境...")
import Pkg
Pkg.activate(".")

println("📦 项目环境已激活: $(Pkg.project().path)")
println("=" ^ 70)

# 菜单函数
function show_menu()
    println("🎯 Rotation_PNJL 项目功能菜单")
    println("=" ^ 50)
    println("📊 测试功能:")
    println("  1. 运行贝叶斯优化测试")
    println("  2. 运行可视化功能测试")
    println("  3. 运行气液相变温度查找测试")
    println("  4. 运行ForwardDiff温度扫描测试")
    println()
    println("🎨 演示功能:")
    println("  5. 气液相变温度查找示例")
    println("  6. 绘制温度扫描结果")
    println()
    println("🔧 维护功能:")
    println("  7. 重新安装项目环境")
    println("  8. 检查包状态")
    println("  9. 运行所有测试")
    println()
    println("  0. 退出")
    println("=" ^ 50)
    print("请选择功能 (0-9): ")
end

# 执行选择的功能
function run_selection(choice)
    try
        if choice == "1"
            println("\n🧠 运行贝叶斯优化测试...")
            include("test/test_bayesian_fixed.jl")
        elseif choice == "2"
            println("\n🎨 运行可视化功能测试...")
            include("test/run_visualization_test.jl")
        elseif choice == "3"
            println("\n🌡️ 运行气液相变温度查找测试...")
            include("test/Gas_Liquid/test_find_temperature.jl")
        elseif choice == "4"
            println("\n📈 运行ForwardDiff温度扫描测试...")
            include("test/Gas_Liquid/test_forwarddiff_temperature_scan.jl")
        elseif choice == "5"
            println("\n💡 运行气液相变温度查找示例...")
            include("scripts/Gas_Liquid/example_find_temperature.jl")
        elseif choice == "6"
            println("\n📊 绘制温度扫描结果...")
            include("scripts/Gas_Liquid/plot_temperature_scan.jl")
        elseif choice == "7"
            println("\n🔧 重新安装项目环境...")
            include("install.jl")
        elseif choice == "8"
            println("\n📋 检查包状态...")
            Pkg.status()
        elseif choice == "9"
            println("\n🧪 运行所有测试...")
            test_files = [
                "test/test_bayesian_fixed.jl",
                "test/run_visualization_test.jl", 
                "test/Gas_Liquid/test_find_temperature.jl",
                "test/Gas_Liquid/test_forwarddiff_temperature_scan.jl"
            ]
            for test_file in test_files
                if isfile(test_file)
                    println("\n" * "="^50)
                    println("运行: $test_file")
                    println("="^50)
                    try
                        include(test_file)
                        println("✅ $test_file 完成")
                    catch e
                        println("❌ $test_file 失败: $e")
                    end
                else
                    println("⚠️  文件不存在: $test_file")
                end
            end
        elseif choice == "0"
            println("\n👋 感谢使用Rotation_PNJL项目!")
            return false
        else
            println("\n❌ 无效选择，请输入0-9之间的数字")
        end
    catch e
        println("❌ 执行失败: $e")
        println("💡 提示: 确保已运行 'julia install.jl' 安装所有依赖")
    end
    return true
end

# 主程序循环
function main()
    println("🎉 欢迎使用Rotation_PNJL项目!")
    println("📁 当前目录: $(pwd())")
    
    # 检查项目环境
    if !isfile("Project.toml")
        println("❌ 错误: 请在项目根目录运行此脚本")
        return
    end
    
    while true
        println("\n" ^ 70)
        show_menu()
        
        choice = strip(readline())
        
        if !run_selection(choice)
            break
        end
        
        if choice != "0"
            println("\n⏸️  按Enter键继续...")
            readline()
        end
    end
end

# 如果直接运行此脚本，启动主程序
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end