"""
绘制温度扫描结果
使用CSV文件中的数据绘制kappa3/kappa1和kappa4/kappa2随温度T_MeV的变化
"""

# 激活项目环境
import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using Plots, CSV, DataFrames

# 设置绘图后端，避免中文字符问题
gr() # 使用GR后端

function read_csv_with_metadata(csv_file::String)
    """
    读取带有元数据头部的CSV文件，跳过以#开头的注释行
    返回DataFrame和元数据字典
    """
    metadata = Dict{String, String}()
    
    # 读取所有行并提取元数据
    lines = readlines(csv_file)
    for line in lines
        if startswith(line, "#") && occursin("=", line)
            parts = split(replace(line, "#" => ""), "=")
            if length(parts) == 2
                key = strip(parts[1])
                value = strip(parts[2])
                metadata[key] = value
            end
        end
    end
    
    # 创建临时文件，只包含数据部分（非注释行）
    temp_lines = filter(line -> !startswith(line, "#") && !isempty(strip(line)), lines)
    
    if isempty(temp_lines)
        error("CSV文件中没有找到数据行")
    end
    
    # 写入临时文件
    temp_file = tempname() * ".csv"
    open(temp_file, "w") do io
        for line in temp_lines
            println(io, line)
        end
    end
    
    # 读取临时文件
    df = CSV.read(temp_file, DataFrame)
    
    # 删除临时文件
    rm(temp_file)
    
    return df, metadata
end

function plot_temperature_scan_results(csv_file::String, save_path::String="")
    """
    从CSV文件读取数据并绘制kappa比值随温度的变化
    
    参数:
    - csv_file: CSV文件路径
    - save_path: 图像保存路径（可选），如果提供则保存图像
    """
    
    # 读取CSV数据和元数据
    println("正在读取数据文件: $csv_file")
    df, metadata = read_csv_with_metadata(csv_file)
    
    # 显示找到的元数据
    if !isempty(metadata)
        println("发现运行参数元数据:")
        for (key, value) in metadata
            println("  $key = $value")
        end
    end
    
    # 检查数据列
    required_cols = ["T_MeV", "kappa3_over_kappa1", "kappa4_over_kappa2"]
    for col in required_cols
        if !(col in names(df))
            error("CSV文件中缺少必要的列: $col")
        end
    end
    
    # 过滤掉NaN值
    valid_idx = .!isnan.(df.kappa3_over_kappa1) .& .!isnan.(df.kappa4_over_kappa2)
    df_filtered = df[valid_idx, :]
    
    println("有效数据点数: $(nrow(df_filtered))")
    
    # 创建绘图标题，包含参数信息
    title_text = "Kappa Ratios vs Temperature"
    if haskey(metadata, "fs") && haskey(metadata, "mu_B")
        title_text *= "\n(fs=$(metadata["fs"]), μB=$(metadata["mu_B"]))"
    end
    
    # 创建绘图 - 使用英文标签避免字体问题
    p = plot(
        title=title_text,
        xlabel="Temperature T (MeV)",
        ylabel="Kappa Ratios",
        size=(800, 600),
        legend=:topright,
        grid=true,
        gridwidth=1,
        gridcolor=:lightgray,
        dpi=300  # 高分辨率
    )
    
    # 绘制kappa3/kappa1
    plot!(p, df_filtered.T_MeV, df_filtered.kappa3_over_kappa1,
          label="κ₃/κ₁",
          linewidth=2,
          color=:blue,
          marker=:circle,
          markersize=3)
    
    # 绘制kappa4/kappa2
    plot!(p, df_filtered.T_MeV, df_filtered.kappa4_over_kappa2,
          label="κ₄/κ₂",
          linewidth=2,
          color=:red,
          marker=:square,
          markersize=3)
    
    # 添加水平参考线(比值=1)
    hline!(p, [1.0], 
           label="Ratio = 1",
           color=:black,
           linestyle=:dash,
           linewidth=1)
    
    # 保存图像到文件（如果指定了路径）
    if !isempty(save_path)
        savefig(p, save_path)
        println("图像已保存到: $save_path")
    end
    
    # 在SSH环境下不显示图形，只保存
    if isempty(save_path)
        display(p)
    end
    
    # 打印一些统计信息
    println("\n数据统计:")
    println("温度范围: $(minimum(df_filtered.T_MeV)) - $(maximum(df_filtered.T_MeV)) MeV")
    println("κ₃/κ₁ 范围: $(minimum(df_filtered.kappa3_over_kappa1)) - $(maximum(df_filtered.kappa3_over_kappa1))")
    println("κ₄/κ₂ 范围: $(minimum(df_filtered.kappa4_over_kappa2)) - $(maximum(df_filtered.kappa4_over_kappa2))")
    
    return p
end

function plot_individual_kappas(csv_file::String, save_path::String="")
    """
    分别绘制各个kappa值随温度的变化（可选功能）
    """
    
    df, metadata = read_csv_with_metadata(csv_file)
    
    println("原始数据点数: $(nrow(df))")
    
    # 检查每个kappa列的有效数据
    κ1_valid = sum(.!isnan.(df.kappa1))
    κ2_valid = sum(.!isnan.(df.kappa2))
    κ3_valid = sum(.!isnan.(df.kappa3))
    κ4_valid = sum(.!isnan.(df.kappa4))
    
    println("各kappa有效数据点:")
    println("  κ₁: $κ1_valid 个")
    println("  κ₂: $κ2_valid 个")
    println("  κ₃: $κ3_valid 个")
    println("  κ₄: $κ4_valid 个")
    
    # 分别为每个kappa过滤有效数据
    df_κ1 = df[.!isnan.(df.kappa1), :]
    df_κ2 = df[.!isnan.(df.kappa2), :]
    df_κ3 = df[.!isnan.(df.kappa3), :]
    df_κ4 = df[.!isnan.(df.kappa4), :]
    
    # 检查是否有足够的数据绘图
    if nrow(df_κ1) == 0 && nrow(df_κ2) == 0 && nrow(df_κ3) == 0 && nrow(df_κ4) == 0
        println("警告：所有kappa值都是NaN，无法绘制图形")
        return nothing
    end
    
    p = plot(
        title="Individual Kappa Values vs Temperature",
        xlabel="Temperature T (MeV)",
        ylabel="Kappa Values",
        size=(800, 600),
        legend=:topright,
        grid=true,
        yscale=:log10,  # 使用对数坐标，因为kappa值可能变化范围很大
        dpi=300
    )
    
    # 分别绘制每个有效的kappa
    if nrow(df_κ1) > 0
        plot!(p, df_κ1.T_MeV, abs.(df_κ1.kappa1), label="κ₁", linewidth=2, color=:blue, marker=:circle, markersize=2)
    end
    
    if nrow(df_κ2) > 0
        plot!(p, df_κ2.T_MeV, abs.(df_κ2.kappa2), label="κ₂", linewidth=2, color=:red, marker=:square, markersize=2)
    end
    
    if nrow(df_κ3) > 0
        plot!(p, df_κ3.T_MeV, abs.(df_κ3.kappa3), label="κ₃", linewidth=2, color=:green, marker=:diamond, markersize=2)
    end
    
    if nrow(df_κ4) > 0
        plot!(p, df_κ4.T_MeV, abs.(df_κ4.kappa4), label="κ₄", linewidth=2, color=:orange, marker=:pentagon, markersize=2)
    end
    
    # 保存图像到文件（如果指定了路径）
    if !isempty(save_path)
        savefig(p, save_path)
        println("图像已保存到: $save_path")
    end
    
    # 在SSH环境下不显示图形，只保存
    if isempty(save_path)
        display(p)
    end
    
    return p
end

# 主执行部分
function main(csv_file::String="")
    """
    主函数：绘制温度扫描结果图
    
    参数:
    - csv_file: CSV文件路径（可选）
        - 如果提供，使用指定的文件路径
        - 如果未提供或为空字符串，则：
          1. 优先使用命令行参数 ARGS[1]
          2. 其次使用默认路径 "../../output/Gas_Liquid/forwarddiff_temperature_scan.csv"
    """
    
    # 获取脚本所在目录
    script_dir = dirname(@__FILE__)
    
    # 确定CSV文件路径的优先级：
    # 1. 函数参数 csv_file（如果非空）
    # 2. 命令行参数 ARGS[1]（如果存在）
    # 3. 默认路径
    if !isempty(csv_file)
        # 使用函数参数
        target_csv = csv_file
        println("使用函数参数指定的CSV文件: $target_csv")
    elseif length(ARGS) > 0
        # 使用命令行参数
        target_csv = ARGS[1]
        println("使用命令行参数指定的CSV文件: $target_csv")
    else
        # 使用默认路径
        target_csv = joinpath(script_dir, "../../output/Gas_Liquid/forwarddiff_temperature_scan.csv")
        println("使用默认CSV文件路径: $target_csv")
    end
    
    # 如果提供的是相对路径，则相对于脚本目录
    if !isabspath(target_csv)
        target_csv = joinpath(script_dir, target_csv)
    end
    
    csv_file = target_csv
    
    # 输出图像路径
    output_dir = joinpath(script_dir, "../../output/Gas_Liquid")
    ratio_plot_path = joinpath(output_dir, "kappa_ratios_temperature_scan.png")
    individual_plot_path = joinpath(output_dir, "individual_kappas_temperature_scan.png")
    
    # 检查输出目录是否存在，如果不存在则创建
    if !isdir(output_dir)
        mkpath(output_dir)
        println("创建输出目录: $output_dir")
    end
    
    # 检查文件是否存在
    if !isfile(csv_file)
        println("错误：找不到CSV文件: $csv_file")
        println("当前工作目录: $(pwd())")
        println("脚本所在目录: $script_dir")
        println("请检查文件路径是否正确")
        error("找不到CSV文件: $csv_file")
    end
    
    println("开始绘制温度扫描结果...")
    println("注意：在SSH环境下，图像将保存为PNG文件而不是直接显示")
    
    try
        # 绘制主要的比值图并保存
        println("\n1. 绘制Kappa比值图...")
        p1 = plot_temperature_scan_results(csv_file, ratio_plot_path)
        
        # 可选：绘制各个kappa值并保存
        println("\n2. 绘制各个Kappa值图...")
        p2 = plot_individual_kappas(csv_file, individual_plot_path)
        
        println("\n绘图完成！")
        println("图像文件保存位置:")
        println("- Kappa比值图: $ratio_plot_path")
        println("- 各个Kappa值图: $individual_plot_path")
        println("\n您可以通过以下方式查看图像:")
        println("1. 使用SCP/SFTP下载图像文件到本地查看")
        println("2. 在本地使用: scp user@server:$(ratio_plot_path) .")
        println("3. 或者使用VS Code的文件浏览器直接查看PNG文件")
        
    catch e
        println("绘图过程中出现错误: $e")
        println("请检查CSV文件格式和数据完整性")
    end
end

# 自动执行main函数的逻辑
# 方法1：检查是否是直接运行脚本（命令行方式）
if abspath(PROGRAM_FILE) == @__FILE__
    println("检测到命令行运行，自动执行main()...")
    main()
# 方法2：检查是否在REPL中include此文件且没有命令行参数
elseif isinteractive() && length(ARGS) == 0
    println("检测到REPL环境，自动执行main()...")
    println("如果不希望自动执行，请使用: include(\"plot_temperature_scan.jl\"); # 然后手动调用main()")
    #main("../../output/Gas_Liquid/forwarddiff_optimization_params_scan.csv")
    main("../../output/Gas_Liquid/forwarddiff_temperature_scan.csv")
    
else
    println("脚本已加载。请手动调用 main() 或 main(\"csv文件路径\") 来绘制图形。")
    println("示例:")
    println("  main()  # 使用默认路径")
    println("  main(\"../../output/Gas_Liquid/your_file.csv\")  # 使用指定路径")
end
