# PNJL Python贝叶斯优化文件结构说明

## 📁 文件组织结构

### 🚀 **演示脚本** - `scripts/Gas_Liquid/`
```
scripts/Gas_Liquid/
└── demo_python_bayesian.py    # Python贝叶斯优化演示脚本
```

**运行方式**：
```bash
cd d:/Desktop/Julia/Rotation_PNJL
python scripts/Gas_Liquid/demo_python_bayesian.py
```

### 🧪 **测试脚本** - `test/Gas_Liquid/`
```
test/Gas_Liquid/
└── test_python_bayesian.py    # Python贝叶斯优化功能测试
```

**运行方式**：
```bash
cd d:/Desktop/Julia/Rotation_PNJL
python test/Gas_Liquid/test_python_bayesian.py
```

### 📂 **输出目录** - `output/Gas_Liquid/`
```
output/Gas_Liquid/
├── demo_optimization_result.csv      # 演示结果CSV文件
├── demo_optimization_result.png      # 演示收敛图
├── test_optimization_result.csv      # 测试结果CSV文件
├── optimization_history.png          # 优化历史图
└── *.config.json                     # 配置文件
```

### 📚 **核心模块** - `src/Gas_Liquid/`
```
src/Gas_Liquid/
├── Advanced_Bayesian.py              # 完整的Python贝叶斯优化库
├── Advanced_BayesianOptimization.jl  # Julia版本
└── Advanced_FindTforDiff.jl           # Julia目标函数
```

## 🔧 主要改进

### ✅ **修复的问题**
1. **文件组织**：按功能分类到不同目录
2. **输出路径**：统一输出到 `output/Gas_Liquid/`
3. **字体问题**：修复matplotlib中文显示
4. **路径自适应**：自动检测项目根目录

### 🎯 **字体显示优化**
- 自动检测并使用系统中可用的中文字体
- 支持的字体列表：
  - SimHei (黑体)
  - Microsoft YaHei (微软雅黑)
  - Arial Unicode MS
- 如果没有中文字体，自动切换到英文标签

### 📊 **输出目录自动管理**
- 所有CSV文件自动保存到 `output/Gas_Liquid/`
- 所有PNG图片自动保存到 `output/Gas_Liquid/`
- 自动创建不存在的目录
- 支持绝对路径和相对路径

## 🚀 使用指南

### **快速开始**
```bash
# 1. 运行演示
cd d:/Desktop/Julia/Rotation_PNJL
python scripts/Gas_Liquid/demo_python_bayesian.py

# 2. 查看结果
# CSV文件: output/Gas_Liquid/demo_optimization_result.csv
# 图片: output/Gas_Liquid/demo_optimization_result.png
```

### **功能测试**
```bash
# 运行完整测试套件
python test/Gas_Liquid/test_python_bayesian.py
```

### **使用完整版本**
```python
from src.Gas_Liquid.Advanced_Bayesian import PNJLBayesianOptimizer

# 创建优化器
optimizer = PNJLBayesianOptimizer()

# 运行优化
result = optimizer.optimize_with_warmup(...)

# 结果自动保存到 output/Gas_Liquid/
```

## 📈 **性能与稳定性**

### **字体兼容性**
- ✅ Windows：支持微软雅黑、黑体
- ✅ macOS：支持Arial Unicode MS
- ✅ Linux：支持文泉驿字体
- ✅ 备用方案：英文标签

### **路径兼容性**
- ✅ 自动检测项目根目录
- ✅ 支持不同操作系统路径分隔符
- ✅ 自动创建输出目录
- ✅ 处理相对路径和绝对路径

### **内存管理**
- ✅ 自动关闭matplotlib图形
- ✅ 使用非交互式后端
- ✅ 优化大文件输出

## 🔄 **迁移说明**

### **从旧版本迁移**
1. 原来的 `scripts/demo_python_bayesian.py` → `scripts/Gas_Liquid/demo_python_bayesian.py`
2. 原来的 `scripts/test_python_bayesian.py` → `test/Gas_Liquid/test_python_bayesian.py`
3. 所有输出文件 → `output/Gas_Liquid/`

### **兼容性保证**
- 所有API保持向后兼容
- 原有的命令行参数不变
- 配置文件格式不变

## 🎉 **新功能特性**

1. **智能字体检测**：自动选择最佳中文字体
2. **统一输出管理**：所有结果集中管理
3. **增强错误处理**：更好的错误信息和恢复
4. **模块化设计**：清晰的文件组织结构
5. **自动路径处理**：无需手动配置路径

---

**升级完成！** 🎊 现在您可以享受更好的文件组织和无字体问题的可视化体验。