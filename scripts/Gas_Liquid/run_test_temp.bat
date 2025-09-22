@echo off
set "JULIA_DEPOT_PATH=C:\Users\82306\.julia"
cd /d "C:\Users\82306\Desktop\Julia\Julia\Rotation_PNJL\test\Gas_Liquid"
"C:\ProgramData\chocolatey\bin\julia.exe" "C:\Users\82306\Desktop\Julia\Julia\Rotation_PNJL\test\Gas_Liquid\test_temp.jl" > "C:\Users\82306\Desktop\Julia\Julia\Rotation_PNJL\test\Gas_Liquid\test_temp.log" 2>&1


#运行任务：
#schtasks /Run /TN RunJulia_test_temp