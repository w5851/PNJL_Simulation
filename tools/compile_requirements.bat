@echo off
REM Batch wrapper to run the PowerShell compile_requirements.ps1 script
REM Usage: compile_requirements.bat

n:: Resolve script directory
nsetlocal enabledelayedexpansion
set "SCRIPT_DIR=%~dp0"
:: Navigate to project root (parent of tools)
cd /d "%SCRIPT_DIR%.."
:: Call PowerShell script using pwsh if available, otherwise fallback to powershell
where pwsh >nul 2>&1
if %ERRORLEVEL%==0 (
  pwsh -NoProfile -ExecutionPolicy Bypass -File "%SCRIPT_DIR%compile_requirements.ps1"
  set "RC=%ERRORLEVEL%"
) else (
  powershell -NoProfile -ExecutionPolicy Bypass -File "%SCRIPT_DIR%compile_requirements.ps1"
  set "RC=%ERRORLEVEL%"
)
endlocal & exit /b %RC%
