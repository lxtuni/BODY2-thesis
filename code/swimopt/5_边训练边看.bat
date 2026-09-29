@echo off
chcp 65001 >nul
set PYTHONIOENCODING=utf-8
cd /d "%~dp0"
title Optimize with live view (every candidate)
echo ============================================================
echo   Optimize with live view (every candidate)
echo ============================================================
echo.
if exist "..\mjenv\Scripts\activate.bat" (
  call "..\mjenv\Scripts\activate.bat"
) else if exist "mjenv\Scripts\activate.bat" (
  call "mjenv\Scripts\activate.bat"
) else (
  echo [!] Virtual env not found -- using system python
)
python -c "import mujoco" 2>nul
if errorlevel 1 (echo [X] mujoco not installed. Run the "0_" setup .bat first. & pause & exit /b 1)
python optimize_view.py config.json 250
echo.
echo   Finished. Press any key to close.
pause >nul
