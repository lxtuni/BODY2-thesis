@echo off
chcp 65001 >nul
set PYTHONIOENCODING=utf-8
cd /d "%~dp0"
title swimopt Control Panel
if exist "..\mjenv\Scripts\activate.bat" (
  call "..\mjenv\Scripts\activate.bat"
) else if exist "mjenv\Scripts\activate.bat" (
  call "mjenv\Scripts\activate.bat"
)
python -c "import mujoco" 2>nul
if errorlevel 1 (echo [X] mujoco not installed. Run the "0_" setup .bat first. & pause & exit /b 1)
start "" pythonw ui.py
if errorlevel 1 python ui.py
exit
