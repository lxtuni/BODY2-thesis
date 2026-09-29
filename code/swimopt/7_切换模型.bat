@echo off
chcp 65001 >nul
set PYTHONIOENCODING=utf-8
cd /d "%~dp0"
title Switch model
echo ============================================================
echo   Switch which robot model to use
echo ============================================================
echo.
if exist "..\mjenv\Scripts\activate.bat" call "..\mjenv\Scripts\activate.bat"
python set_model.py
echo.
pause >nul
