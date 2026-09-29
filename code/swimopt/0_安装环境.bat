@echo off
chcp 65001 >nul
cd /d "%~dp0"
rem If a venv already exists one level up (original SA workspace layout), reuse it.
rem Otherwise create a self-contained venv inside this folder (standalone git clone).
if exist "..\mjenv\Scripts\activate.bat" cd /d "%~dp0.."
title 0. Setup Python environment
echo ============================================================
echo   Setting up MuJoCo environment in:  %CD%
echo ============================================================
echo.
python --version
if errorlevel 1 (
  echo [X] Python not found. Install Python 3.10-3.12 from python.org
  echo     and CHECK "Add Python to PATH" during installation.
  pause
  exit /b 1
)
if not exist "mjenv" (
  echo Creating virtual environment mjenv ...
  python -m venv mjenv
)
call "mjenv\Scripts\activate.bat"
echo Installing mujoco / cma / numpy / matplotlib ...
python -m pip install --upgrade pip
python -m pip install mujoco cma numpy matplotlib
echo.
python -c "import mujoco,cma;print('OK  mujoco',mujoco.__version__,' cma',cma.__version__)"
echo.
echo ============================================================
echo   Done. Now double-click  swimopt\1_看示例步态.bat
echo ============================================================
pause >nul
