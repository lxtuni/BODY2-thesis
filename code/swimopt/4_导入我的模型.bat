@echo off
chcp 65001 >nul
set PYTHONIOENCODING=utf-8
cd /d "%~dp0"
title Import your URDF model
echo ============================================================
echo   Import URDF into the pipeline
echo ============================================================
echo.
if exist "..\mjenv\Scripts\activate.bat" call "..\mjenv\Scripts\activate.bat"
set /p URDF=Drag your .urdf file here then press Enter: 
set URDF=%URDF:"=%
if "%URDF%"=="" (echo No file given. & pause & exit /b 1)
set /p NAME=Output name (e.g. body3) [default: mymodel]: 
if "%NAME%"=="" set NAME=mymodel
echo.
python import_model.py "%URDF%" "robots/%NAME%.xml" --drive "2.1,1.1" --armature 3e-4 --damping 0.02 --kp 15 --force 3
echo.
echo ------------------------------------------------------------
echo  Next: open config.json, set  "model": "robots/%NAME%.xml"
echo ------------------------------------------------------------
pause >nul
