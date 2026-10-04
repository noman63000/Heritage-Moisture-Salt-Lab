@echo off
setlocal
cd /d "%~dp0"
where py >nul 2>nul
if %errorlevel%==0 (
    py -3 analysis_bridge.py
) else (
    python analysis_bridge.py
)
if errorlevel 1 pause
