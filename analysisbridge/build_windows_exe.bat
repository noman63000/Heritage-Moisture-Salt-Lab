@echo off
setlocal
cd /d "%~dp0"
where py >nul 2>nul
if %errorlevel%==0 (
    set PY=py -3
) else (
    set PY=python
)
%PY% -m pip install --upgrade pyinstaller
if errorlevel 1 goto :fail
%PY% -m PyInstaller --noconfirm --clean --onefile --windowed --name AnalysisBridge analysis_bridge.py
if errorlevel 1 goto :fail
echo.
echo Build complete: %CD%\dist\AnalysisBridge.exe
pause
exit /b 0
:fail
echo.
echo Build failed. Review the messages above.
pause
exit /b 1
