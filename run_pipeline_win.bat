@echo off
setlocal

echo =================================================
echo    Bulk RNA-seq Analysis Pipeline Launcher
echo =================================================

REM 1. Check for Python (Try python, py, then python3)
set PYTHON_CMD=python

python --version >nul 2>&1
if %errorlevel% equ 0 goto :FOUND_PYTHON

set PYTHON_CMD=py
py --version >nul 2>&1
if %errorlevel% equ 0 goto :FOUND_PYTHON

set PYTHON_CMD=python3
python3 --version >nul 2>&1
if %errorlevel% equ 0 goto :FOUND_PYTHON

echo [ERROR] Python is not detected.
echo Please install Python 3 from https://www.python.org/downloads/
echo IMPORTANT: Check "Add Python to PATH" during installation.
pause
exit /b

:FOUND_PYTHON
echo [INFO] Using Python command: %PYTHON_CMD%


REM 2. Create Virtual Environment if not exists
if not exist ".venv" (
    echo [INFO] Creating virtual environment (.venv)...
    %PYTHON_CMD% -m venv .venv
    if %errorlevel% neq 0 (
        echo [ERROR] Failed to create virtual environment.
        pause
        exit /b
    )
) else (
    echo [INFO] Found existing virtual environment.
)

REM 3. Activate Virtual Environment
call .venv\Scripts\activate

REM 4. Install/Update Requirements
echo [INFO] Checking and installing requirements...
pip install --upgrade pip
pip install -r requirements.txt
if %errorlevel% neq 0 (
    echo [ERROR] Failed to install requirements.
    pause
    exit /b
)

REM 5. Run Streamlit App
echo [INFO] Launching GUI...
echo The browser should open automatically. Close this window to stop the app.
echo -------------------------------------------------

streamlit run app.py

pause
