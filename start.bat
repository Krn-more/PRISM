@echo off
echo ==========================================
echo Starting Cheminformatics Engine...
echo ==========================================

:: Check if Python is installed
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Python is not installed or not in your PATH. Please install Python 3.10+.
    pause
    exit /b
)

:: Create Virtual Environment if it doesn't exist
if not exist "venv\Scripts\activate.bat" (
    echo [INFO] Creating Python virtual environment...
    python -m venv venv
)

:: Activate the virtual environment
echo [INFO] Activating virtual environment...
call venv\Scripts\activate.bat

:: Install Requirements
echo [INFO] Installing required dependencies...
python -m pip install --upgrade pip >nul
pip install -r requirements.txt

:: Start the Application
echo [INFO] Starting the server...
python main.py

pause
