#!/bin/bash
echo "=========================================="
echo "Starting Cheminformatics Engine..."
echo "=========================================="

# Check if Python is installed
if ! command -v python3 &> /dev/null
then
    echo "[ERROR] Python3 is not installed. Please install Python 3.10+."
    exit
fi

# Create Virtual Environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo "[INFO] Creating Python virtual environment..."
    python3 -m venv venv
fi

# Activate the virtual environment
echo "[INFO] Activating virtual environment..."
source venv/bin/activate

# Install Requirements
echo "[INFO] Installing required dependencies..."
python -m pip install --upgrade pip > /dev/null
pip install -r requirements.txt

# Start the Application
echo "[INFO] Starting the server..."
python main.py
