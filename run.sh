#!/bin/bash

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Path to your virtual environment's activate script
VENV_ACTIVATE_PATH="$SCRIPT_DIR/venv/bin/activate"

# Check if the virtual environment exists
if [ ! -f "$VENV_ACTIVATE_PATH" ]; then
    echo "Error: Virtual environment not found at $VENV_ACTIVATE_PATH"
    exit 1
fi

# Activate the virtual environment
source "$VENV_ACTIVATE_PATH"

# Run the Python file
python3 "$SCRIPT_DIR/main.py"

# Deactivate the virtual environment after running the Python file
deactivate
