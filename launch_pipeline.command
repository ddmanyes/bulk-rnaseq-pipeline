#!/bin/bash

# 1. Initialize Conda (Adjust path if Anaconda is installed elsewhere)
# Common paths: /opt/anaconda3, /usr/local/anaconda3, ~/anaconda3, ~/opt/anaconda3
CONDA_PATH="/opt/anaconda3"

if [ -f "$CONDA_PATH/etc/profile.d/conda.sh" ]; then
    source "$CONDA_PATH/etc/profile.d/conda.sh"
else
    echo "Conda not found at $CONDA_PATH. Please edit this script with the correct path."
    exit 1
fi

# 2. Activate the environment
conda activate omicverse_env

# 3. Navigate to the script's directory
cd "$(dirname "$0")"

# 4. Run the Streamlit app
echo "Starting RNA-seq Pipeline..."
streamlit run app.py
