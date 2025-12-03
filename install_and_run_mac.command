#!/bin/bash

# Get the directory where this script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

echo "========================================================"
echo "   RNA-seq Pipeline - Auto Installer & Launcher"
echo "========================================================"

# Define paths
MAMBA_EXE="$DIR/bin/micromamba"
ENV_DIR="$DIR/env"

# 1. Check if Micromamba is installed locally
if [ ! -f "$MAMBA_EXE" ]; then
    echo "📦 Micromamba not found. Downloading..."
    mkdir -p "$DIR/bin"
    
    # Detect architecture (Intel vs Apple Silicon)
    ARCH=$(uname -m)
    if [ "$ARCH" == "arm64" ]; then
        URL="https://micro.mamba.pm/api/micromamba/osx-arm64/latest"
    else
        URL="https://micro.mamba.pm/api/micromamba/osx-64/latest"
    fi
    
    curl -Ls "$URL" | tar -xvj -C "$DIR/bin" bin/micromamba
    mv "$DIR/bin/bin/micromamba" "$DIR/bin/micromamba"
    rmdir "$DIR/bin/bin"
    
    echo "✅ Micromamba downloaded."
else
    echo "✅ Micromamba found."
fi

# 2. Check if environment exists
if [ ! -d "$ENV_DIR" ]; then
    echo "⚙️  Creating Python environment (This may take a few minutes)..."
    "$MAMBA_EXE" create -p "$ENV_DIR" -f environment.yml -y
    
    if [ $? -eq 0 ]; then
        echo "✅ Environment created successfully."
    else
        echo "❌ Failed to create environment."
        exit 1
    fi
else
    echo "✅ Environment found."
fi

# 3. Launch App
echo "🚀 Starting App..."
"$ENV_DIR/bin/streamlit" run app.py
