#!/bin/bash

# 取得腳本所在的目錄，確保在正確的資料夾執行
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

echo "================================================="
echo "   Bulk RNA-seq Analysis Pipeline 啟動器"
echo "================================================="

# 1. 檢查 Python 是否安裝
if ! command -v python3 &> /dev/null; then
    echo "❌ 錯誤: 未偵測到 Python 3。"
    echo "請先前往 https://www.python.org/downloads/ 安裝 Python 3。"
    read -p "按 Enter 鍵退出..."
    exit 1
fi

# 2. 建立虛擬環境 (如果不存在)
if [ ! -d ".venv" ]; then
    echo "🔨 正在建立虛擬環境 (.venv)..."
    python3 -m venv .venv
    if [ $? -ne 0 ]; then
        echo "❌ 建立虛擬環境失敗。"
        read -p "按 Enter 鍵退出..."
        exit 1
    fi
else
    echo "✅ 偵測到現有虛擬環境。"
fi

# 3. 啟動虛擬環境
source .venv/bin/activate

# 4. 安裝/更新 套件
echo "📦 正在檢查並安裝必要套件..."
pip install --upgrade pip
pip install -r requirements.txt
if [ $? -ne 0 ]; then
    echo "❌ 套件安裝失敗。"
    read -p "按 Enter 鍵退出..."
    exit 1
fi

# 5. 執行 Streamlit 應用程式
echo "🚀 正在啟動 GUI..."
echo "網頁瀏覽器將自動開啟。若要關閉，請關閉此視窗或按 Ctrl+C。"
echo "-------------------------------------------------"

streamlit run app.py

# 保持視窗開啟 (如果程式崩潰)
read -p "程式已結束，按 Enter 鍵關閉視窗..."
