# Bulk RNA-seq 分析流程與互動式瀏覽器

本專案提供了一個全面且自動化的 Bulk RNA-seq 分析流程，並搭配一個互動式網頁應用程式來視覺化分析結果。

此流程從原始的 `featureCounts` 矩陣開始，執行差異基因表現分析 (DGE)、富集分析 (Enrichment Analysis) 和時間序列分群 (Time-Series Clustering)，並生成一套完整的發表級圖表和表格。整個工作流程都可以透過單一個 `config.yml` 檔案輕鬆控制。

## 主要功能

- **使用者工作階段隔離 (Session Isolation)**：為每個使用者建立獨立的暫存目錄，確保在共享主機（如 Hugging Face Spaces）上的數據隱私。
- **可自訂元數據 (Customizable Metadata)**：支援自訂時間欄位名稱（如 "Stage", "Day"）和非數字的時間點。
- **自動同步元數據 (Auto-Sync Metadata)**：一鍵從計數矩陣同步 Sample IDs 到元數據表格。
- **端到端分析**：只需一個指令，即可從原始計數數據產生最終圖表。
- **主成分分析 (PCA)**：生成 PCA 圖，根據整體基因表現視覺化樣本間的關係和批次效應。
- **批次效應校正 (Batch Correction)**：整合 ComBat 批次校正功能，可移除不必要的技術變異，並自動生成「校正前後」的比較圖。
- **豐富的視覺化**：自動生成火山圖 (Volcano Plots)、聚類熱圖 (Clustered Heatmaps，針對顯著基因和高變異基因) 以及時間序列趨勢圖。
- **富集分析**：針對多個資料庫 (如 GO, KEGG, Reactome) 執行過度表現分析 (ORA)，並生成結果的點圖 (Dot Plots)。
- **時間序列分群**：使用 K-means 分群法識別時間序列數據中共同表現的基因模組，並視覺化其趨勢。
- **完全可配置**：在一個中央 `config.yml` 檔案中輕鬆管理所有檔案路徑、分析參數和閾值。
- **互動式瀏覽器**：包含一個 Streamlit 網頁應用程式，可透過箱形圖 (Boxplots) 動態探索基因表現。
- **數據與圖片導出**：可直接從應用程式下載處理後的數據、差異表現結果、通路分析表格，以及**高品質的圖表圖片 (PNG)**。
- **Mac 一鍵安裝**：內建 macOS 專用的一鍵安裝檔 (`install_and_run_mac.command`)，無需手動設定 Python 環境。
- **可重現性**：附帶 `requirements.txt` 檔案，確保環境的一致性。
- **Kallisto 整合**：包含一個輔助腳本，可將 Kallisto 的轉錄本豐度 (transcript abundance) 聚合為與此流程相容的基因級計數矩陣 (gene-level count matrix)。

## 專案結構

```
├── input/
│   ├── featureCount files                   # 您的原始 featureCounts 數據
│   └── metadata.csv                         # 樣本元數據 (Condition, Time, Type)
├── genesets/
│   └── pair_GRCm39.tsv                      # 基因註釋檔案
├── output/                                  # 所有生成的結果都儲存在這裡
│   ├── PCA_plot.png                         # 主成分分析圖
│   ├── PCA_Batch_Correction_Comparison.png  # 批次校正前後的 PCA 比較圖
│   ├── counts.csv                           # 預處理並映射基因 ID 後的計數
│   ├── counts_batch_corrected.csv           # 經批次校正後的計數矩陣
│   ├── log2cpm_normalized_counts.csv        # Log2CPM 標準化計數
│   ├── DEG_*.csv                            # 每個比較的 DGE 結果
│   ├── Volcano_*.png                        # 每個比較的火山圖
│   ├── Heatmap_Significant_Genes.png        # 所有顯著基因的熱圖
│   ├── Heatmap_Top50_Variable_Genes.png     # 前 N 個高變異基因的熱圖
│   ├── Master_Results_Table.csv             # 結合 DGE、分群和表現數據的總表
│   ├── Master_Results_Table.xlsx            # 結合 DGE、分群和表現數據的總表
│   ├── TimeSeries_Clustering_k4.png         # 時間序列分群圖
│   ├── k_evaluation_plot.png                # 選擇 k 值的輔助圖 (可選)
│   └── enrichment/                          # GO/KEGG/Reactome 結果 (CSV 和圖表)
├── time_series_enrichment.ipynb             # 用於時間序列分群與富集分析的獨立 Notebook
├── analysis_pipeline.py                     # 主要的可執行分析腳本
├── kallisto_to_matrix.py                    # 將 Kallisto 輸出轉換為基因計數的輔助腳本
├── app.py                                   # 互動式 Streamlit 網頁應用程式
├── config.yml                               # 所有使用者可配置的參數
├── requirements.txt                         # 需要的 Python 套件
└── README.md                                # 本文件
```

## 下載專案

您可以在透過以下兩種方式下載此專案：

### 選項 1：Git Clone (推薦)

如果您已安裝 Git，請在終端機中執行以下指令：

```bash
git clone https://github.com/ddmanyes/bulk-rnaseq-pipeline.git
cd bulk-rnaseq-pipeline
```

### 選項 2：下載 ZIP 檔

1. 前往 GitHub 專案頁面：[https://github.com/ddmanyes/bulk-rnaseq-pipeline](https://github.com/ddmanyes/bulk-rnaseq-pipeline)
2. 點擊綠色的 **Code** 按鈕。
3. 選擇 **Download ZIP**。
4. 將下載的檔案解壓縮到您想要的資料夾。

## 🚀 Mac 一鍵安裝 (無需寫程式)

如果您使用的是 macOS，且不想手動處理 Python 環境，請使用內附的 **一鍵安裝檔**：

1. 下載此專案 (ZIP 或 Git)。
2. 雙擊 **`install_and_run_mac.command`**。
3. 等待自動設定完成 (僅第一次需要)。它會：
    - 下載輕量級的環境管理器 (Micromamba)。
    - 在專案資料夾內安裝所有必要的套件。
    - 自動啟動應用程式。

## 手動設定與安裝

如果您偏好手動設定環境：

1. **複製 (Clone) 此儲存庫或下載專案檔案。**

2. **確保您已安裝 Python 3.9+。**

3. **建立並啟用虛擬環境 (建議)：**

    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    # 在 Windows 上使用: .venv\Scripts\activate
    ```

4. **安裝所有需要的套件：**

    ```bash
    pip install torch
    pip install -r requirements.txt
    ```

### Windows 相容性說明

大多數套件在 Windows 上都能順利運作。然而，`omicverse` 和 `torch` (PyTorch) 有時會因為系統依賴性問題而較難安裝。

- **如果自動安裝失敗**：請嘗試先手動安裝 PyTorch（請參考 [pytorch.org](https://pytorch.org/get-started/locally/) 的說明），然後再執行 `pip install -r requirements.txt`。
- **CPU 與 GPU**：預設安裝通常包含 PyTorch 的 CPU 版本，這對於此分析流程已經足夠。

## 輸入數據來源說明

本流程專注於 **下游分析 (Downstream Analysis)**，意即從序列比對與定量完成後的數據開始。

預期的輸入格式為 **基因計數矩陣 (Gene Count Matrix)**（例如來自 **featureCounts** 的輸出），這通常是由以下標準上游流程產生的：

1. **品質控制 (QC)**：使用 FastQC, Trimmomatic 清理原始序列。
2. **序列比對 (Alignment)**：使用 STAR 或 HISAT2 將序列比對回參考基因組 (產生 `.bam` 檔)。
3. **定量 (Quantification)**：使用 **featureCounts** 或 htseq-count 計算每個基因的 Read 數 (產生 `.txt` 檔)。

您可以直接將這些 `featureCounts.txt` 檔案上傳到 App 中進行分析！

## ☁️ Hugging Face Spaces 使用說明

**🔗 App 連結:** [https://huggingface.co/spaces/chanchiru/bulk-rnaseq-pipeline](https://huggingface.co/spaces/chanchiru/bulk-rnaseq-pipeline)

本工具已部署於 Hugging Face Spaces 供您直接使用。但請注意以下限制與使用規範：

### 1. 硬體資源與效能限制

Hugging Face 免費版提供 **2 vCPU 和 16GB RAM**。

- **適用範圍**：Count Matrix 下游分析 (DESeq2, Enrichment)，建議樣本數 **< 100**。
- **不適用**：序列比對 (Alignment, 如 STAR/HISAT2) 或超大數據集 (數千個樣本)，這會導致記憶體不足 (OOM) 而崩潰。
- **建議**：請僅上傳 **Gene Count Matrix**。請勿嘗試在此進行 Raw FASTQ 的比對工作。

### 2. 冷啟動 (Cold Start) 與休眠機制

免費版 Spaces 在閒置一段時間後會自動 **休眠 (Sleep)**。

- **冷啟動**：如果您是近期第一位使用者，可能會看到 "Building" 或 "Starting" 畫面。請耐心 **等待 1-3 分鐘** 讓伺服器喚醒。

### 3. 輸入格式與範例數據

為避免錯誤，請確保您的輸入檔案符合格式。

- **範例數據**：點擊側邊欄的 **"📥 Load Demo Data"** 按鈕，可以載入一組範例數據。**這將會自動載入預先計算好的分析結果，讓您可以直接體驗互動式圖表功能。** 您也可以下載生成的 `merged_counts.csv` 作為模板參考。
- **格式**：程式支援 **CSV** 或 **Tab 分隔** 的文字檔。

### 4. 數據隱私聲明
>
> ⚠️ **免責聲明**：本工具運行於 **公開雲端環境**。
>
> - **工作階段隔離**：我們已實作工作階段隔離。每位使用者都會獲得一個獨立的暫存目錄，確保您的數據 **不會** 被其他同時在線的使用者看到。
> - **暫存性**：檔案僅暫存於記憶體中，工作階段結束後即刪除。
> - **請勿上傳**：病患個資 (PII) 或極度機密的未發表數據。
> - **用途**：僅供學術交流、測試與教學使用。

## 如何執行流程

您可以使用 **互動式 GUI** (推薦) 或 **命令列** 來執行流程。

### 選項 1：使用互動式 GUI (推薦)

Streamlit 應用程式為整個工作流程提供了友善的使用者介面。

1. **啟動應用程式**：

    ```bash
    streamlit run app.py
    ```

2. **瀏覽分頁**：
    - **📂 Data Preprocessing (資料前處理)**：選擇包含多個計數檔案 (CSV, Excel, featureCounts) 的資料夾，並將它們合併為單一矩陣。
    - **📝 Metadata Creation (建立元數據)**：互動式地建立或編輯您的 `metadata.csv` 檔案。
        - **上傳元數據**：**[新增]** 直接上傳您現有的 `metadata.csv` 或 `.xlsx` 檔案來填入表格。
        - **同步 Sample IDs**：自動從您的合併計數矩陣中抓取 Sample IDs。
        - **選擇樣本**：使用 **「Include」** 勾選框來選擇要儲存的樣本。未勾選的樣本將被排除。
        - **可選的時間欄位**：切換 **「Include Time Column」** 來啟用/停用時間序列功能。如果停用，將跳過時間序列分析。
        - **自訂時間欄位**：您可以重新命名「時間」欄位（例如改為 "Stage"）並使用文字內容。
    - **⚙️ Run Analysis (執行分析)**：設定所有分析參數。
        - **基因註釋**：從下拉選單中選擇適合您物種的註釋檔案。
        - 點擊 **"Run Analysis Pipeline"** 開始分析。您可以即時查看執行日誌。
    - **📊 Static Plots (靜態圖表)**：查看分析產生的所有預生成圖表 (PCA, 火山圖, 熱圖)。
    - **🌋 Volcano Plot (火山圖)**：查看靜態火山圖，並透過互動式表格探索底層數據（依顯著性篩選、下載結果）。
    - **Interactive Visualization (互動式視覺化)**：
        - **Boxplot (箱形圖)**：視覺化不同組別間的基因表現分佈。**自動 Log 轉換**：如果檢測到原始計數 (數值 > 100)，數據將自動進行 log 轉換 (`log2(x+1)`) 以獲得更好的視覺效果。
    - **🧬 Pathway Analysis (通路分析)**：探索富集分析結果。選擇特定的比較和通路圖表（例如 GO Biological Process, KEGG）並下載結果表格。

### 選項 2：命令列工作流程

1. **準備輸入**：將 `featureCounts` 檔案放入 `input/` 並確保 `metadata.csv` 已準備好。
2. **設定**：手動編輯 `config.yml`。
3. **執行**：

    ```bash
    python analysis_pipeline.py
    ```

### 第一步：資料前處理 (選用)

如果您有多個計數檔案（例如 featureCounts 的輸出），請前往 **「📂 Data Preprocessing」** 分頁：

1. **上傳檔案**：將您的計數檔案 (`.txt`, `.csv`, `.xlsx`) 拖曳到上傳區塊中。
2. **執行合併**：點擊 "Merge Files"。程式會將它們合併成一個 `merged_counts.csv` 並存放在 `input` 資料夾中。

## 批次效應校正

1. **元數據 (Metadata)**：在您的元數據中新增 "Batch" 欄位 (透過 GUI 或 CSV)。
2. **設定 (Config)**：在 **"Run Analysis"** 分頁或 `config.yml` 中啟用批次校正。
3. **結果**：檢查 `output/PCA_Batch_Correction_Comparison.png`。

## 輸出檔案

流程將在 `output/` 目錄中生成以下檔案：

- **圖表**：`PCA_plot.png`, `Volcano_*.png`, `Heatmap_*.png`, `TimeSeries_*.png`。
- **數據**：`counts.csv`, `counts_batch_corrected.csv`, `DEG_*.csv`, `Master_Results_Table.csv`。
- **富集分析**：`enrichment/` 資料夾，包含 GO/KEGG 分析的 dotplots 和 CSV。

## 時間序列富集分析 Notebook

對於想要進行**深入時間序列分析**（分群 + 功能富集）的使用者，我們提供了一個獨立的 Jupyter Notebook：`time_series_enrichment.ipynb`。

### 功能特色

- **自訂分群**：可自由調整 K 值與視覺化參數。
- **數據處理**：自動處理重複樣本 (replicates) 與 0h 對照組時間點。
- **富集分析**：針對**每個識別出的 Cluster** 執行 GO Biological Process 和 KEGG pathway 富集分析。
- **互動性**：直接在 Notebook 中修改基因集或分群邏輯。

### 如何使用

1. 至少執行一次主流程以生成 `output/normalized_counts_for_plotting.csv`。
2. 在 Jupyter Lab 或 Notebook 中開啟 `time_series_enrichment.ipynb`。
3. 執行所有單元格 (Cells) 以生成分群熱圖和富集分析點圖。

## 參考文獻與引用 (References)

如果您在研究中使用了此流程，請考慮引用以下核心工具與套件：

- **OmicVerse**: [https://github.com/StarlitOnes/OmicVerse](https://github.com/StarlitOnes/OmicVerse)
- **GSEApy (Enrichr)**: Fang, Z., Liu, X., & Peltz, G. (2023). GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics*.
- **ComBat (Batch Correction)**: Johnson, W. E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*.
- **Streamlit**: [https://streamlit.io/](https://streamlit.io/)

## 授權 (License)

本專案採用 MIT License 授權 - 詳情請見 [LICENSE](LICENSE) 檔案。
