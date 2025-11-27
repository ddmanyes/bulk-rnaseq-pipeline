# Bulk RNA-seq 分析流程與互動式瀏覽器

本專案提供了一個全面且自動化的 Bulk RNA-seq 分析流程，並搭配一個互動式網頁應用程式來視覺化分析結果。

此流程從原始的 `featureCounts` 矩陣開始，執行差異基因表現分析 (DGE)、富集分析 (Enrichment Analysis) 和時間序列分群 (Time-Series Clustering)，並生成一套完整的發表級圖表和表格。整個工作流程都可以透過單一個 `config.yml` 檔案輕鬆控制。

## 主要功能

- **端到端分析**：只需一個指令，即可從原始計數數據產生最終圖表。
- **主成分分析 (PCA)**：生成 PCA 圖，根據整體基因表現視覺化樣本間的關係和批次效應。
- **批次效應校正 (Batch Correction)**：整合 ComBat 批次校正功能，可移除不必要的技術變異，並自動生成「校正前後」的比較圖。
- **豐富的視覺化**：自動生成火山圖 (Volcano Plots)、聚類熱圖 (Clustered Heatmaps，針對顯著基因和高變異基因) 以及時間序列趨勢圖。
- **富集分析**：針對多個資料庫 (如 GO, KEGG, Reactome) 執行過度表現分析 (ORA)，並生成結果的點圖 (Dot Plots)。
- **時間序列分群**：使用 K-means 分群法識別時間序列數據中共同表現的基因模組，並視覺化其趨勢。
- **完全可配置**：在一個中央 `config.yml` 檔案中輕鬆管理所有檔案路徑、分析參數和閾值。
- **互動式瀏覽器**：包含一個 Streamlit 網頁應用程式，可透過箱形圖 (Boxplots) 和互動式熱圖 (Interactive Heatmaps) 動態探索基因表現。
- **數據導出**：可直接從應用程式下載處理後的數據、差異表現結果和通路分析表格。
- **可重現性**：附帶 `requirements.txt` 檔案，確保環境的一致性。
- **Kallisto 整合**：包含一個輔助腳本，可將 Kallisto 的轉錄本豐度 (transcript abundance) 聚合為與此流程相容的基因級計數矩陣 (gene-level count matrix)。

## 專案結構

```
├── input/
│   ├── KwtfWz_TS241219002_featureCount.txt  # 您的原始 featureCounts 數據
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
├── analysis_pipeline.py                     # 主要的可執行分析腳本
├── kallisto_to_matrix.py                    # 將 Kallisto 輸出轉換為基因計數的輔助腳本
├── app.py                                   # 互動式 Streamlit 網頁應用程式
├── config.yml                               # 所有使用者可配置的參數
├── requirements.txt                         # 需要的 Python 套件
└── README.md                                # 本文件
```

## 下載專案

您可以透過以下兩種方式下載此專案：

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

## 如何執行流程

您可以使用 **互動式 GUI** (推薦) 或 **命令列** 來執行流程。

### 選項 1：使用互動式 GUI (推薦)

Streamlit 應用程式為整個工作流程提供了友善的使用者介面。

1. **啟動應用程式**：

    ```bash
    streamlit run app.py
    ```

2. **瀏覽分頁**：
    - **📂 Merge Counts (合併計數)**：選擇包含多個計數檔案 (CSV, Excel, featureCounts) 的資料夾，並將它們合併為單一矩陣。
    - **📝 Metadata Creation (建立元數據)**：互動式地建立或編輯您的 `metadata.csv` 檔案。確保包含 `Condition`, `Time`, 和 `Type` 欄位。
    - **⚙️ Run Analysis (執行分析)**：設定所有分析參數 (檔案、閾值、批次校正等) 並點擊 **"Run Analysis Pipeline"**。您可以即時查看執行日誌。
    - **📊 Interactive Visualization (互動式視覺化)**：
        - **Boxplot (箱形圖)**：視覺化不同組別間的基因表現分佈。
        - **Heatmap (熱圖)**：透過 Z-score 標準化和可自訂的色階探索表現模式。
    - **🌋 Volcano Plot (火山圖)**：查看靜態火山圖，並透過互動式表格探索底層數據（依顯著性篩選、下載結果）。
    - **📈 Static Plots (靜態圖表)**：查看分析產生的所有預生成圖表 (PCA, 火山圖, 熱圖)。
    - **🧬 Pathway Analysis (通路分析)**：探索富集分析結果。選擇特定的比較和通路圖表（例如 GO Biological Process, KEGG）並下載結果表格。

### 選項 2：命令列工作流程

1. **準備輸入**：將 `featureCounts` 檔案放入 `input/` 並確保 `metadata.csv` 已準備好。
2. **設定**：手動編輯 `config.yml`。
3. **執行**：

    ```bash
    python analysis_pipeline.py
    ```

## 合併多個數據集

如果您有多個計數檔案，請使用 GUI 中的 **"Merge Counts"** 分頁或執行輔助腳本：

```bash
python merge_counts.py input/batch1.csv input/batch2.xlsx -o input/merged_counts.csv
```

## 批次效應校正

1. **元數據 (Metadata)**：在您的元數據中新增 "Batch" 欄位 (透過 GUI 或 CSV)。
2. **設定 (Config)**：在 **"Run Analysis"** 分頁或 `config.yml` 中啟用批次校正。
3. **結果**：檢查 `output/PCA_Batch_Correction_Comparison.png`。

## 輸出檔案

流程將在 `output/` 目錄中生成以下檔案：

- **圖表**：`PCA_plot.png`, `Volcano_*.png`, `Heatmap_*.png`, `TimeSeries_*.png`。
- **數據**：`counts.csv`, `counts_batch_corrected.csv`, `DEG_*.csv`, `Master_Results_Table.csv`。
- **富集分析**：`enrichment/` 資料夾，包含 GO/KEGG 分析的 dotplots 和 CSV。

## 參考文獻與引用 (References)

如果您在研究中使用了此流程，請考慮引用以下核心工具與套件：

- **OmicVerse**: [https://github.com/StarlitOnes/OmicVerse](https://github.com/StarlitOnes/OmicVerse)
- **GSEApy (Enrichr)**: Fang, Z., Liu, X., & Peltz, G. (2023). GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics*.
- **ComBat (Batch Correction)**: Johnson, W. E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*.
- **Streamlit**: [https://streamlit.io/](https://streamlit.io/)

## 授權 (License)

本專案採用 MIT License 授權 - 詳情請見 [LICENSE](LICENSE) 檔案。
