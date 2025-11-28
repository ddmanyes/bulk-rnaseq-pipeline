---
title: Bulk RNAseq Pipeline
emoji: 🧬
colorFrom: blue
colorTo: purple
sdk: streamlit
sdk_version: 1.39.0
app_file: app.py
pinned: false
license: mit
---

# Bulk RNA-seq Analysis Pipeline & Interactive Explorer

This project provides a comprehensive and automated pipeline for bulk RNA-seq analysis, paired with an interactive web application for visualizing the results.

The pipeline starts from a raw `featureCounts` matrix, performs differential gene expression, enrichment analysis, and time-series clustering, and generates a full suite of publication-ready plots and tables. The entire workflow is easily controlled through a single `config.yml` file.

## Key Features

- **End-to-End Analysis**: From raw counts to final plots with one command.
- **Principal Component Analysis (PCA)**: Generates a PCA plot to visualize sample-level relationships and batch effects based on overall gene expression.
- **Batch Correction**: Integrated ComBat batch correction to remove unwanted technical variation, with automatic "Before vs After" visualization.
- **Rich Visualizations**: Automatically generates Volcano Plots, Clustered Heatmaps (for significant and top-variable genes), and Time-Series trend plots.
- **Enrichment Analysis**: Performs Over-Representation Analysis (ORA) against multiple databases (e.g., GO, KEGG, Reactome) and generates dot plots for results.
- **Time-Series Clustering**: Identifies co-expressed gene modules in time-course data using K-means clustering and visualizes their trends.
- **Fully Configurable**: Easily manage all file paths, analysis parameters, and thresholds in one central `config.yml` file.
- **Interactive Explorer**: Includes a Streamlit web app to dynamically explore gene expression via Boxplots and Interactive Heatmaps.
- **Data Export**: Easily download processed data, differential expression results, and pathway analysis tables directly from the app.
- **Reproducible**: Bundled with a `requirements.txt` file to ensure a consistent environment.
- **Kallisto Integration**: Includes a helper script to aggregate Kallisto transcript abundance into a gene-level count matrix compatible with this pipeline.

## Project Structure

```
├── input/
│   ├── featureCount files                  # Your raw featureCounts data
│   └── metadata.csv                        # Sample metadata (Condition, Time, Type)
├── genesets/
│   └── pair_GRCm39.tsv                     # Gene annotation files
├── output/                                 # All generated results are saved here
│   ├── PCA_plot.png                        # Principal Component Analysis plot
│   ├── PCA_Batch_Correction_Comparison.png  # PCA showing before vs after batch correction
│   ├── counts.csv                           # Preprocessed and gene-ID-mapped counts
│   ├── counts_batch_corrected.csv           # Batch-corrected count matrix
│   ├── log2cpm_normalized_counts.csv        # Log2CPM normalized counts
│   ├── DEG_*.csv                            # DGE results for each comparison
│   ├── Volcano_*.png                        # Volcano plots for each comparison
│   ├── Heatmap_Significant_Genes.png        # Heatmap of all significant genes
│   ├── Heatmap_Top50_Variable_Genes.png     # Heatmap of top N variable genes
│   ├── Master_Results_Table.csv             # Combined DGE, cluster, and expression data
│   ├── Master_Results_Table.xlsx            # Combined DGE, cluster, and expression data
│   ├── TimeSeries_Clustering_k4.png         # Time-series cluster plot
│   ├── k_evaluation_plot.png                # Optional plot to help choose k
│   └── enrichment/                          # GO/KEGG/Reactome results (CSVs & plots)
├── analysis_pipeline.py                     # The main executable analysis script
├── kallisto_to_matrix.py                    # Helper script to convert Kallisto output to gene counts
├── app.py                                   # The interactive Streamlit web application
├── config.yml                               # All user-configurable parameters
├── requirements.txt                         # Required Python libraries
└── README.md                                # This documentation file
```

## Download the Project

You can download this project in two ways:

### Option 1: Git Clone (Recommended)

If you have Git installed, run the following command in your terminal:

```bash
git clone https://github.com/ddmanyes/bulk-rnaseq-pipeline.git
cd bulk-rnaseq-pipeline
```

### Option 2: Download ZIP

1. Go to the GitHub repository page: [https://github.com/ddmanyes/bulk-rnaseq-pipeline](https://github.com/ddmanyes/bulk-rnaseq-pipeline)
2. Click the green **Code** button.
3. Select **Download ZIP**.
4. Unzip the downloaded file to your desired folder.

## Manual Setup and Installation

If you prefer to set up the environment manually:

1. **Clone the repository or download the project files.**

2. **Ensure you have Python 3.9+ installed.**

3. **Create and activate a virtual environment (recommended):**

    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    # On Windows, use: .venv\Scripts\activate
    ```

4. **Install all required libraries:**

    ```bash
    pip install torch
    pip install -r requirements.txt
    ```

### Windows Compatibility Note

Most packages work seamlessly on Windows. However, `omicverse` and `torch` (PyTorch) can sometimes be tricky to install due to system-specific dependencies.

- **If the automatic installer fails**: Try installing PyTorch manually first, following the instructions at [pytorch.org](https://pytorch.org/get-started/locally/), and then run `pip install -r requirements.txt`.
- **CPU vs GPU**: The default installation typically includes the CPU version of PyTorch, which is sufficient for this pipeline.

## How to Run the Pipeline

You can run the pipeline using the **Interactive GUI** (Recommended) or the **Command Line**.

### Option 1: Using the Interactive GUI (Recommended)

The Streamlit app provides a user-friendly interface for the entire workflow.

1. **Launch the App**:

    ```bash
    streamlit run app.py
    ```

2. **Navigate the Tabs**:
    - **📂 Merge Counts**: Select a folder containing multiple count files (CSV, Excel, featureCounts) and merge them into a single matrix.
    - **📝 Metadata Creation**: Create or edit your `metadata.csv` file interactively. Ensure you have `Condition`, `Time`, and `Type` columns.
    - **⚙️ Run Analysis**: Configure all analysis parameters (files, thresholds, batch correction, etc.) and click **"Run Analysis Pipeline"**. You can see real-time logs of the execution.
    - **📊 Interactive Visualization**:
        - **Boxplot**: Visualize gene expression distribution across groups.
        - **Heatmap**: Explore expression patterns with Z-score normalization and customizable color scales.
    - **🌋 Volcano Plot**: View static Volcano Plots and explore the underlying data with an interactive table (filter by significance, download results).
    - **📈 Static Plots**: View all the pre-generated plots (PCA, Volcano, Heatmaps) from the analysis.
    - **🧬 Pathway Analysis**: Explore enrichment analysis results. Select specific comparisons and pathway plots (e.g., GO Biological Process, KEGG) and download result tables.

### Option 2: Command Line Workflow

1. **Prepare Input**: Place `featureCounts` file in `input/` and ensure `metadata.csv` is ready.
2. **Configure**: Edit `config.yml` manually.
3. **Run**:

    ```bash
    python analysis_pipeline.py
    ```

### Step 1: Data Preprocessing (Optional)

If you have multiple count files (e.g., from featureCounts), go to the **"📂 Data Preprocessing"** tab:

1. **Upload Files**: Drag and drop your count files (`.txt`, `.csv`, `.xlsx`) into the file uploader.
2. **Merge**: Click "Merge Files". The app will combine them into a single `merged_counts.csv` in the `input` folder.

## Batch Correction

1. **Metadata**: Add a "Batch" column to your metadata (via GUI or CSV).
2. **Config**: Enable batch correction in the **"Run Analysis"** tab or `config.yml`.
3. **Result**: Check `output/PCA_Batch_Correction_Comparison.png`.

## Output Files

The pipeline generates the following in the `output/` directory:

- **Plots**: `PCA_plot.png`, `Volcano_*.png`, `Heatmap_*.png`, `TimeSeries_*.png`.
- **Data**: `counts.csv`, `counts_batch_corrected.csv`, `DEG_*.csv`, `Master_Results_Table.csv`.
- **Enrichment**: `enrichment/` folder containing dotplots and CSVs for GO/KEGG analysis.

## References & Citations

If you use this pipeline in your research, please consider citing the following tools and packages that made this work possible:

- **OmicVerse**: [https://github.com/StarlitOnes/OmicVerse](https://github.com/StarlitOnes/OmicVerse)
- **GSEApy (Enrichr)**: Fang, Z., Liu, X., & Peltz, G. (2023). GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics*.
- **ComBat (Batch Correction)**: Johnson, W. E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*.
- **Streamlit**: [https://streamlit.io/](https://streamlit.io/)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
