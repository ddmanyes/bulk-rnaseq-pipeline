import os
import pandas as pd
import omicverse as ov
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import warnings
import io
import gseapy as gp
from contextlib import redirect_stdout
import yaml
from adjustText import adjust_text
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
from scipy import stats
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import math
import anndata

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# --- Helper Functions ---

def load_config(path='config.yml'):
    """Loads the YAML config file."""
    with open(path, 'r') as f:
        return yaml.safe_load(f)

# --- Step 1: Data Preprocessing ---
def preprocess_data(feature_counts_file, gene_annotation_file, output_dir):
    """
    Reads raw featureCounts, cleans it, converts gene IDs, and saves the count matrix.
    """
    print("--- Step 1: Data Preprocessing ---")
    
    # 1.1 Read file
    print(f"Reading counts from '{feature_counts_file}'...")
    
    # Determine file type based on extension
    if feature_counts_file.endswith('.csv'):
        df = pd.read_csv(feature_counts_file, index_col=0)
        sep = ','
    elif feature_counts_file.endswith('.xlsx'):
        df = pd.read_excel(feature_counts_file, index_col=0)
        sep = 'excel'
    else:
        df = pd.read_csv(feature_counts_file, sep='\t', comment='#', index_col=0)
        sep = '\t'

    # Check if it's a featureCounts file (has 'Chr', 'Start', etc.)
    if 'Chr' in df.columns and 'Start' in df.columns:
        print("  - Detected featureCounts format. Cleaning columns...")
        counts_matrix = df.iloc[:, 5:]
        counts_matrix.columns = [i.split('/')[-1].replace('.sorted.bam', '') for i in counts_matrix.columns]
    else:
        print("  - Detected simple count matrix format. Using as is.")
        counts_matrix = df

    # 1.2 Gene ID Conversion
    print("Downloading gene ID annotation table (if not present)...")
    ov.utils.download_geneid_annotation_pair()
    
    print("Converting gene IDs to gene symbols...")
    
    # Resolve annotation file path
    if not gene_annotation_file:
        gene_annotation_file = "genesets/pair_GRCm39.tsv" # Default fallback
        
    if not os.path.isabs(gene_annotation_file):
        # Allow resolving relative to the script directory if not found in CWD
        candidate_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), gene_annotation_file)
        if os.path.exists(candidate_path):
            gene_annotation_file = candidate_path
    
    if os.path.exists(gene_annotation_file):
        try:
            original_counts = counts_matrix.copy()
            counts_matrix = ov.bulk.Matrix_ID_mapping(counts_matrix, gene_annotation_file)
            
            if counts_matrix.empty:
                print("Warning: Gene ID mapping resulted in 0 genes (likely wrong species or IDs already symbols). Reverting to original matrix.")
                counts_matrix = original_counts
        except Exception as e:
            print(f"Warning: Gene ID mapping failed: {e}. Proceeding with original IDs.")
            counts_matrix = original_counts
    else:
        print(f"Warning: Annotation file '{gene_annotation_file}' not found. Skipping ID mapping.")

    
    # 1.3 Save processed counts
    processed_counts_path_csv = os.path.join(output_dir, "counts.csv")
    processed_counts_path_xlsx = os.path.join(output_dir, "counts.xlsx")
    counts_matrix.to_csv(processed_counts_path_csv)
    counts_matrix.to_excel(processed_counts_path_xlsx)
    
    print(f"Processed count matrix saved to '{processed_counts_path_csv}' and '{processed_counts_path_xlsx}'")
    print(f"Final data dimension: {counts_matrix.shape} (Genes x Samples)")
    print("-" * 40)
    return counts_matrix

# --- Step 2: Metadata Creation/Loading ---
def load_metadata(metadata_file, output_dir, cfg):
    """
    Loads metadata from a file or creates a default one if missing.
    """
    print("--- Step 2: Loading Metadata ---")
    
    if os.path.exists(metadata_file):
        print(f"Loading metadata from '{metadata_file}'...")
        if metadata_file.endswith('.xlsx'):
            metadata = pd.read_excel(metadata_file, index_col=0)
        else:
            metadata = pd.read_csv(metadata_file, index_col=0)
    else:
        print(f"⚠️ Metadata file '{metadata_file}' not found. Creating default metadata...")
        sample_ids = [
            'C_1', 'C_2', 'T_6_1', 'T_6_2',
            'T_24_1', 'T_24_2', 'T_72_1', 'T_72_2'
        ]
        conditions = [
            'Control', 'Control', 'Treat_6h', 'Treat_6h',
            'Treat_24h', 'Treat_24h', 'Treat_72h', 'Treat_72h'
        ]
        times = [0, 0, 6, 6, 24, 24, 72, 72]
        types = ['Control', 'Control', 'Treat', 'Treat', 'Treat', 'Treat', 'Treat', 'Treat']
        
        metadata = pd.DataFrame({
            'Condition': conditions,
            'Time': times,
            'Type': types
        }, index=sample_ids)
        metadata.index.name = 'SampleID'
        metadata.to_csv(metadata_file)
        print(f"  - Default metadata saved to '{metadata_file}'")

    # Ensure required columns exist
    required_cols = ['Condition', 'Type']
    
    # Check if Time Series is enabled and what the column name is
    ts_cfg = cfg.get('analysis', {}).get('time_series', {})
    if ts_cfg.get('enabled', True):
        time_col = ts_cfg.get('time_column', 'Time')
        # Only require it if it's supposed to be present
        if ts_cfg.get('time_column_present', True):
            required_cols.append(time_col)

    missing = [c for c in required_cols if c not in metadata.columns]
    if missing:
        # If Time is missing but we expected it, try to be helpful
        if 'Time' in missing and 'Time' not in metadata.columns:
             # Check if maybe the user renamed it in the file but didn't update config?
             # For now, just raise error as before but with better logic
             pass
        raise ValueError(f"Metadata file is missing required columns: {missing}")

    # Save a copy to output for reference
    metadata_path_xlsx = os.path.join(output_dir, "metadata.xlsx")
    metadata.to_excel(metadata_path_xlsx)

    print(f"Metadata loaded with {len(metadata)} samples.")
    print("-" * 40)
    return metadata

# --- Step 2.5: Batch Correction (Optional) ---
def perform_batch_correction(counts, metadata, batch_key, output_dir):
    """
    Performs ComBat batch correction using omicverse.
    Returns the corrected count matrix (DataFrame).
    """
    print("\n--- Step 2.5: Batch Correction ---")
    
    # 0. Check and Log-Transform Data if needed
    # ComBat expects data to be roughly normal (log-transformed counts), not raw integer counts.
    # If max value is large (> 100), assume it's raw counts and log-transform.
    is_log_transformed = False
    if counts.max().max() > 100:
        print("  - Data appears to be raw counts (Max > 100). applying log2(x+1) before batch correction...")
        # Add 1 pseudo-count to handle zeros
        counts_for_correction = np.log2(counts + 1)
        is_log_transformed = True
    else:
        print("  - Data appears to be already log-transformed (Max <= 100). Using as is.")
        counts_for_correction = counts

    # 1. Prepare AnnData
    # AnnData expects X to be (n_obs, n_vars), so we transpose counts (n_vars, n_obs)
    adata = anndata.AnnData(counts_for_correction.T)
    
    # Add metadata
    # Ensure metadata index matches adata.obs_names
    adata.obs = metadata.loc[adata.obs_names].copy()
    
    # Check if batch key exists
    if batch_key not in adata.obs.columns:
        raise ValueError(f"Batch key '{batch_key}' not found in metadata columns: {adata.obs.columns.tolist()}")
    
    print(f"  - Found batch key: '{batch_key}'")
    print(f"  - Batches detected: {adata.obs[batch_key].unique().tolist()}")
    
    # 2. Run Batch Correction
    # Note: omicverse's batch_correction modifies adata in-place or adds a layer
    print("  - Running ComBat batch correction...")
    try:
        ov.bulk.batch_correction(adata, batch_key=batch_key)
    except Exception as e:
        print(f"  ❌ Batch correction failed: {e}")
        raise e

    # 3. Extract Corrected Data
    # The result is stored in adata.layers['batch_correction']
    if 'batch_correction' in adata.layers:
        corrected_data = adata.to_df(layer='batch_correction').T # Transpose back to (Genes x Samples)
    else:
        # Fallback if specific layer not found (though ov usually adds it)
        print("  ⚠️ Warning: 'batch_correction' layer not found, checking X...")
        corrected_data = adata.to_df().T

    # 4. Save Results
    corrected_path = os.path.join(output_dir, "counts_batch_corrected.csv")
    corrected_data.to_csv(corrected_path)
    print(f"  - Corrected counts saved to '{corrected_path}'")

    # 5. Visualization: PCA Comparison (Before vs After)
    print("  - Generating PCA comparison plot...")
    
    # Helper to calc PCA
    def get_pca_df(data_df, meta_df, label):
        # Log-transform for PCA (pseudo-count + 1) -> standard practice for visualization
        # Note: ComBat output might already be normalized/transformed depending on usage, 
        # but usually it returns corrected expression values. 
        # If values are large (counts), we log. If they are small (log-like), we don't.
        # Simple check: max value
        if data_df.max().max() > 100:
            mat = np.log2(data_df + 1)
        else:
            mat = data_df
            
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(mat.T)
        df = pd.DataFrame(pcs, columns=['PC1', 'PC2'], index=mat.columns)
        df = df.merge(meta_df, left_index=True, right_index=True)
        df['Type'] = label
        return df

    try:
        pca_before = get_pca_df(counts, metadata, 'Before Correction')
        pca_after = get_pca_df(corrected_data, metadata, 'After Correction')
        
        combined_pca = pd.concat([pca_before, pca_after])
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Plot Before
        sns.scatterplot(data=pca_before, x='PC1', y='PC2', hue=batch_key, style='Condition', s=100, ax=axes[0])
        axes[0].set_title(f"Before Correction (Color={batch_key})")
        
        # Plot After
        sns.scatterplot(data=pca_after, x='PC1', y='PC2', hue=batch_key, style='Condition', s=100, ax=axes[1])
        axes[1].set_title(f"After Correction (Color={batch_key})")
        
        plt.tight_layout()
        comp_plot_path = os.path.join(output_dir, "PCA_Batch_Correction_Comparison.png")
        plt.savefig(comp_plot_path, dpi=300)
        plt.close()
        print(f"  - PCA comparison saved to '{comp_plot_path}'")
        
    except Exception as plot_e:
        print(f"  ⚠️ PCA plotting failed: {plot_e}")

    print("-" * 40)
    return corrected_data

# --- Step 3: Batch Differential Expression Analysis ---
def run_batch_deg_analysis(dds, metadata, cfg, gene_annotation_file):
    """
    Automatically runs DESeq2 and enrichment analysis for all treatment groups against the control.
    """
    print("--- Step 3: Batch Differential Expression and Enrichment ---")
    
    # Extract config parameters
    analysis_cfg = cfg['analysis']
    control_group = analysis_cfg['diff_expression']['control_group']
    filter_low = analysis_cfg['diff_expression']['filter_low_expression']
    enrichment_cfg = analysis_cfg['enrichment']
    sig_cfg = analysis_cfg['significance']
    output_dir = cfg['files']['output_dir']

    ctrl_samples = metadata[metadata['Condition'] == control_group].index.tolist()
    all_groups = metadata['Condition'].unique()
    treat_groups = [g for g in all_groups if g != control_group]
    
    print(f"Control group: {control_group} (n={len(ctrl_samples)})")
    print(f"Detected {len(treat_groups)} treatment groups: {treat_groups}")
    
    results_dict = {}
    
    for treat in treat_groups:
        print(f"\nAnalyzing: {treat} vs {control_group}...")
        treat_samples = metadata[metadata['Condition'] == treat].index.tolist()
        
        try:
            res = dds.deg_analysis(treat_samples, ctrl_samples, method='DEseq2')

            if filter_low:
                print(f"  - Filtering low expression genes (log2(BaseMean) > 1)...")
                res = res.loc[res['log2(BaseMean)'] > 1]

            if 'log2FoldChange' in res.columns: res['log2FC'] = res['log2FoldChange']
            if 'padj' in res.columns: res['qvalue'] = res['padj']
            
            comparison_name = f"{treat}_vs_{control_group}"
            results_dict[comparison_name] = res
            
            res_path_csv = os.path.join(output_dir, f"DEG_{comparison_name}.csv")
            res.to_csv(res_path_csv)
            print(f"  - Result saved to '{res_path_csv}'")
            
            if not res.empty:
                # Volcano plot
                try:
                    plot_title = f'{treat} vs {control_group}'
                    
                    # --- 0. 欄位名稱防呆 (新增) ---
                    fc_col = 'log2FC' if 'log2FC' in res.columns else 'log2FoldChange'
                    p_col = 'qvalue' if 'qvalue' in res.columns else 'padj'

                    # Get thresholds from config
                    lfc_thresh = sig_cfg['log2fc_threshold']
                    padj_thresh = sig_cfg['padj_threshold']

                    # --- 1. 數據預處理 (處理 P=0 問題) ---
                    # 避免 log10(0) = inf
                    if (res[p_col] == 0).any():
                        min_nonzero = res.loc[res[p_col] > 0, p_col].min()
                        res[p_col] = res[p_col].replace(0, min_nonzero / 10)
                        print("  ⚠️ Warning: Some p-values were 0, replaced with min_nonzero/10 for plotting.")

                    # --- 2. Auto-calculate symmetric X-axis range ---
                    max_fc = res[fc_col].abs().max() * 1.1
                    limit = max(max_fc, lfc_thresh * 1.2) 

                    # --- 3. Create canvas and assign colors ---
                    fig, ax = plt.subplots(figsize=(7, 7))
                    
                    # 設定顏色
                    res['color'] = 'grey'
                    res.loc[(res[fc_col] > lfc_thresh) & (res[p_col] < padj_thresh), 'color'] = '#E64B35' # Red (Up)
                    res.loc[(res[fc_col] < -lfc_thresh) & (res[p_col] < padj_thresh), 'color'] = '#3182BD' # Blue (Down)
                    
                    # 排序：把顯著的點放到 DataFrame 後面，這樣畫圖時會蓋在灰色點上面 (Z-order trick)
                    res = res.sort_values(by='color', key=lambda x: x != 'grey')

                    # --- 4. Draw points and lines ---
                    ax.scatter(res[fc_col], -np.log10(res[p_col]), 
                               c=res['color'], s=15, alpha=0.6, edgecolors='none')
                    
                    # 輔助線
                    ax.axhline(y=-np.log10(padj_thresh), color='black', linestyle='--', lw=1, alpha=0.5)
                    ax.axvline(x=lfc_thresh, color='black', linestyle='--', lw=1, alpha=0.5)
                    ax.axvline(x=-lfc_thresh, color='black', linestyle='--', lw=1, alpha=0.5)

                    # --- 5. Set labels and limits ---
                    ax.set_xlim(-limit, limit)
                    ax.set_title(plot_title)
                    ax.set_xlabel(f"{fc_col} (Log2 Fold Change)")
                    ax.set_ylabel(f"-log10({p_col})")
                    
                    # --- 6. Add smart gene labels ---
                    # 只標示顯著的基因
                    sig_genes = res[res['color'] != 'grey']
                    # 挑選 P-value 最小的前 10 個
                    top_genes = sig_genes.sort_values(p_col).head(10)
                    
                    texts = []
                    if not top_genes.empty:
                        for gene_name, row in top_genes.iterrows():
                            # gene_name 是 index，確保這是你要的 label
                            texts.append(ax.text(row[fc_col], -np.log10(row[p_col]), gene_name, fontsize=9))
                        
                        # 自動避讓標籤
                        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))
                    
                    # --- 7. Save file ---
                    volcano_path = os.path.join(output_dir, f"Volcano_{comparison_name}.png")
                    plt.savefig(volcano_path, dpi=300, bbox_inches='tight')
                    plt.close(fig)
                    print(f"  - Successfully saved new volcano plot to '{volcano_path}'")

                except Exception as plot_e:
                    print(f"  - Error creating/saving new volcano plot: {plot_e}")
                    import traceback
                    traceback.print_exc() # 印出詳細錯誤訊息以便除錯

                # Enrichment analysis
                run_enrichment_analysis(
                    deg_results=res,
                    comparison_name=comparison_name,
                    enrichment_dir=os.path.join(output_dir, "enrichment"),
                    cfg=enrichment_cfg,
                    sig_cfg=sig_cfg,
                    gene_annotation_file=gene_annotation_file
                )
        except Exception as e:
            print(f"  - Analysis failed for {treat}: {e}")
            
    print("-" * 40)
    return results_dict

# --- ORA Analysis Function (from user's notebook) ---
def run_ora_analysis(deg_df, prefix, species, library, padj_thresh, log2fc_thresh, top_term, save_dir, figsize=None):
    species = species.capitalize()
    library_name = library
    
    # Use specific library names from config
    print(f"📊 分析設定: {prefix} | {species} | {library_name}")

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    fc_col = 'log2FoldChange' if 'log2FoldChange' in deg_df.columns else 'log2FC'
    p_col = 'padj' if 'padj' in deg_df.columns else 'qvalue'

    up_genes = deg_df[(deg_df[fc_col] > log2fc_thresh) & (deg_df[p_col] < padj_thresh)].index.tolist()
    down_genes = deg_df[(deg_df[fc_col] < -log2fc_thresh) & (deg_df[p_col] < padj_thresh)].index.tolist()

    print(f"   - Up: {len(up_genes)} | Down: {len(down_genes)}")

    def process_group(gene_list, direction):
        if len(gene_list) < 5:
            print(f"   ⚠️ {direction}: Not enough genes to analyze (found {len(gene_list)}), skipping.")
            return

        try:
            print(f"🚀 執行 {direction} 分析...")
            enr = gp.enrichr(gene_list=gene_list, gene_sets=library_name, organism=species, outdir=None)
            sig_res = enr.results[enr.results['Adjusted P-value'] < 0.05].copy()

            if sig_res.empty:
                print(f"   ❌ {direction}: No significant results found.")
                return

            clean_lib = prefix.split('_')[0]
            filename = f"{prefix}_{direction}_{clean_lib}.csv"
            sig_res.to_csv(os.path.join(save_dir, filename))
            print(f"   💾 {direction}: Results saved to {os.path.join(save_dir, filename)}")

            from gseapy.plot import dotplot
            actual_terms = min(top_term, len(sig_res))
            final_figsize = figsize if figsize is not None else (6, max(4, actual_terms * 0.5))
            
            fig, ax = plt.subplots(figsize=final_figsize)
            dotplot(enr.res2d, title=f'{prefix} {direction}\n{library_name}', cmap='viridis_r', size=10, top_term=actual_terms, ax=ax)
            plot_filename = f"Dotplot_{prefix}_{direction}_{clean_lib}.png"
            plt.savefig(os.path.join(save_dir, plot_filename), bbox_inches='tight', dpi=300)
            plt.close(fig)
            print(f"   🖼️ {direction}: Plot saved to {os.path.join(save_dir, plot_filename)}")
        except Exception as e:
            print(f"   ❌ Error during ORA for {direction} genes: {e}")

    print("-" * 40)
    process_group(up_genes, "Up")
    print("-" * 40)
    process_group(down_genes, "Down")

# --- Enrichment Analysis Wrapper ---
def run_enrichment_analysis(deg_results, comparison_name, enrichment_dir, cfg, sig_cfg, gene_annotation_file):
    """
    Performs ORA for a list of libraries defined in the config.
    """
    print(f"--- Running ORA Analysis Wrapper for {comparison_name} ---")
    
    species = cfg['species']
    libraries_to_run = cfg.get('libraries', {}) # Use .get for safety
    
    if not libraries_to_run:
        print("  - No libraries defined in config['analysis']['enrichment']['libraries']. Skipping ORA.")
        return

    print(f"  - Using species: {species}")
    print(f"  - Libraries to analyze: {list(libraries_to_run.keys())}")

    for lib_prefix, lib_name in libraries_to_run.items():
        print(f"\n  - Preparing {lib_prefix} analysis...")
        run_ora_analysis(
            deg_df=deg_results, 
            prefix=f"{lib_prefix}_{comparison_name}", 
            species=species,
            library=lib_name, 
            padj_thresh=sig_cfg['padj_threshold'],
            log2fc_thresh=sig_cfg['log2fc_threshold'], 
            top_term=cfg['top_terms_to_plot'],
            save_dir=os.path.join(enrichment_dir, comparison_name)
        )

    print(f"--- ORA Wrapper for {comparison_name} Complete ---")

# --- Step 4: Create Master Table and Heatmap ---
def get_sig_genes_data(results_dict, norm_counts, padj_thresh, log2fc_thresh):
    print("--- Step 4: Identifying Significant Genes for Heatmap ---")
    all_sig_genes = set()
    for name, res in results_dict.items():
        if 'qvalue' in res.columns and 'log2FC' in res.columns:
            sig_mask = (res['qvalue'] < padj_thresh) & (res['log2FC'].abs() > log2fc_thresh)
            sig_genes = res[sig_mask].index.tolist()
            print(f"  - {name}: Found {len(sig_genes)} significant genes")
            all_sig_genes.update(sig_genes)
    print(f"\nTotal unique significant genes: {len(all_sig_genes)}")
    if not all_sig_genes: return pd.DataFrame()
    
    valid_genes = [g for g in all_sig_genes if g in norm_counts.index]
    plot_data = norm_counts.loc[valid_genes]
    return plot_data.loc[plot_data.std(axis=1) > 0].dropna()

def create_master_and_heatmap(results_dict, norm_counts, metadata, output_dir, padj_thresh, log2fc_thresh, cluster_results=None):
    plot_data = get_sig_genes_data(results_dict, norm_counts, padj_thresh, log2fc_thresh)
    if plot_data.empty: return

    print("\nCreating and saving heatmap...")
    unique_conditions = sorted(metadata['Condition'].unique())
    palette = plt.cm.viridis(np.linspace(0, 1, len(unique_conditions)))
    color_map = {cond: color for cond, color in zip(unique_conditions, palette)}
    col_colors = metadata['Condition'].map(color_map)
    g = sns.clustermap(plot_data, z_score=0, cmap="vlag", center=0, col_cluster=False, row_cluster=True, col_colors=col_colors, yticklabels=False, figsize=(8, 10))
    for condition, color in color_map.items():
        g.ax_heatmap.bar(0, 0, color=color, label=condition, linewidth=0)
    g.ax_heatmap.legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left')
    heatmap_path = os.path.join(output_dir, "Heatmap_Significant_Genes.png")
    plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Heatmap saved to '{heatmap_path}'")
    
    print("\nCreating and saving master results table...")
    master_df = pd.DataFrame(index=norm_counts.index)
    for name, res_df in results_dict.items():
        subset = res_df[['log2FC', 'qvalue']].copy()
        subset.rename(columns={'log2FC': f'log2FC_{name}', 'qvalue': f'padj_{name}'}, inplace=True)
        master_df = master_df.merge(subset, left_index=True, right_index=True, how='left')
    
    # Add cluster information if available
    if cluster_results is not None:
        print("Adding time-series cluster information to master table...")
        # The cluster info is for a subset of genes, so we merge with 'how=left'
        cluster_info = cluster_results[['Cluster']].copy()
        cluster_info['Cluster'] = cluster_info['Cluster'] + 1 # Make clusters 1-based
        master_df = master_df.merge(cluster_info, left_index=True, right_index=True, how='left')

    master_df = pd.concat([master_df, norm_counts], axis=1)
    master_path_csv = os.path.join(output_dir, "Master_Results_Table.csv")
    master_df.to_csv(master_path_csv)
    print(f"Master results table saved to '{master_path_csv}'")

    master_path_xlsx = os.path.join(output_dir, "Master_Results_Table.xlsx")
    master_df.to_excel(master_path_xlsx)
    print(f"Master results table saved to '{master_path_xlsx}'")
    print("-" * 40)

def create_top_variable_heatmap(norm_counts, metadata, output_dir, top_n):
    print(f"\n--- Creating Heatmap for Top {top_n} Variable Genes ---")
    if norm_counts.empty: return
    top_genes = norm_counts.var(axis=1).nlargest(top_n).index
    plot_data = norm_counts.loc[top_genes].dropna()
    print(f"Selected {plot_data.shape[0]} genes for the heatmap.")
    
    unique_conditions = sorted(metadata['Condition'].unique())
    palette = plt.cm.viridis(np.linspace(0, 1, len(unique_conditions)))
    color_map = {cond: color for cond, color in zip(unique_conditions, palette)}
    col_colors = metadata['Condition'].map(color_map)
    
    g = sns.clustermap(plot_data, z_score=0, cmap="vlag", center=0, col_cluster=False, row_cluster=True, col_colors=col_colors, yticklabels=True, figsize=(8, 12))
    heatmap_path = os.path.join(output_dir, "Heatmap_Top50_Variable_Genes.png")
    plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Top {top_n} variable genes heatmap saved to '{heatmap_path}'")
    print("-" * 40)

def generate_pca_plot(norm_counts, metadata, output_dir, pca_cfg={}):
    """
    Generates a PCA plot from normalized counts.
    Supports filtering for high variance genes based on config.
    """
    print("\n--- Step 6: Generating PCA Plot ---")
    
    if norm_counts.empty:
        print("Warning: Normalized counts matrix is empty. Skipping PCA.")
        return

    # Check for filtering
    filter_hvg = pca_cfg.get('filter_high_variance', False)
    top_n = pca_cfg.get('top_n_genes', 500)
    
    if filter_hvg:
        print(f"  - Filtering top {top_n} highest variance genes for PCA...")
        # Calculate variance
        gene_vars = norm_counts.var(axis=1)
        top_genes = gene_vars.nlargest(top_n).index
        data_subset = norm_counts.loc[top_genes]
        print(f"  - Selected {data_subset.shape[0]} genes.")
        plot_title = f'Principal Component Analysis (PCA)\nTop {top_n} Variable Genes'
    else:
        print("  - Using ALL genes for PCA (No filtering).")
        data_subset = norm_counts
        plot_title = 'Principal Component Analysis (PCA)\nAll Genes'

    if data_subset.empty:
        print("Warning: Data for PCA is empty after filtering. Skipping PCA.")
        return


    # PCA works on samples x features, so we transpose the counts matrix
    data = data_subset.T
    
    # Perform PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(data)
    
    # Create a new DataFrame for plotting
    pca_df = pd.DataFrame(data=principal_components, 
                          columns=['PC1', 'PC2'], 
                          index=data.index)
    
    # Merge with metadata to get conditions
    pca_df = pca_df.merge(metadata, left_index=True, right_index=True)
    
    # Get variance explained
    explained_variance = pca.explained_variance_ratio_
    
    # Plotting
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='PC1', y='PC2', hue='Condition', data=pca_df, s=100, alpha=0.8)
    
    # Add labels to points
    texts = []
    for i, row in pca_df.iterrows():
        texts.append(plt.text(row['PC1'] + 0.05, row['PC2'], row.name))
        
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

    plt.title(plot_title)
    plt.xlabel(f'PC1 ({explained_variance[0]:.1%})')
    plt.ylabel(f'PC2 ({explained_variance[1]:.1%})')
    plt.legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)
    
    # Save the plot
    filename_suffix = "TopHVG" if filter_hvg else "AllGenes"
    pca_plot_path = os.path.join(output_dir, f"PCA_plot_{filename_suffix}.png")
    plt.savefig(pca_plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"PCA plot saved to '{pca_plot_path}'")
    print("-" * 40)

def _get_timeseries_data_for_clustering(norm_counts, metadata, cfg):
    """Internal helper to prepare data for clustering and evaluation."""
    ts_cfg = cfg.get('analysis', {}).get('time_series', {})
    top_n = ts_cfg.get('top_n_genes', 200)
    time_col = ts_cfg.get('time_column', 'Time') # Default to 'Time'

    # Select only treatment samples
    treat_samples = metadata[metadata['Type'] == 'Treat'].index
    treat_counts = norm_counts[treat_samples]
    treat_meta = metadata.loc[treat_samples]

    # Average replicates based on time column
    if time_col not in treat_meta.columns:
        print(f"⚠️ Warning: Time column '{time_col}' not found in metadata. Using 'Time' if available.")
        time_col = 'Time' if 'Time' in treat_meta.columns else treat_meta.columns[1] # Fallback
    
    time_avg_counts = treat_counts.T.groupby(treat_meta[time_col]).mean().T
    
    # Check if we have enough time points
    if time_avg_counts.shape[1] < 2:
        print(f"⚠️ Warning: Not enough time points for time-series analysis. Found {time_avg_counts.shape[1]}, need at least 2.")
        return None
    
    # Select top N variable genes
    top_genes = time_avg_counts.var(axis=1).nlargest(top_n).index
    plot_data = time_avg_counts.loc[top_genes]

    # Z-score standardization
    # Use nan_policy='omit' or handle NaNs after.
    # If a gene has constant expression across time points, std is 0 -> zscore is NaN.
    # We fill these with 0 to avoid clustering errors.
    zscore_data = pd.DataFrame(stats.zscore(plot_data, axis=1, nan_policy='omit'), index=plot_data.index, columns=plot_data.columns)
    zscore_data = zscore_data.fillna(0)
    return zscore_data

# --- Time-Series K-Value Evaluation ---
def evaluate_clustering_k(norm_counts, metadata, cfg, max_k):
    """
    Calculates and plots Elbow Method and Silhouette Score to help choose the optimal k.
    """
    print("\n--- Evaluating Optimal K for Time-Series Clustering ---")
    data_for_clustering = _get_timeseries_data_for_clustering(norm_counts, metadata, cfg)
    
    if data_for_clustering is None:
        print("Skipping K-evaluation due to insufficient data.")
        return

    output_dir = cfg['files']['output_dir']
    
    k_range = range(2, max_k + 1)
    inertias = []
    silhouette_scores = []

    print(f"Testing k from 2 to {max_k}...")
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        clusters = kmeans.fit_predict(data_for_clustering)
        
        inertias.append(kmeans.inertia_)
        score = silhouette_score(data_for_clustering, clusters)
        silhouette_scores.append(score)
        print(f"  k={k}, Inertia={kmeans.inertia_:.2f}, Silhouette Score={score:.4f}")

    # Plotting
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Elbow Method Plot
    ax1.plot(k_range, inertias, 'bo-')
    ax1.set_xlabel('Number of Clusters (k)')
    ax1.set_ylabel('Inertia (WCSS)')
    ax1.set_title('Elbow Method for Optimal k')
    ax1.grid(True)

    # Silhouette Score Plot
    ax2.plot(k_range, silhouette_scores, 'go-')
    ax2.set_xlabel('Number of Clusters (k)')
    ax2.set_ylabel('Average Silhouette Score')
    ax2.set_title('Silhouette Score for Optimal k')
    ax2.grid(True)
    
    plt.tight_layout()
    eval_plot_path = os.path.join(output_dir, "k_evaluation_plot.png")
    plt.savefig(eval_plot_path, dpi=300)
    plt.close(fig)
    print(f"Cluster evaluation plot saved to '{eval_plot_path}'")
    print("-" * 40)

# --- Step 8: Time-Series Analysis ---
def run_timeseries_analysis(norm_counts, metadata, cfg, output_dir):
    """
    Performs time-series clustering and generates a combined trend and heatmap plot.
    Optimized for dynamic layout and better visualization.
    """
    print("\n--- Step 8: Time-Series Analysis ---")
    ts_cfg = cfg.get('analysis', {}).get('time_series', {})
    if not ts_cfg.get('enabled', False):
        print("Time-series analysis is disabled in the config. Skipping.")
        return None

    # --- 1. Get Parameters ---
    n_clusters = ts_cfg.get('num_clusters', 4)
    # output_dir is passed as argument
    output_file = os.path.join(output_dir, f"TimeSeries_Clustering_k{n_clusters}.png")

    print(f"Running with {n_clusters} clusters...")

    # --- 2. Prepare Data (With Robust Sorting) ---
    zscore_data = _get_timeseries_data_for_clustering(norm_counts, metadata, cfg)
    
    if zscore_data is None:
        print("Skipping Time-Series Analysis due to insufficient data.")
        return None
    
    # Try to sort numerically if possible, otherwise string sort
    try:
        sorted_cols = sorted(zscore_data.columns, key=lambda x: float(x))
        is_numeric_time = True
    except ValueError:
        sorted_cols = sorted(zscore_data.columns) # Fallback to string sort
        is_numeric_time = False
        
    zscore_data = zscore_data[sorted_cols] # Reorder dataframe columns
    
    # Add suffix only if numeric, otherwise use value as is
    if is_numeric_time:
        time_points = [f"{t}h" for t in sorted_cols]
    else:
        time_points = [str(t) for t in sorted_cols]
        
    zscore_data.columns = time_points

    # --- 3. K-Means Clustering ---
    print(f"Performing K-Means clustering...")
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(zscore_data)
    zscore_data['Cluster'] = clusters
    
    # --- 4. Calculate Membership ---
    print("Calculating cluster membership...")
    membership = []
    for i in range(n_clusters):
        # 僅取數值欄位計算
        cluster_indices = zscore_data[zscore_data['Cluster'] == i].index
        cluster_data_numeric = zscore_data.loc[cluster_indices, time_points]
        
        if not cluster_data_numeric.empty:
            centroid = kmeans.cluster_centers_[i]
            distances = np.linalg.norm(cluster_data_numeric.values - centroid, axis=1)
            
            min_dist, max_dist = distances.min(), distances.max()
            # 避免除以零
            if max_dist == min_dist:
                normalized_membership = np.ones_like(distances)
            else:
                # 距離越小，Membership 越高
                normalized_membership = 1 - (distances - min_dist) / (max_dist - min_dist)
            
            for gene_idx, score in zip(cluster_indices, normalized_membership):
                membership.append({'gene': gene_idx, 'membership': score})

    if membership:
        membership_df = pd.DataFrame(membership).set_index('gene')
        zscore_data = zscore_data.merge(membership_df, left_index=True, right_index=True, how='left')
        zscore_data['membership'] = zscore_data['membership'].fillna(0) # Fill NaN just in case
    else:
        zscore_data['membership'] = 0.5 

    # 排序：先按 Cluster，再按 Membership (為了 Heatmap 美觀)
    zscore_data = zscore_data.sort_values(by=['Cluster', 'membership'], ascending=[True, False])
    
    # --- 5. Plotting (Dynamic Layout Fix) ---
    print("Generating plots...")
    plt.style.use('default')
    
    # 動態計算高度：如果有 8 個 clusters，圖需要高一點
    n_cols_plot = 2
    n_rows_plot = math.ceil(n_clusters / n_cols_plot)
    fig_height = max(8, 3 * n_rows_plot) 
    
    fig = plt.figure(figsize=(15, fig_height))
    fig.suptitle(f'Time-Series Gene Expression Clustering (k={n_clusters})', fontsize=16, fontweight='bold')
    
    # 主網格：左邊 Heatmap，右邊 Trend Plots
    gs_main = fig.add_gridspec(1, 2, width_ratios=[1, 1.3], wspace=0.25)

    # --- Part A: Heatmap ---
    ax_heatmap = fig.add_subplot(gs_main[0, 0])
    heatmap_vals = zscore_data[time_points]
    cluster_ids = zscore_data['Cluster']
    
    sns.heatmap(heatmap_vals, ax=ax_heatmap, cmap='vlag', center=0, 
                cbar_kws={'label': 'Z-score', 'orientation': 'horizontal', 'pad': 0.05, 'shrink': 0.6},
                yticklabels=False)
    ax_heatmap.set_title("Clustered Heatmap (Sorted by Membership)", fontsize=14, fontweight='bold')
    ax_heatmap.set_xlabel("Time Points")
    ax_heatmap.set_ylabel(f"Genes (n={len(heatmap_vals)})")

    # Heatmap 標註 Cluster 範圍
    y_pos = 0
    for i in sorted(cluster_ids.unique()):
        cluster_gene_count = (cluster_ids == i).sum()
        if cluster_gene_count > 0:
            # 右側文字標籤
            ax_heatmap.text(heatmap_vals.shape[1] * 1.02, y_pos + cluster_gene_count / 2, f"C{i+1}", 
                            verticalalignment='center', fontsize=12, fontweight='bold')
            # 白色分隔線
            ax_heatmap.hlines(y_pos + cluster_gene_count, *ax_heatmap.get_xlim(), colors='white', linewidth=2)
            y_pos += cluster_gene_count

    # --- Part B: Trend Plots (Dynamic) ---
    gs_trends = gs_main[0, 1].subgridspec(n_rows_plot, n_cols_plot, hspace=0.4, wspace=0.25)
    
    trend_axes = []
    # 依序建立子圖，並處理最後一行可能空白的情況
    for r in range(n_rows_plot):
        for c in range(n_cols_plot):
            if len(trend_axes) < n_clusters:
                trend_axes.append(fig.add_subplot(gs_trends[r, c]))
            else:
                # 隱藏多餘的格子
                fig.add_subplot(gs_trends[r, c]).axis('off')

    cmap = plt.get_cmap('viridis')
    
    for i, ax in enumerate(trend_axes):
        cluster_subset = zscore_data[zscore_data['Cluster'] == i]
        
        if not cluster_subset.empty:
            # 1. 繪製所有基因的線條 (根據 membership 上色)
            # 使用 alpha=0.4 讓重疊處更明顯
            for _, row in cluster_subset.iterrows():
                ax.plot(time_points, row[time_points], color=cmap(row['membership']), alpha=0.3, linewidth=1)
            
            # 2. 【新增】繪製平均趨勢線 (Mean Trend) - 黑色粗虛線
            mean_trend = cluster_subset[time_points].mean()
            ax.plot(time_points, mean_trend, color='black', linewidth=2.5, linestyle='--', label='Mean')

        ax.set_title(f"Cluster {i+1} (n={len(cluster_subset)})", fontsize=12, fontweight='bold')
        ax.set_ylim(-3, 3) # 固定 Y 軸範圍以便比較
        
        # 只在最後一列顯示 X 軸標籤，保持整潔 (可選)
        # if i < n_clusters - 2: ax.set_xticklabels([]) 
        
        ax.tick_params(axis='x', rotation=45)
        if i % n_cols_plot == 0:
            ax.set_ylabel("Z-score")
            
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # --- Part C: Membership Colorbar ---
    # 調整 Colorbar 位置，放在 Heatmap 下方或右側空位
    cax = fig.add_axes([0.92, 0.3, 0.015, 0.4]) 
    norm = mcolors.Normalize(vmin=0, vmax=1)
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cbar.set_label('Membership Score (Centrality)')

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Time-series analysis plot saved to '{output_file}'")
    print("-" * 40)
    return zscore_data


# --- Main Execution Logic ---
if __name__ == "__main__":
    import sys
    
    # Load configuration
    config_path = 'config.yml'
    if len(sys.argv) > 1:
        config_path = sys.argv[1]
        print(f"Using config file: {config_path}")
    
    config = load_config(config_path)
    files_cfg = config['files']
    analysis_cfg = config['analysis']
    
    # Ensure output directories exist
    output_dir = files_cfg['output_dir']
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs("genesets", exist_ok=True)
    os.makedirs(os.path.join(output_dir, "enrichment"), exist_ok=True)
    print(f"--- Outputs will be saved to '{output_dir}' ---")

    # 1. Preprocess
    counts_matrix = preprocess_data(
        feature_counts_file=files_cfg.get('input', files_cfg.get('feature_counts')),
        gene_annotation_file=files_cfg['gene_annotation'],
        output_dir=output_dir
    )
    
    # 2. Load metadata
    metadata = load_metadata(
        metadata_file=files_cfg.get('metadata', 'input/metadata.csv'),
        output_dir=output_dir,
        cfg=config
    )
    
    # Check for sample consistency
    missing_samples = [s for s in metadata.index if s not in counts_matrix.columns]
    if missing_samples:
        print(f"⚠️ Warning: The following samples in metadata are missing from counts matrix: {missing_samples}")
        # Optional: Filter metadata to match counts
        metadata = metadata.loc[counts_matrix.columns.intersection(metadata.index)]
    
    # Align counts to metadata
    counts_matrix = counts_matrix[metadata.index]
    
    if metadata.empty or counts_matrix.empty:
        raise ValueError(
            "CRITICAL ERROR: No matching samples found between Count Matrix and Metadata.\n"
            "Please check your Metadata 'SampleID' column. They MUST match the column names in your count matrix.\n"
            f"Count Matrix Columns: {list(counts_matrix.columns)[:5]} ...\n"
            f"Metadata SampleIDs: {list(metadata.index)[:5]} ..."
        )

    
    # 3. Initialize DESeq2 object
    # Ensure counts are integers (raw counts are required for DESeq2/pyDEG)
    print("Initializing DESeq2 object (ensuring integer counts)...")
    counts_matrix_int = counts_matrix.round().astype(int)
    dds = ov.bulk.pyDEG(counts_matrix_int)
    dds.drop_duplicates_index()
    
    # 4. Run batch analysis (DEG + Enrichment)
    # NOTE: We use the ORIGINAL raw counts (dds) for DE analysis as per best practices.
    all_results = run_batch_deg_analysis(
        dds, metadata, 
        cfg=config, # Pass the whole config object
        gene_annotation_file=files_cfg['gene_annotation']
    )
    
    # 5. Normalize data for plotting
    print("\n--- Step 5: Normalizing Data for Plots ---")
    
    # Check if batch correction is enabled
    batch_cfg = files_cfg.get('batch_correction', {}) # Check in files first (legacy) or root? 
    # Actually config structure is:
    # files: ...
    # batch_correction: ... (root level based on my edit)
    # Let's check safely
    batch_cfg = config.get('batch_correction', {})
    
    plotting_counts = None
    
    if batch_cfg.get('enabled', False):
        batch_key = batch_cfg.get('batch_key', 'Batch')
        print(f"Batch correction enabled. Key: {batch_key}")
        
        # Perform correction
        # We pass the raw counts_matrix
        corrected_counts = perform_batch_correction(
            counts_matrix, metadata, batch_key, output_dir
        )
        
        # For plotting (PCA, Heatmap), we use the corrected data.
        # We assume ComBat returns expression-like values. 
        # If they are still count-like, we might need log transform.
        # Usually ComBat on counts (if using specific methods) or log-counts.
        # omicverse's batch_correction documentation implies it standardizes data.
        # Let's inspect the output range in the function, but here we assume 
        # we can use it for visualization directly or with log2.
        
        # Let's treat the corrected_counts as our "norm_counts" basis.
        # If max value is large, we log2 it.
        if corrected_counts.max().max() > 100:
             norm_counts = np.log2(corrected_counts + 1)
        else:
             norm_counts = corrected_counts
             
        print("✅ Using BATCH-CORRECTED data for PCA and Heatmaps.")
        
    else:
        # Standard Normalization (Log2CPM)
        raw_counts = dds.data 
        library_size = raw_counts.sum(axis=0)
        cpm = raw_counts.div(library_size, axis=1) * 1e6
        norm_counts = np.log2(cpm + 1)
        print('✅ Normalization complete (Log2CPM) on RAW data.')

    norm_path = os.path.join(output_dir, "normalized_counts_for_plotting.csv")
    norm_counts.to_csv(norm_path)
    print(f"Normalized plotting data saved to '{norm_path}'")
    print("-" * 40)

    # 6. Generate PCA plot
    generate_pca_plot(norm_counts, metadata, output_dir, analysis_cfg.get('pca', {}))
    
    # --- Optional: Evaluate optimal k for time-series clustering ---
    ts_cfg = analysis_cfg.get('time_series', {})
    k_eval_cfg = ts_cfg.get('k_evaluation', {})
    
    # Only run if Time Series is enabled AND K-eval is enabled AND Time column is present
    if ts_cfg.get('enabled', True) and ts_cfg.get('time_column_present', True) and k_eval_cfg.get('enabled', False):
        print("\n--- Running K-Value Evaluation for Time-Series Clustering ---")
        evaluate_clustering_k(
            norm_counts, 
            metadata, 
            config, 
            max_k=k_eval_cfg.get('max_k', 10)
        )

    # 7. Time-Series Analysis (Clustering)
    cluster_results = None # Initialize to None
    if ts_cfg.get('enabled', True) and ts_cfg.get('time_column_present', True):
        cluster_results = run_timeseries_analysis(norm_counts, metadata, config, output_dir)
    else:
        print("\n--- Time-Series Analysis Skipped (Disabled or Time Column Missing) ---")

    # 8. Create heatmap and master table (now with cluster info)
    if all_results:
        create_master_and_heatmap(
            results_dict=all_results, norm_counts=norm_counts, metadata=metadata,
            output_dir=output_dir, 
            padj_thresh=analysis_cfg['significance']['padj_threshold'],
            log2fc_thresh=analysis_cfg['significance']['log2fc_threshold'],
            cluster_results=cluster_results
        )

    # 7. Create top variable heatmap
    create_top_variable_heatmap(
        norm_counts=norm_counts, metadata=metadata, output_dir=output_dir,
        top_n=analysis_cfg['variable_heatmap']['top_n_genes']
    )

    print("\n🎉 Analysis Pipeline Completed! 🎉")