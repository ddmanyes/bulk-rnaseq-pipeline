import shutil
import sys
import streamlit as st
import pandas as pd
import numpy as np
import omicverse as ov
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import yaml
import subprocess
import time
import tempfile

# Set Pandas Styler limit to avoid errors with large tables
pd.set_option("styler.render.max_elements", 2000000)

# --- Global Constants ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# --- Session Isolation ---
if 'temp_dir' not in st.session_state:
    st.session_state.temp_dir = tempfile.mkdtemp(prefix="rnaseq_session_")
    # Initialize default directories
    os.makedirs(os.path.join(st.session_state.temp_dir, "input"), exist_ok=True)
    os.makedirs(os.path.join(st.session_state.temp_dir, "output"), exist_ok=True)
    
    # Copy default config to session dir
    
    # Copy default config to session dir
    # Copy default config to session dir
    default_config_path = os.path.join(BASE_DIR, "config.yml")
    if os.path.exists(default_config_path):
        shutil.copy(default_config_path, os.path.join(st.session_state.temp_dir, "config.yml"))
    else:
        # Fallback: Create a minimal config if file is missing to avoid crash
        minimal_config = {
            'files': {'input': '', 'metadata': '', 'output_dir': 'output'},
            'analysis': {
                'diff_expression': {}, 'significance': {}, 'enrichment': {}, 'time_series': {}
            }
        }
        with open(os.path.join(st.session_state.temp_dir, "config.yml"), 'w') as f:
            yaml.dump(minimal_config, f)

    # NEW: Sync project 'input' directory to session 'input' directory
    # Moved to sync_input_files() function to run on every script execution

# Helper to get session paths
def get_session_path(rel_path):
    return os.path.join(st.session_state.temp_dir, rel_path)

# --- Sessions & Sync ---
def sync_input_files():
    """Syncs files from project input dir to session input dir."""
    project_input_dir = os.path.join(BASE_DIR, "input")
    session_input_dir = get_session_path("input")
    
    if not os.path.exists(session_input_dir):
        os.makedirs(session_input_dir)
    
    if os.path.exists(project_input_dir):
        for item in os.listdir(project_input_dir):
            if item.startswith('.'): continue # Skip hidden files
            s = os.path.join(project_input_dir, item)
            d = os.path.join(session_input_dir, item)
            if os.path.isfile(s):
                # Update only if destination doesn't exist (don't overwrite user session data)
                if not os.path.exists(d):
                    try:
                        shutil.copy2(s, d)
                    except Exception as e:
                        print(f"Error syncing {item}: {e}")

# Execute Sync on every rerun
sync_input_files()

# --- Configuration & Helper Functions ---
CONFIG_FILE = get_session_path("config.yml")

def load_config():
    """Loads the YAML config file from session dir."""
    # ... (existing load_config code) ...
    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, 'r') as f:
            return yaml.safe_load(f)
    return {}

def save_config(config_data):
    """Saves the dictionary to the YAML config file in session dir."""
    with open(CONFIG_FILE, 'w') as f:
        yaml.dump(config_data, f, default_flow_style=False, sort_keys=False)

def run_pipeline():
    """Executes the analysis pipeline script with session config."""
    try:
        # Run the script and capture output
        # Pass the session config path as argument
        process = subprocess.Popen(
            [sys.executable, os.path.join(BASE_DIR, "analysis_pipeline.py"), CONFIG_FILE],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        return process
    except Exception as e:
        st.error(f"Failed to start pipeline: {e}")
        return None

# --- Main App Layout ---
st.set_page_config(layout="wide", page_title="RNA-seq Pipeline GUI")
st.title("🧬 Bulk RNA-seq Analysis Pipeline & Interactive Explorer v1.2")

# --- Sidebar ---
with st.sidebar:
    st.header("🔧 Utility")
    
    # ... (Demo Data section) ...
    st.warning("⚠️ **Privacy Note**: This app runs on a public cloud. Do NOT upload PII or confidential data.")
    if os.path.exists("demo_data"):
        if st.button("📥 Load Demo Data"):
            st.caption("Loads pre-computed results for immediate visualization.")
            try:
                # Copy input to session dir
                session_input = get_session_path("input")
                if not os.path.exists(session_input): os.makedirs(session_input)
                
                if os.path.exists("demo_data/input/metadata.csv"):
                    shutil.copy("demo_data/input/metadata.csv", os.path.join(session_input, "metadata.csv"))
                
                # Copy output to session dir
                session_output = get_session_path("output")
                if os.path.exists(session_output): shutil.rmtree(session_output)
                shutil.copytree("demo_data/demo_results", session_output)
                
                # Update config to point to session dirs
                cfg = load_config()
                cfg['files']['input'] = os.path.join(session_input, "merged_counts.csv") # Placeholder
                cfg['files']['metadata'] = os.path.join(session_input, "metadata.csv")
                cfg['files']['output_dir'] = session_output
                save_config(cfg)
                
                st.success("Demo data loaded! Please refresh or go to Visualization tabs.")
                time.sleep(1)
                st.rerun()
            except Exception as e:
                st.error(f"Failed to load demo data: {e}")
        
        st.markdown("---")
        st.markdown("""
        **Acknowledgement**
        
        Demo data courtesy of **Prof. Sung-Jan Lin's Lab**, 
        *Institute of Biomedical Engineering, NTU*.
        """)

# Load initial config
config = load_config()

# Ensure config has necessary keys (Robustness Check)
if 'files' not in config: config['files'] = {}
if 'analysis' not in config: config['analysis'] = {}
if 'input' not in config['files']: config['files']['input'] = ""
if 'metadata' not in config['files']: config['files']['metadata'] = ""
if 'output_dir' not in config['files']: config['files']['output_dir'] = "output"

# Create Tabs
tab_merge, tab0, tab1, tab3, tab_volcano, tab2, tab4 = st.tabs(["📂 Data Preprocessing", "📝 Metadata Creation", "⚙️ Run Analysis", "📈 Static Plots", "🌋 Volcano Plot", "📊 Interactive Visualization", "🧬 Pathway Analysis"])

# ==========================================
# TAB MERGE: Data Preprocessing
# ==========================================
with tab_merge:
    st.header("Data Preprocessing (Merge & Clean)")
    st.markdown("Upload multiple count files (CSV, Excel, featureCounts) to clean and merge them into a single matrix for analysis.")
    
    st.markdown("Combine multiple count files (CSV, Excel, featureCounts) into a single matrix for analysis.")
    
    # File Uploader
    uploaded_files = st.file_uploader("Upload Count Files", type=['csv', 'xlsx', 'txt', 'tsv'], accept_multiple_files=True)
    
    if uploaded_files:
        st.success(f"Uploaded {len(uploaded_files)} files.")
        
        output_filename = st.text_input("Output Filename", value="merged_counts.csv")

        # Show success message if persistent
        if 'merge_success' in st.session_state:
            st.success(st.session_state['merge_success'])
            # Provide Download Button (Re-read file)
            merged_path = st.session_state.get('merge_file')
            if merged_path and os.path.exists(merged_path):
                 with open(merged_path, 'rb') as f:
                    st.download_button(
                        label="📥 Download Merged Matrix (CSV)",
                        data=f,
                        file_name=os.path.basename(merged_path),
                        mime="text/csv",
                        type="primary"
                    )
            st.caption("Note: Go to 'Metadata Creation' or 'Run Analysis' tab to use this file.")
        
        if st.button("Merge Files"):
            try:
                merged_df = pd.DataFrame()
                
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                for i, uploaded_file in enumerate(uploaded_files):
                    status_text.text(f"Processing {uploaded_file.name}...")
                    
                    # Determine file type and read
                    if uploaded_file.name.endswith('.csv'):
                        df = pd.read_csv(uploaded_file, index_col=0)
                    elif uploaded_file.name.endswith('.xlsx'):
                        df = pd.read_excel(uploaded_file, index_col=0)
                    elif uploaded_file.name.endswith(('.txt', '.tsv')):
                        df = pd.read_csv(uploaded_file, sep='\t', comment='#')
                        # Handle featureCounts format
                        if 'Geneid' in df.columns:
                            df = df.set_index('Geneid')
                        elif df.index.name != 'Geneid' and 'Geneid' not in df.columns and df.shape[1] > 1: 
                                # Assume first column is index if numeric data follows
                                if df.iloc[:,0].dtype.kind in 'biufc': 
                                    pass 
                                else:
                                    pass

                        # Clean featureCounts columns
                        if 'Chr' in df.columns and 'Start' in df.columns:
                            df = df.iloc[:, 5:] 
                            df.columns = [c.split('/')[-1].split('\\')[-1].replace('.sorted.bam', '').replace('.bam', '') for c in df.columns]

                    # Handle duplicate columns
                    new_columns = []
                    for col in df.columns:
                        if not merged_df.empty and col in merged_df.columns:
                            new_columns.append(f"{col}_{uploaded_file.name}")
                        else:
                            new_columns.append(col)
                    df.columns = new_columns
                    
                    # Merge
                    if merged_df.empty:
                        merged_df = df
                    else:
                        merged_df = merged_df.join(df, how='outer')
                    
                    progress_bar.progress((i + 1) / len(uploaded_files))
                
                # Fill NAs and Save
                merged_df = merged_df.fillna(0)
                
                session_input = get_session_path("input")
                if not os.path.exists(session_input):
                    os.makedirs(session_input)
                output_path_session = os.path.join(session_input, output_filename)
                merged_df.to_csv(output_path_session)

                st.session_state['merge_success'] = f"Successfully merged {len(uploaded_files)} files to `{output_path_session}`"
                st.session_state['merge_file'] = output_path_session
                st.rerun()

            except Exception as e:
                # Indentation fixed v1.2
                st.error(f"An error occurred during merging: {e}")

    else:
        st.info("Please upload files to begin.")

# ==========================================
# TAB 0: Metadata Creation
# ==========================================
with tab0:
    st.header("Metadata Editor")
    st.markdown("Create or edit your sample metadata file here. This file is required for the analysis.")

    # --- 1. File Uploader ---
    uploaded_meta = st.file_uploader("📂 Import Metadata File (Optional)", type=['csv', 'xlsx'], help="Upload a CSV or Excel file to populate the table below. It should have SampleID as the first column (or index).")

    # Path to metadata file (from config or default)
    # Ensure we use the Session Path for robust read/write permissions
    config_meta_path = config.get('files', {}).get('metadata', '')
    if config_meta_path and not os.path.isabs(config_meta_path):
        # If relative, anchor to session_input_dir or just session_dir
        # Check if it starts with 'input/'
        if config_meta_path.startswith("input/"):
             base_name = os.path.basename(config_meta_path)
             meta_file_path = get_session_path(f"input/{base_name}")
        else:
             meta_file_path = get_session_path(config_meta_path)
    elif os.path.isabs(config_meta_path):
        meta_file_path = config_meta_path
    else:
        meta_file_path = get_session_path('input/metadata.csv')

    # Load logic
    existing_df = None
    
    # Priority 1: Uploaded File (if selected/present)
    if uploaded_meta:
        try:
            if uploaded_meta.name.endswith('.xlsx'):
                existing_df = pd.read_excel(uploaded_meta, index_col=0)
            else:
                existing_df = pd.read_csv(uploaded_meta, index_col=0)
            st.success(f"Previewing uploaded metadata: `{uploaded_meta.name}`")
        except Exception as e:
            st.error(f"Error reading uploaded file: {e}")
            existing_df = None

    # Priority 2: Existing File on Disk (only if no upload or upload failed)
    if existing_df is None:
        if os.path.exists(meta_file_path):
            st.info(f"Loading existing metadata from: `{meta_file_path}`")
            try:
                existing_df = pd.read_csv(meta_file_path, index_col=0)
            except Exception as e:
                st.error(f"Error loading metadata: {e}")
                existing_df = None

    # Priority 3: Default Template
    if existing_df is None:
        if not uploaded_meta: # Only show warning if not trying to upload
             st.warning(f"File `{meta_file_path}` not found. Starting with a template.")
        existing_df = pd.DataFrame({
            'Condition': ['Control', 'Treat'],
            'Time': [0, 6],
            'Type': ['Control', 'Treat']
        }, index=['Sample1', 'Sample2'])
        existing_df.index.name = 'SampleID'

    # --- Sync with Count Matrix ---
    st.subheader("Sync with Count Data")
    
    # Dynamic discovery of count files
    session_input_dir = get_session_path("input")
    count_candidates = [f for f in os.listdir(session_input_dir) if f.endswith(('.csv', '.xlsx', '.txt', '.tsv'))] if os.path.exists(session_input_dir) else []
    
    # Determine default selection
    default_idx = 0
    # Prioritize config if valid
    current_config_input = config.get('files', {}).get('input', '')
    current_config_basename = os.path.basename(current_config_input) if current_config_input else ""
    
    if current_config_basename in count_candidates:
        default_idx = count_candidates.index(current_config_basename)
    elif "kallisto_gene_counts.csv" in count_candidates:
         default_idx = count_candidates.index("kallisto_gene_counts.csv")
    elif "merged_counts.csv" in count_candidates:
         default_idx = count_candidates.index("merged_counts.csv")
    
    selected_count_file = st.selectbox(
        "Select Count Matrix to Sync From",
        options=count_candidates,
        index=default_idx if count_candidates else 0,
        help="Select the count matrix file (CSV/Excel) containing Sample IDs."
    ) if count_candidates else None
    
    if selected_count_file:
         count_file_path = os.path.join(session_input_dir, selected_count_file)
    else:
         count_file_path = None

    if st.button("🔄 Sync Sample IDs", help="Read column names from the selected matrix and update Metadata rows."):
        if count_file_path and os.path.exists(count_file_path):
            try:
                # Read only the header to get columns
                if count_file_path.endswith('.xlsx'):
                     counts_df = pd.read_excel(count_file_path, index_col=0, nrows=0)
                else:
                     # Detect separator
                     try:
                        counts_df = pd.read_csv(count_file_path, index_col=0, nrows=0)
                        if len(counts_df.columns) < 1: # Failed to parse proper cols, try tab
                            counts_df = pd.read_csv(count_file_path, index_col=0, sep='\t', nrows=0)
                     except:
                        counts_df = pd.read_csv(count_file_path, index_col=0, sep='\t', nrows=0)

                sample_ids = counts_df.columns.tolist()
                
                if sample_ids:
                    # Reindex existing metadata to match count matrix samples
                    existing_df = existing_df.reindex(sample_ids)
                    existing_df.index.name = 'SampleID'
                    existing_df = existing_df.fillna("")
                    
                    # Save the updated metadata to file so it persists after rerun
                    # Ensure directory exists
                    os.makedirs(os.path.dirname(meta_file_path), exist_ok=True)
                    existing_df.to_csv(meta_file_path)
                    
                    # Update config input path if sync is successful
                    config['files']['input'] = count_file_path
                    save_config(config)
                    
                    st.success(f"Synced {len(sample_ids)} samples from `{selected_count_file}`.")
                    # Rerun to update the table editor below
                    st.rerun()
                else:
                    st.warning("Count matrix seems to have no sample columns.")
            except Exception as e:
                st.error(f"Failed to read count matrix: {e}")
        else:
            st.error("No valid count file selected or file not found.")

    # Editable Dataframe
    df_to_edit = existing_df.reset_index()
    
    # Add 'Include' column for row selection (default True)
    if 'Include' not in df_to_edit.columns:
        df_to_edit.insert(0, 'Include', True)
    
    # Get Time Column Name from config
    time_col = config.get('analysis', {}).get('time_series', {}).get('time_column', 'Time')
    
    # Toggle for Time Column
    use_time_col = st.checkbox(f"Include Time Column ('{time_col}')", value=config.get('analysis', {}).get('time_series', {}).get('time_column_present', True))
    
    st.markdown("### Edit Metadata Table")
    req_cols = "**Condition**, **Type**" + (f", **{time_col}**" if use_time_col else "")
    st.markdown(f"Check **Include** to save specific rows. Required columns: {req_cols}.")
    
    # Ensure object columns are strings to avoid Streamlit TextColumn error with NaNs (floats)
    if 'Condition' in df_to_edit.columns:
        df_to_edit['Condition'] = df_to_edit['Condition'].fillna('').astype(str)
    else:
        df_to_edit['Condition'] = ""

    if 'Type' in df_to_edit.columns:
         df_to_edit['Type'] = df_to_edit['Type'].fillna('Control').astype(str) # Default to Control if missing
    else:
         df_to_edit['Type'] = "Control"

    # Dynamic column config
    cols_cfg = {
        "Include": st.column_config.CheckboxColumn("Include", help="Check to save this sample.", default=True),
        "SampleID": st.column_config.TextColumn("Sample ID (Must match count matrix columns)", required=True),
        "Condition": st.column_config.TextColumn("Condition (e.g. Control, Treat)", required=True),
        "Type": st.column_config.SelectboxColumn("Type", options=["Control", "Treat"], required=True)
    }
    
    # Handle Time Column
    if use_time_col:
        if time_col in df_to_edit.columns:
            df_to_edit[time_col] = df_to_edit[time_col].astype(str)
            cols_cfg[time_col] = st.column_config.TextColumn(f"{time_col} (e.g. 0, 6h, Day1)", required=True)
        elif 'Time' in df_to_edit.columns:
            df_to_edit = df_to_edit.rename(columns={'Time': time_col})
            df_to_edit[time_col] = df_to_edit[time_col].astype(str)
            cols_cfg[time_col] = st.column_config.TextColumn(f"{time_col} (e.g. 0, 6h, Day1)", required=True)
        else:
            df_to_edit[time_col] = ""
            cols_cfg[time_col] = st.column_config.TextColumn(f"{time_col} (e.g. 0, 6h, Day1)", required=True)
    else:
        # If disabled, remove from view if present
        if time_col in df_to_edit.columns:
            df_to_edit = df_to_edit.drop(columns=[time_col])
        if 'Time' in df_to_edit.columns:
            df_to_edit = df_to_edit.drop(columns=['Time'])

    edited_df = st.data_editor(
        df_to_edit, 
        num_rows="dynamic", 
        use_container_width=True,
        column_config=cols_cfg,
        key="metadata_editor"
    )

    if st.button("💾 Save Metadata"):
        try:
            # Filter rows where Include is True
            final_df = edited_df[edited_df['Include'] == True].copy()
            final_df = final_df.drop(columns=['Include'])
            
            # Check for empty SampleIDs
            if final_df['SampleID'].isnull().any() or (final_df['SampleID'] == '').any():
                st.error("Error: SampleID cannot be empty.")
            elif final_df['SampleID'].duplicated().any():
                st.error("Error: Duplicate SampleIDs found.")
            else:
                # Set index back to SampleID
                final_df = final_df.set_index('SampleID')
                
                # Ensure directory exists
                os.makedirs(os.path.dirname(meta_file_path), exist_ok=True)
                
                final_df.to_csv(meta_file_path)
                
                # Update config for time column presence
                config['analysis']['time_series']['time_column_present'] = use_time_col
                save_config(config)
                
                st.success(f"Metadata saved successfully to `{meta_file_path}`! ({len(final_df)} samples)")
                
                # Show full table
                st.dataframe(final_df, use_container_width=True)
                
                # Provide Download Button
                meta_csv = final_df.to_csv().encode('utf-8')
                st.download_button(
                    label="📥 Download Metadata (CSV)",
                    data=meta_csv,
                    file_name="metadata.csv",
                    mime="text/csv"
                )
                
        except Exception as e:
            st.error(f"Failed to save metadata: {e}")

# ==========================================
# TAB 1: Run Analysis (Configuration & Execution)
# ==========================================
with tab1:
    st.header("Pipeline Configuration")
    
    # --- Configuration Section (No Form) ---
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("1. File Paths")
    
        # Default to merged_counts.csv
        # Input File Selection
        session_input_dir = get_session_path("input")
        input_candidates = [f for f in os.listdir(session_input_dir) if f.endswith(('.csv', '.xlsx', '.txt', '.tsv'))] if os.path.exists(session_input_dir) else []
        
        default_candidate_index = 0
        current_input_basename = os.path.basename(config.get('files', {}).get('input', ''))
        
        if current_input_basename in input_candidates:
             default_candidate_index = input_candidates.index(current_input_basename)
        elif "kallisto_gene_counts.csv" in input_candidates:
            default_candidate_index = input_candidates.index("kallisto_gene_counts.csv")
        elif "merged_counts.csv" in input_candidates:
            default_candidate_index = input_candidates.index("merged_counts.csv")
            
        if input_candidates:
            selected_input_file = st.selectbox(
                "Select Input Count Matrix", 
                options=input_candidates, 
                index=default_candidate_index,
                key="run_analysis_input_select"
            )
            default_input_path = os.path.join(session_input_dir, selected_input_file)
            st.success(f"Selected input file: `{default_input_path}`")
        else:
            default_input_path = get_session_path("input/merged_counts.csv")
            st.warning("No input files found in directory. Please upload/merge files first.")
        
        config['files']['input'] = default_input_path

        config['files']['metadata'] = st.text_input("Metadata File", value=config.get('files', {}).get('metadata', get_session_path('input/metadata.csv')))

        
        # Gene Annotation Selection
        genesets_dir = "genesets"
        if os.path.exists(genesets_dir):
            available_genesets = [f for f in os.listdir(genesets_dir) if f.endswith('.tsv')]
            available_genesets.sort()
        else:
            available_genesets = []
            
        default_geneset = config.get('files', {}).get('gene_annotation', 'genesets/pair_GRCm39.tsv')
        # Extract filename from path if needed
        default_geneset_name = os.path.basename(default_geneset)
        
        if default_geneset_name not in available_genesets and available_genesets:
            default_index = 0
        elif default_geneset_name in available_genesets:
            default_index = available_genesets.index(default_geneset_name)
        else:
            default_index = 0

        selected_geneset = st.selectbox(
            "Gene Annotation File", 
            options=available_genesets,
            index=default_index,
            help="Select the appropriate gene annotation file for your species."
        )
        
        gene_annotation_path = os.path.join(genesets_dir, selected_geneset) if selected_geneset else ""
        
        st.info("""
        **Species Guide:**
        *   🐭 **Mouse**: `pair_GRCm39.tsv`
        *   👱 **Human**: `pair_GRCh38.tsv` (or `GRCh37`, `T2TCHM13`)
        *   🐟 **Zebrafish**: `pair_danRer11.tsv` (or `danRer7`)
        """)
        output_dir_path = st.text_input("Output Directory", value=config.get('files', {}).get('output_dir', get_session_path('output')))

        st.subheader("🧪 Batch Correction")
        # Batch Correction
        batch_cfg = config.get('files', {}).get('batch_correction', {})
        batch_enabled = st.checkbox("Enable Batch Correction", value=batch_cfg.get('enabled', False))
        batch_key = st.text_input("Batch Key (Column in Metadata)", value=batch_cfg.get('batch_key', 'Batch'))
        st.info("ℹ️ **Method**: Uses **ComBat** (via `omicverse` / `combat`) to adjust for batch effects. It normalizes the data to remove technical variation while preserving biological differences.")

    with col2:
        st.subheader("🔬 Analysis Parameters")
        # Analysis
        analysis_cfg = config.get('analysis', {})
        
        # Diff Expression
        diff_cfg = analysis_cfg.get('diff_expression', {})
        control_group = st.text_input("Control Group Name", value=diff_cfg.get('control_group', 'Control'))
        filter_low = st.checkbox(
            "Filter Low Expression Genes", 
            value=diff_cfg.get('filter_low_expression', True),
            help="If enabled, filters out genes with log2(BaseMean) <= 1. This helps remove noisy, lowly expressed genes."
        )
        
        # Significance
        sig_cfg = analysis_cfg.get('significance', {})
        padj_threshold = st.number_input("Adjusted P-value Threshold", value=sig_cfg.get('padj_threshold', 0.05), format="%.3f")
        log2fc_threshold = st.number_input("Log2 Fold Change Threshold", value=sig_cfg.get('log2fc_threshold', 1.0), format="%.2f")
        
        # PCA Settings
        st.markdown("---")
        st.write("📊 **PCA Settings**")
        pca_cfg = analysis_cfg.get('pca', {})
        pca_filter_hvg = st.checkbox("Filter High Variance Genes for PCA", value=pca_cfg.get('filter_high_variance', True), help="If enabled, PCA will be performed only on the most variable genes for better separation.")
        pca_top_n = st.number_input("Top N Genes for PCA", value=pca_cfg.get('top_n_genes', 2000), min_value=100, step=100)

        # Enrichment
        enrich_cfg = analysis_cfg.get('enrichment', {})
        species = st.selectbox("Species", ["Mouse", "Human", "Yeast", "Fly"], index=["Mouse", "Human", "Yeast", "Fly"].index(enrich_cfg.get('species', 'Mouse')))
        
        # Define common libraries per species
        common_libraries = {
            "Mouse": [
                "GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
                "KEGG_2019_Mouse", "WikiPathways_2019_Mouse", "Reactome_2022", "BioPlanet_2019"
            ],
            "Human": [
                "GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
                "KEGG_2021_Human", "WikiPathways_2019_Human", "Reactome_2022", "BioPlanet_2019"
            ],
            "Yeast": ["GO_Biological_Process_2018", "KEGG_2019_Yeast"],
            "Fly": ["GO_Biological_Process_2018", "KEGG_2019_Fly"]
        }
        
        # Get current libraries from config
        current_libs_dict = enrich_cfg.get('libraries', {})
        # If it's a dict, get values (the actual library names)
        current_lib_names = list(current_libs_dict.values()) if isinstance(current_libs_dict, dict) else []
        
        # Available options: Defaults for species + any currently configured ones (in case they are custom)
        available_options = sorted(list(set(common_libraries.get(species, []) + current_lib_names)))
        
        selected_libraries = st.multiselect(
            "Enrichr Libraries",
            options=available_options,
            default=[lib for lib in current_lib_names if lib in available_options],
            help="Gene sets source: https://maayanlab.cloud/Enrichr/#libraries"
        )
        
        custom_libraries_input = st.text_area(
            "Custom Enrichr Libraries (Optional)",
            help="Enter additional library names (one per line or comma-separated). e.g., TRRUST_Transcription_Factors_2019",
            height=100
        )
        
        # Time Series
        ts_cfg = analysis_cfg.get('time_series', {})
        time_col_present = ts_cfg.get('time_column_present', True)
        
        if not time_col_present:
            ts_enabled = st.checkbox("Enable Time-Series Analysis", value=False, disabled=True, help="Disabled because 'Include Time Column' was unchecked in Metadata Creation.")
            st.caption("⚠️ Time column disabled in Metadata.")
        else:
            ts_enabled = st.checkbox("Enable Time-Series Analysis", value=ts_cfg.get('enabled', True))
        
        time_col_name = st.text_input(
            "Time Column Name", 
            value=ts_cfg.get('time_column', 'Time'),
            help="The column in your metadata that represents time points (e.g., 'Time', 'Day', 'Stage')."
        )
        
        num_clusters = st.number_input("Number of Clusters (k)", value=ts_cfg.get('num_clusters', 4), min_value=2)
        ts_top_n_genes = st.number_input(
            "Top N Variable Genes for Clustering", 
            value=ts_cfg.get('top_n_genes', 2000), 
            min_value=100, 
            step=100,
            help="Recommended range: 1000 - 2000. \n\nSelects the top N genes with the highest variance across time points. \n- < 500: Very strict, only most drastic changes. \n- > 3000: May include too much noise.",
            key="ts_top_n"
        )

    # --- Construct Config Object Dynamically ---
    new_config = config.copy()
    # new_config['files']['feature_counts'] = feature_counts_path # Deprecated
    new_config['files']['gene_annotation'] = gene_annotation_path
    new_config['files']['metadata'] = config['files']['metadata']
    new_config['files']['output_dir'] = output_dir_path
    new_config['files']['batch_correction'] = {'enabled': batch_enabled, 'batch_key': batch_key}
    
    new_config['analysis']['diff_expression']['control_group'] = control_group
    new_config['analysis']['diff_expression']['filter_low_expression'] = filter_low
    new_config['analysis']['significance']['padj_threshold'] = padj_threshold
    new_config['analysis']['significance']['log2fc_threshold'] = log2fc_threshold
    new_config['analysis']['enrichment']['species'] = species
    
    # Convert selected list back to dict format for config
    # Combine selected and custom libraries
    custom_libs_list = [lib.strip() for lib in custom_libraries_input.replace('\n', ',').split(',') if lib.strip()]
    all_selected_libs = list(set(selected_libraries + custom_libs_list))

    new_libs = {}
    for lib in all_selected_libs:
        # Create a short key if possible, or just use the name
        if "GO_" in lib:
            key = "GO_" + lib.split("GO_")[1].split("_")[0] # e.g. GO_Biological
        elif "KEGG" in lib:
            key = "KEGG"
        elif "Reactome" in lib:
            key = "Reactome"
        elif "WikiPathways" in lib:
            key = "WikiPathways"
        else:
            key = lib
        
        # Handle duplicates (e.g. multiple GO libs)
        if key in new_libs:
            key = lib # Fallback to full name if collision
        
        new_libs[key] = lib
    
    new_config['analysis']['enrichment']['libraries'] = new_libs
    
    new_config['analysis']['time_series']['enabled'] = ts_enabled
    new_config['analysis']['time_series']['num_clusters'] = num_clusters
    new_config['analysis']['time_series']['top_n_genes'] = ts_top_n_genes
    new_config['analysis']['time_series']['time_column'] = time_col_name

    # Check/Create PCA section if missing
    if 'pca' not in new_config['analysis']: new_config['analysis']['pca'] = {}
    new_config['analysis']['pca']['filter_high_variance'] = pca_filter_hvg
    new_config['analysis']['pca']['top_n_genes'] = pca_top_n

    # Save Button (Manual)
    if st.button("💾 Save Configuration"):
        save_config(new_config)
        st.success("Configuration saved successfully!")
        # Reload config to ensure consistency
        config = load_config()

    st.markdown("---")
    st.header("🚀 Execute Pipeline")
    
    if st.button("Run Analysis Pipeline", type="primary"):
        # Auto-save before running
        save_config(new_config)
        st.info("Configuration auto-saved. Starting pipeline... This may take a while.")
        
        # Create a placeholder for logs
        log_container = st.empty()
        full_logs = ""
        
        process = run_pipeline()
        
        if process:
            # Stream output
            while True:
                line = process.stdout.readline()
                if not line and process.poll() is not None:
                    break
                if line:
                    full_logs += line
                    # Update the log container (showing last 20 lines for readability)
                    log_container.code("".join(full_logs.splitlines(keepends=True)[-20:]))
            
            if process.returncode == 0:
                st.success("✅ Pipeline completed successfully!")
                st.balloons()
                st.cache_data.clear() # Clear cache to force reload of new data
            else:
                st.error("❌ Pipeline failed. Check the logs above.")
                st.code(full_logs) # Show full logs on error

# ==========================================
# TAB 2: Interactive Visualization
# ==========================================
with tab2:
    st.header("Interactive Results Viewer")
    
    # --- Helper to load data ---
    OUTPUT_DIR = config.get('files', {}).get('output_dir', get_session_path('output'))
    
    # Check for batch corrected data first
    BATCH_CORRECTED_FILE = os.path.join(OUTPUT_DIR, "counts_batch_corrected.csv")
    RAW_COUNTS_FILE = os.path.join(OUTPUT_DIR, "counts.csv")
    METADATA_FILE = config.get('files', {}).get('metadata', get_session_path('input/metadata.csv'))
    
    if os.path.exists(BATCH_CORRECTED_FILE):
        COUNTS_FILE = BATCH_CORRECTED_FILE
        # RAW_COUNTS_FILE is already defined above
        st.info(f"Using Batch Corrected Data for Visualizations: `{os.path.basename(COUNTS_FILE)}`")
    else:
        COUNTS_FILE = RAW_COUNTS_FILE
        # RAW_COUNTS_FILE = RAW_COUNTS_FILE

    @st.cache_data
    def load_data(counts_path, metadata_path):
        if not os.path.exists(counts_path) or not os.path.exists(metadata_path):
            return None, None
        counts = pd.read_csv(counts_path, index_col=0)
        metadata = pd.read_csv(metadata_path, index_col=0)
        return counts, metadata
    
    # Load RAW counts for DE analysis (DESeq2 requires raw integer counts)
    @st.cache_data
    def load_raw_counts(raw_path):
        if os.path.exists(raw_path):
            return pd.read_csv(raw_path, index_col=0)
        return None

    @st.cache_resource
    def get_dds_object(_counts):
        # Ensure counts are integers (common issue if loaded from non-raw)
        try:
             # Round to nearest integer if not already
            _counts = _counts.round().astype(int)
        except:
            pass # If fails, let omicverse handle or fail
        dds = ov.bulk.pyDEG(_counts)
        dds.drop_duplicates_index()
        return dds

    counts, metadata = load_data(COUNTS_FILE, METADATA_FILE)
    raw_counts_df = load_raw_counts(RAW_COUNTS_FILE)

    if counts is None or metadata is None:
        st.warning(f"Data not found. Please run the pipeline first or check paths: \n- {COUNTS_FILE}\n- {METADATA_FILE}")
    else:
        # --- Interactive Visualization (Enhanced) ---
        st.subheader("1. Gene Expression Visualization")
        
        # Use RAW counts for DE Object
        if raw_counts_df is not None:
             dds = get_dds_object(raw_counts_df)
        else:
             st.warning("Raw counts file not found. Using displayed counts for DE (may be incorrect if normalized).")
             dds = get_dds_object(counts)
        all_genes = counts.index.tolist()
        all_conditions = metadata['Condition'].unique().tolist()
        
        col1, col2 = st.columns([1, 3])
        
        with col1:
            # viz_type = st.radio("Visualization Type", ["Boxplot", "Heatmap"], horizontal=True)
            # st.markdown("---")
            st.markdown("### Visualization Settings")
            
            selected_genes = st.multiselect("Select Gene(s)", options=all_genes, default=all_genes[:2] if len(all_genes)>1 else all_genes[:1])
            
            control_group_plot = st.selectbox("Control Group", all_conditions, index=0, key="viz_ctrl")
            treat_options = [c for c in all_conditions if c != control_group_plot]
            treat_groups_plot = st.multiselect("Treatment Group(s)", treat_options, default=treat_options, key="viz_treat")
            
            st.markdown("---")
            st.write("📊 **Plot Settings**")
            
            plot_w = st.slider("Width", 400, 1200, 800, step=50)
            plot_h = st.slider("Height", 300, 1000, 500, step=50)

            btn_plot = st.button("Generate Plot")

        with col2:
            if btn_plot and selected_genes and treat_groups_plot:
                # Prepare Data
                selected_conditions = [control_group_plot] + treat_groups_plot
                target_samples = metadata[metadata['Condition'].isin(selected_conditions)].index.tolist()
                
                # --- BOXPLOT LOGIC ---
                expr_data = counts.loc[selected_genes, target_samples].T
                
                # Check if data needs log transformation for visualization
                y_label = "Normalized Counts"
                if expr_data.max().max() > 100:
                     st.info("Data values > 100 detected. Applying Log2(Expression + 1) for better visualization.")
                     expr_data = np.log2(expr_data + 1)
                     y_label = "Log2(Normalized Counts + 1)"

                plot_df = expr_data.merge(metadata[['Condition']], left_index=True, right_index=True)
                df_long = plot_df.melt(id_vars='Condition', value_vars=selected_genes, var_name='Gene', value_name='Expression')
                df_long['Condition'] = pd.Categorical(df_long['Condition'], categories=selected_conditions, ordered=True)
                df_long = df_long.sort_values('Condition')
                
                try:
                    import seaborn as sns
                    import matplotlib.pyplot as plt
                    
                    # Convert pixels to inches for matplotlib (approx)
                    fig, ax = plt.subplots(figsize=(plot_w/100, plot_h/100))
                    
                    sns.boxplot(data=df_long, x='Gene', y='Expression', hue='Condition', ax=ax, palette='Set2', showfliers=False)
                    sns.stripplot(data=df_long, x='Gene', y='Expression', hue='Condition', ax=ax, dodge=True, color='black', alpha=0.6, size=4)
                    
                    handles, labels = ax.get_legend_handles_labels()
                    n_cond = len(selected_conditions)
                    ax.legend(handles[:n_cond], labels[:n_cond], title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left')
                    
                    ax.set_title("Gene Expression Levels")
                    ax.set_ylabel(y_label)
                    plt.tight_layout()
                    
                    st.pyplot(fig)
                    
                    # Save plot to buffer for download
                    import io
                    buf = io.BytesIO()
                    fig.savefig(buf, format="png", dpi=300, bbox_inches='tight')
                    buf.seek(0)
                    
                    st.download_button(
                        label="Download Plot (PNG)",
                        data=buf,
                        file_name="Gene_Expression_Boxplot.png",
                        mime="image/png"
                    )
                    
                except Exception as e:
                    st.error(f"Plotting error: {e}")

                # --- Statistics Table (Shared) ---
                st.markdown("---")
                st.subheader("Differential Expression Statistics")
                
                with st.spinner("Calculating statistics..."):
                    all_stats = []
                    control_samples = metadata[metadata['Condition'] == control_group_plot].index.tolist()

                    for treat_grp in treat_groups_plot:
                        treat_samples = metadata[metadata['Condition'] == treat_grp].index.tolist()
                        
                        # Run DEG analysis
                        deg_res = dds.deg_analysis(treat_samples, control_samples, method='DEseq2')
                        
                        if not deg_res.empty:
                            # Filter for selected genes
                            sub_res = deg_res.loc[deg_res.index.intersection(selected_genes)].copy()
                            sub_res['Comparison'] = f"{treat_grp} vs {control_group_plot}"
                            sub_res['Gene'] = sub_res.index
                            all_stats.append(sub_res)
                    
                    if all_stats:
                        final_stats = pd.concat(all_stats).reset_index(drop=True)
                        # Reorder columns
                        cols = ['Gene', 'Comparison'] + [c for c in final_stats.columns if c not in ['Gene', 'Comparison']]
                        st.dataframe(final_stats[cols], use_container_width=True)
                    else:
                        st.warning("No significant results found for selected genes.")

            elif btn_plot:
                st.info("Please select genes and at least one treatment group.")



# ==========================================
# TAB: Volcano Plot
# ==========================================
with tab_volcano:
    st.header("🌋 Volcano Plot Viewer")
    
    OUTPUT_DIR = config.get('files', {}).get('output_dir', get_session_path('output'))
    
    if os.path.exists(OUTPUT_DIR):
        # Find all Volcano Plot files in the output directory
        # Pattern: Volcano_{Comparison}.png
        volcano_files = [f for f in os.listdir(OUTPUT_DIR) if f.startswith("Volcano_") and f.endswith(".png")]
        
        if volcano_files:
            # Extract comparison names
            comparisons = set()
            file_map = {} # Map comparison -> png path
            
            for f in volcano_files:
                # Remove extension
                base = os.path.splitext(f)[0]
                # Remove prefix
                comp_name = base.replace("Volcano_", "", 1)
                comparisons.add(comp_name)
                file_map[comp_name] = os.path.join(OUTPUT_DIR, f)
            
            sorted_comparisons = sorted(list(comparisons))
            
            col_sel, col_slider = st.columns([3, 1])
            with col_sel:
                selected_comp_volcano = st.selectbox("Select Comparison", options=sorted_comparisons, key="volcano_comp_select")
            with col_slider:
                volcano_width = st.slider("Image Width", min_value=300, max_value=2000, value=800, step=50, key="volcano_width")
            
            if selected_comp_volcano:
                png_path = file_map[selected_comp_volcano]
                st.image(png_path, caption=os.path.basename(png_path), width=volcano_width)
                
                with open(png_path, "rb") as f:
                    st.download_button(
                        label="Download Plot (PNG)",
                        data=f,
                        file_name=os.path.basename(png_path),
                        mime="image/png"
                    )
                
                # --- Data Table ---
                st.markdown("---")
                st.subheader("📋 Gene Data Table")
                
                # Construct CSV path
                csv_filename = f"DEG_{selected_comp_volcano}.csv"
                csv_path = os.path.join(OUTPUT_DIR, csv_filename)
                
                if os.path.exists(csv_path):
                    try:
                        df = pd.read_csv(csv_path, index_col=0)
                        
                        # Add Status column if missing
                        if 'Status' not in df.columns and 'log2FoldChange' in df.columns and 'padj' in df.columns:
                            log2fc_thresh = config['analysis']['significance']['log2fc_threshold']
                            padj_thresh = config['analysis']['significance']['padj_threshold']
                            
                            conditions = [
                                (df['padj'] < padj_thresh) & (df['log2FoldChange'] > log2fc_thresh),
                                (df['padj'] < padj_thresh) & (df['log2FoldChange'] < -log2fc_thresh)
                            ]
                            choices = ['Up', 'Down']
                            df['Status'] = np.select(conditions, choices, default='Not Sig')
                        
                        # Filter options
                        filter_opt = st.radio("Show:", ["All Genes", "Significant Only (Up/Down)"], horizontal=True, key="volcano_table_filter")
                        
                        if filter_opt == "Significant Only (Up/Down)":
                            if 'Status' in df.columns:
                                table_df = df[df['Status'] != 'Not Sig']
                            else:
                                st.warning("Status column not found, showing all genes.")
                                table_df = df
                        else:
                            table_df = df
                        
                        # Display with formatting
                        format_dict = {
                            'padj': '{:.2e}',
                            'pvalue': '{:.2e}',
                            'log2FoldChange': '{:.2f}',
                            'baseMean': '{:.2f}',
                            'lfcSE': '{:.3f}'
                        }
                        # Only format columns that exist
                        final_format = {k: v for k, v in format_dict.items() if k in table_df.columns}
                        
                        st.dataframe(table_df.style.format(final_format), use_container_width=True)
                        
                        # Download button
                        csv = table_df.to_csv().encode('utf-8')
                        st.download_button(
                            label="Download Data (CSV)",
                            data=csv,
                            file_name=f"Volcano_Data_{selected_comp_volcano}.csv",
                            mime="text/csv"
                        )
                        
                    except Exception as e:
                        st.error(f"Error loading data table: {e}")
                else:
                    st.info(f"Data file `{csv_filename}` not found.")
        else:
            st.warning("No Volcano Plots (PNG) found in output directory.")
    else:
        st.warning("Analysis output not found. Please run the pipeline first.")

# ==========================================
# TAB 3: Static Plots Gallery
# ==========================================
with tab3:
    st.header("📈 Static Plots Gallery")
    
    OUTPUT_DIR = config.get('files', {}).get('output_dir', get_session_path('output'))
    
    # Find all PNGs in output directory
    if os.path.exists(OUTPUT_DIR):
        png_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('.png') and not f.startswith('Volcano_')]
        
        if png_files:
            col_sel, col_slider = st.columns([3, 1])
            with col_sel:
                selected_plot = st.selectbox("Select a plot to view", options=png_files)
            with col_slider:
                img_width = st.slider("Image Width", min_value=300, max_value=2000, value=800, step=50, key="static_width")
            
            if selected_plot:
                plot_path = os.path.join(OUTPUT_DIR, selected_plot)
                st.image(plot_path, caption=selected_plot, width=img_width)
                
                with open(plot_path, "rb") as f:
                    st.download_button(
                        label="Download Plot (PNG)",
                        data=f,
                        file_name=selected_plot,
                        mime="image/png"
                    )
        else:
            st.info("No plot images found in output directory.")
    else:
        st.warning("Output directory does not exist.")

# ==========================================
# TAB 4: Pathway Analysis (Enrichment)
# ==========================================
with tab4:
    st.header("🧬 Pathway Enrichment Analysis")
    
    ENRICHMENT_DIR = os.path.join(config.get('files', {}).get('output_dir', get_session_path('output')), "enrichment")
    
    if os.path.exists(ENRICHMENT_DIR):
        # 1. Select Comparison
        comparisons = [d for d in os.listdir(ENRICHMENT_DIR) if os.path.isdir(os.path.join(ENRICHMENT_DIR, d))]
        
        if comparisons:
            selected_comparison = st.selectbox("Select Comparison Group", options=comparisons)
            
            if selected_comparison:
                comp_dir = os.path.join(ENRICHMENT_DIR, selected_comparison)
                
                # 2. List Dotplots
                # Filter for Dotplot pngs
                plot_files = [f for f in os.listdir(comp_dir) if f.startswith("Dotplot_") and f.endswith(".png")]
                
                if plot_files:
                    # Create a mapping for readable names
                    plot_map = {}
                    for f in plot_files:
                        name = f.replace("Dotplot_", "").replace(".png", "")
                        plot_map[name] = f
                    
                    sorted_plots = sorted(list(plot_map.keys()))
                    
                    col_p_sel, col_p_slider = st.columns([3, 1])
                    with col_p_sel:
                        selected_plot_name = st.selectbox("Select Pathway Result", options=sorted_plots)
                    with col_p_slider:
                        enrich_width = st.slider("Image Width", min_value=400, max_value=2000, value=800, step=50, key="enrich_width")
                    
                    if selected_plot_name:
                        plot_file = plot_map[selected_plot_name]
                        plot_path = os.path.join(comp_dir, plot_file)
                        st.image(plot_path, width=enrich_width)
                        
                        with open(plot_path, "rb") as f:
                            st.download_button(
                                label="Download Plot (PNG)",
                                data=f,
                                file_name=plot_file,
                                mime="image/png"
                            )
                        
                        # Optional: Check for corresponding CSV
                        # Pattern: {name}.csv (without Dotplot_)
                        # The name variable above is e.g. GO_Treat_24h_vs_Control_Down_GO
                        # The file is likely GO_Treat_24h_vs_Control_Down_GO.csv in the same dir
                        csv_file = f"{selected_plot_name}.csv"
                        csv_path = os.path.join(comp_dir, csv_file)
                        
                        if os.path.exists(csv_path):
                            with open(csv_path, "rb") as f:
                                st.download_button(
                                    label="Download Result Table (CSV)",
                                    data=f,
                                    file_name=csv_file,
                                    mime="text/csv"
                                )
                else:
                    st.warning(f"No dotplots found in {comp_dir}. This might mean no significant pathways were found.")
        else:
            st.warning("No comparison results found in the enrichment directory.")
    else:
        st.warning(f"Enrichment directory not found: {ENRICHMENT_DIR}. Please run the pipeline first.")

