import streamlit as st
import pandas as pd
import omicverse as ov
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import yaml
import subprocess
import time

# --- Configuration & Helper Functions ---
CONFIG_FILE = "config.yml"

def load_config():
    """Loads the YAML config file."""
    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, 'r') as f:
            return yaml.safe_load(f)
    return {}

def save_config(config_data):
    """Saves the dictionary to the YAML config file."""
    with open(CONFIG_FILE, 'w') as f:
        yaml.dump(config_data, f, default_flow_style=False, sort_keys=False)

def run_pipeline():
    """Executes the analysis pipeline script."""
    try:
        # Run the script and capture output
        process = subprocess.Popen(
            ["python", "analysis_pipeline.py"],
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
st.title("🧬 Bulk RNA-seq Analysis Pipeline")

# Load initial config
config = load_config()

# Create Tabs
tab_merge, tab0, tab1, tab2, tab3, tab4 = st.tabs(["📂 Merge Counts", "📝 Metadata Creation", "⚙️ Run Analysis", "📊 Interactive Visualization", "📈 Static Plots", "🧬 Pathway Analysis"])

# ==========================================
# TAB MERGE: Merge Count Files
# ==========================================
with tab_merge:
    st.header("Merge Count Matrices")
    st.markdown("Combine multiple count files (CSV, Excel, featureCounts) into a single matrix for analysis.")
    
    col1, col2 = st.columns([3, 1])
    with col1:
        input_folder = st.text_input("Input Folder Path", value="input", help="Folder containing your count files")
    with col2:
        if st.button("Scan Folder"):
            st.rerun()

    if os.path.exists(input_folder):
        # Find potential count files
        all_files = [f for f in os.listdir(input_folder) if f.endswith(('.csv', '.xlsx', '.txt', '.tsv')) and not f.startswith('merged_')]
        
        if all_files:
            st.success(f"Found {len(all_files)} potential count files.")
            
            selected_files = st.multiselect("Select files to merge", options=all_files, default=all_files)
            
            output_filename = st.text_input("Output Filename", value="merged_counts.csv")
            
            if st.button("Merge Selected Files", type="primary"):
                if not selected_files:
                    st.warning("Please select at least one file.")
                else:
                    with st.spinner("Merging files..."):
                        try:
                            # Import the merge logic dynamically or reimplement it here for simplicity
                            # Reimplementing core logic to avoid external script dependency issues in Streamlit
                            merged_df = pd.DataFrame()
                            
                            for filename in selected_files:
                                filepath = os.path.join(input_folder, filename)
                                st.text(f"Processing {filename}...")
                                
                                # Read file
                                if filename.endswith('.csv'):
                                    df = pd.read_csv(filepath, index_col=0)
                                elif filename.endswith('.xlsx'):
                                    df = pd.read_excel(filepath, index_col=0)
                                else:
                                    # Try tab then comma
                                    try:
                                        df = pd.read_csv(filepath, sep='\t', comment='#', index_col=0)
                                    except:
                                        df = pd.read_csv(filepath, index_col=0)

                                # Clean featureCounts
                                if 'Chr' in df.columns and 'Start' in df.columns:
                                    df = df.iloc[:, 5:]
                                    df.columns = [c.split('/')[-1].replace('.sorted.bam', '').replace('.bam', '') for c in df.columns]
                                
                                # Merge
                                if merged_df.empty:
                                    merged_df = df
                                else:
                                    # Handle duplicates
                                    common_cols = set(merged_df.columns) & set(df.columns)
                                    if common_cols:
                                        df.rename(columns={c: f"{c}_{filename}" for c in common_cols}, inplace=True)
                                    
                                    merged_df = merged_df.merge(df, left_index=True, right_index=True, how='outer')
                            
                            # Fill NaNs
                            merged_df.fillna(0, inplace=True)
                            
                            # Save
                            output_path = os.path.join(input_folder, output_filename)
                            merged_df.to_csv(output_path)
                            
                            st.success(f"✅ Successfully merged {len(selected_files)} files into `{output_path}`")
                            st.dataframe(merged_df.head())
                            
                        except Exception as e:
                            st.error(f"Error during merge: {e}")
        else:
            st.warning("No compatible files (.csv, .xlsx, .txt) found in directory.")
    else:
        st.error("Directory does not exist.")

# ==========================================
# TAB 0: Metadata Creation
# ==========================================
with tab0:
    st.header("Metadata Editor")
    st.markdown("Create or edit your sample metadata file here. This file is required for the analysis.")

    # Path to metadata file (from config or default)
    meta_file_path = config.get('files', {}).get('metadata', 'input/metadata.csv')
    
    # Load existing or create new
    if os.path.exists(meta_file_path):
        st.info(f"Loading existing metadata from: `{meta_file_path}`")
        try:
            existing_df = pd.read_csv(meta_file_path, index_col=0)
        except Exception as e:
            st.error(f"Error loading metadata: {e}")
            existing_df = pd.DataFrame(columns=['Condition', 'Time', 'Type'])
    else:
        st.warning(f"File `{meta_file_path}` not found. Starting with a template.")
        existing_df = pd.DataFrame({
            'Condition': ['Control', 'Treat'],
            'Time': [0, 6],
            'Type': ['Control', 'Treat']
        }, index=['Sample1', 'Sample2'])
        existing_df.index.name = 'SampleID'

    # Editable Dataframe
    # We reset index to make SampleID editable as a column, then set it back before saving
    df_to_edit = existing_df.reset_index()
    
    st.markdown("### Edit Metadata Table")
    st.markdown("Make sure to have unique **SampleID**s and required columns: **Condition**, **Time**, **Type**.")
    
    edited_df = st.data_editor(
        df_to_edit, 
        num_rows="dynamic", 
        use_container_width=True,
        column_config={
            "SampleID": st.column_config.TextColumn("Sample ID (Must match count matrix columns)", required=True),
            "Condition": st.column_config.TextColumn("Condition (e.g. Control, Treat)", required=True),
            "Time": st.column_config.NumberColumn("Time (e.g. 0, 6, 24)", required=True),
            "Type": st.column_config.SelectboxColumn("Type", options=["Control", "Treat"], required=True)
        }
    )

    if st.button("💾 Save Metadata"):
        try:
            # Check for empty SampleIDs
            if edited_df['SampleID'].isnull().any() or (edited_df['SampleID'] == '').any():
                st.error("Error: SampleID cannot be empty.")
            elif edited_df['SampleID'].duplicated().any():
                st.error("Error: Duplicate SampleIDs found.")
            else:
                # Set index back to SampleID
                final_df = edited_df.set_index('SampleID')
                
                # Ensure directory exists
                os.makedirs(os.path.dirname(meta_file_path), exist_ok=True)
                
                final_df.to_csv(meta_file_path)
                st.success(f"Metadata saved successfully to `{meta_file_path}`!")
                
                # Update config if needed (though we used the path from config)
                # If the user wants to change the path, they should do it in the Config tab.
                
                # Show preview
                st.dataframe(final_df.head())
                
        except Exception as e:
            st.error(f"Failed to save metadata: {e}")

# ==========================================
# TAB 1: Run Analysis (Configuration & Execution)
# ==========================================
with tab1:
    st.header("Pipeline Configuration")
    
    with st.form("config_form"):
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("📁 File Paths")
            # Files
            files_cfg = config.get('files', {})
            feature_counts_path = st.text_input("Feature Counts File", value=files_cfg.get('feature_counts', 'input/featureCounts.txt'))
            gene_annotation_path = st.text_input("Gene Annotation File", value=files_cfg.get('gene_annotation', 'genesets/pair_GRCm39.tsv'))
            # Metadata path input - sync with Tab 0 implicitly by reading config
            metadata_path = st.text_input("Metadata File", value=files_cfg.get('metadata', 'input/metadata.csv'), help="Path where metadata is saved (see 'Metadata Creation' tab)")
            output_dir_path = st.text_input("Output Directory", value=files_cfg.get('output_dir', 'output'))

            st.subheader("🧪 Batch Correction")
            # Batch Correction
            batch_cfg = files_cfg.get('batch_correction', {})
            batch_enabled = st.checkbox("Enable Batch Correction", value=batch_cfg.get('enabled', False))
            batch_key = st.text_input("Batch Key (Column in Metadata)", value=batch_cfg.get('batch_key', 'Batch'))

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
                default=[lib for lib in current_lib_names if lib in available_options]
            )
            
            # Time Series
            ts_cfg = analysis_cfg.get('time_series', {})
            ts_enabled = st.checkbox("Enable Time-Series Analysis", value=ts_cfg.get('enabled', True))
            num_clusters = st.number_input("Number of Clusters (k)", value=ts_cfg.get('num_clusters', 4), min_value=2)
            top_n_genes = st.number_input("Top N Variable Genes for Clustering", value=ts_cfg.get('top_n_genes', 2000), min_value=100, step=100)

        # Save Button
        submitted = st.form_submit_button("💾 Save Configuration")
        
        if submitted:
            # Update config dictionary
            new_config = config.copy()
            new_config['files']['feature_counts'] = feature_counts_path
            new_config['files']['gene_annotation'] = gene_annotation_path
            new_config['files']['metadata'] = metadata_path
            new_config['files']['output_dir'] = output_dir_path
            new_config['files']['batch_correction'] = {'enabled': batch_enabled, 'batch_key': batch_key}
            
            new_config['analysis']['diff_expression']['control_group'] = control_group
            new_config['analysis']['diff_expression']['filter_low_expression'] = filter_low
            new_config['analysis']['significance']['padj_threshold'] = padj_threshold
            new_config['analysis']['significance']['log2fc_threshold'] = log2fc_threshold
            new_config['analysis']['enrichment']['species'] = species
            
            # Convert selected list back to dict format for config
            new_libs = {}
            for lib in selected_libraries:
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
            new_config['analysis']['time_series']['top_n_genes'] = top_n_genes
            
            save_config(new_config)
            st.success("Configuration saved successfully!")
            # Reload config to ensure consistency
            config = load_config()

    st.markdown("---")
    st.header("🚀 Execute Pipeline")
    
    if st.button("Run Analysis Pipeline", type="primary"):
        st.info("Starting pipeline... This may take a while.")
        
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
            else:
                st.error("❌ Pipeline failed. Check the logs above.")
                st.code(full_logs) # Show full logs on error

# ==========================================
# TAB 2: Interactive Visualization
# ==========================================
with tab2:
    st.header("Interactive Results Viewer")
    
    # --- Helper to load data ---
    COUNTS_FILE = os.path.join(config.get('files', {}).get('output_dir', 'output'), "counts.csv")
    METADATA_FILE = config.get('files', {}).get('metadata', 'input/metadata.csv')
    OUTPUT_DIR = config.get('files', {}).get('output_dir', 'output')

    @st.cache_data
    def load_data(counts_path, metadata_path):
        if not os.path.exists(counts_path) or not os.path.exists(metadata_path):
            return None, None
        counts = pd.read_csv(counts_path, index_col=0)
        metadata = pd.read_csv(metadata_path, index_col=0)
        return counts, metadata

    @st.cache_resource
    def get_dds_object(_counts):
        dds = ov.bulk.pyDEG(_counts)
        dds.drop_duplicates_index()
        return dds

    counts, metadata = load_data(COUNTS_FILE, METADATA_FILE)

    if counts is None or metadata is None:
        st.warning(f"Data not found. Please run the pipeline first or check paths: \n- {COUNTS_FILE}\n- {METADATA_FILE}")
    else:
        # --- Interactive Boxplot (Original Logic) ---
        st.subheader("1. Gene Expression Boxplot")
        
        dds = get_dds_object(counts)
        all_genes = counts.index.tolist()
        all_conditions = metadata['Condition'].unique().tolist()
        
        col1, col2 = st.columns([1, 3])
        
        with col1:
            selected_genes = st.multiselect("Select Gene(s)", options=all_genes, default=all_genes[:2] if len(all_genes)>1 else all_genes[:1])
            
            control_group_plot = st.selectbox("Control Group", all_conditions, index=0, key="viz_ctrl")
            treat_options = [c for c in all_conditions if c != control_group_plot]
            treat_group_plot = st.selectbox("Treatment Group", treat_options, index=0 if treat_options else -1, key="viz_treat")
            
            btn_plot = st.button("Generate Boxplot")

        with col2:
            if btn_plot and selected_genes and treat_group_plot:
                control_samples = metadata[metadata['Condition'] == control_group_plot].index.tolist()
                treat_samples = metadata[metadata['Condition'] == treat_group_plot].index.tolist()
                
                try:
                    fig, ax = dds.plot_boxplot(
                        genes=selected_genes,
                        treatment_groups=treat_samples,
                        control_groups=control_samples,
                        figsize=(max(4, 2 * len(selected_genes)), 5),
                        fontsize=12,
                        legend_bbox=(1.05, 1)
                    )
                    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
                    st.pyplot(fig)
                    
                    # Display data table
                    st.markdown("---")
                    st.subheader("Differential Expression Statistics")
                    
                    with st.spinner("Calculating statistics..."):
                        # Run DEG analysis for the selected pair to show the table
                        # Note: This might take a moment depending on dataset size
                        deg_res = dds.deg_analysis(treat_samples, control_samples, method='DEseq2')
                        
                        # Filter for selected genes and display
                        if not deg_res.empty:
                            # Add gene symbol as a column if it's the index
                            display_df = deg_res.loc[selected_genes].copy()
                            st.dataframe(display_df, use_container_width=True)
                        else:
                            st.warning("No result from DEG analysis.")

                except Exception as e:
                    st.error(f"Plotting error: {e}")
            elif btn_plot:
                st.info("Please select genes and groups.")

# ==========================================
# TAB 3: Static Plots Gallery
# ==========================================
with tab3:
    st.header("📈 Static Plots Gallery")
    
    OUTPUT_DIR = config.get('files', {}).get('output_dir', 'output')
    
    # Find all PNGs in output directory
    if os.path.exists(OUTPUT_DIR):
        png_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('.png')]
        
        if png_files:
            selected_plot = st.selectbox("Select a plot to view", options=png_files)
            if selected_plot:
                st.image(os.path.join(OUTPUT_DIR, selected_plot), caption=selected_plot, use_container_width=True)
        else:
            st.info("No plot images found in output directory.")
    else:
        st.warning("Output directory does not exist.")

# ==========================================
# TAB 4: Pathway Analysis (Enrichment)
# ==========================================
with tab4:
    st.header("🧬 Pathway Enrichment Analysis")
    
    ENRICHMENT_DIR = os.path.join(config.get('files', {}).get('output_dir', 'output'), "enrichment")
    
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
                    st.info(f"Found {len(plot_files)} enrichment plots for {selected_comparison}.")
                    
                    # Display in a grid or list
                    for plot_file in plot_files:
                        st.markdown(f"### {plot_file.replace('Dotplot_', '').replace('.png', '')}")
                        st.image(os.path.join(comp_dir, plot_file), use_container_width=True)
                        st.markdown("---")
                else:
                    st.warning(f"No dotplots found in {comp_dir}. This might mean no significant pathways were found.")
        else:
            st.warning("No comparison results found in the enrichment directory.")
    else:
        st.warning(f"Enrichment directory not found: {ENRICHMENT_DIR}. Please run the pipeline first.")

