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

# Set Pandas Styler limit to avoid errors with large tables
pd.set_option("styler.render.max_elements", 2000000)

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
tab_merge, tab0, tab1, tab2, tab_volcano, tab3, tab4 = st.tabs(["📂 Merge Counts", "📝 Metadata Creation", "⚙️ Run Analysis", "📊 Interactive Visualization", "🌋 Volcano Plot", "📈 Static Plots", "🧬 Pathway Analysis"])

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
    
    # --- Configuration Section (No Form) ---
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
        ts_enabled = st.checkbox("Enable Time-Series Analysis", value=ts_cfg.get('enabled', True))
        num_clusters = st.number_input("Number of Clusters (k)", value=ts_cfg.get('num_clusters', 4), min_value=2)
        top_n_genes = st.number_input(
            "Top N Variable Genes for Clustering", 
            value=ts_cfg.get('top_n_genes', 2000), 
            min_value=100, 
            step=100,
            help="Recommended range: 1000 - 2000. \n\nSelects the top N genes with the highest variance across time points. \n- < 500: Very strict, only most drastic changes. \n- > 3000: May include too much noise."
        )

    # --- Construct Config Object Dynamically ---
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
    new_config['analysis']['time_series']['top_n_genes'] = top_n_genes

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
        # --- Interactive Visualization (Enhanced) ---
        st.subheader("1. Gene Expression Visualization")
        
        dds = get_dds_object(counts)
        all_genes = counts.index.tolist()
        all_conditions = metadata['Condition'].unique().tolist()
        
        col1, col2 = st.columns([1, 3])
        
        with col1:
            viz_type = st.radio("Visualization Type", ["Boxplot", "Heatmap"], horizontal=True)
            st.markdown("---")
            
            selected_genes = st.multiselect("Select Gene(s)", options=all_genes, default=all_genes[:2] if len(all_genes)>1 else all_genes[:1])
            
            control_group_plot = st.selectbox("Control Group", all_conditions, index=0, key="viz_ctrl")
            treat_options = [c for c in all_conditions if c != control_group_plot]
            treat_groups_plot = st.multiselect("Treatment Group(s)", treat_options, default=treat_options[:1] if treat_options else [], key="viz_treat")
            
            st.markdown("---")
            st.write("📊 **Plot Settings**")
            
            if viz_type == "Heatmap":
                use_zscore = st.checkbox("Apply Z-score (Row-wise)", value=True, help="Standardize expression per gene to highlight relative changes.")
                color_scale = st.selectbox("Color Scale", ["RdBu_r", "Viridis", "Magma", "Plasma"], index=0)
            
            plot_w = st.slider("Width", 400, 1200, 800, step=50)
            plot_h = st.slider("Height", 300, 1000, 500, step=50)

            btn_plot = st.button("Generate Plot")

        with col2:
            if btn_plot and selected_genes and treat_groups_plot:
                # Prepare Data
                selected_conditions = [control_group_plot] + treat_groups_plot
                target_samples = metadata[metadata['Condition'].isin(selected_conditions)].index.tolist()
                
                if viz_type == "Boxplot":
                    # --- BOXPLOT LOGIC ---
                    expr_data = counts.loc[selected_genes, target_samples].T
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
                        ax.set_ylabel("Normalized Counts")
                        plt.tight_layout()
                        
                        st.pyplot(fig)
                        
                    except Exception as e:
                        st.error(f"Plotting error: {e}")
                
                elif viz_type == "Heatmap":
                    # --- HEATMAP LOGIC ---
                    try:
                        import plotly.express as px
                        from scipy.stats import zscore
                        
                        # Get data: Genes x Samples
                        heatmap_data = counts.loc[selected_genes, target_samples]
                        
                        # Sort samples by condition
                        sample_conditions = metadata.loc[target_samples, 'Condition']
                        # Create a sorter based on selected_conditions order
                        condition_sorter = {c: i for i, c in enumerate(selected_conditions)}
                        sorted_samples = sorted(target_samples, key=lambda x: (condition_sorter.get(sample_conditions[x], 999), x))
                        
                        heatmap_data = heatmap_data[sorted_samples]
                        
                        if use_zscore:
                            # Apply Z-score row-wise (axis=1)
                            # zscore returns numpy array, need to rebuild DF
                            heatmap_data = pd.DataFrame(zscore(heatmap_data, axis=1), index=heatmap_data.index, columns=heatmap_data.columns)
                            title_suffix = "(Z-score)"
                            color_mid = 0
                        else:
                            title_suffix = "(Normalized Counts)"
                            color_mid = None
                        
                        # Create Heatmap
                        fig = px.imshow(
                            heatmap_data,
                            labels=dict(x="Sample", y="Gene", color="Expression"),
                            x=heatmap_data.columns,
                            y=heatmap_data.index,
                            color_continuous_scale=color_scale,
                            color_continuous_midpoint=color_mid,
                            aspect="auto",
                            title=f"Gene Expression Heatmap {title_suffix}"
                        )
                        
                        fig.update_layout(
                            width=plot_w,
                            height=plot_h,
                            xaxis_title="Samples",
                            yaxis_title="Genes"
                        )
                        
                        st.plotly_chart(fig)
                        
                    except Exception as e:
                        st.error(f"Heatmap error: {e}")

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
    
    OUTPUT_DIR = config.get('files', {}).get('output_dir', 'output')
    
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
    
    OUTPUT_DIR = config.get('files', {}).get('output_dir', 'output')
    
    # Find all PNGs in output directory
    if os.path.exists(OUTPUT_DIR):
        png_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('.png')]
        
        if png_files:
            col_sel, col_slider = st.columns([3, 1])
            with col_sel:
                selected_plot = st.selectbox("Select a plot to view", options=png_files)
            with col_slider:
                img_width = st.slider("Image Width", min_value=300, max_value=2000, value=800, step=50, key="static_width")
            
            if selected_plot:
                st.image(os.path.join(OUTPUT_DIR, selected_plot), caption=selected_plot, width=img_width)
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
                        st.image(os.path.join(comp_dir, plot_file), width=enrich_width)
                        
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

