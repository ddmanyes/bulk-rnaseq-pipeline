import pandas as pd
import os
import argparse
import glob

def read_count_file(filepath):
    """
    Reads a count file and returns a DataFrame.
    Supports .csv, .xlsx, .txt (featureCounts).
    """
    print(f"Reading '{filepath}'...")
    filename = os.path.basename(filepath)
    
    if filepath.endswith('.csv'):
        df = pd.read_csv(filepath, index_col=0)
    elif filepath.endswith('.xlsx'):
        df = pd.read_excel(filepath, index_col=0)
    else:
        # Assume tab-separated or featureCounts
        try:
            df = pd.read_csv(filepath, sep='\t', comment='#', index_col=0)
        except:
            print(f"  ❌ Failed to read '{filename}' as tab-separated. Trying comma...")
            df = pd.read_csv(filepath, index_col=0)

    # Handle featureCounts format (remove metadata columns)
    if 'Chr' in df.columns and 'Start' in df.columns:
        print("  - Detected featureCounts format. Cleaning columns...")
        df = df.iloc[:, 5:]
        # Clean column names (remove path and extension)
        df.columns = [c.split('/')[-1].replace('.sorted.bam', '').replace('.bam', '') for c in df.columns]
    
    return df

def merge_files(input_paths, output_file):
    """
    Merges a list of file paths into a single CSV.
    """
    merged_df = pd.DataFrame()
    
    for filepath in input_paths:
        if not os.path.exists(filepath):
            print(f"⚠️ File not found: {filepath}")
            continue
            
        df = read_count_file(filepath)
        
        if merged_df.empty:
            merged_df = df
        else:
            # Check for duplicate columns
            common_cols = set(merged_df.columns) & set(df.columns)
            if common_cols:
                print(f"  ⚠️ Warning: Duplicate sample names found: {common_cols}")
                print("  - Renaming duplicates in new file with suffix '_new'")
                df.rename(columns={c: f"{c}_new" for c in common_cols}, inplace=True)
            
            # Merge on index (GeneID)
            print(f"  - Merging {len(df.columns)} samples...")
            merged_df = merged_df.merge(df, left_index=True, right_index=True, how='outer')
    
    # Fill NaN with 0 (assuming missing means 0 counts)
    print("Filling missing values with 0...")
    merged_df.fillna(0, inplace=True)
    
    # Save
    print(f"Saving merged counts to '{output_file}'...")
    merged_df.to_csv(output_file)
    print("✅ Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge multiple count matrices into one.")
    parser.add_argument('inputs', nargs='+', help="Input files (e.g., batch1.csv batch2.xlsx) or a glob pattern")
    parser.add_argument('-o', '--output', default='merged_counts.csv', help="Output filename (default: merged_counts.csv)")
    
    args = parser.parse_args()
    
    # Expand globs if shell didn't do it (e.g. on Windows sometimes)
    all_files = []
    for inp in args.inputs:
        expanded = glob.glob(inp)
        if expanded:
            all_files.extend(expanded)
        else:
            all_files.append(inp) # It might be a specific file that doesn't exist yet but let read_count_file handle it
            
    if not all_files:
        print("❌ No input files found.")
    else:
        merge_files(all_files, args.output)
