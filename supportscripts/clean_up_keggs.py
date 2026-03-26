import pandas as pd
import sys

def process_humann_matrix(input_file, output_file):
    # Load the HUMAnN output
    df = pd.read_csv(input_file, sep='\t', index_offset=0)
    
    # Rename the first column for clarity
    df.rename(columns={df.columns[0]: 'KO_Description'}, inplace=True)
    
    # 1. Filter for UNSTRATIFIED results
    # HUMAnN marks the total abundance of a KO without a pipe "|" symbol.
    # Stratified lines look like "KO | s__Species"
    df_unstratified = df[~df['KO_Description'].str.contains('\|')].copy()
    
    # 2. Split KO ID and Description
    # Format is usually "K00001: Description text"
    # We want to keep both or separate them.
    df_unstratified[['KO_ID', 'Description']] = df_unstratified['KO_Description'].str.split(': ', n=1, expand=True)
    
    # 3. Clean up the sample names 
    # HUMAnN often appends "_Abundance-RPKs" or similar to column headers
    new_cols = {col: col.replace('_Abundance-CPM', '') for col in df_unstratified.columns}
    df_unstratified.rename(columns=new_cols, inplace=True)
    
    # Reorder columns to have ID and Description first
    cols = ['KO_ID', 'Description'] + [c for c in df_unstratified.columns if c not in ['KO_ID', 'Description', 'KO_Description']]
    df_final = df_unstratified[cols]
    
    # Save the cleaned matrix
    df_final.to_csv(output_file, sep='\t', index=False)
    print(f"Cleaned matrix saved to {output_file}")

if __name__ == "__main__":
    process_humann_matrix("ko_cpm_named.tsv", "KO_Sample_Matrix_Final.tsv")