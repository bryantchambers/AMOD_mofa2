#!/bin/bash

# --- Configuration ---
INPUT_DIR="./fastq_input"
OUTPUT_DIR="./humann_results"
THREADS=40
# Path to your humann utility mapping files (usually in your conda env)
# You need the uniref90_to_ko mapping file.
MAPPING_FILE="uniref90_to_ko.tsv.gz" 

mkdir -p $OUTPUT_DIR/main_outputs
mkdir -p $OUTPUT_DIR/ko_tables

# --- Step 1: Run HUMAnN on individual samples ---
for f in $INPUT_DIR/*.fastq.gz; do
    sample_name=$(basename "$f" .fastq.gz)
    echo "Processing sample: $sample_name"
    
    humann --input "$f" \
           --output "$OUTPUT_DIR/main_outputs" \
           --threads $THREADS \
           --search-mode rapsearch \
           --taxonomic-profile "$OUTPUT_DIR/main_outputs/${sample_name}_humann_temp/${sample_name}_metaphlan_bugs_list.tsv"
done

# --- Step 2: Join individual gene family tables ---
# HUMAnN defaults to UniRef90 gene families.
humann_join_tables --input "$OUTPUT_DIR/main_outputs" \
                   --output "$OUTPUT_DIR/genefamilies_merged.tsv" \
                   --file_name "genefamilies"

# --- Step 3: Regroup UniRef90 to KEGG Orthologs (KO) ---
humann_regroup_table --input "$OUTPUT_DIR/genefamilies_merged.tsv" \
                     --output "$OUTPUT_DIR/ko_merged.tsv" \
                     --groups uniref90_to_ko

# --- Step 4: Normalize to CPM (Copies Per Million) ---
humann_renorm_table --input "$OUTPUT_DIR/ko_merged.tsv" \
                    --output "$OUTPUT_DIR/ko_cpm.tsv" \
                    --units cpm \
                    --update-snames

# --- Step 5: Attach Human-Readable Names ---
humann_rename_table --input "$OUTPUT_DIR/ko_cpm.tsv" \
                    --output "$OUTPUT_DIR/ko_cpm_named.tsv" \
                    --names kegg-orthology


echo "HUMAnN KEGG mapping and normalization complete. Results are in $OUTPUT_DIR/ko_cpm_named.tsv"