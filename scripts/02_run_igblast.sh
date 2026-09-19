#!/bin/bash
# Step 2 of SangerBCR-Flow: annotate the cleaned reads with IgBLAST.
#
#   1. download IgBLAST 1.22.0 into igblast_tool/ (skipped when already present)
#   2. download the IMGT human V/D/J and the NCBI human C-region references, pass them
#      through IgBLAST's edit_imgt_file.pl and build the BLAST databases
#   3. run igblastn in AIRR rearrangement format (outfmt 19)
#
# Input : ../results/01_cleaned_sequences.fasta   produced by 01_ab1_to_fasta.py
# Output: ../results/02_igblast_results.tsv
#
# Requires bash, wget, tar, perl and internet access on the first run.
# Run from this directory:  bash 02_run_igblast.sh

set -euo pipefail

PROJECT_DIR=~/bcr
IG_DIR="$PROJECT_DIR/igblast_tool"
REF_DIR="$PROJECT_DIR/references"
RES_DIR="$PROJECT_DIR/results"
THREADS=4

export IGDATA="$IG_DIR"

mkdir -p "$IG_DIR" "$REF_DIR" "$RES_DIR"

# Download $url to $out, unless $out is already there.
download_file() {
    local url="$1"
    local out="$2"
    if [ ! -f "$out" ]; then
        wget -c --tries=3 --timeout=30 -O "$out" "$url"
    fi
}

# Build a BLAST nucleotide database next to $fasta, unless it is already there.
build_db_if_needed() {
    local fasta="$1"
    if [ ! -f "${fasta}.nsq" ]; then
        "$IG_DIR/bin/makeblastdb" -dbtype nucl -parse_seqids -in "$fasta"
    fi
}

echo "[1/5] IgBLAST"
if [ ! -x "$IG_DIR/bin/igblastn" ]; then
    cd "$IG_DIR"
    download_file \
        "https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/1.22.0/ncbi-igblast-1.22.0-x64-linux.tar.gz" \
        "ncbi-igblast-1.22.0-x64-linux.tar.gz"
    tar -xzf ncbi-igblast-1.22.0-x64-linux.tar.gz --strip-components=1
    rm -f ncbi-igblast-1.22.0-x64-linux.tar.gz
else
    echo "IgBLAST already present, skipping the download."
fi

echo "[2/5] Germline references"
cd "$REF_DIR"

if [ ! -f "$REF_DIR/imgt_ig_v_human.fa" ] || \
   [ ! -f "$REF_DIR/imgt_ig_d_human.fa" ] || \
   [ ! -f "$REF_DIR/imgt_ig_j_human.fa" ]; then

    for gene in IGHV IGKV IGLV IGHD IGHJ IGKJ IGLJ; do
        download_file \
            "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/${gene}.fasta" \
            "${gene}.fasta"
    done

    cat IGHV.fasta IGKV.fasta IGLV.fasta > imgt_ig_v_human.fasta
    cat IGHD.fasta > imgt_ig_d_human.fasta
    cat IGHJ.fasta IGKJ.fasta IGLJ.fasta > imgt_ig_j_human.fasta

    # IgBLAST needs the IMGT numbering gaps that edit_imgt_file.pl inserts
    "$IG_DIR/bin/edit_imgt_file.pl" imgt_ig_v_human.fasta > imgt_ig_v_human.fa
    "$IG_DIR/bin/edit_imgt_file.pl" imgt_ig_d_human.fasta > imgt_ig_d_human.fa
    "$IG_DIR/bin/edit_imgt_file.pl" imgt_ig_j_human.fasta > imgt_ig_j_human.fa

    rm -f IGHV.fasta IGKV.fasta IGLV.fasta IGHD.fasta IGHJ.fasta IGKJ.fasta IGLJ.fasta
else
    echo "V/D/J references already present, skipping the download."
fi

if [ ! -f "$REF_DIR/human_C_genes.fasta" ]; then
    download_file \
        "https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/database/ncbi_human_c_genes.tar" \
        "ncbi_human_c_genes.tar"
    tar -xf ncbi_human_c_genes.tar
    rm -f ncbi_human_c_genes.tar

    "$IG_DIR/bin/blastdbcmd" -db ncbi_human_c_genes -entry all -out human_C_genes.fasta
else
    echo "C-region reference already present, skipping the download."
fi

echo "[3/5] BLAST databases"
build_db_if_needed "$REF_DIR/imgt_ig_v_human.fa"
build_db_if_needed "$REF_DIR/imgt_ig_d_human.fa"
build_db_if_needed "$REF_DIR/imgt_ig_j_human.fa"

if [ -f "$REF_DIR/human_C_genes.fasta" ]; then
    build_db_if_needed "$REF_DIR/human_C_genes.fasta"
fi

echo "[4/5] Input check"
QUERY_FILE="$RES_DIR/01_cleaned_sequences.fasta"
OUT_FILE="$RES_DIR/02_igblast_results.tsv"

if [ ! -s "$QUERY_FILE" ]; then
    echo "ERROR: input FASTA missing or empty: $QUERY_FILE"
    exit 1
fi

if [ ! -f "$IG_DIR/optional_file/human_gl.aux" ]; then
    echo "ERROR: auxiliary data file missing: $IG_DIR/optional_file/human_gl.aux"
    exit 1
fi

echo "[5/5] igblastn"
CMD=(
    "$IG_DIR/bin/igblastn"
    -germline_db_V "$REF_DIR/imgt_ig_v_human.fa"
    -germline_db_D "$REF_DIR/imgt_ig_d_human.fa"
    -germline_db_J "$REF_DIR/imgt_ig_j_human.fa"
    -auxiliary_data "$IG_DIR/optional_file/human_gl.aux"
    -organism human
    -domain_system imgt
    -ig_seqtype Ig
    -query "$QUERY_FILE"
    -outfmt 19
    -num_threads "$THREADS"
    -out "$OUT_FILE"
)

if [ -f "$REF_DIR/human_C_genes.fasta" ]; then
    CMD+=(-c_region_db "$REF_DIR/human_C_genes.fasta")
fi

"${CMD[@]}"

echo "IgBLAST finished."
echo "Output: $OUT_FILE"
