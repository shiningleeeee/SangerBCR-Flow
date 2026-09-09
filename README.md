# SangerBCR-Flow

SangerBCR-Flow is a simple and reproducible pipeline for the analysis of small-scale BCR sequencing datasets derived from Sanger sequencing. Starting from raw .ab1 files, the workflow performs sequence extraction and trimming, IgBLAST-based V(D)J annotation, and heavy-chain-only clonotype clustering. Clonotypes are defined using a general framework based on identical heavy-chain V gene, J gene, junction nucleotide length, and junction nucleotide similarity. The pipeline further provides standard descriptive visualizations, including gene usage, somatic hypermutation profiles, CDR3 length distributions, and clonal composition plots. SangerBCR-Flow is intended for low-throughput BCR studies, particularly in settings such as single-cell antibody discovery, targeted validation, and exploratory repertoire analysis.

## Repository layout

```
SangerBCR-Flow/
├── README.md                 This file
├── requirements.txt          Python dependencies (Biopython, pandas, seaborn, matplotlib, numpy)
├── .gitignore
├── docs/
│   └── walkthrough.md        Step-by-step setup & usage tutorial (WSL / Linux, 中文)
└── scripts/
    ├── 01_ab1_to_fasta.py    .ab1 → quality-trimmed FASTA + per-file summary TSV
    ├── 02_run_igblast.sh     IgBLAST download, IMGT reference preparation & V(D)J annotation (AIRR TSV)
    └── 03_analyze_clones.py  Gene usage / SHM / CDR3 / heavy-chain clonotyping plots & tables
```

## Quick start (WSL / Linux)

1. Set up WSL, Miniconda and the `bcr` Python environment — see [docs/walkthrough.md](docs/walkthrough.md).
2. Put your `.ab1` files into `raw_data/` and copy the three scripts into `scripts/` next to it.
3. Run the pipeline in order:

```bash
cd ~/bcr/scripts
python3 01_ab1_to_fasta.py    # raw_data/*.ab1 -> results/01_cleaned_sequences.fasta
bash 02_run_igblast.sh        # IgBLAST 1.22.0 + IMGT human references -> results/02_igblast_results.tsv
python3 03_analyze_clones.py  # plots + clonotype tables -> results/analysis_output/
```

## Requirements & notes

- Linux (WSL recommended). Python 3.9+; Biopython **≥ 1.71** (the `abi-trim` SeqIO format used in step 01 requires it); pandas, seaborn, matplotlib, numpy (see `requirements.txt`).
- IgBLAST 1.22.0 and the IMGT human germline references are downloaded automatically by `02_run_igblast.sh` on first run (requires internet access and `wget`, `tar`, `perl`).
- Human IGH/IGK/IGL only (`-organism human` fixed in the pipeline).
- Move `results/` outputs away between runs so consecutive runs do not mix.
- Clone sizes count retained sequences/wells sharing a clonotype — a sampling-based metric, not an absolute B-cell frequency.
- The three scripts were extracted verbatim from the original single-file `code`; no code changes were made in this reorganization.
