# SangerBCR-Flow

A small, readable pipeline that turns Sanger-sequenced B-cell receptor amplicons into
clonotype calls and standard repertoire plots.

Starting from raw `.ab1` files, the workflow performs quality trimming, IgBLAST V(D)J
annotation and heavy-chain clonotype clustering, and then produces descriptive figures
(V gene usage, somatic hypermutation, CDR3 length, clonal composition) together with
clonotype tables. It is aimed at low-throughput work: single-cell antibody discovery,
targeted validation of a few dozen wells, and exploratory repertoire analysis.

A step-by-step setup guide (WSL / Linux, with a troubleshooting table) is in
[docs/walkthrough.md](docs/walkthrough.md).

## Requirements

- Linux, or Windows with WSL (the walkthrough uses WSL + Ubuntu).
- Python 3.9 or newer.
- Biopython >= 1.71 (the `abi-trim` SeqIO format used in step 1 was added in 1.71),
  pandas, seaborn, matplotlib, numpy — see `requirements.txt`.
- `wget`, `tar` and `perl`, plus internet access on the first run: step 2 downloads
  IgBLAST 1.22.0 and the IMGT human germline references by itself.

Tested with Python 3.9, Biopython 1.71+, pandas 1.5+/2.x, seaborn 0.13+, matplotlib 3.x.

## Directory layout

The scripts are run from `scripts/` and reach their neighbours through relative paths
(`../raw_data`, `../results`), so keep this layout:

```
bcr/
├── raw_data/       your Sanger .ab1 files, one file per well / per chain
├── scripts/        the three scripts from this repository
├── results/        all output, created by the scripts
├── references/     germline references, created by step 2
└── igblast_tool/   IgBLAST, created by step 2
```

## Install

```bash
conda create -n bcr python=3.9
conda activate bcr
pip install -r requirements.txt
```

## Run

```bash
cd ~/bcr/scripts
python3 01_ab1_to_fasta.py     # .ab1  -> quality-trimmed FASTA
bash 02_run_igblast.sh         # FASTA -> IgBLAST annotation (AIRR TSV)
python3 03_analyze_clones.py   # TSV   -> plots and clonotype tables
```

The three steps run in order: step 1 must finish before step 2, and step 2 before step 3.

## Input

- `raw_data/*.ab1`: one Sanger read per file. The file name becomes the record ID
  (`A1_VH.ab1` -> `A1_VH`), so name files in a way that identifies the well and the
  chain.
- Reads are trimmed with Biopython's `abi-trim` (Mott algorithm, default parameters).
  Reads shorter than `MIN_LEN` nucleotides after trimming are dropped, and every file is
  listed in the step 1 summary table with its raw and trimmed length.

## Output

| File | Content |
| --- | --- |
| `results/01_cleaned_sequences.fasta` | quality-trimmed reads that passed the length filter |
| `results/01_raw_sequences.fasta` | untrimmed reads, kept for reference |
| `results/01_sequence_summary.tsv` | per file: raw length, trimmed length, status, note |
| `results/02_igblast_results.tsv` | IgBLAST annotation in AIRR rearrangement format (`-outfmt 19`) |
| `results/analysis_output/Heavy_Clonotypes.tsv` | every heavy-chain record with its clonotype ID and clonotype size |
| `results/analysis_output/Clone_Summary.tsv` | one row per clonotype |
| `results/analysis_output/*.png` | V gene usage (heavy and light), SHM distribution, CDR3 length distribution, clone pie chart |

## Parameters

Tuning happens at the top of `scripts/03_analyze_clones.py`, and `MIN_LEN` at the top of
`scripts/01_ab1_to_fasta.py`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `MIN_LEN` | 100 | minimum trimmed read length (nt) kept by step 1 |
| `MAX_DISTANCE` | 0.15 | maximum normalised Hamming distance between junction nucleotide sequences for two heavy chains to be called clonal |
| `DROP_DUPLICATE_SEQUENCE_ID` | False | drop records sharing a `sequence_id` before clonotyping |
| `TOP_N_CLONES_IN_PIE` | 19 | clones drawn individually in the pie chart; smaller clones are pooled into "Others" |

## How clonotypes are defined

Clonotyping is done on the heavy chain only. Two records belong to the same clonotype
when they share the V gene, the J gene and the junction nucleotide length, and their
junction nucleotide sequences differ by at most `MAX_DISTANCE` (mismatches divided by
length). Links are merged transitively (single linkage), so a chain of pairwise-similar
junctions ends up in one clonotype. At the default threshold this allows up to 3
mismatches in a 24 nt junction and up to 9 mismatches in a 60 nt junction.

Records that are not productive, have no V or J call, or have no usable junction
sequence are excluded from clonotyping. The console output reports how many records
were retained at each step.

## Notes and limitations

- **Human immunoglobulin only.** `-organism human` and the human IMGT references are
  fixed in step 2; another species needs its own germline references.
- **Heavy chain only for clonotyping.** Light chains are annotated and plotted, but they
  are not clustered.
- **Somatic hypermutation is a quick descriptive measure.** It is read off IgBLAST's
  `v_identity`, the percentage of identical positions in the V gene alignment
  (`mutation_rate = 100 - v_identity`). It is not a substitution-model estimate.
- **Clone sizes are a sampling metric.** They count retained sequences (wells) sharing a
  clonotype, not absolute B-cell frequencies.
- **Rename or move `results/` between runs** so that two experiments do not mix. Step 1
  recreates `results/` on its own.

## Troubleshooting

| Symptom | Cause / fix |
| --- | --- |
| `no .ab1 files found in raw_data/` | run step 1 from `scripts/`, with the `.ab1` files directly in `raw_data/` |
| error mentioning `abi-trim` | Biopython older than 1.71: `pip install -U "biopython>=1.71"` |
| step 2 stops while downloading | check the network, delete the incomplete file in `igblast_tool/` or `references/`, then run step 2 again |
| `ERROR: auxiliary data file missing` | the IgBLAST download did not complete: delete `igblast_tool/` and run step 2 again |
| `conda: command not found` | `source ~/.bashrc` first |
| figures are empty or nearly empty | check `results/01_sequence_summary.tsv`: reads may have been dropped as too short, or no record may have been called productive |
