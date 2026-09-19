# Walkthrough (WSL / Linux)

Step-by-step setup for SangerBCR-Flow on Windows with WSL, or on any Linux machine.
If you already have a Linux shell with Python 3.9 or newer, jump to
[section 4](#4-get-the-scripts).

The pipeline itself is three commands; sections 1-3 are one-time setup.

## 1. Install WSL (Windows only)

```powershell
# Install Ubuntu into D:\wsl
wsl --install --web-download --location D:\wsl
# List distributions with their state
wsl --list --verbose
# Stop one distribution
wsl --terminate Ubuntu
# Stop everything
wsl --shutdown
```

Then open Ubuntu from the Start menu and continue in the shell it opens.

## 2. System packages

```bash
sudo apt update
sudo apt install -y libgomp1 perl
```

`libgomp1` is a runtime library IgBLAST needs. `perl` runs IgBLAST's
`edit_imgt_file.pl` on the first run of step 2 (Ubuntu normally ships it; installing it
again is harmless).

## 3. Python environment

```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc          # reload the shell so that conda is on PATH
conda --version           # prints a version once the installation worked

conda create -n bcr python=3.9
conda activate bcr
```

## 4. Get the scripts

```bash
mkdir -p ~/bcr
cd ~/bcr
git clone https://github.com/shiningleeeee/SangerBCR-Flow.git .
mkdir -p raw_data references results      # igblast_tool/ is created by step 2
pip install -r requirements.txt
```

The scripts are run from `scripts/` and reach their neighbours through relative paths
(`../raw_data`, `../results`), so this layout matters:

```
~/bcr/
├── raw_data/       your .ab1 files
├── scripts/        01_ab1_to_fasta.py, 02_run_igblast.sh, 03_analyze_clones.py
├── results/        output, created by the scripts
├── references/     germline references, created by step 2
└── igblast_tool/   IgBLAST, created by step 2
```

## 5. Put your data in place

Copy your Sanger `.ab1` files into `raw_data/`. The file name becomes the record ID, so
name each file after the well, and after the chain if you sequence heavy and light
chains separately, for example `A1_VH.ab1` and `A1_VL.ab1`.

## 6. Run the pipeline

```bash
conda activate bcr
cd ~/bcr/scripts
python3 01_ab1_to_fasta.py     # .ab1  -> quality-trimmed FASTA
bash 02_run_igblast.sh         # FASTA -> IgBLAST annotation (downloads IgBLAST once)
python3 03_analyze_clones.py   # TSV   -> plots and clonotype tables
```

| Step | Output |
| --- | --- |
| 1 | `results/01_cleaned_sequences.fasta`, `results/01_raw_sequences.fasta`, `results/01_sequence_summary.tsv` |
| 2 | `results/02_igblast_results.tsv` |
| 3 | `results/analysis_output/` with four PNG figures, `Heavy_Clonotypes.tsv` and `Clone_Summary.tsv` |

The first run of step 2 takes a few minutes: it downloads IgBLAST and the germline
references. Later runs reuse `igblast_tool/` and `references/` and go straight to the
search.

## 7. Between runs

Rename or move `results/` before starting a new experiment, so that two experiments do
not end up in the same folder:

```bash
mv ~/bcr/results ~/bcr/results_2026-01-31
```

The next run of step 1 recreates `results/` on its own.

## 8. Notes

- Parameters such as `MIN_LEN`, `MAX_DISTANCE` and `TOP_N_CLONES_IN_PIE` are set at the
  top of the corresponding script.
- Step 2 is human-specific: `-organism human` and the human IMGT references are fixed.
- Clonotyping uses the heavy chain only; light chains are annotated and plotted.

## 9. Common problems

| Symptom | What to do |
| --- | --- |
| `no .ab1 files found in raw_data/` | run step 1 from `scripts/`, and make sure the files are directly in `raw_data/`, not in a subfolder |
| error mentioning `abi-trim` | the environment has Biopython < 1.71: `pip install -U "biopython>=1.71"` |
| step 2 fails while downloading | check the network, delete the incomplete file under `igblast_tool/` or `references/`, then run step 2 again |
| `ERROR: auxiliary data file missing` | the IgBLAST download did not complete: delete `igblast_tool/` and run step 2 again |
| `conda: command not found` | run `source ~/.bashrc` first |
| figures are empty or nearly empty | check `results/01_sequence_summary.tsv`: reads may have been dropped as too short, or no record may have been called productive |
