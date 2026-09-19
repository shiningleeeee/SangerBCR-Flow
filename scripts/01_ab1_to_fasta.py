#!/usr/bin/env python3
"""Step 1 of SangerBCR-Flow: turn Sanger .ab1 reads into quality-trimmed FASTA.

Input : ../raw_data/*.ab1
Output: ../results/01_cleaned_sequences.fasta   reads that passed the length filter
        ../results/01_raw_sequences.fasta       untrimmed reads, kept for reference
        ../results/01_sequence_summary.tsv      per file: raw length, trimmed length, status

Run from this directory:  python3 01_ab1_to_fasta.py
"""

import os
import glob
import csv
from Bio import SeqIO

INPUT_DIR = "../raw_data"                                # folder holding the .ab1 files
TRIMMED_FASTA = "../results/01_cleaned_sequences.fasta"  # quality-trimmed reads
RAW_FASTA = "../results/01_raw_sequences.fasta"          # untrimmed reads
SUMMARY_FILE = "../results/01_sequence_summary.tsv"      # one row per input file
MIN_LEN = 100                                            # minimum trimmed length (nt) to keep


def process_ab1():
    ab1_files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.ab1")))
    if not ab1_files:
        print("ERROR: no .ab1 files found in raw_data/")
        return

    os.makedirs(os.path.dirname(TRIMMED_FASTA), exist_ok=True)

    print(f"Found {len(ab1_files)} .ab1 file(s); trimming...")

    kept = 0
    skipped = 0
    failed = 0

    with open(TRIMMED_FASTA, "w") as trim_handle, \
         open(RAW_FASTA, "w") as raw_handle, \
         open(SUMMARY_FILE, "w", newline="", encoding="utf-8") as summary_handle:

        writer = csv.writer(summary_handle, delimiter="\t")
        writer.writerow(["seq_id", "raw_len", "trimmed_len", "status", "note"])

        for f in ab1_files:
            file_name = os.path.basename(f)
            seq_id = file_name.rsplit(".", 1)[0]

            try:
                # "abi" keeps the trace as stored, "abi-trim" applies Mott quality trimming
                raw_record = SeqIO.read(f, "abi")
                trim_record = SeqIO.read(f, "abi-trim")

                raw_seq = str(raw_record.seq).upper()
                trim_seq = str(trim_record.seq).upper()

                raw_handle.write(f">{seq_id}\n{raw_seq}\n")

                raw_len = len(raw_seq)
                trimmed_len = len(trim_seq)

                if trimmed_len >= MIN_LEN:
                    trim_handle.write(f">{seq_id}\n{trim_seq}\n")
                    writer.writerow([seq_id, raw_len, trimmed_len, "kept", ""])
                    kept += 1
                else:
                    writer.writerow([seq_id, raw_len, trimmed_len, "skipped", f"trimmed_len < {MIN_LEN}"])
                    print(f"Skipped short read: {seq_id} (raw={raw_len}, trimmed={trimmed_len})")
                    skipped += 1

            except Exception as e:
                writer.writerow([seq_id, "", "", "failed", str(e)])
                print(f"Could not read {file_name}: {e}")
                failed += 1

    print(f"Done. kept={kept}, skipped={skipped}, failed={failed}")
    print(f"Trimmed FASTA: {TRIMMED_FASTA}")
    print(f"Raw FASTA: {RAW_FASTA}")
    print(f"Summary: {SUMMARY_FILE}")

if __name__ == "__main__":
    process_ab1()
