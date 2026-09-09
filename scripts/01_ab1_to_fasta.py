#!/usr/bin/env python3
import os
import glob
import csv
from Bio import SeqIO

INPUT_DIR = "../raw_data"
TRIMMED_FASTA = "../results/01_cleaned_sequences.fasta"
RAW_FASTA = "../results/01_raw_sequences.fasta"
SUMMARY_FILE = "../results/01_sequence_summary.tsv"
MIN_LEN = 100

def process_ab1():
    ab1_files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.ab1")))
    if not ab1_files:
        print("❌ 错误: raw_data 文件夹里没有找到 .ab1 文件！")
        return

    os.makedirs(os.path.dirname(TRIMMED_FASTA), exist_ok=True)

    print(f"📂 发现 {len(ab1_files)} 个 .ab1 文件，开始处理...")

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
                    print(f"⚠️ 跳过短序列: {seq_id} (raw={raw_len}, trimmed={trimmed_len})")
                    skipped += 1

            except Exception as e:
                writer.writerow([seq_id, "", "", "failed", str(e)])
                print(f"❌ 读取错误 {file_name}: {e}")
                failed += 1

    print(f"✅ 完成！保留 {kept} 条，跳过 {skipped} 条，失败 {failed} 条")
    print(f"📄 Trimmed fasta: {TRIMMED_FASTA}")
    print(f"📄 Raw fasta: {RAW_FASTA}")
    print(f"📄 Summary: {SUMMARY_FILE}")

if __name__ == "__main__":
    process_ab1()
