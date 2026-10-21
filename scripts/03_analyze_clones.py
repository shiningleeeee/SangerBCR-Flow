#!/usr/bin/env python3
"""Step 3 of SangerBCR-Flow: summary plots and heavy-chain clonotype calling.

Input : ../results/02_igblast_results.tsv   AIRR rearrangement TSV from step 2
Output: ../results/analysis_output/
          Heavy_V_Usage.png, Light_V_Usage.png
          SHM_Distribution.png
          CDR3_Length_Distribution.png
          Clone_Pie_Chart.png
          Heavy_Clonotypes.tsv, Clone_Summary.tsv
          stage_accounting.tsv   how many records survived each filter
          run_info.json          parameters, library versions and the input checksum

Clonotypes are called on the heavy chain only: records sharing the V gene, the J gene
and the junction length are linked when their junction nucleotide sequences differ by at
most MAX_DISTANCE (normalised Hamming distance), and linked records are merged
transitively (single linkage). This is the clonotype definition implemented by Change-O
(Gupta et al., Bioinformatics 2015, doi:10.1093/bioinformatics/btv359) and used by the
AIRR Community to populate clone_id (Vander Heiden et al., Front Immunol 2018,
doi:10.3389/fimmu.2018.02206).

Clonotype IDs are assigned deterministically (largest clone first, ties broken on the
junction sequence), so the same input always yields the same IDs.

Run from this directory:  python3 03_analyze_clones.py
Run with --help to see the options; the defaults below are only the fallbacks.
"""


import argparse
import hashlib
import json
import os
import platform
import re
import sys
from datetime import datetime

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt


# =========================
# 1. Paths and parameters
# =========================
INPUT_FILE = "../results/02_igblast_results.tsv"
OUTPUT_DIR = "../results/analysis_output"

# Maximum normalised Hamming distance between junction nucleotide sequences
MAX_DISTANCE = 0.15

# Drop records that share a sequence_id before clonotyping
DROP_DUPLICATE_SEQUENCE_ID = False

# Clones drawn individually in the pie chart; smaller clones are pooled into "Others"
TOP_N_CLONES_IN_PIE = 19


# =========================
# 2. Shared constants
# =========================
# Values that count as true in the AIRR TSV (IgBLAST writes T / F)
TRUTHY = {"T", "TRUE", "1", "YES", "Y"}
# Strings that stand for "no value" in the AIRR TSV
FAKE_NULL = {"", "NA", "NAN", "NONE", "NULL"}
# Gap characters that IMGT-gapped sequences may contain
GAP_CHARS = (".", "-")


# =========================
# 3. Helper functions
# =========================
def is_truthy(x):
    """Return True when an AIRR boolean field (T / F, True / False, 1 / 0) is true."""
    if pd.isna(x):
        return False
    return str(x).strip().upper() in TRUTHY


def get_chain_type(v_call):
    """Return Heavy, Light or Unknown from a V gene call."""
    if pd.isna(v_call):
        return "Unknown"
    s = str(v_call).upper()
    if "IGH" in s:
        return "Heavy"
    if "IGK" in s or "IGL" in s:
        return "Light"
    return "Unknown"


def is_heavy_chain(v_call):
    """Return True when the V gene call belongs to a heavy chain."""
    return get_chain_type(v_call) == "Heavy"


def first_call(x):
    """Take the first of the possibly comma-separated gene calls reported by IgBLAST."""
    if pd.isna(x):
        return ""
    return str(x).split(",")[0].strip()


def strip_allele(x):
    """Drop the allele suffix of a gene call: IGHV3-23*01 becomes IGHV3-23."""
    x = first_call(x)
    if not x:
        return ""
    return x.split("*")[0].strip()


def strip_gaps(x):
    """Remove IMGT gap characters from a nucleotide or amino acid string."""
    for gap in GAP_CHARS:
        x = x.replace(gap, "")
    return x


def normalize_cdr3_aa(x):
    """Return an uppercase CDR3 amino acid string, or None when the value is unusable.

    Gaps are removed and values containing a stop codon are rejected, so that neither
    inflates the length distribution. Empty fields stay empty rather than entering the
    plot as "nan".
    """
    if pd.isna(x):
        return None

    s = strip_gaps(str(x).strip().upper())

    if s in FAKE_NULL:
        return None

    if not re.fullmatch(r"[A-Z]+", s):
        return None

    return s


def normalize_junction_nt(x):
    """Return the junction nucleotide sequence in upper case, or None when unusable.

    Gaps are removed first: AIRR allows an IMGT-gapped junction, and a gap would
    otherwise fail the character check and silently drop every record. Empty fields stay
    empty rather than entering clonotyping as "nan".
    """
    if pd.isna(x):
        return None

    s = strip_gaps(str(x).strip().upper())

    if s in FAKE_NULL:
        return None

    if not re.fullmatch(r"[A-Z]+", s):
        return None

    return s


def normalized_hamming_distance(seq1, seq2):
    """Return the normalised Hamming distance (mismatches / length) of two sequences.

    Both sequences must have the same length.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Hamming distance requires equal-length sequences.")
    if len(seq1) == 0:
        return 0.0

    mismatches = sum(a != b for a, b in zip(seq1, seq2))
    return mismatches / len(seq1)


def file_sha256(path):
    """Return the SHA-256 checksum of a file, or None when it cannot be read."""
    try:
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError:
        return None


# =========================
# 4. Union-Find
# =========================
class UnionFind:
    """Disjoint-set structure used to merge junction sequences that are linked."""

    def __init__(self, items):
        self.parent = {x: x for x in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra = self.find(a)
        rb = self.find(b)
        if ra != rb:
            self.parent[rb] = ra


# =========================
# 5. Plot functions
# =========================
def plot_gene_usage(df, output_dir):
    """Write the heavy and light chain V gene usage bar charts."""
    print("\nGene usage plots ...")

    # Heavy chain V genes
    heavy_df = df[df["chain_type"] == "Heavy"].copy()
    heavy_v = heavy_df["v_gene"].dropna().value_counts().head(20)

    if not heavy_v.empty:
        heavy_plot = (
            heavy_v
            .rename_axis("v_gene")
            .reset_index(name="count")
        )

        plt.figure(figsize=(12, 8))
        sns.barplot(
            data=heavy_plot,
            x="count",
            y="v_gene",
            hue="v_gene",
            palette="viridis",
            legend=False
        )
        plt.title(f"Top 20 Heavy Chain V-Gene Usage (n = {len(heavy_df)} records)")
        plt.xlabel("Count")
        plt.ylabel("V Gene")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "Heavy_V_Usage.png"), dpi=300)
        plt.close()

    # Light chain V genes
    light_df = df[df["chain_type"] == "Light"].copy()
    light_v = light_df["v_gene"].dropna().value_counts().head(20)

    if not light_v.empty:
        light_plot = (
            light_v
            .rename_axis("v_gene")
            .reset_index(name="count")
        )

        plt.figure(figsize=(12, 8))
        sns.barplot(
            data=light_plot,
            x="count",
            y="v_gene",
            hue="v_gene",
            palette="plasma",
            legend=False
        )
        plt.title(f"Top 20 Light Chain V-Gene Usage (n = {len(light_df)} records)")
        plt.xlabel("Count")
        plt.ylabel("V Gene")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "Light_V_Usage.png"), dpi=300)
        plt.close()


def plot_shm(df, output_dir):
    """Write the somatic hypermutation histogram for both chain types."""
    print("\nSHM analysis ...")

    if "v_identity" not in df.columns:
        print("No v_identity column found; skipping the SHM plot.")
        return

    plot_df = df.copy()
    plot_df["v_identity_num"] = pd.to_numeric(plot_df["v_identity"], errors="coerce")
    plot_df = plot_df.dropna(subset=["v_identity_num"]).copy()

    if plot_df.empty:
        print("No usable v_identity values; skipping the SHM plot.")
        return

    # IgBLAST reports v_identity as a percentage (0-100), not as a fraction. Note that
    # AIRR defines v_identity as a fraction, so the value is converted here and the
    # column in the TSV is left as IgBLAST wrote it.
    plot_df["mutation_rate"] = 100 - plot_df["v_identity_num"]

    plt.figure(figsize=(8, 6))
    sns.histplot(data=plot_df, x="mutation_rate", hue="chain_type", kde=True, bins=30)
    plt.title(f"Somatic Hypermutation (n = {len(plot_df)} records)")
    plt.xlabel("Mutation Rate (%) (100 - v_identity)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "SHM_Distribution.png"), dpi=300)
    plt.close()


def plot_cdr3_length(df, output_dir):
    """Write the CDR3 amino acid length distribution."""
    print("\nCDR3 length distribution ...")

    # IgBLAST writes cdr3_aa whenever the domain system is IMGT. With -domain_system
    # imgt the junction *is* the CDR3, so junction_aa is an equivalent fallback for
    # older IgBLAST builds that do not emit cdr3_aa.
    column = next((c for c in ("cdr3_aa", "junction_aa") if c in df.columns), None)

    if column is None:
        print("Neither cdr3_aa nor junction_aa found; skipping the CDR3 length plot.")
        return

    plot_df = df.copy()
    plot_df["cdr3_aa_clean"] = plot_df[column].apply(normalize_cdr3_aa)
    plot_df = plot_df.dropna(subset=["cdr3_aa_clean"]).copy()

    if plot_df.empty:
        print(f"No usable {column} values; skipping the CDR3 length plot.")
        return

    print(f"Using the {column} column.")

    plot_df["cdr3_length"] = plot_df["cdr3_aa_clean"].str.len()

    plt.figure(figsize=(8, 6))
    sns.histplot(data=plot_df, x="cdr3_length", hue="chain_type", kde=True, element="step")
    plt.title(f"CDR3 Amino Acid Length Distribution (n = {len(plot_df)} records)")
    plt.xlabel("CDR3 Length (aa)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "CDR3_Length_Distribution.png"), dpi=300)
    plt.close()


def plot_clone_donut(clone_counts, total_cells, output_dir, top_n=19):
    """Write the clone size pie chart, showing the top_n largest clones."""
    print("\nClone pie chart ...")

    if clone_counts.empty:
        print("No clone counts available; skipping the pie chart.")
        return

    top_clones = clone_counts.head(top_n)
    others_count = (
        int(clone_counts.iloc[top_n:]["clone_size"].sum())
        if len(clone_counts) > top_n else 0
    )

    labels = [f"Clone {int(cid)}" for cid in top_clones["clone_id"]]
    sizes = list(top_clones["clone_size"])

    if others_count > 0:
        labels.append("Others")
        sizes.append(others_count)

    plt.figure(figsize=(10, 10))
    plt.pie(
        sizes,
        labels=labels,
        autopct="%1.1f%%",
        startangle=140,
        pctdistance=0.82,
        textprops={"fontsize": 12}
    )

    centre_circle = plt.Circle(
        (0, 0),
        0.70,
        fc="white",
        edgecolor="black",
        linewidth=1.2
    )
    fig = plt.gcf()
    fig.gca().add_artist(centre_circle)

    plt.text(
        0, 0,
        f"N = {total_cells}",
        ha="center",
        va="center",
        fontsize=18,
        fontweight="bold"
    )

    plt.title(f"Clonal Distribution (Top {top_n})", fontsize=18, pad=20)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Clone_Pie_Chart.png"), dpi=300)
    plt.close()


# =========================
# 6. Input loading and filtering
# =========================
def load_raw_data(input_file):
    """Read the AIRR TSV and keep the productive records."""
    print(f"Loading data from: {input_file}")
    df = pd.read_csv(input_file, sep="\t", dtype=str)

    required_cols = ["sequence_id", "productive", "v_call", "j_call"]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Input file is missing required column(s): {missing_cols}\n"
            f"Required at least: {required_cols}"
        )

    # Productive records only
    df["productive_bool"] = df["productive"].apply(is_truthy)
    df = df[df["productive_bool"]].copy()

    # Chain type from the V gene call
    df["chain_type"] = df["v_call"].apply(get_chain_type)

    # V and J gene without the allele suffix
    df["v_gene"] = df["v_call"].apply(strip_allele)
    df["j_gene"] = df["j_call"].apply(strip_allele)

    return df.reset_index(drop=True)


def prepare_heavy_clonotyping_input(df, drop_duplicate_sequence_id=False, accounting=None):
    """Keep the productive heavy-chain records that can be clonotyped."""
    if "junction" not in df.columns:
        raise ValueError("Input file has no 'junction' column; clonotyping is not possible.")

    def record(stage, value, note=""):
        if accounting is not None:
            accounting.append({"stage": stage, "records": value, "note": note})

    heavy_df = df[df["chain_type"] == "Heavy"].copy()
    record("productive heavy-chain records", len(heavy_df))

    heavy_df["junction_nt"] = heavy_df["junction"].apply(normalize_junction_nt)
    no_v = int(heavy_df["v_gene"].eq("").sum())
    no_j = int(heavy_df["j_gene"].eq("").sum())
    no_junction = int(heavy_df["junction_nt"].isna().sum())
    record("without a V gene call", no_v)
    record("without a J gene call", no_j)
    record("without a usable junction", no_junction)

    heavy_df = heavy_df[
        heavy_df["v_gene"].ne("")
        & heavy_df["j_gene"].ne("")
        & heavy_df["junction_nt"].notna()
    ].copy()

    heavy_df["junction_length"] = heavy_df["junction_nt"].str.len()

    if drop_duplicate_sequence_id:
        before = len(heavy_df)
        heavy_df = heavy_df.drop_duplicates(subset=["sequence_id"]).copy()
        after = len(heavy_df)
        print(f"Dropped duplicated sequence_id: {before - after}")
        record("duplicated sequence_id dropped", before - after)

    heavy_df = heavy_df.reset_index(drop=True)

    record("retained for clonotyping", len(heavy_df))
    print(f"\nRetained productive heavy-chain records for clonotyping: {len(heavy_df)}")
    return heavy_df


# =========================
# 7. Heavy-chain clonotyping
# =========================
def run_heavy_clonotyping(df, max_distance=0.15):
    """Assign a clonotype ID to every heavy-chain record.

    Records are grouped by V gene, J gene and junction length. Inside a group, pairs
    whose normalised Hamming distance is at most max_distance are linked, and linked
    records are merged transitively (single linkage).

    The union-find roots are then renumbered so that the IDs are a stable property of
    the data rather than of the input order: clones are ranked by size (largest first)
    and ties are broken on V gene, J gene and the representative junction, so that
    re-running on the same reads always yields the same clone_id.
    """
    df = df.copy().reset_index(drop=True)

    df["group_key"] = (
        df["v_gene"] + "|" +
        df["j_gene"] + "|" +
        df["junction_length"].astype(str)
    )

    df["clone_root"] = -1
    root_counter = 1

    for _, group in df.groupby("group_key", sort=False):
        indices = group.index.tolist()

        if len(indices) == 1:
            df.at[indices[0], "clone_root"] = root_counter
            root_counter += 1
            continue

        uf = UnionFind(indices)

        for i in range(len(indices)):
            idx_i = indices[i]
            seq_i = df.at[idx_i, "junction_nt"]

            for j in range(i + 1, len(indices)):
                idx_j = indices[j]
                seq_j = df.at[idx_j, "junction_nt"]

                dist = normalized_hamming_distance(seq_i, seq_j)
                if dist <= max_distance:
                    uf.union(idx_i, idx_j)

        root_to_root = {}
        for idx in indices:
            root = uf.find(idx)
            if root not in root_to_root:
                root_to_root[root] = root_counter
                root_counter += 1
            df.at[idx, "clone_root"] = root_to_root[root]

    # Deterministic renumbering: size first, then the junction sequence itself.
    sizes = df.groupby("clone_root").size()
    representatives = (
        df.sort_values(["junction_nt", "sequence_id"])
        .groupby("clone_root")
        .first()[["v_gene", "j_gene", "junction_nt"]]
    )

    order = sorted(
        sizes.index,
        key=lambda root: (
            -int(sizes[root]),
            representatives.at[root, "v_gene"],
            representatives.at[root, "j_gene"],
            representatives.at[root, "junction_nt"],
        )
    )
    root_to_clone = {root: rank + 1 for rank, root in enumerate(order)}

    df["clone_id"] = df["clone_root"].map(root_to_clone)
    df = df.drop(columns=["clone_root"])

    return df.copy()


# =========================
# 8. Output tables
# =========================
def make_output_tables(df):
    """Return the per-record table, the per-clonotype summary and the clone sizes."""
    # One copy keeps the many assignments above from fragmenting the frame
    df = df.copy()

    clone_counts = (
        df
        .groupby("clone_id", as_index=False)
        .size()
        .rename(columns={"size": "clone_size"})
        .sort_values(["clone_size", "clone_id"], ascending=[False, True])
        .reset_index(drop=True)
    )

    clone_size_map = clone_counts.set_index("clone_id")["clone_size"]
    # clone_rank is the position in the size ranking (largest clone = 1), kept for the
    # figure labels. It follows the order of clone_counts, so it is known before the
    # summary table is built.
    clone_rank_map = pd.Series(
        range(1, len(clone_counts) + 1), index=clone_counts["clone_id"]
    )

    # Both columns are assigned at once: a frame that is extended one column at a time
    # becomes fragmented and pandas warns about it further down.
    result_df = df.assign(
        clone_size=df["clone_id"].map(clone_size_map),
        clone_rank=df["clone_id"].map(clone_rank_map),
    ).copy()

    # Representative record of each clonotype: the most frequent junction in the clone,
    # with ties broken on sequence_id so that the choice is deterministic. The number of
    # distinct junctions and the largest within-clone distance are reported alongside, so
    # that a "clone" made of several distinct junctions can be told apart from one made
    # of repeated reads of the same molecule.
    junctions = (
        result_df
        .sort_values(["clone_id", "junction_nt", "sequence_id"])
        .groupby(["clone_id", "junction_nt"], as_index=False)
        .agg(n_records=("sequence_id", "size"), sequence_id=("sequence_id", "first"))
    )
    representatives = (
        junctions
        .sort_values(["clone_id", "n_records", "junction_nt"], ascending=[True, False, True])
        .groupby("clone_id", as_index=False)
        .first()[["clone_id", "junction_nt", "sequence_id"]]
        .rename(columns={
            "junction_nt": "representative_junction_nt",
            "sequence_id": "representative_sequence_id"
        })
    )

    distinct = (
        result_df
        .groupby("clone_id")["junction_nt"]
        .nunique()
        .reset_index()
        .rename(columns={"junction_nt": "n_distinct_junctions"})
    )

    spread = []
    for clone_id, group in result_df.groupby("clone_id"):
        sequences = group["junction_nt"].tolist()
        worst = 0.0
        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                if len(sequences[i]) == len(sequences[j]):
                    worst = max(worst, normalized_hamming_distance(sequences[i], sequences[j]))
        spread.append({"clone_id": clone_id, "max_within_clone_distance": round(worst, 4)})
    spread = pd.DataFrame(spread, columns=["clone_id", "max_within_clone_distance"])

    summary_df = (
        representatives
        .merge(clone_counts, on="clone_id", how="left")
        .merge(distinct, on="clone_id", how="left")
        .merge(spread, on="clone_id", how="left")
        .merge(
            result_df.groupby("clone_id", as_index=False)
            .first()[["clone_id", "v_gene", "j_gene", "junction_length"]],
            on="clone_id",
            how="left"
        )
        [["clone_id", "clone_size", "n_distinct_junctions", "max_within_clone_distance",
          "v_gene", "j_gene", "junction_length",
          "representative_junction_nt", "representative_sequence_id"]]
        .sort_values(["clone_size", "clone_id"], ascending=[False, True])
        .reset_index(drop=True)
    )

    # clone_rank is the position in the size ranking, kept for display
    summary_df.insert(1, "clone_rank", range(1, len(summary_df) + 1))

    return result_df, summary_df, clone_counts


# =========================
# 9. Provenance
# =========================
def write_run_info(path, args, input_file, n_records, n_clones):
    """Write the parameters, library versions and input checksum used by this run."""
    info = {
        "tool": "SangerBCR-Flow step 3 (03_analyze_clones.py)",
        "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
        "parameters": {
            "input": input_file,
            "max_distance": args.max_distance,
            "drop_duplicate_sequence_id": bool(args.drop_duplicate_sequence_id),
            "top_n_clones_in_pie": args.top_n_clones,
        },
        "input_sha256": file_sha256(input_file),
        "results": {"records_clonotyped": n_records, "clones": n_clones},
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "seaborn": sns.__version__,
            "matplotlib": matplotlib.__version__,
        },
    }

    with open(path, "w", encoding="utf-8") as handle:
        json.dump(info, handle, indent=2, sort_keys=False)
        handle.write("\n")


def write_stage_accounting(path, accounting):
    """Write the stage-accounting table (records surviving each filter)."""
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("stage\trecords\tnote\n")
        for row in accounting:
            handle.write(f"{row['stage']}\t{row['records']}\t{row['note']}\n")


# =========================
# 10. Command line
# =========================
def parse_args(argv=None):
    """Return the command-line arguments, with the module constants as defaults."""
    parser = argparse.ArgumentParser(
        description="Summary plots and heavy-chain clonotype calling for SangerBCR-Flow.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", default=INPUT_FILE,
                        help="AIRR rearrangement TSV produced by step 2")
    parser.add_argument("--output-dir", default=OUTPUT_DIR,
                        help="directory for the figures and tables")
    parser.add_argument("--max-distance", type=float, default=MAX_DISTANCE,
                        help="maximum normalised Hamming distance between junctions")
    parser.add_argument("--top-n-clones", type=int, default=TOP_N_CLONES_IN_PIE,
                        help="clones drawn individually in the pie chart")
    parser.add_argument("--drop-duplicate-sequence-id", action="store_true",
                        default=DROP_DUPLICATE_SEQUENCE_ID,
                        help="drop records that share a sequence_id before clonotyping")
    return parser.parse_args(argv)


# =========================
# 11. Main
# =========================
def main(argv=None):
    args = parse_args(argv)
    output_dir = args.output_dir

    os.makedirs(output_dir, exist_ok=True)

    accounting = []

    # Load and filter the productive records
    raw_df = load_raw_data(args.input)
    print(f"Total productive sequences: {len(raw_df)}")
    print(raw_df["chain_type"].value_counts())
    accounting.append({"stage": "input records", "records": len(raw_df), "note": "productive only"})
    accounting.append({"stage": "  of which heavy chain",
                       "records": int((raw_df["chain_type"] == "Heavy").sum()), "note": ""})
    accounting.append({"stage": "  of which light chain",
                       "records": int((raw_df["chain_type"] == "Light").sum()), "note": ""})

    # Descriptive plots, on all productive records
    plot_gene_usage(raw_df, output_dir)
    plot_shm(raw_df, output_dir)
    plot_cdr3_length(raw_df, output_dir)

    # Heavy-chain input for clonotyping
    heavy_df = prepare_heavy_clonotyping_input(
        raw_df,
        drop_duplicate_sequence_id=args.drop_duplicate_sequence_id,
        accounting=accounting,
    )

    if len(heavy_df) == 0:
        raise ValueError("No heavy-chain records left after filtering; nothing to clonotype.")

    print("\nHeavy-chain clonotyping ...")
    cloned_df = run_heavy_clonotyping(heavy_df, max_distance=args.max_distance)

    print(f"Total Clones Identified: {cloned_df['clone_id'].nunique()}")

    # Tables
    result_df, summary_df, clone_counts = make_output_tables(cloned_df)

    result_file = os.path.join(output_dir, "Heavy_Clonotypes.tsv")
    summary_file = os.path.join(output_dir, "Clone_Summary.tsv")

    result_df.to_csv(result_file, sep="\t", index=False)
    summary_df.to_csv(summary_file, sep="\t", index=False)

    # Clone composition
    plot_clone_donut(
        clone_counts=clone_counts,
        total_cells=len(cloned_df),
        output_dir=output_dir,
        top_n=args.top_n_clones
    )

    # Provenance files
    write_run_info(
        os.path.join(output_dir, "run_info.json"),
        args=args,
        input_file=args.input,
        n_records=len(cloned_df),
        n_clones=int(cloned_df["clone_id"].nunique()),
    )
    write_stage_accounting(os.path.join(output_dir, "stage_accounting.tsv"), accounting)

    # Console summary
    n_seq = len(result_df)
    n_clone = result_df["clone_id"].nunique()
    largest_clone = result_df["clone_size"].max()
    n_singleton = int((summary_df["clone_size"] == 1).sum())

    print("\n===== Summary =====")
    print(f"Sequences retained for clonotyping : {n_seq}")
    print(f"Clones identified                  : {n_clone}")
    print(f"  of which singletons              : {n_singleton}")
    print(f"Largest clone size                 : {largest_clone}")
    if n_singleton / max(n_clone, 1) > 0.9:
        print("NOTE: more than 90% of the clonotypes are singletons. That is expected for")
        print("      independent single-molecule amplicons, and it means the clone-size")
        print("      figures describe sampling, not clonal expansion.")
    print(f"Result file                        : {result_file}")
    print(f"Summary file                       : {summary_file}")
    print("Done.")


if __name__ == "__main__":
    sys.exit(main())
