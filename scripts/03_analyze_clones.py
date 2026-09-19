#!/usr/bin/env python3
"""Step 3 of SangerBCR-Flow: summary plots and heavy-chain clonotype calling.

Input : ../results/02_igblast_results.tsv   AIRR rearrangement TSV from step 2
Output: ../results/analysis_output/
          Heavy_V_Usage.png, Light_V_Usage.png
          SHM_Distribution.png
          CDR3_Length_Distribution.png
          Clone_Pie_Chart.png
          Heavy_Clonotypes.tsv, Clone_Summary.tsv

Clonotypes are called on the heavy chain only: records sharing the V gene, the J gene
and the junction length are linked when their junction nucleotide sequences differ by at
most MAX_DISTANCE (normalised Hamming distance), and linked records are merged
transitively (single linkage).

Run from this directory:  python3 03_analyze_clones.py
"""


import os
import re
import pandas as pd
import seaborn as sns
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


def normalize_cdr3_aa(x):
    """Return an uppercase CDR3 amino acid string, or None when the value is unusable.

    Keeps empty fields, so that a missing CDR3 never enters the length plot as "nan".
    """
    if pd.isna(x):
        return None
    s = str(x).strip().upper()
    if s in FAKE_NULL:
        return None
    if not re.fullmatch(r"[A-Z\*\-]+", s):
        return None
    return s


def normalize_junction_nt(x):
    """Return the junction nucleotide sequence in upper case, or None when unusable.

    Keeps empty fields, so that a missing junction never enters clonotyping as "nan".
    """
    if pd.isna(x):
        return None

    s = str(x).strip().upper()

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
        plt.title("Top 20 Heavy Chain V-Gene Usage")
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
        plt.title("Top 20 Light Chain V-Gene Usage")
        plt.xlabel("Count")
        plt.ylabel("V Gene")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "Light_V_Usage.png"), dpi=300)
        plt.close()


def plot_shm(df, output_dir):
    """Write the somatic hypermutation histogram for both chain types."""
    print("\nSHM analysis ...")

    plot_df = df.copy()
    plot_df["v_identity_num"] = pd.to_numeric(plot_df["v_identity"], errors="coerce")
    plot_df = plot_df.dropna(subset=["v_identity_num"]).copy()

    if plot_df.empty:
        print("No usable v_identity values; skipping the SHM plot.")
        return

    # IgBLAST reports v_identity as a percentage (0-100), not as a fraction
    plot_df["mutation_rate"] = 100 - plot_df["v_identity_num"]

    plt.figure(figsize=(8, 6))
    sns.histplot(data=plot_df, x="mutation_rate", hue="chain_type", kde=True, bins=30)
    plt.title("Somatic Hypermutation (SHM) Distribution")
    plt.xlabel("Mutation Rate (%) (100 - v_identity)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "SHM_Distribution.png"), dpi=300)
    plt.close()


def plot_cdr3_length(df, output_dir):
    """Write the CDR3 amino acid length distribution."""
    print("\nCDR3 length distribution ...")

    plot_df = df.copy()
    plot_df["cdr3_aa_clean"] = plot_df["cdr3_aa"].apply(normalize_cdr3_aa) if "cdr3_aa" in plot_df.columns else None

    if "cdr3_aa_clean" not in plot_df.columns:
        print("No cdr3_aa column found; skipping the CDR3 length plot.")
        return

    plot_df = plot_df.dropna(subset=["cdr3_aa_clean"]).copy()

    if plot_df.empty:
        print("No usable cdr3_aa values; skipping the CDR3 length plot.")
        return

    plot_df["cdr3_length"] = plot_df["cdr3_aa_clean"].str.len()

    plt.figure(figsize=(8, 6))
    sns.histplot(data=plot_df, x="cdr3_length", hue="chain_type", kde=True, element="step")
    plt.title("CDR3 Amino Acid Length Distribution")
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
    others_count = clone_counts.iloc[top_n:]["size"].sum() if len(clone_counts) > top_n else 0

    labels = [f"Clone {cid}" for cid in top_clones["clone_id"]]
    sizes = list(top_clones["size"])

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


def prepare_heavy_clonotyping_input(df):
    """Keep the productive heavy-chain records that can be clonotyped."""
    if "junction" not in df.columns:
        raise ValueError("Input file has no 'junction' column; clonotyping is not possible.")

    heavy_df = df[df["chain_type"] == "Heavy"].copy()

    heavy_df["junction_nt"] = heavy_df["junction"].apply(normalize_junction_nt)

    heavy_df = heavy_df[
        heavy_df["v_gene"].ne("")
        & heavy_df["j_gene"].ne("")
        & heavy_df["junction_nt"].notna()
    ].copy()

    heavy_df["junction_length"] = heavy_df["junction_nt"].str.len()

    if DROP_DUPLICATE_SEQUENCE_ID:
        before = len(heavy_df)
        heavy_df = heavy_df.drop_duplicates(subset=["sequence_id"]).copy()
        after = len(heavy_df)
        print(f"Dropped duplicated sequence_id: {before - after}")

    heavy_df = heavy_df.reset_index(drop=True)

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
    """
    df = df.copy().reset_index(drop=True)

    df["group_key"] = (
        df["v_gene"] + "|" +
        df["j_gene"] + "|" +
        df["junction_length"].astype(str)
    )

    df["clone_id"] = -1
    clone_counter = 1

    for _, group in df.groupby("group_key", sort=False):
        indices = group.index.tolist()

        if len(indices) == 1:
            df.at[indices[0], "clone_id"] = clone_counter
            clone_counter += 1
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

        root_to_clone = {}
        for idx in indices:
            root = uf.find(idx)
            if root not in root_to_clone:
                root_to_clone[root] = clone_counter
                clone_counter += 1
            df.at[idx, "clone_id"] = root_to_clone[root]

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
        .rename(columns={"size": "size"})
        .sort_values(["size", "clone_id"], ascending=[False, True])
        .reset_index(drop=True)
    )

    clone_size_map = clone_counts.set_index("clone_id")["size"]

    result_df = df.assign(
        clone_size=df["clone_id"].map(clone_size_map)
    ).copy()

    summary_df = (
        result_df
        .sort_values(["clone_size", "clone_id"], ascending=[False, True])
        .groupby("clone_id", as_index=False)
        .first()[[
            "clone_id",
            "clone_size",
            "v_gene",
            "j_gene",
            "junction_length",
            "junction_nt",
            "sequence_id"
        ]]
        .rename(columns={
            "junction_nt": "representative_junction_nt",
            "sequence_id": "representative_sequence_id"
        })
        .sort_values(["clone_size", "clone_id"], ascending=[False, True])
        .reset_index(drop=True)
        .copy()
    )

    return result_df, summary_df, clone_counts


# =========================
# 9. Main
# =========================
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load and filter the productive records
    raw_df = load_raw_data(INPUT_FILE)
    print(f"Total productive sequences: {len(raw_df)}")
    print(raw_df["chain_type"].value_counts())

    # Descriptive plots, on all productive records
    plot_gene_usage(raw_df, OUTPUT_DIR)
    plot_shm(raw_df, OUTPUT_DIR)
    plot_cdr3_length(raw_df, OUTPUT_DIR)

    # Heavy-chain input for clonotyping
    heavy_df = prepare_heavy_clonotyping_input(raw_df)

    if len(heavy_df) == 0:
        raise ValueError("No heavy-chain records left after filtering; nothing to clonotype.")

    print("\nHeavy-chain clonotyping ...")
    cloned_df = run_heavy_clonotyping(heavy_df, max_distance=MAX_DISTANCE)

    print(f"Total Clones Identified: {cloned_df['clone_id'].nunique()}")

    # Tables
    result_df, summary_df, clone_counts = make_output_tables(cloned_df)

    result_file = os.path.join(OUTPUT_DIR, "Heavy_Clonotypes.tsv")
    summary_file = os.path.join(OUTPUT_DIR, "Clone_Summary.tsv")

    result_df.to_csv(result_file, sep="\t", index=False)
    summary_df.to_csv(summary_file, sep="\t", index=False)

    # Clone composition
    plot_clone_donut(
        clone_counts=clone_counts,
        total_cells=len(cloned_df),
        output_dir=OUTPUT_DIR,
        top_n=TOP_N_CLONES_IN_PIE
    )

    # Console summary
    n_seq = len(result_df)
    n_clone = result_df["clone_id"].nunique()
    largest_clone = result_df["clone_size"].max()

    print("\n===== Summary =====")
    print(f"Sequences retained for clonotyping : {n_seq}")
    print(f"Clones identified                  : {n_clone}")
    print(f"Largest clone size                 : {largest_clone}")
    print(f"Result file                        : {result_file}")
    print(f"Summary file                       : {summary_file}")
    print("Done.")


if __name__ == "__main__":
    main()
