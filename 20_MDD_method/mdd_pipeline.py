#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# ==============================
# CONFIG
# ==============================
DATA_DIR = "/bio2/PMD+DMR/data/"
OUT_DIR = "/bio2/PMD_result/0/code/result/"
BIN_SIZES = [10000]  # main pipeline uses 10kb bins


# ==============================
# STEP 1. BUILD BINS
# ==============================
GENOME = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468
}


def make_bins(bin_size, outdir):
    os.makedirs(outdir, exist_ok=True)
    rows = []
    for chrom, length in GENOME.items():
        for start in range(0, length, bin_size):
            end = min(start + bin_size, length)
            rows.append([chrom, start, end])
    df = pd.DataFrame(rows, columns=["chr", "start", "end"])
    df.to_csv(f"{outdir}/hg38_bins_{bin_size}.bed", sep="\t", index=False)


# ==============================
# STEP 2. MAP CpG METH TO BINS
# ==============================
def map_to_bins(sample_file, bin_df, bin_size):
    df = pd.read_csv(sample_file, sep="\t", header=None,
                     names=["chr", "pos", "meth", "mc", "uc"])
    df["bin"] = df["pos"] // bin_size
    bin_df["bin"] = bin_df["start"] // bin_size
    merged = df.merge(bin_df, on=["chr", "bin"])
    out = merged.groupby(["chr", "start", "end"])["meth"].mean().reset_index()
    return out


# ==============================
# STEP 3. FIND BIMODAL THRESHOLD
# ==============================
def find_threshold(values):
    X = np.array(values).reshape(-1, 1)
    gm = GaussianMixture(n_components=2).fit(X)
    m = np.sort(gm.means_.flatten())
    thr = m.mean()
    return float(thr)


# ==============================
# STEP 4. MERGE LOW-METH BINS → MDDs
# ==============================
def merge_regions(df):
    merged = []
    current = None
    for _, r in df.iterrows():
        if current is None:
            current = r
        else:
            if r["start"] == current["end"]:
                current["end"] = r["end"]
            else:
                merged.append(current.copy())
                current = r
    if current is not None:
        merged.append(current)
    return pd.DataFrame(merged)


# ==============================
# STEP 5. FIND EAC / ESCC SPECIFIC
# ==============================
def load_regions(files):
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        dfs.append(df)
    return pd.concat(dfs)


def intersect(a, b):
    results = []
    for _, r in a.iterrows():
        overlap = b[
            (b["chr"] == r["chr"]) &
            (b["start"] < r["end"]) &
            (b["end"] > r["start"])
        ]
        if len(overlap) > 0:
            results.append(r)
    return pd.DataFrame(results)


# ==============================
# MAIN PIPELINE
# ==============================
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1. Build bins (10kb)
    for bs in BIN_SIZES:
        make_bins(bs, f"{OUT_DIR}/bins")

    # 2. Map samples to bins
    samples = glob.glob(f"{DATA_DIR}/*.bed")
    for bin_size in BIN_SIZES:
        bin_df = pd.read_csv(f"{OUT_DIR}/bins/hg38_bins_{bin_size}.bed", sep="\t")
        os.makedirs(f"{OUT_DIR}/bin_meth", exist_ok=True)

        for f in tqdm(samples, desc=f"Mapping to {bin_size}bp bins"):
            name = os.path.basename(f).replace(".bed", "")
            out = map_to_bins(f, bin_df.copy(), bin_size)
            out.to_csv(f"{OUT_DIR}/bin_meth/{name}_bin{bin_size}.tsv", sep="\t", index=False)

    # 3. Find threshold from all samples' methylation distribution
    for bin_size in BIN_SIZES:
        vals = []
        for f in glob.glob(f"{OUT_DIR}/bin_meth/*bin{bin_size}.tsv"):
            df = pd.read_csv(f, sep="\t")
            vals.extend(df["meth"].dropna().tolist())
        thr = find_threshold(vals)
        with open(f"{OUT_DIR}/threshold_bin{bin_size}.txt", "w") as fh:
            fh.write(str(thr))

    # 4. Call MDDs
    for f in glob.glob(f"{OUT_DIR}/bin_meth/*.tsv"):
        name = os.path.basename(f).split("_bin")[0]
        bin_size = int(os.path.basename(f).split("bin")[1].split(".")[0])
        thr = float(open(f"{OUT_DIR}/threshold_bin{bin_size}.txt").read())
        df = pd.read_csv(f, sep="\t")
        low = df[df["meth"] < thr]
        merged = low.groupby("chr").apply(merge_regions).reset_index(drop=True)
        os.makedirs(f"{OUT_DIR}/MDDs", exist_ok=True)
        merged.to_csv(f"{OUT_DIR}/MDDs/{name}_MDDs.bed", sep="\t", index=False)

    # 5. Compare ESCC vs EAC
    indir = f"{OUT_DIR}/MDDs/"
    ESCC = [f for f in glob.glob(f"{indir}/ESCC*.bed")]
    EAC  = [f for f in glob.glob(f"{indir}/EAC*.bed")] + \
           [f for f in glob.glob(f"{indir}/GEJ*.bed")]

    escc_all = load_regions(ESCC)
    eac_all  = load_regions(EAC)

    shared = intersect(escc_all, eac_all)
    escc_specific = escc_all[~escc_all.index.isin(shared.index)]
    eac_specific =  eac_all[~eac_all.index.isin(shared.index)]

    os.makedirs(f"{OUT_DIR}/final_MDDs", exist_ok=True)
    shared.to_csv(f"{OUT_DIR}/final_MDDs/Shared_MDDs.bed", sep="\t", index=False)
    escc_specific.to_csv(f"{OUT_DIR}/final_MDDs/ESCC_specific_MDDs.bed", sep="\t", index=False)
    eac_specific.to_csv(f"{OUT_DIR}/final_MDDs/EAC_specific_MDDs.bed", sep="\t", index=False)

    print("MDD pipeline completed!")


if __name__ == "__main__":
    main()
