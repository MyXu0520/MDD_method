#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parameter_sweep_full_with_metrics.py

"""

import os
import glob
import math
import random
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

from sklearn.mixture import GaussianMixture
from scipy.signal import argrelextrema
from scipy.stats import ttest_ind
from sklearn.metrics import jaccard_score

# ------------------------
# USER CONFIGURATION
# ------------------------
DATA_DIR = "/bio2/MDD+DMR/data/"        # directory containing ESCC_*.bed files
OUT_DIR = "/bio2/MDD_result/0/code/result/sweep/"

BIN_SIZES = [2000, 5000, 10000, 20000, 30000, 50000]
THRESHOLDS = np.arange(0.2, 0.61, 0.05)

MAX_SAMPLES_PER_BIN_SIZE = None   # set to int to limit samples for speed, or None for all
N_BOOTSTRAP = 50                  # bootstrap iterations for Jaccard stability
BOOTSTRAP_SAMPLE_FRAC = 0.8
N_PERMUTATIONS = 200              # permutation iterations for low_frac significance

VERBOSE = True
RANDOM_SEED = 2025

# ------------------------
# Safety / seeds
# ------------------------
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

# ------------------------
# Helpers: create bins
# ------------------------
def create_bins(bin_size):
    # Use chr1-12 lengths (modify if you need full genome)
    genome = {
        "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
        "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
        "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309
    }
    rows = []
    for chrom, length in genome.items():
        # create bins starting at 0
        for s in range(0, length, bin_size):
            rows.append([chrom, s, min(s + bin_size, length)])
    return pd.DataFrame(rows, columns=["chr", "start", "end"])

# ------------------------
# Map a sample to full aligned bin vector
# ------------------------
def map_sample_aligned(file, bins_df, bin_size):
    """
    Returns numpy array aligned to bins_df rows (same order), with avg methylation per bin or np.nan.
    Expects file with columns: chr, pos, meth, mc, uc (tab separated).
    """
    try:
        df = pd.read_csv(file, sep="\t", header=None, names=["chr", "pos", "meth", "mc", "uc"], dtype={"chr": str})
    except Exception as e:
        if VERBOSE: print(f"  ERROR reading {file}: {e}")
        return np.array([np.nan] * len(bins_df))

    # compute bin index (start pos)
    df["bin"] = (df["pos"] // bin_size) * bin_size
    # group by chr and bin
    bin_meth = df.groupby(["chr", "bin"])["meth"].mean()
    # create arrays aligned to bins_df
    avg_list = []
    # optimization: create dict for quick lookup
    bin_dict = { (chrom, int(bin_start)): float(m) for (chrom, bin_start), m in bin_meth.items() }
    for idx, row in bins_df.iterrows():
        key = (row["chr"], int(row["start"]))
        avg_list.append(bin_dict.get(key, np.nan))
    return np.array(avg_list, dtype=float)

# ------------------------
# GMM valley finder
# ------------------------
def find_gmm_valley(values, n_components=2, random_state=RANDOM_SEED, nbins=200):
    vals = np.asarray(values)
    vals = vals[~np.isnan(vals)]
    if len(vals) < 100:
        return np.nan
    # histogram for smoothing
    counts, edges = np.histogram(vals, bins=nbins, range=(0,1), density=True)
    mids = (edges[:-1] + edges[1:]) / 2.0
    try:
        gmm = GaussianMixture(n_components=n_components, random_state=random_state)
        gmm.fit(vals.reshape(-1,1))
        logprob = gmm.score_samples(mids.reshape(-1,1))
        pdf = np.exp(logprob)
        minima_idx = argrelextrema(pdf, np.less)[0]
        if len(minima_idx) == 0:
            means = np.sort(gmm.means_.ravel())
            return float(np.mean(means))
        valley_mids = mids[minima_idx]
        valley = float(valley_mids[np.argmin(pdf[minima_idx])])
        return valley
    except Exception:
        # fallback: use histogram valley
        try:
            minima_idx = argrelextrema(counts, np.less)[0]
            if len(minima_idx)==0:
                return np.nan
            valley = float(mids[minima_idx[np.argmin(counts[minima_idx])]])
            return valley
        except:
            return np.nan

# ------------------------
# Cohen's d
# ------------------------
def cohens_d(x, y):
    x = np.asarray(x); y = np.asarray(y)
    if len(x) < 2 or len(y) < 2:
        return np.nan
    nx, ny = len(x), len(y)
    sx = np.nanstd(x, ddof=1)
    sy = np.nanstd(y, ddof=1)
    pooled = np.sqrt(((nx-1)*(sx**2) + (ny-1)*(sy**2)) / (nx+ny-2))
    if pooled == 0:
        return np.nan
    return (np.nanmean(x) - np.nanmean(y)) / pooled

# ------------------------
# bootstrap jaccard stability
# ------------------------
def bootstrap_jaccard(all_bins_by_sample, threshold, n_boot=N_BOOTSTRAP, sample_frac=BOOTSTRAP_SAMPLE_FRAC):
    """
    all_bins_by_sample: list of 1D numpy arrays (aligned to same bins), values = avg meth or nan
    returns average pairwise jaccard across bootstrap resamples
    """
    n_samples = len(all_bins_by_sample)
    if n_samples < 2:
        return np.nan
    jaccs = []
    # prepare arrays
    arrs = [np.array(a) for a in all_bins_by_sample]
    for _ in range(n_boot):
        chosen_idx = random.sample(range(n_samples), max(2, int(math.ceil(n_samples * sample_frac))))
        masks = []
        for i in chosen_idx:
            arr = arrs[i]
            mask = (~np.isnan(arr)) & (arr < threshold)
            masks.append(mask.astype(int))
        # compute pairwise jaccard among masks
        for i in range(len(masks)):
            for j in range(i+1, len(masks)):
                a = masks[i]; b = masks[j]
                # need at least one positive in union to compute jaccard
                if (a.sum() + b.sum()) == 0:
                    continue
                try:
                    j = jaccard_score(a, b)
                    jaccs.append(j)
                except Exception:
                    continue
    return float(np.nanmean(jaccs)) if len(jaccs) > 0 else np.nan

# ------------------------
# permutation test for low_frac
# ------------------------
def permutation_test_lowfrac(flat_all_bin_avg_meth, thr, n_perm=N_PERMUTATIONS):
    arr = np.array(flat_all_bin_avg_meth)
    arr = arr[~np.isnan(arr)]
    if len(arr) == 0:
        return np.nan, np.nan
    obs = float(np.mean(arr < thr))
    perms = []
    for _ in range(n_perm):
        perm = np.random.permutation(arr)
        perms.append(np.mean(perm < thr))
    perms = np.array(perms)
    p_value = (np.sum(perms >= obs) + 1) / (n_perm + 1)
    return obs, float(p_value)

# ------------------------
# compute metrics for given bin_size
# ------------------------
def compute_metrics_for_bin_size(bin_size, samples, max_samples=None, thresholds=THRESHOLDS):
    if VERBOSE: print(f"\n=== Bin size {bin_size} ===")
    bins_df = create_bins(bin_size)
    # load per-sample aligned arrays
    sample_paths = samples.copy()
    if max_samples is not None:
        sample_paths = sample_paths[:max_samples]
    if len(sample_paths) == 0:
        return None

    all_by_sample = []
    for p in tqdm(sample_paths, desc=f"mapping samples (bin={bin_size})", leave=False):
        arr = map_sample_aligned(p, bins_df, bin_size)
        all_by_sample.append(arr)

    # convert to 2D array (n_samples x n_bins)
    arr_stack = np.vstack(all_by_sample)  # shape (n_samples, n_bins)
    n_samples, n_bins = arr_stack.shape
    if VERBOSE: print(f"  samples={n_samples}, bins={n_bins}")

    # compute per-bin mean across samples
    mean_per_bin = np.nanmean(arr_stack, axis=0)

    # flatten all non-nan values for GMM and global low_frac/permutation
    flat_all = mean_per_bin[~np.isnan(mean_per_bin)]

    # find gmm valley for this bin size
    gmm_valley = find_gmm_valley(flat_all, n_components=2, nbins=300)

    # metrics per threshold
    rows = []
    for thr in thresholds:
        # low_frac computed across all samples/bins combined: proportion of bins (mean_per_bin) < thr
        low_frac = float(np.sum((~np.isnan(mean_per_bin)) & (mean_per_bin < thr)) / np.sum(~np.isnan(mean_per_bin))) if np.sum(~np.isnan(mean_per_bin))>0 else np.nan

        # bootstrap jaccard (uses per-sample masks)
        stab = bootstrap_jaccard(all_by_sample, thr, n_boot=N_BOOTSTRAP, sample_frac=BOOTSTRAP_SAMPLE_FRAC)

        # effect size: compare mean values of bins labelled MDD vs non-MDD (by mean_per_bin)
        pmask = (~np.isnan(mean_per_bin)) & (mean_per_bin < thr)
        npmask = (~np.isnan(mean_per_bin)) & (mean_per_bin >= thr)
        if pmask.sum() >= 5 and npmask.sum() >= 5:
            d = cohens_d(mean_per_bin[pmask], mean_per_bin[npmask])
            # t-test as complementary stat
            try:
                tstat, pval = ttest_ind(mean_per_bin[pmask], mean_per_bin[npmask], equal_var=False, nan_policy='omit')
            except Exception:
                tstat, pval = np.nan, np.nan
        else:
            d, tstat, pval = np.nan, np.nan, np.nan

        # permutation test for low_frac (on flattened mean_per_bin)
        obs_low_frac, perm_p = permutation_test_lowfrac(flat_all, thr, n_perm=N_PERMUTATIONS)

        rows.append({
            "bin_size": int(bin_size),
            "threshold": float(thr),
            "low_frac": float(low_frac),
            "bootstrap_jaccard": float(stab) if not np.isnan(stab) else np.nan,
            "cohens_d": float(d) if not np.isnan(d) else np.nan,
            "t_stat": float(tstat) if not np.isnan(tstat) else np.nan,
            "t_pval": float(pval) if not np.isnan(pval) else np.nan,
            "perm_low_frac": float(obs_low_frac) if not np.isnan(obs_low_frac) else np.nan,
            "perm_pval": float(perm_p) if not np.isnan(perm_p) else np.nan,
            "gmm_valley": float(gmm_valley) if not np.isnan(gmm_valley) else np.nan,
            "n_samples": int(n_samples),
            "n_bins_non_na": int(np.sum(~np.isnan(mean_per_bin))),
            "mean_cpg_per_bin_est": float(np.nanmean(np.sum(~np.isnan(arr_stack), axis=0)))  # approx count of samples with data per bin
        })
    return pd.DataFrame(rows), bins_df, all_by_sample

# ------------------------
# Main pipeline
# ------------------------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    samples = sorted(glob.glob(os.path.join(DATA_DIR, "ESCC_*.bed")))
    if len(samples) == 0:
        print("No samples found - check DATA_DIR and file patterns.")
        return

    # optionally limit samples for testing
    if MAX_SAMPLES_PER_BIN_SIZE is not None:
        if VERBOSE: print(f"Limiting to first {MAX_SAMPLES_PER_BIN_SIZE} samples for each bin size.")
    else:
        if VERBOSE: print(f"Using {len(samples)} samples for each bin size.")

    all_metrics = []
    summary_by_bin = []

    # store per-bin_size all_by_sample for optional downstream analyses (be mindful of memory)
    per_bin_all_by_sample = {}

    for bin_size in BIN_SIZES:
        df_metrics, bins_df, all_by_sample = compute_metrics_for_bin_size(
            bin_size, samples, max_samples=MAX_SAMPLES_PER_BIN_SIZE, thresholds=THRESHOLDS)
        if df_metrics is None:
            continue
        # save
        all_metrics.append(df_metrics)
        per_bin_all_by_sample[bin_size] = {
            "bins_df": bins_df,
            "all_by_sample": all_by_sample
        }
        # summary per bin size: pick threshold that maximizes stability*effect (simple composite)
        # compute simple composite: stability_mean * effect_mean (ignoring nans)
        stab_mean = df_metrics["bootstrap_jaccard"].mean()
        eff_mean = df_metrics["cohens_d"].mean()
        summary_by_bin.append({
            "bin_size": int(bin_size),
            "gmm_valley_overall": float(df_metrics["gmm_valley"].iloc[0]),
            "stab_mean": float(stab_mean) if not np.isnan(stab_mean) else np.nan,
            "effect_mean": float(eff_mean) if not np.isnan(eff_mean) else np.nan,
            "n_bins_non_na": int(df_metrics["n_bins_non_na"].max()),
            "n_samples": int(df_metrics["n_samples"].max())
        })

    if len(all_metrics) == 0:
        print("No metrics computed.")
        return

    df_all = pd.concat(all_metrics, ignore_index=True)
    df_all.to_csv(os.path.join(OUT_DIR, "metrics_per_bin_threshold.tsv"), sep="\t", index=False)
    pd.concat([pd.DataFrame(summary_by_bin)], ignore_index=True).to_csv(os.path.join(OUT_DIR, "summary_by_bin_size.tsv"), sep="\t", index=False)
    # also save parameter sweep (bin_size vs threshold vs low_frac)
    df_all[["bin_size", "threshold", "low_frac"]].to_csv(os.path.join(OUT_DIR, "parameter_sweep.tsv"), sep="\t", index=False)

    # ============================
    # PLOTTING - compile PDF
    # ============================
    pdf_path = os.path.join(OUT_DIR, "parameter_sweep_plots_with_metrics.pdf")
    pdf = PdfPages(pdf_path)

    sns.set(style="whitegrid")

    # Plot A: threshold sensitivity curves (low_frac) per bin_size
    plt.figure(figsize=(10,6))
    for b in BIN_SIZES:
        sub = df_all[df_all["bin_size"] == b]
        if not sub.empty:
            plt.plot(sub["threshold"], sub["low_frac"], marker="o", label=f"{b}bp", linewidth=2)
    plt.legend(title="Bin Size")
    plt.xlabel("Threshold")
    plt.ylabel("Fraction low-methyl bins")
    plt.title("Threshold sensitivity: low_frac by bin size")
    plt.tight_layout()
    pdf.savefig(); plt.close()

    # Plot B: bin size effect (low_frac) per threshold
    plt.figure(figsize=(10,6))
    for thr in THRESHOLDS:
        sub = df_all[df_all["threshold"] == thr]
        if not sub.empty:
            plt.plot(sub["bin_size"], sub["low_frac"], marker="s", label=f"thr={thr:.2f}", linewidth=2)
    plt.xscale("log")
    plt.xticks(BIN_SIZES, BIN_SIZES)
    plt.xlabel("Bin size (bp)")
    plt.ylabel("Fraction low-methyl bins")
    plt.title("Effect of bin size on low_frac")
    plt.legend(title="Threshold", bbox_to_anchor=(1.05,1), loc='upper left')
    plt.tight_layout()
    pdf.savefig(); plt.close()

    # Heatmap (bin_size x threshold) for low_frac
    heat = df_all.pivot(index="bin_size", columns="threshold", values="low_frac")
    plt.figure(figsize=(8,6))
    sns.heatmap(heat, cmap="viridis", annot=True, fmt=".2f")
    plt.title("Heatmap: low_frac (bin_size x threshold)")
    plt.tight_layout()
    pdf.savefig(); plt.close()

    # Heatmap: stability (bootstrap_jaccard)
    heat_stab = df_all.pivot(index="bin_size", columns="threshold", values="bootstrap_jaccard")
    plt.figure(figsize=(8,6))
    sns.heatmap(heat_stab, cmap="magma", annot=True, fmt=".2f")
    plt.title("Heatmap: bootstrap_jaccard (stability)")
    plt.tight_layout()
    pdf.savefig(); plt.close()

    # Heatmap: effect (Cohen's d)
    heat_eff = df_all.pivot(index="bin_size", columns="threshold", values="cohens_d")
    plt.figure(figsize=(8,6))
    sns.heatmap(heat_eff, cmap="coolwarm", center=0, annot=True, fmt=".2f")
    plt.title("Heatmap: Cohen's d (MDD vs non-MDD separation)")
    plt.tight_layout()
    pdf.savefig(); plt.close()

    # Plot: stability vs threshold for each bin_size
    plt.figure(figsize=(10,6))
    for b in BIN_SIZES:
        sub = df_all[df_all["bin_size"] == b]
        if not sub.empty:
            plt.plot(sub["threshold"], sub["bootstrap_jaccard"], marker="o", label=f"{b}bp")
    plt.xlabel("Threshold"); plt.ylabel("Bootstrap Jaccard (stability)")
    plt.title("Stability across thresholds")
    plt.legend(title="Bin Size", bbox_to_anchor=(1.05,1), loc='upper left')
    plt.tight_layout()
    pdf.savefig(); plt.close()

    # Plot: Cohen's d vs threshold for each bin_size
    plt.figure(figsize=(10,6))
    for b in BIN_SIZES:
        sub = df_all[df_all["bin_size"] == b]
        if not sub.empty:
            plt.plot(sub["threshold"], sub["cohens_d"], marker="o", label=f"{b}bp")
    plt.xlabel("Threshold"); plt.ylabel("Cohen's d")
    plt.title("Effect size across thresholds")
    plt.legend(title="Bin Size", bbox_to_anchor=(1.05,1), loc='upper left')
    plt.tight_layout()
    pdf.savefig(); plt.close()

    # Boxplot: distribution of low_frac per bin_size
    plt.figure(figsize=(8,6))
    sns.boxplot(data=df_all, x="bin_size", y="low_frac")
    plt.xlabel("Bin size (bp)"); plt.ylabel("low_frac")
    plt.title("Distribution of low_frac per bin_size")
    plt.tight_layout()
    pdf.savefig(); plt.close()

    # Boxplot: distribution of bootstrap_jaccard per bin_size
    plt.figure(figsize=(8,6))
    sns.boxplot(data=df_all, x="bin_size", y="bootstrap_jaccard")
    plt.xlabel("Bin size (bp)"); plt.ylabel("bootstrap_jaccard")
    plt.title("Distribution of stability (bootstrap_jaccard) per bin_size")
    plt.tight_layout()
    pdf.savefig(); plt.close()

    # Scatter: effect_mean vs stability_mean from summary_by_bin_size
    summary_df = pd.read_csv(os.path.join(OUT_DIR, "summary_by_bin_size.tsv"), sep="\t")
    plt.figure(figsize=(8,6))
    sns.scatterplot(data=summary_df, x="stab_mean", y="effect_mean", size="n_bins_non_na", hue="bin_size", palette="tab10", legend="brief", sizes=(50,300))
    plt.xlabel("Mean stability (bootstrap_jaccard)"); plt.ylabel("Mean effect (Cohen's d)")
    plt.title("Bin size comparison: stability vs effect")
    plt.tight_layout()
    pdf.savefig(); plt.close()

    # Add GMM valley per bin size plot (valley vs bin_size)
    valley_vals = summary_df[["bin_size","gmm_valley_overall"]].dropna()
    if not valley_vals.empty:
        plt.figure(figsize=(8,5))
        plt.plot(valley_vals["bin_size"], valley_vals["gmm_valley_overall"], marker="o")
        plt.xscale("log")
        plt.xticks(BIN_SIZES, BIN_SIZES)
        plt.xlabel("Bin size (bp)"); plt.ylabel("GMM valley (beta)")
        plt.title("GMM valley (bimodal minimum) by bin size")
        plt.tight_layout()
        pdf.savefig(); plt.close()

    # Save a few CSV summary files too
    df_all.to_csv(os.path.join(OUT_DIR, "metrics_per_bin_threshold_full.tsv"), sep="\t", index=False)
    summary_df.to_csv(os.path.join(OUT_DIR, "summary_by_bin_size_full.tsv"), sep="\t", index=False)

    pdf.close()
    print(f"All outputs saved into: {OUT_DIR}")
    print(f"PDF: {pdf_path}")
    print("Done.")

# ------------------------
# run
# ------------------------
if __name__ == "__main__":
    main()
