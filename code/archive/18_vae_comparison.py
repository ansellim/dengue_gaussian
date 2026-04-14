#!/usr/bin/env python3
"""
18_vae_comparison.py

VAE Latent Space Comparison with GP Residual (f_residual) from Bayesian
renewal equation model for dengue Rt in Singapore.

Train a sliding-window VAE on log-transformed weekly case counts, then compare
the learned latent representation with the GP residual and serotype dynamics.

Prerequisites:
  Rscript code/17_export_for_vae.R   # exports GP draws and case data

Usage:
  python code/18_vae_comparison.py          # from project root
  cd code && python 18_vae_comparison.py    # from code/

Outputs (results/figures/):
  vae_latent_space.png            - 2x2 panel: PCA of latent encodings
  vae_latent_trajectory.png       - latent PCA over time vs f_residual
  vae_reconstruction_error.png    - anomaly detection via reconstruction error
  vae_uncertainty_comparison.png  - GP vs VAE uncertainty over time
  vae_training_loss.png           - training curves

Outputs (results/):
  vae_comparison_summary.csv      - correlation & anomaly statistics
"""

import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
from sklearn.decomposition import PCA
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

warnings.filterwarnings("ignore", category=FutureWarning)

# =============================================================================
# Resolve project root
# =============================================================================

_this_dir = Path(__file__).resolve().parent
# Walk up until we find a parent containing both code/ and results/
_p = _this_dir
for _ in range(4):
    if (_p / "code").exists() and (_p / "results").exists():
        PROJECT_ROOT = _p
        break
    _p = _p.parent
else:
    PROJECT_ROOT = _this_dir

RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
VAE_EXPORT = RESULTS_DIR / "vae_export"

FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Reproducibility
SEED = 42
np.random.seed(SEED)
torch.manual_seed(SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)

# =============================================================================
# SECTION A: Data Loading
# =============================================================================

print("=" * 70)
print("VAE LATENT SPACE COMPARISON WITH GP RESIDUAL")
print("=" * 70)

# Check that export files exist
required_files = [
    VAE_EXPORT / "cases_weekly.csv",
    VAE_EXPORT / "f_residual_summary.csv",
    VAE_EXPORT / "serotype_info.csv",
]
missing = [str(f) for f in required_files if not f.exists()]
if missing:
    print("\nERROR: Missing required files:")
    for m in missing:
        print(f"  {m}")
    print("\nPlease run 17_export_for_vae.R first:")
    print("  Rscript code/17_export_for_vae.R")
    sys.exit(1)

print("\nA. Loading data...")

cases_df = pd.read_csv(VAE_EXPORT / "cases_weekly.csv", parse_dates=["date"])
f_resid_df = pd.read_csv(VAE_EXPORT / "f_residual_summary.csv", parse_dates=["date"])
serotype_df = pd.read_csv(VAE_EXPORT / "serotype_info.csv", parse_dates=["year_month"])

# Drop duplicate dates from f_residual (pre-existing artefact of postprocess
# regeneration). Keep the first occurrence per date.
f_resid_df = f_resid_df.drop_duplicates(subset="date", keep="first").reset_index(drop=True)

# Serotype switch timing (original file)
switch_file = RESULTS_DIR / "serotype_switch_timing.csv"
if switch_file.exists():
    switch_df = pd.read_csv(switch_file, parse_dates=["switch_date"])
else:
    switch_df = pd.read_csv(VAE_EXPORT / "serotype_switch_timing.csv",
                            parse_dates=["switch_date"])

switch_dates = switch_df["switch_date"].values

# -------------------------------------------------------------------------
# Restrict analysis to 2013-01-01 onwards to match serotype data coverage.
# The serotype surveillance data starts January 2013; running the VAE on the
# 2012 weeks (which have no serotype labels) produced "Unknown" points in the
# latent-space serotype panel. Period-matching the inputs eliminates that
# confusion and makes the sanity check comparable across panels.
# -------------------------------------------------------------------------
ANALYSIS_START = pd.Timestamp("2013-01-01")
cases_df = cases_df[cases_df["date"] >= ANALYSIS_START].reset_index(drop=True)
f_resid_df = f_resid_df[f_resid_df["date"] >= ANALYSIS_START].reset_index(drop=True)

N_total = len(cases_df)
N_model = len(f_resid_df)

print(f"  Analysis window: {ANALYSIS_START.date()} onwards (period-matched to serotype data)")
print(f"  Case weeks in window: {N_total}")
print(f"  f_residual weeks in window: {N_model}")
print(f"  Cases range: {cases_df['date'].min().date()} to {cases_df['date'].max().date()}")
print(f"  Serotype switches: {len(switch_dates)}")

# =============================================================================
# SECTION B: Sliding Window Dataset
# =============================================================================

print("\nB. Building sliding window dataset...")

W = 26       # window width in weeks
STRIDE = 1

# Log-transform cases
log_cases = np.log1p(cases_df["cases"].values).astype(np.float32)

# Min-max normalize to [0, 1]
log_min = log_cases.min()
log_max = log_cases.max()
norm_cases = (log_cases - log_min) / (log_max - log_min + 1e-8)

# Create sliding windows
windows = []
window_start_indices = []
for i in range(0, N_total - W + 1, STRIDE):
    windows.append(norm_cases[i : i + W])
    window_start_indices.append(i)

windows = np.array(windows, dtype=np.float32)
window_start_indices = np.array(window_start_indices)
N_windows = len(windows)

# Midpoint index (in the full time series) for each window
window_mid_indices = window_start_indices + W // 2  # week i+13

# Dates for each window midpoint
all_dates = cases_df["date"].values
window_mid_dates = all_dates[window_mid_indices]

print(f"  Window size: {W}, stride: {STRIDE}")
print(f"  Number of windows: {N_windows}")

# 80/20 chronological split
split_idx = int(0.8 * N_windows)
X_train = torch.tensor(windows[:split_idx])
X_val = torch.tensor(windows[split_idx:])

train_loader = DataLoader(TensorDataset(X_train), batch_size=32, shuffle=True)
val_loader = DataLoader(TensorDataset(X_val), batch_size=64, shuffle=False)

print(f"  Train: {len(X_train)} windows, Val: {len(X_val)} windows")

# =============================================================================
# SECTION C: VAE Architecture
# =============================================================================

print("\nC. Setting up VAE...")

INPUT_DIM = W   # 26
LATENT_DIM = 4
HIDDEN = 32


class Encoder(nn.Module):
    def __init__(self):
        super().__init__()
        self.shared = nn.Sequential(
            nn.Linear(INPUT_DIM, HIDDEN),
            nn.ReLU(),
            nn.Linear(HIDDEN, HIDDEN // 2),
            nn.ReLU(),
        )
        self.fc_mu = nn.Linear(HIDDEN // 2, LATENT_DIM)
        self.fc_logvar = nn.Linear(HIDDEN // 2, LATENT_DIM)

    def forward(self, x):
        h = self.shared(x)
        return self.fc_mu(h), self.fc_logvar(h)


class VAE(nn.Module):
    def __init__(self):
        super().__init__()
        self.encoder = Encoder()
        self.decoder = nn.Sequential(
            nn.Linear(LATENT_DIM, HIDDEN // 2),
            nn.ReLU(),
            nn.Linear(HIDDEN // 2, HIDDEN),
            nn.ReLU(),
            nn.Linear(HIDDEN, INPUT_DIM),
            nn.Sigmoid(),
        )

    def forward(self, x):
        mu, logvar = self.encoder(x)
        std = torch.exp(0.5 * logvar)
        z = mu + std * torch.randn_like(std)
        x_hat = self.decoder(z)
        return x_hat, mu, logvar

    def encode_mu(self, x):
        """Deterministic encoding (posterior mean)."""
        mu, _ = self.encoder(x)
        return mu


def elbo_loss(x_hat, x, mu, logvar, beta=1.0):
    recon = F.mse_loss(x_hat, x, reduction="sum") / x.size(0)
    kl = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp()) / x.size(0)
    return recon + beta * kl, recon.item(), kl.item()


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = VAE().to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

print(f"  Device: {device}")
print(f"  Architecture: {INPUT_DIM} -> {HIDDEN} -> {HIDDEN // 2} -> {LATENT_DIM} (latent)")
print(f"  Parameters: {sum(p.numel() for p in model.parameters()):,}")

# =============================================================================
# SECTION D: Training
# =============================================================================

print("\nD. Training VAE...")

N_EPOCHS = 200
BETA_WARMUP = 50

history = {"epoch": [], "train_loss": [], "train_recon": [], "train_kl": [],
           "val_loss": [], "val_recon": [], "val_kl": []}

for epoch in range(1, N_EPOCHS + 1):
    # Beta annealing: linearly from 0.01 to 1.0 over first BETA_WARMUP epochs
    if epoch <= BETA_WARMUP:
        beta = 0.01 + (1.0 - 0.01) * (epoch - 1) / (BETA_WARMUP - 1)
    else:
        beta = 1.0

    # --- Train ---
    model.train()
    train_loss_sum, train_recon_sum, train_kl_sum, n_train = 0, 0, 0, 0
    for (batch_x,) in train_loader:
        batch_x = batch_x.to(device)
        x_hat, mu, logvar = model(batch_x)
        loss, recon_val, kl_val = elbo_loss(x_hat, batch_x, mu, logvar, beta=beta)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        train_loss_sum += loss.item() * batch_x.size(0)
        train_recon_sum += recon_val * batch_x.size(0)
        train_kl_sum += kl_val * batch_x.size(0)
        n_train += batch_x.size(0)

    # --- Validate ---
    model.eval()
    val_loss_sum, val_recon_sum, val_kl_sum, n_val = 0, 0, 0, 0
    with torch.no_grad():
        for (batch_x,) in val_loader:
            batch_x = batch_x.to(device)
            x_hat, mu, logvar = model(batch_x)
            loss, recon_val, kl_val = elbo_loss(x_hat, batch_x, mu, logvar, beta=beta)
            val_loss_sum += loss.item() * batch_x.size(0)
            val_recon_sum += recon_val * batch_x.size(0)
            val_kl_sum += kl_val * batch_x.size(0)
            n_val += batch_x.size(0)

    history["epoch"].append(epoch)
    history["train_loss"].append(train_loss_sum / n_train)
    history["train_recon"].append(train_recon_sum / n_train)
    history["train_kl"].append(train_kl_sum / n_train)
    history["val_loss"].append(val_loss_sum / n_val)
    history["val_recon"].append(val_recon_sum / n_val)
    history["val_kl"].append(val_kl_sum / n_val)

    if epoch % 50 == 0 or epoch == 1:
        print(f"  Epoch {epoch:3d}/{N_EPOCHS}  "
              f"beta={beta:.3f}  "
              f"train={history['train_loss'][-1]:.4f}  "
              f"val={history['val_loss'][-1]:.4f}  "
              f"(recon={history['val_recon'][-1]:.4f}, kl={history['val_kl'][-1]:.4f})")

# Save model weights
weights_path = RESULTS_DIR / "vae_weights.pt"
torch.save(model.state_dict(), weights_path)
print(f"  Saved model weights: {weights_path}")

# =============================================================================
# SECTION E: Analysis & Figures
# =============================================================================

print("\nE. Generating analysis figures...")

# --- Encode all windows ---
model.eval()
X_all = torch.tensor(windows).to(device)
with torch.no_grad():
    z_mu_all = model.encode_mu(X_all).cpu().numpy()  # (N_windows, LATENT_DIM)

# Sanity check latent encodings
n_nan = np.isnan(z_mu_all).sum()
n_inf = np.isinf(z_mu_all).sum()
print(f"  Latent matrix: shape={z_mu_all.shape}, NaN={n_nan}, Inf={n_inf}, "
      f"mean={z_mu_all.mean():.4f}, std={z_mu_all.std():.4f}")
if n_nan > 0 or n_inf > 0:
    print("  WARNING: non-finite values in latent matrix — replacing with 0")
    z_mu_all = np.nan_to_num(z_mu_all, nan=0.0, posinf=0.0, neginf=0.0)

# PCA on latent encodings via covariance eigendecomposition. sklearn's
# SVD-based solvers have repeatedly thrown LAPACK convergence errors on this
# (497, 4) matrix — a hand-rolled eigendecomposition on the 4x4 covariance
# matrix is trivially stable for this size.
_zc = z_mu_all - z_mu_all.mean(axis=0, keepdims=True)
_cov = (_zc.T @ _zc) / max(_zc.shape[0] - 1, 1)
_evals, _evecs = np.linalg.eigh(_cov)  # ascending
# Take the top 2 components (last 2 in ascending order), sign-flip if needed
_order = np.argsort(_evals)[::-1][:2]
_pc_axes = _evecs[:, _order]
z_pca = _zc @ _pc_axes  # (N_windows, 2)
print(f"  PCA explained variance ratio (top 2): "
      f"{_evals[_order[0]] / _evals.sum():.3f}, "
      f"{_evals[_order[1]] / _evals.sum():.3f}")

# ---------------------------------------------------------------------------
# Helper: map window midpoint dates to f_residual values
# ---------------------------------------------------------------------------

# f_residual is defined for dates_model (starting S weeks after the first case week).
# Window midpoint is at week_index = start + W//2 (0-based in full series).
# The corresponding f_residual index is (midpoint_week_index - S), if >= 0.

f_resid_dates = f_resid_df["date"].values
f_resid_median = f_resid_df["median"].values
f_resid_q025 = f_resid_df["q025"].values
f_resid_q975 = f_resid_df["q975"].values

# Build a date -> index lookup for f_residual
f_date_to_idx = {d: i for i, d in enumerate(f_resid_dates)}

# For each window, look up f_residual at the midpoint date
f_at_mid = np.full(N_windows, np.nan)
f_q025_at_mid = np.full(N_windows, np.nan)
f_q975_at_mid = np.full(N_windows, np.nan)
for w in range(N_windows):
    mid_date = window_mid_dates[w]
    if mid_date in f_date_to_idx:
        idx = f_date_to_idx[mid_date]
        f_at_mid[w] = f_resid_median[idx]
        f_q025_at_mid[w] = f_resid_q025[idx]
        f_q975_at_mid[w] = f_resid_q975[idx]

# ---------------------------------------------------------------------------
# Helper: map window midpoint dates to serotype info (monthly)
# ---------------------------------------------------------------------------

# Serotype info is monthly; map each window midpoint to its month
sero_date_to_row = {}
for i, row in serotype_df.iterrows():
    sero_date_to_row[pd.Timestamp(row["year_month"]).to_period("M")] = i

dominant_at_mid = []
entropy_at_mid = np.full(N_windows, np.nan)
n_unmatched = 0
for w in range(N_windows):
    mid_period = pd.Timestamp(window_mid_dates[w]).to_period("M")
    if mid_period in sero_date_to_row:
        row_idx = sero_date_to_row[mid_period]
        dominant_at_mid.append(serotype_df.iloc[row_idx]["dominant"])
        entropy_at_mid[w] = serotype_df.iloc[row_idx]["entropy"]
    else:
        dominant_at_mid.append("Unmatched")
        n_unmatched += 1

dominant_at_mid = np.array(dominant_at_mid)
if n_unmatched > 0:
    print(f"  WARNING: {n_unmatched} windows had no matching serotype month — "
          f"check analysis window alignment")

# Year (fractional) for time coloring
mid_dates_ts = pd.to_datetime(window_mid_dates)
year_frac = mid_dates_ts.year + mid_dates_ts.day_of_year / 365.25

# =========================================================================
# Figure 1: vae_latent_space.png  (2x2 panel)
# =========================================================================

print("  Figure 1: vae_latent_space.png")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# (a) Colored by time
sc_a = axes[0, 0].scatter(z_pca[:, 0], z_pca[:, 1], c=year_frac, cmap="viridis",
                           s=8, alpha=0.7)
axes[0, 0].set_title("(a) Colored by year")
axes[0, 0].set_xlabel("PC1")
axes[0, 0].set_ylabel("PC2")
plt.colorbar(sc_a, ax=axes[0, 0], label="Year")

# (b) Colored by f_residual median
mask_b = ~np.isnan(f_at_mid)
sc_b = axes[0, 1].scatter(z_pca[mask_b, 0], z_pca[mask_b, 1], c=f_at_mid[mask_b],
                            cmap="RdBu_r", s=8, alpha=0.7,
                            vmin=-1.5, vmax=1.5)
axes[0, 1].set_title("(b) Colored by f_residual (GP)")
axes[0, 1].set_xlabel("PC1")
axes[0, 1].set_ylabel("PC2")
plt.colorbar(sc_b, ax=axes[0, 1], label="f_residual median")

# (c) Colored by dominant serotype.
# Note: DENV-4 was never the argmax of GAM-smoothed proportions in Singapore
# 2013-2022 (max smoothed ~23%), so no DENV-4 points appear on this panel.
# This is a descriptive artefact of collapsing proportions to a single label,
# not an absence of DENV-4 in circulation.
serotype_colors = {"DENV-1": "#e41a1c", "DENV-2": "#377eb8",
                   "DENV-3": "#4daf4a", "DENV-4": "#984ea3"}
for sero in ["DENV-1", "DENV-2", "DENV-3", "DENV-4"]:
    mask_s = dominant_at_mid == sero
    if mask_s.sum() > 0:
        axes[1, 0].scatter(z_pca[mask_s, 0], z_pca[mask_s, 1],
                           c=serotype_colors.get(sero, "#999999"),
                           label=sero, s=8, alpha=0.7)
axes[1, 0].set_title("(c) Colored by dominant serotype")
axes[1, 0].set_xlabel("PC1")
axes[1, 0].set_ylabel("PC2")
axes[1, 0].legend(fontsize=8, markerscale=2, title="DENV-4: never dominant\n(max smoothed prop. 0.23)")

# (d) Colored by Shannon entropy
mask_d = ~np.isnan(entropy_at_mid)
sc_d = axes[1, 1].scatter(z_pca[mask_d, 0], z_pca[mask_d, 1],
                            c=entropy_at_mid[mask_d], cmap="magma",
                            s=8, alpha=0.7)
axes[1, 1].set_title("(d) Colored by Shannon entropy")
axes[1, 1].set_xlabel("PC1")
axes[1, 1].set_ylabel("PC2")
plt.colorbar(sc_d, ax=axes[1, 1], label="Entropy")

fig.suptitle("VAE Latent Space (PCA of encoder means)", fontsize=14, y=1.01)
fig.tight_layout()
fig.savefig(FIGURES_DIR / "vae_latent_space.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# =========================================================================
# Figure 2: vae_latent_trajectory.png
# =========================================================================

print("  Figure 2: vae_latent_trajectory.png")

fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True)

# Valid f_residual mask for windows
valid_f = ~np.isnan(f_at_mid)

# --- Panel 0: PC1 vs f_residual ---
ax1 = axes[0]
ax1.plot(mid_dates_ts, z_pca[:, 0], color="#377eb8", alpha=0.8, linewidth=0.8,
         label="VAE PC1")
ax1.set_ylabel("VAE PC1", color="#377eb8")
ax1.tick_params(axis="y", labelcolor="#377eb8")

ax1b = ax1.twinx()
ax1b.plot(mid_dates_ts[valid_f], f_at_mid[valid_f], color="#e41a1c", alpha=0.6,
          linewidth=0.8, label="f_residual")
ax1b.set_ylabel("f_residual (GP)", color="#e41a1c")
ax1b.tick_params(axis="y", labelcolor="#e41a1c")

# Switch dates
for sd in switch_dates:
    ax1.axvline(pd.Timestamp(sd), color="grey", linestyle="--", alpha=0.5, linewidth=0.7)

# Spearman correlation
valid_both = valid_f
rho1, p1 = stats.spearmanr(z_pca[valid_both, 0], f_at_mid[valid_both])
ax1.set_title(f"PC1 vs f_residual  (Spearman r = {rho1:.3f}, p = {p1:.2e})")

# --- Panel 1: PC2 vs f_residual ---
ax2 = axes[1]
ax2.plot(mid_dates_ts, z_pca[:, 1], color="#377eb8", alpha=0.8, linewidth=0.8,
         label="VAE PC2")
ax2.set_ylabel("VAE PC2", color="#377eb8")
ax2.tick_params(axis="y", labelcolor="#377eb8")

ax2b = ax2.twinx()
ax2b.plot(mid_dates_ts[valid_f], f_at_mid[valid_f], color="#e41a1c", alpha=0.6,
          linewidth=0.8, label="f_residual")
ax2b.set_ylabel("f_residual (GP)", color="#e41a1c")
ax2b.tick_params(axis="y", labelcolor="#e41a1c")

for sd in switch_dates:
    ax2.axvline(pd.Timestamp(sd), color="grey", linestyle="--", alpha=0.5, linewidth=0.7)

rho2, p2 = stats.spearmanr(z_pca[valid_both, 1], f_at_mid[valid_both])
ax2.set_title(f"PC2 vs f_residual  (Spearman r = {rho2:.3f}, p = {p2:.2e})")

ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
ax2.xaxis.set_major_locator(mdates.YearLocator())

fig.suptitle("VAE Latent Trajectory vs GP Residual", fontsize=14)
fig.tight_layout()
fig.savefig(FIGURES_DIR / "vae_latent_trajectory.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# =========================================================================
# Figure 3: vae_reconstruction_error.png  (Anomaly Detection)
# =========================================================================

print("  Figure 3: vae_reconstruction_error.png")

with torch.no_grad():
    x_hat_all = model(X_all)[0].cpu().numpy()

recon_error = np.mean((windows - x_hat_all) ** 2, axis=1)  # per-window MSE

# Top 10% anomalous
threshold_90 = np.percentile(recon_error, 90)
is_anomaly = recon_error >= threshold_90

# Check if serotype switch windows have higher reconstruction error
# A window is a "switch window" if any switch date falls within it
is_switch_window = np.zeros(N_windows, dtype=bool)
for sd in switch_dates:
    sd_ts = pd.Timestamp(sd)
    for w in range(N_windows):
        win_start = pd.Timestamp(all_dates[window_start_indices[w]])
        win_end = pd.Timestamp(all_dates[window_start_indices[w] + W - 1])
        if win_start <= sd_ts <= win_end:
            is_switch_window[w] = True

# Mann-Whitney test: switch vs non-switch
if is_switch_window.sum() > 0 and (~is_switch_window).sum() > 0:
    mw_stat, mw_p = stats.mannwhitneyu(
        recon_error[is_switch_window],
        recon_error[~is_switch_window],
        alternative="greater"
    )
    switch_median = np.median(recon_error[is_switch_window])
    nonswitch_median = np.median(recon_error[~is_switch_window])
else:
    mw_stat, mw_p = np.nan, np.nan
    switch_median, nonswitch_median = np.nan, np.nan

fig, ax = plt.subplots(figsize=(14, 4))
ax.plot(mid_dates_ts, recon_error, color="#377eb8", linewidth=0.7, alpha=0.8)
ax.fill_between(mid_dates_ts, 0, recon_error,
                where=is_anomaly, color="#e41a1c", alpha=0.3, label="Top 10% anomalies")
ax.axhline(threshold_90, color="#e41a1c", linestyle=":", alpha=0.5, linewidth=0.7)

for sd in switch_dates:
    ax.axvline(pd.Timestamp(sd), color="green", linestyle="--", alpha=0.6, linewidth=1,
               label="Serotype switch")

# De-duplicate legend
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), fontsize=9)

ax.set_ylabel("Reconstruction error (MSE)")
ax.set_title(f"VAE Reconstruction Error as Anomaly Detector\n"
             f"Switch vs non-switch median MSE: {switch_median:.5f} vs {nonswitch_median:.5f} "
             f"(Mann-Whitney p = {mw_p:.3e})")
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
ax.xaxis.set_major_locator(mdates.YearLocator())

fig.tight_layout()
fig.savefig(FIGURES_DIR / "vae_reconstruction_error.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# =========================================================================
# Figure 4: vae_uncertainty_comparison.png
# =========================================================================

print("  Figure 4: vae_uncertainty_comparison.png")

# GP uncertainty: 95% CrI width for f_residual per week
gp_ci_width = f_resid_q975 - f_resid_q025  # length N_model
gp_dates = pd.to_datetime(f_resid_dates)

# VAE uncertainty: for each window, sample z 50 times, decode, compute std
N_SAMPLES = 50
model.eval()
vae_uncertainty = np.zeros(N_windows)

with torch.no_grad():
    for w in range(N_windows):
        x_in = X_all[w : w + 1]
        mu_w, logvar_w = model.encoder(x_in)
        std_w = torch.exp(0.5 * logvar_w)
        decoded_samples = []
        for _ in range(N_SAMPLES):
            z_sample = mu_w + std_w * torch.randn_like(std_w)
            x_decoded = model.decoder(z_sample)
            decoded_samples.append(x_decoded.cpu().numpy())
        decoded_samples = np.concatenate(decoded_samples, axis=0)  # (N_SAMPLES, W)
        # Mean std across the W positions
        vae_uncertainty[w] = np.mean(np.std(decoded_samples, axis=0))

fig, ax1 = plt.subplots(figsize=(14, 4))

ax1.plot(gp_dates, gp_ci_width, color="#e41a1c", alpha=0.7, linewidth=0.8,
         label="GP 95% CrI width")
ax1.set_ylabel("GP f_residual 95% CrI width", color="#e41a1c")
ax1.tick_params(axis="y", labelcolor="#e41a1c")

ax2 = ax1.twinx()
ax2.plot(mid_dates_ts, vae_uncertainty, color="#377eb8", alpha=0.7, linewidth=0.8,
         label="VAE decode std")
ax2.set_ylabel("VAE mean decode std (50 samples)", color="#377eb8")
ax2.tick_params(axis="y", labelcolor="#377eb8")

for sd in switch_dates:
    ax1.axvline(pd.Timestamp(sd), color="grey", linestyle="--", alpha=0.4, linewidth=0.7)

ax1.set_title("Uncertainty Comparison: GP vs VAE")
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
ax1.xaxis.set_major_locator(mdates.YearLocator())

# Combined legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left", fontsize=9)

fig.tight_layout()
fig.savefig(FIGURES_DIR / "vae_uncertainty_comparison.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# =========================================================================
# Figure 5: vae_training_loss.png
# =========================================================================

print("  Figure 5: vae_training_loss.png")

hist_df = pd.DataFrame(history)
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Total loss
axes[0].plot(hist_df["epoch"], hist_df["train_loss"], label="Train", linewidth=0.8)
axes[0].plot(hist_df["epoch"], hist_df["val_loss"], label="Val", linewidth=0.8)
axes[0].set_title("Total ELBO Loss")
axes[0].set_xlabel("Epoch")
axes[0].set_ylabel("Loss")
axes[0].legend()

# Reconstruction loss
axes[1].plot(hist_df["epoch"], hist_df["train_recon"], label="Train", linewidth=0.8)
axes[1].plot(hist_df["epoch"], hist_df["val_recon"], label="Val", linewidth=0.8)
axes[1].set_title("Reconstruction Loss (MSE)")
axes[1].set_xlabel("Epoch")
axes[1].set_ylabel("Loss")
axes[1].legend()

# KL divergence
axes[2].plot(hist_df["epoch"], hist_df["train_kl"], label="Train", linewidth=0.8)
axes[2].plot(hist_df["epoch"], hist_df["val_kl"], label="Val", linewidth=0.8)
axes[2].set_title("KL Divergence")
axes[2].set_xlabel("Epoch")
axes[2].set_ylabel("KL")
axes[2].legend()

fig.suptitle("VAE Training Curves", fontsize=14)
fig.tight_layout()
fig.savefig(FIGURES_DIR / "vae_training_loss.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# =============================================================================
# SECTION F: Summary Statistics
# =============================================================================

print("\nF. Summary statistics...")
print("-" * 50)

# Spearman correlations: each latent dim vs f_residual
summary_rows = []
for d in range(LATENT_DIM):
    rho_d, p_d = stats.spearmanr(z_mu_all[valid_f, d], f_at_mid[valid_f])
    print(f"  Latent dim {d} vs f_residual: Spearman r = {rho_d:.3f}, p = {p_d:.2e}")
    summary_rows.append({
        "metric": f"spearman_latent{d}_vs_fresid",
        "value": rho_d,
        "p_value": p_d,
    })

# PCA components vs f_residual
for pc in range(2):
    rho_pc, p_pc = stats.spearmanr(z_pca[valid_f, pc], f_at_mid[valid_f])
    print(f"  PCA-PC{pc + 1} vs f_residual: Spearman r = {rho_pc:.3f}, p = {p_pc:.2e}")
    summary_rows.append({
        "metric": f"spearman_PC{pc + 1}_vs_fresid",
        "value": rho_pc,
        "p_value": p_pc,
    })

# Fraction of serotype switch windows in top-20% reconstruction error
threshold_80 = np.percentile(recon_error, 80)
in_top20 = recon_error >= threshold_80
if is_switch_window.sum() > 0:
    frac_switch_in_top20 = (is_switch_window & in_top20).sum() / is_switch_window.sum()
else:
    frac_switch_in_top20 = np.nan
print(f"\n  Serotype switch windows in top-20% recon error: "
      f"{frac_switch_in_top20:.1%} "
      f"({(is_switch_window & in_top20).sum()}/{is_switch_window.sum()} windows)")

summary_rows.append({
    "metric": "frac_switch_in_top20_recon",
    "value": frac_switch_in_top20,
    "p_value": np.nan,
})

summary_rows.append({
    "metric": "mannwhitney_switch_recon_p",
    "value": mw_stat,
    "p_value": mw_p,
})

summary_rows.append({
    "metric": "switch_median_mse",
    "value": switch_median,
    "p_value": np.nan,
})

summary_rows.append({
    "metric": "nonswitch_median_mse",
    "value": nonswitch_median,
    "p_value": np.nan,
})

summary_csv = pd.DataFrame(summary_rows)
summary_csv.to_csv(RESULTS_DIR / "vae_comparison_summary.csv", index=False)
print(f"\n  Saved: results/vae_comparison_summary.csv")

print("\n" + "=" * 70)
print("VAE COMPARISON COMPLETE")
print(f"  Figures saved to: {FIGURES_DIR}")
print(f"  Summary saved to: {RESULTS_DIR / 'vae_comparison_summary.csv'}")
print("=" * 70)
