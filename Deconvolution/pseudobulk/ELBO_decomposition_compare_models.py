# ==================================================
# Compare Statescope model runs on pseudobulk data
# ==================================================
#
# This script:
# 1. Loads multiple trained Statescope models
# 2. Loads ground-truth cell type fractions
# 3. Computes ELBO decomposition terms for each model
# 4. Compares predicted vs true cell fractions
# 5. Generates scatterplots colored by cell type
# 6. Summarizes ELBO components across runs
#
# Author: Aryamaan Bose
# Adaptation: Mabel Pronk
# ==================================================

import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt

# --------------------------------------------------
# Import Statescope package from local src directory
# --------------------------------------------------
# Repo src path
import os
import sys
THIS_DIR = os.path.dirname(__file__)
SRC_DIR  = os.path.abspath(os.path.join(THIS_DIR, '..', '..', 'src'))
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from Statescope.Statescope import Statescope

# --------------------------------------------------
# 0) Define input/output paths
# --------------------------------------------------
# Dictionary containing model labels and paths
# to trained Statescope objects
model_paths = {
    "Original": "/net/beegfs/users/P086608/StatescopePro_original/tutorial/pseudobulk_level3/Output_alpha_rem/statescope.pkl",
    "Lamda 0.001": "/net/beegfs/users/P086608/StatescopePro_v2/tutorial/Output_pseudobulk/lam_0_001_rep_10/statescope.pkl",
    "Lamda 0,0001" : "/net/beegfs/users/P086608/StatescopePro_v2/tutorial/Output_pseudobulk/lam_0_0001_rep_10/statescope.pkl",
    'Lamda 0,0001 v2' : '/net/beegfs/users/P086608/StatescopePro_v2/tutorial/Output_pseudobulk/lam_0_0001_rep_10_v2/Output/statescope.pkl'}

gt_path = "/net/beegfs/users/P086608/pseudobulk/level3_celltype/cell_type_fractions.csv"

# --------------------------------------------------
# 1) Helpers
# --------------------------------------------------
def safe_item(x):
    if torch.is_tensor(x):
        return x.detach().cpu().item()
    return float(x)


def safe_pcc(x, y):
    x = np.asarray(x).reshape(-1)
    y = np.asarray(y).reshape(-1)
    if np.allclose(np.std(x), 0) or np.allclose(np.std(y), 0):
        return np.nan
    return np.corrcoef(x, y)[0, 1]


def safe_rmse(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    return np.sqrt(np.mean((x - y) ** 2))


def decompose_parts_using_model_estep(blade, Nu, Beta, Omega):
    PX = safe_item(blade.Estep_PX(Nu, Omega) * (1 / blade.weight))
    PY = safe_item(blade.Estep_PY(Nu, Omega, Beta))
    PF = safe_item(blade.Estep_PF(Beta) * np.sqrt(blade.Ngene / blade.Ncell))
    QX = safe_item(blade.Estep_QX(Omega) * (1 / blade.weight))
    QF = safe_item(blade.Estep_QF(Beta) * np.sqrt(blade.Ngene / blade.Ncell))
    #ELBO = safe_item(blade.E_step(Nu, Beta, Omega))
    #lam = getattr(blade, "lambda_F", "not stored")
    # 2. Check for lambda_F (The Version Switch)
    if hasattr(blade, 'lambda_F'):
        lam = safe_item(blade.lambda_F)
        # New logic: lambda modulates the F terms
        ELBO = PX + PY - QX + lam * (PF - QF)
    else:
        lam = np.nan  # Or np.nan if you want to show it wasn't used
        # Old logic: direct sum/diff
        ELBO = PX + PY + PF - QX - QF

    return {
        "PX": PX,
        "PY": PY,
        "PF": PF,
        "QX": QX,
        "QF": QF,
        "lambda_F": lam,
        "ELBO": ELBO,
    }


def scatter_by_celltype_with_metrics(x_mat, y_mat, title_prefix, xlabel, ylabel, cell_order, n_sample, out_file):
    """
    Generate scatterplot comparing predicted and true
    cell type fractions.

    Parameters
    ----------
    x_mat : tensor
        Predicted fractions (Nsample x Ncell)

    y_mat : tensor
        Ground-truth fractions (Nsample x Ncell)

    Adds:
    - Overall PCC and RMSE in title
    - Per-cell type coloring
    """ 
    x_np = x_mat.detach().cpu().numpy()
    y_np = y_mat.detach().cpu().numpy()

    x = x_np.reshape(-1)
    y = y_np.reshape(-1)

    overall_pcc = safe_pcc(x, y)
    overall_rmse = safe_rmse(x, y)

    ct_labels = np.tile(np.array(cell_order), n_sample)
    cmap = plt.get_cmap("tab20")
    ct_to_color = {ct: cmap(i % 20) for i, ct in enumerate(cell_order)}

    plt.figure(figsize=(8, 6))
    
    # Plot each cell type separately
    for j, ct in enumerate(cell_order):
        xm = x_np[:, j]
        ym = y_np[:, j]
        rmse_ct = safe_rmse(xm, ym)
        pcc_ct = safe_pcc(xm, ym)

        m = (ct_labels == ct)
        plt.scatter(
            x[m], y[m],
            s=18, alpha=0.7,
            label = f'{ct}',
            #label=f"{ct} (RMSE={rmse_ct:.3f}, PCC={pcc_ct:.3f})",
            c=[ct_to_color[ct]]
        )
    # Add diagonal reference line
    lo = min(np.nanmin(x), np.nanmin(y))
    hi = max(np.nanmax(x), np.nanmax(y))
    plt.plot([lo, hi], [lo, hi], linewidth=1)

    plt.title(
        f"{title_prefix}\nOverall PCC = {overall_pcc:.3f}, Total RMSE = {overall_rmse:.3f}"
    )
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.legend(title="Cell types", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.savefig(out_file, dpi=300, bbox_inches="tight")


# --------------------------------------------------
# 2) Load ground-truth data once for alignment
# --------------------------------------------------
first_model = Statescope.load(list(model_paths.values())[0], device="cpu")
cell_order = list(first_model.Celltypes)

gt_df = pd.read_csv(gt_path, index_col=0)

# Ensure GT columns match model cell type order
gt_df = gt_df.reindex(columns=cell_order)

if gt_df.isna().any().any():
    missing_cols = gt_df.columns[gt_df.isna().any()].tolist()
    raise ValueError(
        f"GT has NaNs after reindex (column mismatch). Problem cols: {missing_cols}"
    )

# --------------------------------------------------
# 3) Evaluate all models
# --------------------------------------------------
parts_all = {}

for run_name, model_path in model_paths.items():
    print(f"\nLoading: {run_name}")
    model = Statescope.load(model_path, device="cpu")
    blade = model.BLADE

    Nu = blade.Nu
    Beta = blade.Beta
    Omega = blade.Omega

    Nsample, Ngene, Ncell = blade.Nsample, blade.Ngene, blade.Ncell
    assert gt_df.shape == (Nsample, Ncell), (
        f"GT shape {gt_df.shape} != {(Nsample, Ncell)} for run {run_name}"
    )

    F_model = blade.ExpF(Beta)
    gt = torch.tensor(gt_df.to_numpy(), dtype=Nu.dtype, device=Nu.device)

    parts_all[run_name] = decompose_parts_using_model_estep(blade, Nu, Beta, Omega)
    out_file = f"/net/beegfs/users/P086608/StatescopePro_v2/tutorial/Output_pseudobulk/scatter_{run_name.replace(' ', '_')}.png"
    scatter_by_celltype_with_metrics(
        F_model,
        gt,
        title_prefix=f"Fractions vs GT — {run_name}",
        xlabel="Predicted fractions",
        ylabel="True fractions",
        cell_order=cell_order,
        n_sample=Nsample,
        out_file =out_file
    )

# --------------------------------------------------
# 4) Create combined ELBO decomposition table
# --------------------------------------------------
terms = ["PX", "PY", "PF", "QX", "QF", "lambda_F", "ELBO"]

parts_table = pd.DataFrame(
    {run_name: {term: parts_all[run_name][term] for term in terms}
     for run_name in model_paths.keys()}
)

print("\nCombined decomposition table (ELBO taken directly from blade.E_step)")
print(parts_table)
