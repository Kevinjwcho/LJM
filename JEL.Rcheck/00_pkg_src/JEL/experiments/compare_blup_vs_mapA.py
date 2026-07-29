"""
Compare two equivalent estimators for hat{b}_i in the JEL framework
(longitudinal + prior, no survival contribution):

  (1) BLUP closed form:
      b_hat = c + D Z^T (Omega + Z D Z^T)^{-1} (Y - Z c),  Omega = diag(sigma^2 / w)

  (2) MAP via numerical optimization (= argmax_b [log f(Y|b) + log f(b)]),
      using the *kernel-weighted* form exactly as written in ll_lla.cpp,
      with Delta = 0, l0u = 0  (i.e., survival contribution turned off).
"""

import numpy as np
from scipy.optimize import minimize
np.random.seed(20260501)

# -----------------------------
# Synthetic experiment design
# -----------------------------
n_subject = 50          # number of test subjects
nK        = 2           # number of longitudinal processes
s         = 3.0         # landmark time
h         = 1.0         # kernel bandwidth
ker       = "epan"      # epanechnikov

# Population-level center c (length 2*nK = 4): intercept & slope per process
c_true = np.array([1.0, 0.5,    # process 1: intercept, slope
                   2.0, -0.3])  # process 2: intercept, slope

# Random-effects covariance D (4x4, SPD)
rng = np.random.default_rng(0)
A = rng.normal(size=(4, 4))
D = A @ A.T * 0.2 + np.eye(4) * 0.4
print("D (random-effects covariance):")
print(np.round(D, 3))

sigma2 = np.array([0.25, 0.16])  # residual variance per process

def K_h(t, s, h, ker="epan"):
    u = (t - s) / h
    if ker == "epan":
        return np.where(np.abs(u) <= 1, 0.75 * (1 - u**2) / h, 0.0)
    elif ker == "gaussian":
        return np.exp(-0.5*u*u) / (h*np.sqrt(2*np.pi))
    else:
        raise ValueError(ker)

# -----------------------------
# Generate per-subject data
# -----------------------------
# Each subject has irregular obs times for each process
subjects = []
for i in range(n_subject):
    # draw subject random effects
    b_i = rng.multivariate_normal(c_true, D)
    rec = {"b_true": b_i}

    Y_blocks = []
    Z_blocks = []
    w_blocks = []
    sigma_per_obs = []

    for k in range(nK):
        n_ik = rng.integers(low=4, high=12)         # # obs for process k
        t_ik = np.sort(rng.uniform(0, 6, size=n_ik)) # time grid
        # Z_ik = [1, t - s] -- intercept + slope at landmark
        Z_ik = np.column_stack([np.ones(n_ik), t_ik - s])
        # process k uses b_i[2k:2k+2]
        b_k  = b_i[2*k : 2*k + 2]
        mu_ik = Z_ik @ b_k
        # Gaussian noise scaled by 1/sqrt(K_h(t-s)) -> closed-form Omega = diag(sigma^2/w)
        w_ik = K_h(t_ik, s, h, ker)
        # only keep obs with positive kernel weight (within bandwidth)
        keep = w_ik > 0
        if keep.sum() < 2:
            # ensure at least 2 in-window obs by oversampling near s
            extra_t = rng.uniform(s - 0.8*h, s + 0.8*h, size=4)
            t_ik = np.sort(np.concatenate([t_ik[keep], extra_t]))
            Z_ik = np.column_stack([np.ones(len(t_ik)), t_ik - s])
            mu_ik = Z_ik @ b_k
            w_ik = K_h(t_ik, s, h, ker)
            keep = w_ik > 0

        t_ik, Z_ik, mu_ik, w_ik = t_ik[keep], Z_ik[keep], mu_ik[keep], w_ik[keep]
        # Generate observations with TRUE residual variance sigma2[k]
        # (scaling by kernel weight is a *likelihood* convention, not data-gen)
        eps = rng.normal(scale=np.sqrt(sigma2[k]), size=len(t_ik))
        Y_ik = mu_ik + eps

        Y_blocks.append(Y_ik)
        Z_blocks.append(Z_ik)
        w_blocks.append(w_ik)
        sigma_per_obs.append(np.full(len(t_ik), sigma2[k]))

    rec["Y"]      = np.concatenate(Y_blocks)
    rec["w"]      = np.concatenate(w_blocks)
    rec["sigma2"] = np.concatenate(sigma_per_obs)
    # Block-diagonal Z (mi x q), with q = 2*nK
    from scipy.linalg import block_diag
    rec["Z"]      = block_diag(*Z_blocks)
    subjects.append(rec)

# -----------------------------
# Method 1: BLUP closed form
# -----------------------------
def blup_closed_form(Y, Z, w, sigma2_per_obs, D, c):
    Omega = np.diag(sigma2_per_obs / np.maximum(w, 1e-12))
    A = Omega + Z @ D @ Z.T
    return c + D @ Z.T @ np.linalg.solve(A, Y - Z @ c)

# -----------------------------
# Method 2: MAP from ll_lla form
#   ll_long(b) = -(sum_w/2) log(2 pi)
#                - 1/2 sum_j w_j log(sigma_j^2)
#                - 1/2 sum_j w_j (Y_j - Z_j b)^2 / sigma_j^2
#   ll_b(b)    = -q/2 log(2 pi) - 1/2 log|D| - 1/2 (b-c)^T D^{-1} (b-c)
# Negative of the sum is what ll_lla returns (we MINIMIZE).
# Survival is dropped (Delta = 0, l0u = 0).
# -----------------------------
def neg_ll_lla_no_surv(b, Y, Z, w, sigma2_per_obs, D, c):
    wpos = np.maximum(w, 1e-12)
    resid = Y - Z @ b
    sum_w = wpos.sum()
    quad_long  = float(np.sum(wpos * resid**2 / sigma2_per_obs))
    logdet_long = float(np.sum(wpos * np.log(sigma2_per_obs)))
    ll_long = (-0.5 * sum_w * np.log(2*np.pi)
               - 0.5 * logdet_long
               - 0.5 * quad_long)
    bc = b - c
    sign, logdetD = np.linalg.slogdet(D)
    Dinv_bc = np.linalg.solve(D, bc)
    ll_b = (-0.5 * len(b) * np.log(2*np.pi)
            - 0.5 * logdetD
            - 0.5 * float(bc @ Dinv_bc))
    return -(ll_long + ll_b)

def grad_neg_ll_lla_no_surv(b, Y, Z, w, sigma2_per_obs, D, c):
    wpos = np.maximum(w, 1e-12)
    resid = Y - Z @ b
    grad_long = Z.T @ (wpos * resid / sigma2_per_obs)        # d ll_long / db
    grad_b    = -np.linalg.solve(D, b - c)                    # d ll_b / db
    return -(grad_long + grad_b)

# -----------------------------
# Compute both for each subject
# -----------------------------
b_blup = []
b_mapA = []
b_true_all = []
for rec in subjects:
    b1 = blup_closed_form(rec["Y"], rec["Z"], rec["w"], rec["sigma2"], D, c_true)
    res = minimize(neg_ll_lla_no_surv,
                   x0 = c_true.copy(),
                   args = (rec["Y"], rec["Z"], rec["w"], rec["sigma2"], D, c_true),
                   jac = grad_neg_ll_lla_no_surv,
                   method = "BFGS",
                   options = dict(gtol=1e-8))
    b_blup.append(b1)
    b_mapA.append(res.x)
    b_true_all.append(rec["b_true"])

b_blup    = np.array(b_blup)
b_mapA    = np.array(b_mapA)
b_true_all = np.array(b_true_all)

diff = b_blup - b_mapA
print("\n=== Comparison: BLUP closed form  vs.  MAP (ucminf-style, surv dropped) ===")
print(f"shape: {b_blup.shape} (n_subject x q)")
print(f"max |b_BLUP - b_A|       : {np.max(np.abs(diff)):.3e}")
print(f"mean |b_BLUP - b_A|      : {np.mean(np.abs(diff)):.3e}")
print(f"max relative diff        : {np.max(np.abs(diff) / (np.abs(b_blup)+1e-12)):.3e}")

print("\n=== Sample of first 5 subjects, b_hat (each row is a subject; cols = b1_0,b1_1,b2_0,b2_1) ===")
print("BLUP:")
print(np.round(b_blup[:5], 5))
print("MAP (A):")
print(np.round(b_mapA[:5], 5))

# Sanity: relate to truth
err_blup = b_blup - b_true_all
err_A    = b_mapA - b_true_all
print(f"\nRMSE vs truth, BLUP : {np.sqrt(np.mean(err_blup**2)):.4f}")
print(f"RMSE vs truth, MAP-A: {np.sqrt(np.mean(err_A**2)):.4f}")

# Save numeric outputs
out_path = "/sessions/inspiring-funny-thompson/mnt/outputs/compare_blup_vs_mapA_results.txt"
with open(out_path, "w") as f:
    f.write("Setup\n=====\n")
    f.write(f"n_subject={n_subject}, nK={nK}, s={s}, h={h}, kernel=epanechnikov\n")
    f.write(f"sigma2={sigma2}, c_true={c_true}\nD=\n{np.round(D,4)}\n\n")
    f.write("Method 1 (BLUP closed form) vs Method 2 (MAP via ucminf-style numeric opt, surv dropped)\n\n")
    f.write(f"max |b_BLUP - b_A| = {np.max(np.abs(diff)):.3e}\n")
    f.write(f"mean|b_BLUP - b_A| = {np.mean(np.abs(diff)):.3e}\n")
    f.write(f"RMSE vs truth, BLUP  = {np.sqrt(np.mean(err_blup**2)):.4f}\n")
    f.write(f"RMSE vs truth, MAP-A = {np.sqrt(np.mean(err_A**2)):.4f}\n")
    f.write(f"\nBLUP first 5:\n{np.round(b_blup[:5],6)}\n")
    f.write(f"\nMAP-A first 5:\n{np.round(b_mapA[:5],6)}\n")
print(f"\nSaved: {out_path}")
