#!/usr/bin/env python3
"""Run PGLS with Pagel's lambda estimated by maximum likelihood."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.linalg as la
import scipy.sparse as sp
from Bio import Phylo
from scipy.optimize import minimize_scalar
from scipy.stats import t as tdist


def build_brownian_covariance(tree, names):
    """Brownian phylogenetic covariance as B @ B.T.

    Each branch contributes its length to all descendant-tip covariance pairs.
    """
    name_to_idx = {name: i for i, name in enumerate(names)}
    rows, cols, values = [], [], []
    branch_col = 0

    def visit(clade):
        nonlocal branch_col

        descendants = []
        if clade.is_terminal():
            if clade.name in name_to_idx:
                descendants = [name_to_idx[clade.name]]
        else:
            for child in clade.clades:
                descendants.extend(visit(child))

        branch_length = float(clade.branch_length or 0.0)
        if descendants and branch_length > 0:
            col = branch_col
            branch_col += 1
            value = math.sqrt(branch_length)
            rows.extend(descendants)
            cols.extend([col] * len(descendants))
            values.extend([value] * len(descendants))

        return descendants

    visit(tree.root)

    B = sp.csr_matrix(
        (values, (rows, cols)),
        shape=(len(names), branch_col),
    )
    return (B @ B.T).toarray()


def fit_pgls_pagel_lambda(C, y, x):
    """Fit y ~ 1 + x using PGLS and ML Pagel lambda in [0, 1]."""
    C = np.asarray(C, dtype=float)
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)

    n = len(y)
    X = np.column_stack([np.ones(n), x])

    diag = np.diag(C).copy()
    inv_sqrt_diag = 1.0 / np.sqrt(diag)

    # Convert C to a correlation-like matrix for efficient lambda optimization.
    S = (
        C
        * inv_sqrt_diag[:, None]
        * inv_sqrt_diag[None, :]
    )
    S = (S + S.T) / 2.0

    eigvals, Q = la.eigh(S, check_finite=False)
    eigvals = np.clip(eigvals, 1e-12, None)

    y_rot = Q.T @ (y * inv_sqrt_diag)
    X_rot = Q.T @ (X * inv_sqrt_diag[:, None])
    constant_log_diag = np.log(diag).sum()

    def evaluate(lam, return_details=False):
        # Pagel lambda: off-diagonal phylogenetic covariance is scaled while
        # tip variances remain unchanged.
        transformed_eigvals = 1.0 + lam * (eigvals - 1.0)
        if np.any(transformed_eigvals <= 0):
            return np.inf

        weights = 1.0 / transformed_eigvals

        XtWX = X_rot.T @ (weights[:, None] * X_rot)
        XtWy = X_rot.T @ (weights * y_rot)

        beta = la.solve(XtWX, XtWy, assume_a="sym")

        residuals = y_rot - X_rot @ beta
        rss = np.sum(weights * residuals**2)
        sigma2_ml = rss / n

        logdet = (
            constant_log_diag
            + np.log(transformed_eigvals).sum()
        )
        nll = 0.5 * (
            n * (np.log(2.0 * np.pi * sigma2_ml) + 1.0)
            + logdet
        )

        if not return_details:
            return nll

        covariance_beta = sigma2_ml * la.inv(XtWX)
        se = np.sqrt(np.diag(covariance_beta))
        dfree = n - X.shape[1]
        t_values = beta / se
        p_values = 2.0 * tdist.sf(np.abs(t_values), dfree)
        q975 = tdist.ppf(0.975, dfree)

        return {
            "lambda": float(lam),
            "intercept": float(beta[0]),
            "slope": float(beta[1]),
            "slope_se": float(se[1]),
            "slope_ci_low": float(beta[1] - q975 * se[1]),
            "slope_ci_high": float(beta[1] + q975 * se[1]),
            "slope_p": float(p_values[1]),
            "sigma2_ml": float(sigma2_ml),
            "loglik": float(-nll),
        }

    optimum = minimize_scalar(
        evaluate,
        bounds=(0.0, 1.0),
        method="bounded",
        options={"xatol": 1e-6, "maxiter": 100},
    )

    return evaluate(optimum.x, return_details=True)


def fit_subset(data, full_covariance, mask, label, predictor):
    idx = np.flatnonzero(mask.to_numpy())
    subset_covariance = full_covariance[np.ix_(idx, idx)]

    result = fit_pgls_pagel_lambda(
        subset_covariance,
        data.loc[mask, "pro_over_phe_plus_leutt"].to_numpy(),
        data.loc[mask, predictor].to_numpy(),
    )

    result.update({
        "analysis": label,
        "n": int(mask.sum()),
        "response": "Pro/(Phe+LeuTT)",
        "predictor": predictor,
    })
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = pd.read_csv(args.dataset)
    tree = Phylo.read(str(args.tree), "newick")

    # Dataset is already exact-matched to tree tips by preprocessing.
    covariance = build_brownian_covariance(
        tree,
        data["tree_name"].tolist(),
    )

    all_species = pd.Series(True, index=data.index)
    direct_gl = data["generation_length_direct"].astype(bool)
    zero_ambiguous = data["zero_ambiguous_nt"].astype(bool)

    results = [
        fit_subset(
            data,
            covariance,
            all_species,
            "PRIMARY_all_species_log10_generation_length",
            "log10_generation_length",
        ),
        fit_subset(
            data,
            covariance,
            direct_gl,
            "SENSITIVITY_direct_generation_length_only",
            "log10_generation_length",
        ),
        fit_subset(
            data,
            covariance,
            zero_ambiguous,
            "SENSITIVITY_zero_ambiguous_nucleotides",
            "log10_generation_length",
        ),
        fit_subset(
            data,
            covariance,
            direct_gl & zero_ambiguous,
            "SENSITIVITY_direct_GL_and_zero_ambiguous",
            "log10_generation_length",
        ),
        fit_subset(
            data,
            covariance,
            all_species,
            "SENSITIVITY_untransformed_generation_length",
            "generation_length_d",
        ),
    ]

    output = pd.DataFrame(results)[[
        "analysis",
        "n",
        "response",
        "predictor",
        "lambda",
        "intercept",
        "slope",
        "slope_se",
        "slope_ci_low",
        "slope_ci_high",
        "slope_p",
        "sigma2_ml",
        "loglik",
    ]]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)

    print(output.to_string(index=False))


if __name__ == "__main__":
    main()
