from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
from scipy import stats
from scipy.stats import pearsonr

try:
    from .experiment import whiten_data
except ImportError:
    from experiment import whiten_data


def get_sample_means_vec(X, source, target, site_column,
                         scaled=True, response_column=None,
                         prewhitened=True):
    ns = np.sum(X[site_column] == source)
    nt = np.sum(X[site_column] == target)
    sigma_inv = 1 / np.sqrt((1 / ns + 1 / nt))

    excluded = [site_column]
    if response_column is not None:
        excluded.append(response_column)
    covariate_columns = [c for c in X.columns if c not in excluded]

    if prewhitened:
        X_covariates = X[covariate_columns].copy()
    else:
        X_covariates = whiten_data(X[covariate_columns])

    X_new = X_covariates.copy()
    if response_column is not None:
        response = X[response_column].astype(float)
        X_new[response_column] = (
            response - response.mean()
        ) / response.std(ddof=0)
    X_new[site_column] = X[site_column].to_numpy()

    value_columns = [c for c in X_new.columns if c != site_column]
    source_means = X_new.loc[X_new[site_column] == source, value_columns].mean(axis=0)
    target_means = X_new.loc[X_new[site_column] == target, value_columns].mean(axis=0)
    diff = source_means - target_means

    if scaled:
        return sigma_inv * diff
    return diff


def _qqprep(vec):
    x = np.asarray(vec, dtype=float)
    order = np.argsort(x)
    y = x[order]
    p = (np.arange(1, len(x) + 1) - 0.5) / len(x)
    z = stats.norm.ppf(p)
    q1, q3 = np.quantile(y, [0.25, 0.75])
    z1, z3 = stats.norm.ppf([0.25, 0.75])
    slope = (q3 - q1) / (z3 - z1)
    intercept = q1 - slope * z1
    r2 = float(np.corrcoef(z, y)[0, 1] ** 2)
    return z, y, order, slope, intercept, r2


def get_qqplot(vec, title="", ax=None, response_column=None,
               show_title=True, show_r2=True, show_axis_labels=True):
    z, y, order, slope, intercept, r2 = _qqprep(vec)

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    if response_column is not None:
        if not isinstance(vec, pd.Series) or response_column not in vec.index:
            raise ValueError(f"{response_column} is not in the QQ vector.")
        names = vec.index.to_numpy()[order]
        response = names == response_column
        ax.scatter(z[~response], y[~response], s=18, color="0.1",
                   linewidths=0, alpha=0.95)
        ax.scatter(z[response], y[response], s=95, color="#0072FF",
                   marker="*", linewidths=0.8, zorder=5)
    else:
        ax.scatter(z, y, s=18, color="0.1", linewidths=0, alpha=0.95)

    xx = np.array([z.min(), z.max()])
    ax.plot(xx, intercept + slope * xx, color="#B23A48", linewidth=2)

    if show_title:
        ax.set_title(title, fontsize=18, pad=8, loc="left")
    if show_axis_labels:
        ax.set_xlabel("Theoretical Quantiles", fontsize=11)
        ax.set_ylabel("Sample Quantiles", fontsize=11)
    else:
        ax.set_xlabel("")
        ax.set_ylabel("")

    if show_r2:
        ax.text(
            0.02, 0.98, f"R² = {r2:.3f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=18,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.2),
        )

    ax.set_facecolor("white")
    ax.grid(True, which="major", color="0.92", linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1.0)
    ax.tick_params(axis="both", colors="black", width=1.0, length=4)
    return fig, ax


def plot_target_qq(qq_data, target_id, source_ids, site_column, response,
                   labels=None, target_label=None, title=None, path=None):
    if labels is None:
        labels = [str(s) for s in source_ids]
    if target_label is None:
        target_label = str(target_id)

    n_cols = 5 if len(source_ids) > 6 else 3
    n_rows = int(np.ceil(len(source_ids) / n_cols))
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(5 * n_cols, 4.5 * n_rows),
        facecolor="white",
    )
    axes = np.atleast_1d(axes).ravel()

    for ax, source, label in zip(axes, source_ids, labels):
        vec = get_sample_means_vec(
            qq_data, target_id, source, site_column,
            scaled=True, response_column=response, prewhitened=True,
        )
        get_qqplot(
            vec, f"{target_label} vs {label}", ax=ax,
            response_column=response, show_title=True,
            show_r2=True, show_axis_labels=True,
        )

    for ax in axes[len(source_ids):]:
        ax.axis("off")

    if title is not None:
        fig.suptitle(title, fontsize=18, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_figure(fig, path)
    return fig


def plot_qq_matrix(qq_data, ids, site_column, response,
                   labels=None, title=None, path=None):
    if labels is None:
        labels = [str(i) for i in ids]

    n = len(ids)
    fig, axes = plt.subplots(n, n, figsize=(2.6 * n, 2.6 * n),
                             facecolor="white")
    for i, source in enumerate(ids):
        for j, target in enumerate(ids):
            ax = axes[i, j]
            if source == target:
                ax.set_xticks([])
                ax.set_yticks([])
                ax.text(0.5, 0.5, labels[i], transform=ax.transAxes,
                        ha="center", va="center", fontsize=12)
                for spine in ax.spines.values():
                    spine.set_color("0.6")
            else:
                vec = get_sample_means_vec(
                    qq_data, source, target, site_column,
                    scaled=True, response_column=response, prewhitened=True,
                )
                get_qqplot(
                    vec, ax=ax, response_column=response,
                    show_title=False, show_r2=True, show_axis_labels=False,
                )
                if i < n - 1:
                    ax.set_xticklabels([])
                if j > 0:
                    ax.set_yticklabels([])

    fig.subplots_adjust(wspace=0.08, hspace=0.08)
    if title is not None:
        fig.suptitle(title, fontsize=18, y=0.92)
    save_figure(fig, path)
    return fig


def plot_comparison(summary, title, y_column="avg_rank", y_label=None,
                    highlight_map=None, baseline_summary=None, path=None):
    x_columns = ["duc", "neg_kl", "neg_score_x"]
    x_labels = [
        "Avg Estimated DUC",
        "Avg Estimated Negative KL",
        "Avg Estimated Negative Domain Classifier Score",
    ]
    if y_label is None:
        y_label = y_column
    if baseline_summary is None:
        baseline_summary = summary

    panel_data = [summary, baseline_summary, baseline_summary]
    y_by_label = summary.set_index("label")[y_column]
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True,
                             facecolor="white")

    for ax, df, x_column, x_label in zip(axes, panel_data, x_columns, x_labels):
        x = pd.to_numeric(df[x_column], errors="coerce").to_numpy(float)
        y = pd.to_numeric(
            y_by_label.loc[df.label], errors="coerce"
        ).to_numpy(float)
        mask = np.isfinite(x) & np.isfinite(y)

        x_min = np.nanmin(x[mask])
        x_max = np.nanmax(x[mask])
        x_range = x_max - x_min
        x_pad = 0.1 * (x_range if x_range > 0 else 1)

        ax.scatter(x, y, s=100, color="black", zorder=2)
        if mask.sum() >= 2 and x_range > 0:
            line_x = np.linspace(x_min, x_max, 100)
            slope, intercept = np.polyfit(x[mask], y[mask], 1)
            ax.plot(line_x, intercept + slope * line_x,
                    color="lightgrey", linestyle="--", linewidth=2)

        if highlight_map is None:
            highlights = {df.label.iloc[np.nanargmax(x)]}
        else:
            highlights = highlight_map.get(x_column, set())
        highlight = df.label.isin(highlights).to_numpy()
        ax.scatter(x[highlight], y[highlight], s=140, color="green", zorder=4)

        if mask.sum() >= 3 and np.std(x[mask]) > 0 and np.std(y[mask]) > 0:
            corr_text = f"Correlation = {pearsonr(x[mask], y[mask])[0]:.2f}"
        else:
            corr_text = "Correlation = NA"
        ax.text(0.95, 0.95, corr_text, transform=ax.transAxes,
                ha="right", va="top", fontsize=12, weight="bold")

        dx = 0.02 * (x_range if x_range > 0 else 1)
        for i, label in enumerate(df.label):
            color = "green" if label in highlights else "black"
            ax.text(x[i] + dx, y[i], label, color=color, weight="semibold")

        ax.set_xlim(x_min - x_pad, x_max + x_pad)
        ax.set_xlabel(x_label, fontsize=16)
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.tick_params(axis="both", labelsize=14)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_color("black")
            spine.set_linewidth(1.2)

    axes[0].set_ylabel(y_label, fontsize=16)
    fig.suptitle(title, fontsize=24, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure(fig, path)
    return fig


def save_figure(fig, path):
    if path is not None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white")
