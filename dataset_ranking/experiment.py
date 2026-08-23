import argparse
import hashlib
import importlib.metadata
import os
import pickle
import platform
import random
import re
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import linalg, stats
from scipy.stats import pearsonr, rankdata
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA_DIR = ROOT / "data"
RESULT_DIR = ROOT / "results"
EPS = 1e-6

STATE_MAPPING = {
    "01": "AL", "02": "AK", "04": "AZ", "05": "AR", "06": "CA",
    "08": "CO", "09": "CT", "10": "DE", "11": "DC", "12": "FL",
    "13": "GA", "15": "HI", "16": "ID", "17": "IL", "18": "IN",
    "19": "IA", "20": "KS", "21": "KY", "22": "LA", "23": "ME",
    "24": "MD", "25": "MA", "26": "MI", "27": "MN", "28": "MS",
    "29": "MO", "30": "MT", "31": "NE", "32": "NV", "33": "NH",
    "34": "NJ", "35": "NM", "36": "NY", "37": "NC", "38": "ND",
    "39": "OH", "40": "OK", "41": "OR", "42": "PA", "44": "RI",
    "45": "SC", "46": "SD", "47": "TN", "48": "TX", "49": "UT",
    "50": "VT", "51": "VA", "53": "WA", "54": "WV", "55": "WI",
    "56": "WY", "72": "PR",
}

ACS_STATES = [
    "AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "HI",
    "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI",
    "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC",
    "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT",
    "VT", "VA", "WA", "WV", "WI", "WY", "PR",
]

ACS_COLUMNS = [
    "ST", "AGEP", "COW", "SCHL", "MAR", "OCCP", "POBP",
    "RELP", "WKHP", "SEX", "RAC1P", "PINCP",
]
NHANES_DATASET = "nhanes_cholesterol"
NHANES_YEAR_ORDER = [1999, 2011, 2015, 2001, 2009, 2017, 2013, 2003, 2005, 2007]
NHANES_DEFAULT_SOURCES = ["sex0", "sex1", "race0", "race2", "age2", "age4"]
NHANES_DEFAULT_TARGET = "late_race1"


def create_acs_pickle(data_dir=DATA_DIR):
    data_dir = Path(data_dir)
    path = data_dir / "acs_data_df.p"
    if path.exists():
        return path

    try:
        from folktables import ACSDataSource
    except ImportError as exc:
        raise ImportError(
            "Install folktables to create data/acs_data_df.p automatically."
        ) from exc

    data_dir.mkdir(parents=True, exist_ok=True)
    data_source = ACSDataSource(
        survey_year="2018",
        horizon="1-Year",
        survey="person",
        root_dir=data_dir,
    )
    acs_data = data_source.get_data(states=ACS_STATES, download=True)
    columns = acs_data.columns.intersection(ACS_COLUMNS)
    acs_data_df = acs_data[columns]
    acs_data_df.index = acs_data_df.groupby("ST").cumcount()
    acs_data_df.to_pickle(path)
    return path


def order_nhanes_files(files):
    by_year = {}
    for path in files:
        match = re.search(r".*([0-9]{4})\.XPT", str(path))
        if match is None:
            return files
        year = int(match.group(1))
        by_year.setdefault(year, []).append(path)

    ordered = []
    for year in NHANES_YEAR_ORDER:
        year_files = by_year.pop(year, [])
        if year == 2007:
            first = [p for p in year_files if Path(p).name == "TRIGLY_E2007.XPT"]
            rest = [p for p in year_files if Path(p).name != "TRIGLY_E2007.XPT"]
            year_files = first + rest
        ordered.extend(year_files)

    for year in by_year:
        ordered.extend(by_year[year])
    return ordered


def create_nhanes_pickle(data_dir=DATA_DIR):
    data_dir = Path(data_dir)
    path = data_dir / "nhanes_cholesterol_tableshift.p"
    if path.exists():
        return path

    try:
        from tableshift import get_iid_dataset
        import tableshift.core.data_source as tableshift_data_source
    except ImportError as exc:
        raise ImportError(
            "Install tableshift or add it to PYTHONPATH to create "
            "the NHANES cholesterol pickle automatically."
        ) from exc

    data_dir.mkdir(parents=True, exist_ok=True)
    orig_glob = tableshift_data_source.glob.glob

    def ordered_glob(pattern):
        files = orig_glob(pattern)
        if pattern.endswith("*.XPT"):
            return order_nhanes_files(files)
        return files

    tableshift_data_source.glob.glob = ordered_glob
    try:
        dataset = get_iid_dataset(
            NHANES_DATASET,
            cache_dir=str(data_dir / "tableshift_cache"),
        )
    finally:
        tableshift_data_source.glob.glob = orig_glob

    X_parts = []
    y_parts = []
    for split in ["train", "test", "validation"]:
        X_split, y_split, _, _ = dataset.get_pandas(split)
        X_parts.append(X_split)
        y_parts.append(y_split)

    X = pd.concat(X_parts, axis=0)
    y = pd.concat(y_parts, axis=0)

    pd.concat([X, y], axis=1).to_pickle(path)
    return path


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dataframe_sha256(X):
    digest = hashlib.sha256()
    schema = "\n".join(
        f"{column}:{dtype}" for column, dtype in zip(X.columns, X.dtypes)
    )
    digest.update(schema.encode())
    digest.update(
        pd.util.hash_pandas_object(X, index=False).to_numpy().tobytes()
    )
    return digest.hexdigest()


def environment_versions():
    packages = [
        "joblib", "matplotlib", "numpy", "pandas", "scikit-learn", "scipy"
    ]
    versions = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return {
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": versions,
    }


def whiten_data(X, return_info=False):
    X = X.astype(float)
    if not np.isfinite(X.to_numpy()).all():
        raise ValueError("Whitening input contains non-finite values.")
    X_centered = X - X.mean(axis=0)
    values = X_centered.to_numpy()

    keep = []
    basis = []
    tol = 1e-10 * max(1, np.linalg.norm(values, axis=0).max())
    for j in range(values.shape[1]):
        residual = values[:, j].copy()
        for vector in basis:
            residual -= np.dot(vector, residual) * vector
        # Reorthogonalize so exact dummy-column dependencies stay below tolerance.
        for vector in basis:
            residual -= np.dot(vector, residual) * vector
        norm = np.linalg.norm(residual)
        if norm > tol:
            keep.append(j)
            basis.append(residual / norm)

    X_centered = X_centered.iloc[:, keep]
    values = X_centered.to_numpy()
    cov = values.T @ values / len(values)
    R = np.linalg.cholesky(cov)
    X_white = linalg.solve_triangular(
        R, X_centered.to_numpy().T, lower=True
    ).T

    columns = [f"whitened_{i + 1}" for i in range(len(keep))]
    result = pd.DataFrame(X_white, index=X.index, columns=columns)
    if not return_info:
        return result

    kept_columns = X.columns[keep].tolist()
    info = {
        "input_columns": X.columns.tolist(),
        "kept_columns": kept_columns,
        "dropped_columns": [c for c in X.columns if c not in kept_columns],
        "output_columns": columns,
    }
    return result, info


def scale_data(X):
    return pd.DataFrame(
        StandardScaler(with_mean=False).fit_transform(X),
        index=X.index,
        columns=X.columns,
    )


def get_nhanes_cuts(X, cutoff=2017):
    pre = X.nhanes_year < cutoff
    years = sorted(X.loc[pre, "nhanes_year"].unique())
    if cutoff == 2017:
        periods = [years[:4], years[4:7], years[7:]]
    else:
        periods = np.array_split(years, 3)

    cuts = {
        "early": X.nhanes_year.isin(periods[0]),
        "middle": X.nhanes_year.isin(periods[1]),
        "recent": X.nhanes_year.isin(periods[2]),
    }
    for sex in [0, 1]:
        cuts[f"sex{sex}"] = pre & (X.RIAGENDR == sex)
        cuts[f"recent_sex{sex}"] = cuts["recent"] & (X.RIAGENDR == sex)
    for race in [0, 1, 2]:
        cuts[f"race{race}"] = pre & (X.RIDRETH_merged == race)
        cuts[f"recent_race{race}"] = cuts["recent"] & (X.RIDRETH_merged == race)
    for age in range(5):
        cuts[f"age{age}"] = pre & (X.RIDAGEYR == age)
        cuts[f"recent_age{age}"] = cuts["recent"] & (X.RIDAGEYR == age)
    return cuts


def get_nhanes_target_mask(X, target_name):
    if target_name == "2017":
        return X.nhanes_year == 2017, 2017

    target_mask = X.nhanes_year.isin([2015, 2017])
    cutoff = 2015
    if target_name.startswith("late_sex"):
        target_mask &= X.RIAGENDR == int(target_name[-1])
    elif target_name.startswith("late_race"):
        target_mask &= X.RIDRETH_merged == int(target_name[-1])
    elif target_name.startswith("late_age"):
        target_mask &= X.RIDAGEYR == int(target_name[-1])
    elif target_name != "late":
        raise ValueError(f"Unknown NHANES target: {target_name}")
    return target_mask, cutoff


def prepare_acs(data_dir=DATA_DIR, target_state=6):
    data_path = create_acs_pickle(data_dir)
    X_raw = pd.read_pickle(data_path)
    raw_data_sha256 = dataframe_sha256(X_raw)

    COW_map = {5.0: 99, 7.0: 99, 8.0: 99}
    SCHL_map = {k: 99 for k in [
        23., 17., 24., 14., 15., 13., 1., 12., 11.,
        9., 10., 8., 6., 7., 5., 2., 4., 3.,
    ]}
    MAR_map = {2.0: 99, 4.0: 99}
    RAC1P_map = {9: 99, 3: 99, 5: 99, 7: 99, 4: 99}

    all_states = X_raw.ST.unique()
    all_states = all_states[all_states != 11]
    state_sizes = pd.DataFrame({
        "n": [sum(X_raw.ST == state) for state in all_states],
        "state": all_states,
    }).sort_values(["n", "state"], ascending=[False, True])

    states = state_sizes.iloc[:11].state.astype(int).tolist()
    if target_state not in states:
        raise ValueError(f"State {target_state} is not in the 11-state ACS pool.")
    source_states = sorted(state for state in states if state != target_state)
    qq_source_states = [state for state in states if state != target_state]

    X_raw = X_raw[X_raw.ST.isin(states)].dropna().reset_index(drop=True)
    X_raw = X_raw[X_raw.PINCP >= 0].copy()
    X_raw["PINCP"] = np.log(X_raw.PINCP)
    X_raw["WKHP"] = X_raw.WKHP.map(lambda x: x if x != 40 else 99)

    for column, mapping in [
        ("COW", COW_map), ("MAR", MAR_map),
        ("SCHL", SCHL_map), ("RAC1P", RAC1P_map),
    ]:
        X_raw[column] = X_raw[column].map(
            lambda x: x if x not in mapping else mapping[x]
        )

    categorical = X_raw.drop(
        ["ST", "PINCP", "AGEP", "WKHP", "RELP", "POBP", "OCCP"], axis=1
    )
    X = pd.concat([
        pd.get_dummies(categorical, columns=categorical.columns, dtype=float),
        X_raw[["ST", "PINCP", "AGEP", "WKHP"]],
    ], axis=1).reset_index(drop=True)

    covariates = [c for c in X.columns if c not in ["ST", "PINCP"]]
    model_columns = covariates + ["PINCP"]
    X_model = scale_data(X[model_columns])
    X_model["ST"] = X.ST.to_numpy()

    X_duc, whitening = whiten_data(
        X_model[covariates], return_info=True
    )
    X_duc = X_duc.reset_index(drop=True)
    X_duc["PINCP"] = X_model.PINCP.to_numpy()
    X_duc["ST"] = X.ST.to_numpy()

    white_cov = np.cov(X_duc.drop(columns=["PINCP", "ST"]), rowvar=False, ddof=0)
    whitening_error = np.abs(white_cov - np.eye(white_cov.shape[0])).max()

    duc_columns = [c for c in X_duc.columns if c not in ["ST", "PINCP"]] + ["PINCP"]
    target = X_model[X_model.ST == target_state][model_columns].copy()
    sources = [
        X_model[X_model.ST == state][model_columns].copy()
        for state in source_states
    ]
    target_duc = X_duc[X_duc.ST == target_state][duc_columns].copy()
    sources_duc = [
        X_duc[X_duc.ST == state][duc_columns].copy()
        for state in source_states
    ]
    target_cov = target_duc.drop(columns="PINCP").mean(axis=0).to_numpy()
    target_model_cov = target.drop(columns="PINCP").mean(axis=0).to_numpy()

    labels = [STATE_MAPPING[f"{int(state):02d}"] for state in source_states]
    qq_data = X_duc[[c for c in X_duc.columns if c not in ["ST", "PINCP"]]
                    + ["PINCP", "ST"]].copy()

    return {
        "name": "acs_income",
        "target": target,
        "sources": sources,
        "target_duc": target_duc,
        "sources_duc": sources_duc,
        "target_cov": target_cov,
        "target_model_cov": target_model_cov,
        "labels": labels,
        "source_ids": list(source_states),
        "qq_source_ids": qq_source_states,
        "qq_labels": [STATE_MAPPING[f"{state:02d}"] for state in qq_source_states],
        "target_id": target_state,
        "site_column": "ST",
        "response": "PINCP",
        "qq_data": qq_data,
        "whitening_error": whitening_error,
        "whitening": whitening,
        "data_sha256": dataframe_sha256(
            X[model_columns + ["ST"]].round(12)
        ),
        "raw_data_sha256": raw_data_sha256,
        "data_file_sha256": file_sha256(data_path),
        "model_columns": model_columns,
    }


def prepare_nhanes(data_dir=DATA_DIR, source_names=None, target_name=None):
    data_dir = Path(data_dir)
    if source_names is None:
        source_names = NHANES_DEFAULT_SOURCES
    if target_name is None:
        target_name = NHANES_DEFAULT_TARGET

    data_path = create_nhanes_pickle(data_dir)
    X = pd.read_pickle(data_path).reset_index(drop=True)
    response = "LBDLDL"

    target_mask, cutoff = get_nhanes_target_mask(X, target_name)
    cuts = get_nhanes_cuts(X, cutoff)
    missing = [name for name in source_names if name not in cuts]
    if missing:
        raise ValueError(f"Unknown NHANES source cuts: {missing}")

    drop_columns = ["nhanes_year", "RIAGENDR", "RIDRETH_merged"]
    source_raw = [
        X[cuts[name]].drop(columns=drop_columns).reset_index(drop=True)
        for name in source_names
    ]
    target_raw = (
        X[target_mask]
        .drop(columns=drop_columns)
        .reset_index(drop=True)
    )
    if len(target_raw) < 1030:
        raise ValueError(f"Target needs at least 1030 rows: {len(target_raw)}")
    if min(map(len, source_raw)) < 1000:
        sizes = dict(zip(source_names, map(len, source_raw)))
        raise ValueError(f"Every source needs 1000 rows: {sizes}")

    sources = [scale_data(data) for data in source_raw]
    target = scale_data(target_raw)

    X_all = pd.concat(source_raw + [target_raw], ignore_index=True)
    X_all_standardized = scale_data(X_all)
    datasets = np.concatenate([
        *[np.repeat(i, len(data)) for i, data in enumerate(source_raw)],
        np.repeat(len(source_raw), len(target_raw)),
    ])

    X_duc, whitening = whiten_data(
        X_all_standardized.drop(columns=response), return_info=True
    )
    X_duc = X_duc.reset_index(drop=True)
    X_duc[response] = X_all_standardized[response].to_numpy()
    X_duc["domain"] = datasets

    white_cov = np.cov(
        X_duc.drop(columns=[response, "domain"]), rowvar=False, ddof=0
    )
    whitening_error = np.abs(white_cov - np.eye(white_cov.shape[0])).max()

    duc_columns = [c for c in X_duc.columns if c not in ["domain", response]]
    duc_columns = duc_columns + [response]
    sources_duc = [
        X_duc[X_duc.domain == i][duc_columns].reset_index(drop=True)
        for i in range(len(source_names))
    ]
    target_duc = (
        X_duc[X_duc.domain == len(source_names)][duc_columns]
        .reset_index(drop=True)
    )
    target_cov = target_duc.drop(columns=response).mean(axis=0).to_numpy()
    target_model_cov = target.drop(columns=response).mean(axis=0).to_numpy()

    return {
        "name": "nhanes_cholesterol",
        "target": target,
        "sources": sources,
        "target_duc": target_duc,
        "sources_duc": sources_duc,
        "target_cov": target_cov,
        "target_model_cov": target_model_cov,
        "labels": [f"D{i}" for i in range(1, len(source_names) + 1)],
        "source_ids": list(range(len(source_names))),
        "qq_source_ids": list(range(len(source_names))),
        "qq_labels": [f"D{i}" for i in range(1, len(source_names) + 1)],
        "target_id": len(source_names),
        "site_column": "domain",
        "response": response,
        "source_names": list(source_names),
        "target_name": target_name,
        "qq_data": X_duc.copy(),
        "whitening_error": whitening_error,
        "whitening": whitening,
        "data_sha256": dataframe_sha256(X),
        "raw_data_sha256": dataframe_sha256(X),
        "data_file_sha256": file_sha256(data_path),
        "model_columns": target.columns.tolist(),
    }


def get_weight(y_reg, x_source):
    denom = np.dot(x_source, x_source)
    if denom <= 0:
        return 0, 0

    weight = float(np.dot(x_source, y_reg) / denom)
    weight = np.clip(weight, 0, 1)

    if weight == 0:
        duc = 0
    elif weight == 1:
        duc = np.var(x_source) / np.var(y_reg)
    else:
        duc = pearsonr(y_reg, x_source)[0] ** 2
    return weight, duc


def compute_baseline_scores(target, source, seed):
    pipeline = make_pipeline(
        StandardScaler(),
        PCA(
            n_components=3,
            svd_solver="randomized",
            random_state=seed,
            iterated_power=0,
            power_iteration_normalizer="none",
        ),
    )
    pipeline.fit(pd.concat([target, source]))
    target_pca = pipeline.transform(target)
    source_pca = pipeline.transform(source)

    def fit_density(X):
        grid = GridSearchCV(
            KernelDensity(kernel="gaussian"),
            {"bandwidth": np.logspace(-1, 1, 20)},
            cv=5,
            n_jobs=1,
        )
        grid.fit(X)
        return grid.best_estimator_

    kde_target = fit_density(target_pca)
    kde_source = fit_density(source_pca)
    target_density = np.exp(kde_target.score_samples(target_pca)) + EPS
    source_density = np.exp(kde_source.score_samples(target_pca)) + EPS

    kl = stats.entropy(
        target_density / target_density.sum(),
        source_density / source_density.sum(),
    )
    log_ratio = np.log(target_density) - np.log(source_density)
    kde_log_ratio = np.mean(log_ratio) / np.std(log_ratio)
    return kl, kde_log_ratio


def _run_acs_trial(data, seed, n_test=1000, n_t=30, n_k=1000,
                   n_target_mean=None):
    target = data["target"]
    sources = data["sources"]
    target_duc = data["target_duc"]
    sources_duc = data["sources_duc"]
    response = data["response"]
    target_cov = data["target_cov"]
    target_model_cov = data["target_model_cov"]

    features = [c for c in target.columns if c != response]
    duc_features = [c for c in target_duc.columns if c != response]

    if len(target) < n_test + n_t + (n_target_mean or 0):
        raise ValueError("The target does not contain enough disjoint rows.")
    if any(len(source) < n_k for source in sources):
        raise ValueError("A source does not contain enough rows.")

    py_rng = random.Random(seed)
    np_rng = np.random.RandomState(seed)
    test_index = py_rng.sample(target.index.tolist(), n_test)
    samples_test = target[target.index.isin(test_index)]
    target_train = target.drop(test_index)
    samples_t = target_train.sample(n=n_t, random_state=np_rng)
    target_index = samples_t.index.tolist()
    samples_t_duc = target_duc.loc[samples_t.index]

    target_mean_index = None
    if n_target_mean is not None:
        target_mean_pool = target_train.drop(target_index)
        target_mean_index = py_rng.sample(
            target_mean_pool.index.tolist(), n_target_mean
        )
        target_model_cov = (
            target.loc[target_mean_index, features].mean(axis=0).to_numpy()
        )
        target_cov = (
            target_duc.loc[target_mean_index, duc_features]
            .mean(axis=0).to_numpy()
        )

    Xbar_target = samples_t[features].mean(axis=0).to_numpy()
    y_reg = target_model_cov - Xbar_target
    Xbar_target_duc = samples_t_duc[duc_features].mean(axis=0).to_numpy()
    y_reg_duc = target_cov - Xbar_target_duc

    weights = []
    duc_weights = []
    duc = []
    kl = []
    score_x = []
    pool_mse = []
    weighted_mse = []
    source_hashes = []

    samples_sources = [
        source.sample(n=n_k, random_state=np_rng) for source in sources
    ]

    for i, samples_s in enumerate(samples_sources):
        samples_s_duc = sources_duc[i].loc[samples_s.index]
        x_source = (
            samples_s[features].mean(axis=0).to_numpy()
            - Xbar_target
        )
        weight, _ = get_weight(y_reg, x_source)
        x_source_duc = (
            samples_s_duc[duc_features].mean(axis=0).to_numpy()
            - Xbar_target_duc
        )
        duc_weight, usefulness = get_weight(y_reg_duc, x_source_duc)
        score_kl, score_kde = compute_baseline_scores(
            samples_t[features], samples_s[features], seed
        )

        weights.append(weight)
        duc_weights.append(duc_weight)
        duc.append(usefulness)
        kl.append(score_kl)
        score_x.append(score_kde)
        source_hashes.append(index_hash(samples_s.index))

    fit_target = RandomForestRegressor(random_state=np_rng, n_jobs=1).fit(
        samples_t[features], samples_t[response]
    )
    y_target = fit_target.predict(samples_test[features])

    for samples_s in samples_sources:
        fit_pool = RandomForestRegressor(random_state=np_rng, n_jobs=1).fit(
            pd.concat([samples_s[features], samples_t[features]]),
            pd.concat([samples_s[response], samples_t[response]]),
        )
        y_pool = fit_pool.predict(samples_test[features])
        pool_mse.append(np.mean((samples_test[response] - y_pool) ** 2))

    for weight, samples_s in zip(weights, samples_sources):
        fit_source = RandomForestRegressor(random_state=np_rng, n_jobs=1).fit(
            samples_s[features], samples_s[response]
        )
        y_source = fit_source.predict(samples_test[features])
        y_weighted = weight * y_source + (1 - weight) * y_target
        weighted_mse.append(np.mean((samples_test[response] - y_weighted) ** 2))

    return {
        "target_mse": np.mean((samples_test[response] - y_target) ** 2),
        "pool_mse_list": pool_mse,
        "weighted_mse_list": weighted_mse,
        "weights_list": weights,
        "duc_weights_list": duc_weights,
        "duc_list": duc,
        "kl_list": kl,
        "score_x_list": score_x,
        "seed": seed,
        "test_index_hash": index_hash(samples_test.index),
        "target_index_hash": index_hash(target_index),
        "target_mean_index_hash": (
            None if target_mean_index is None else index_hash(target_mean_index)
        ),
        "source_index_hashes": source_hashes,
    }


def _run_seeded_trial(data, seed, n_test=1000, n_t=30, n_k=1000,
                      n_target_mean=None):
    target = data["target"]
    sources = data["sources"]
    target_duc = data["target_duc"]
    sources_duc = data["sources_duc"]
    response = data["response"]
    target_cov = data["target_cov"]
    target_model_cov = data["target_model_cov"]

    features = [c for c in target.columns if c != response]
    duc_features = [c for c in target_duc.columns if c != response]

    if len(target) < n_test + n_t + (n_target_mean or 0):
        raise ValueError("The target does not contain enough disjoint rows.")
    if any(len(source) < n_k for source in sources):
        raise ValueError("A source does not contain enough rows.")

    rng = random.Random(seed)
    test_index = rng.sample(target.index.tolist(), n_test)
    samples_test = target[target.index.isin(test_index)]
    target_train = target.drop(test_index)
    target_index = rng.sample(target_train.index.tolist(), n_t)
    samples_t = target[target.index.isin(target_index)]
    samples_t_duc = target_duc.loc[samples_t.index]

    target_mean_index = None
    if n_target_mean is not None:
        target_mean_pool = target_train.drop(target_index)
        target_mean_index = rng.sample(
            target_mean_pool.index.tolist(), n_target_mean
        )
        target_model_cov = (
            target.loc[target_mean_index, features].mean(axis=0).to_numpy()
        )
        target_cov = (
            target_duc.loc[target_mean_index, duc_features]
            .mean(axis=0).to_numpy()
        )

    Xbar_target = samples_t[features].mean(axis=0).to_numpy()
    y_reg = target_model_cov - Xbar_target
    Xbar_target_duc = samples_t_duc[duc_features].mean(axis=0).to_numpy()
    y_reg_duc = target_cov - Xbar_target_duc

    weights = []
    duc_weights = []
    duc = []
    kl = []
    score_x = []
    pool_mse = []
    weighted_mse = []
    source_hashes = []

    fit_target = RandomForestRegressor(random_state=seed, n_jobs=1).fit(
        samples_t[features], samples_t[response]
    )
    y_target = fit_target.predict(samples_test[features])

    for i, source in enumerate(sources):
        source_seed = seed + i
        samples_s = source.sample(n=n_k, random_state=source_seed)
        samples_s_duc = sources_duc[i].loc[samples_s.index]
        x_source = (
            samples_s[features].mean(axis=0).to_numpy()
            - Xbar_target
        )
        weight, _ = get_weight(y_reg, x_source)
        x_source_duc = (
            samples_s_duc[duc_features].mean(axis=0).to_numpy()
            - Xbar_target_duc
        )
        duc_weight, usefulness = get_weight(y_reg_duc, x_source_duc)

        score_kl, score_kde = compute_baseline_scores(
            samples_t[features], samples_s[features], source_seed
        )

        fit_pool = RandomForestRegressor(
            random_state=source_seed, n_jobs=1
        ).fit(
            pd.concat([samples_s[features], samples_t[features]]),
            pd.concat([samples_s[response], samples_t[response]]),
        )
        y_pool = fit_pool.predict(samples_test[features])

        fit_source = RandomForestRegressor(
            random_state=source_seed, n_jobs=1
        ).fit(
            samples_s[features], samples_s[response]
        )
        y_source = fit_source.predict(samples_test[features])
        y_weighted = weight * y_source + (1 - weight) * y_target

        weights.append(weight)
        duc_weights.append(duc_weight)
        duc.append(usefulness)
        kl.append(score_kl)
        score_x.append(score_kde)
        source_hashes.append(index_hash(samples_s.index))
        pool_mse.append(np.mean((samples_test[response] - y_pool) ** 2))
        weighted_mse.append(np.mean((samples_test[response] - y_weighted) ** 2))

    return {
        "target_mse": np.mean((samples_test[response] - y_target) ** 2),
        "pool_mse_list": pool_mse,
        "weighted_mse_list": weighted_mse,
        "weights_list": weights,
        "duc_weights_list": duc_weights,
        "duc_list": duc,
        "kl_list": kl,
        "score_x_list": score_x,
        "seed": seed,
        "test_index_hash": index_hash(test_index),
        "target_index_hash": index_hash(target_index),
        "target_mean_index_hash": (
            None if target_mean_index is None else index_hash(target_mean_index)
        ),
        "source_index_hashes": source_hashes,
    }


def run_trial(data, seed, n_test=1000, n_t=30, n_k=1000,
              n_target_mean=None):
    if data["name"] == "acs_income":
        return _run_acs_trial(
            data, seed, n_test, n_t, n_k, n_target_mean
        )
    return _run_seeded_trial(
        data, seed, n_test, n_t, n_k, n_target_mean
    )


def index_hash(index):
    values = np.asarray(index, dtype=np.int64)
    return hashlib.sha256(values.tobytes()).hexdigest()


def experiment_config(data, seed_start, trials, n_test, n_t, n_k,
                      n_target_mean):
    return {
        "code_sha256": file_sha256(__file__),
        "dataset": data["name"],
        "data_sha256": data.get("data_sha256"),
        "raw_data_sha256": data.get("raw_data_sha256"),
        "target_id": data["target_id"],
        "source_ids": list(data["source_ids"]),
        "labels": list(data["labels"]),
        "response": data["response"],
        "model_columns": list(data.get("model_columns", [])),
        "whitening": data.get("whitening"),
        "seed_start": seed_start,
        "trials": trials,
        "n_test": n_test,
        "n_t": n_t,
        "n_k": n_k,
        "n_target_mean": n_target_mean,
        "environment": environment_versions(),
    }


def run_trials(data, trials=1000, seed_start=123, jobs=1,
               n_test=1000, n_t=30, n_k=1000,
               n_target_mean=None,
               checkpoint_path=None, checkpoint_every=4):
    if checkpoint_every < 1:
        raise ValueError("checkpoint_every must be positive.")
    seeds = list(range(seed_start, seed_start + trials))
    config = experiment_config(
        data, seed_start, trials, n_test, n_t, n_k, n_target_mean
    )

    def run_seed(seed):
        return run_trial(
            data, seed, n_test=n_test, n_t=n_t, n_k=n_k,
            n_target_mean=n_target_mean,
        )

    results = []
    if checkpoint_path is not None:
        checkpoint_path = Path(checkpoint_path)
        if checkpoint_path.exists():
            with open(checkpoint_path, "rb") as file:
                checkpoint = pickle.load(file)
            completed_seeds = checkpoint["seeds"]
            results = checkpoint["results"]
            if checkpoint.get("config") != config:
                raise ValueError("Checkpoint settings do not match this run.")
            if completed_seeds != seeds[:len(completed_seeds)]:
                raise ValueError("Checkpoint seeds do not match this run.")
            if len(results) != len(completed_seeds):
                raise ValueError("Checkpoint results are incomplete.")

    for start in range(len(results), len(seeds), checkpoint_every):
        batch_seeds = seeds[start:start + checkpoint_every]
        if jobs == 1:
            batch = [run_seed(seed) for seed in batch_seeds]
        else:
            batch = Parallel(n_jobs=jobs, backend="loky")(
                delayed(run_seed)(seed) for seed in batch_seeds
            )
        results.extend(batch)

        if checkpoint_path is not None:
            checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
            checkpoint = {
                "config": config,
                "seeds": seeds[:len(results)],
                "results": results,
            }
            temporary = checkpoint_path.with_suffix(checkpoint_path.suffix + ".tmp")
            with open(temporary, "wb") as file:
                pickle.dump(checkpoint, file)
            os.replace(temporary, checkpoint_path)
            print(f"Checkpoint: {len(results)}/{len(seeds)} trials", flush=True)
    return results


def summarize_results(results, labels):
    if not results:
        raise ValueError("No trial results were provided.")
    K = len(labels)
    list_keys = [
        "pool_mse_list", "weighted_mse_list", "weights_list", "duc_list",
        "kl_list", "score_x_list",
    ]
    for result in results:
        if any(len(result[key]) != K for key in list_keys):
            raise ValueError("A trial result does not match the source labels.")
        values = np.concatenate([
            np.atleast_1d(np.asarray(result[key], dtype=float))
            for key in ["target_mse"] + list_keys
        ])
        if not np.isfinite(values).all():
            raise ValueError("A trial result contains non-finite values.")

    weighted_mse = np.vstack([
        np.asarray(result["weighted_mse_list"], dtype=float)
        for result in results
    ])
    ranks = np.vstack([rankdata(row, method="average") for row in weighted_mse])

    return pd.DataFrame({
        "label": labels,
        "target_mse": np.mean([r["target_mse"] for r in results]),
        "pool_mse": [
            np.mean([r["pool_mse_list"][i] for r in results])
            for i in range(K)
        ],
        "weighted_mse": weighted_mse.mean(axis=0),
        "weighted_mse_std": weighted_mse.std(axis=0),
        "avg_rank": ranks.mean(axis=0),
        "weights": [
            np.mean([r["weights_list"][i] for r in results])
            for i in range(K)
        ],
        "duc": [
            np.mean([r["duc_list"][i] for r in results])
            for i in range(K)
        ],
        "neg_kl": [
            -np.mean([r["kl_list"][i] for r in results])
            for i in range(K)
        ],
        "neg_score_x": [
            -np.mean([r["score_x_list"][i] for r in results])
            for i in range(K)
        ],
    })


def save_results(data, results, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary = summarize_results(results, data["labels"])
    trials = data.get("trials", len(results))
    seed_start = data.get("seed_start")
    seeds = None if seed_start is None else list(range(seed_start, seed_start + trials))
    payload = {
        "config": experiment_config(
            data, seed_start, trials, data.get("n_test"), data.get("n_t"),
            data.get("n_k"), data.get("n_target_mean"),
        ),
        "dataset": data["name"],
        "labels": data["labels"],
        "target_id": data["target_id"],
        "source_ids": data["source_ids"],
        "site_column": data["site_column"],
        "response": data["response"],
        "whitening_error": data["whitening_error"],
        "prediction_data": "standardized",
        "baseline_data": "standardized",
        "duc_data": "whitened",
        "code_sha256": file_sha256(__file__),
        "data_sha256": data.get("data_sha256"),
        "raw_data_sha256": data.get("raw_data_sha256"),
        "data_file_sha256": data.get("data_file_sha256"),
        "environment": environment_versions(),
        "whitening": data.get("whitening"),
        "trials": trials,
        "seed_start": seed_start,
        "seed_end": None if seed_start is None else seed_start + trials - 1,
        "seeds": seeds,
        "n_test": data.get("n_test"),
        "n_t": data.get("n_t"),
        "n_k": data.get("n_k"),
        "n_target_mean": data.get("n_target_mean"),
        "results": results,
        "summary": summary,
    }
    if "source_names" in data:
        payload["source_names"] = data["source_names"]
    if "target_name" in data:
        payload["target_name"] = data["target_name"]

    with open(output_dir / f"{data['name']}_results.pkl", "wb") as file:
        pickle.dump(payload, file)
    summary.to_csv(output_dir / f"{data['name']}_summary.csv", index=False)
    return summary


def run_dataset(name, data_dir, output_dir, args):
    if name == "acs":
        data = prepare_acs(data_dir, target_state=args.acs_target_state)
    elif name == "nhanes":
        data = prepare_nhanes(
            data_dir,
            source_names=args.nhanes_sources,
            target_name=args.nhanes_target,
        )
    else:
        raise ValueError(name)

    data.update({
        "trials": args.trials,
        "seed_start": args.seed_start,
        "n_test": args.n_test,
        "n_t": args.n_t,
        "n_k": args.n_k,
        "n_target_mean": args.n_target_mean,
    })

    print(f"Running {data['name']} ({args.trials} trials)")
    print(f"Whitening covariance error: {data['whitening_error']:.3e}")
    results = run_trials(
        data,
        trials=args.trials,
        seed_start=args.seed_start,
        jobs=args.jobs,
        n_test=args.n_test,
        n_t=args.n_t,
        n_k=args.n_k,
        n_target_mean=args.n_target_mean,
        checkpoint_path=Path(output_dir) / f"{data['name']}_checkpoint.pkl",
        checkpoint_every=args.checkpoint_every,
    )
    summary = save_results(data, results, output_dir)
    print(summary.to_string(index=False))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", choices=["acs", "nhanes", "all"], default="all")
    parser.add_argument("--data-dir", default=DATA_DIR)
    parser.add_argument("--output-dir", default=RESULT_DIR)
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--seed-start", type=int, default=123)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--n-test", type=int, default=1000)
    parser.add_argument("--n-t", type=int, default=30)
    parser.add_argument("--n-k", type=int, default=1000)
    parser.add_argument("--n-target-mean", type=int, default=None)
    parser.add_argument("--checkpoint-every", type=int, default=4)
    parser.add_argument("--acs-target-state", type=int, default=6)
    parser.add_argument("--nhanes-sources", nargs="+", default=None)
    parser.add_argument("--nhanes-target", default=None)
    args = parser.parse_args()

    if args.dataset in ["acs", "all"]:
        run_dataset("acs", args.data_dir, args.output_dir, args)
    if args.dataset in ["nhanes", "all"]:
        run_dataset("nhanes", args.data_dir, args.output_dir, args)


if __name__ == "__main__":
    main()
