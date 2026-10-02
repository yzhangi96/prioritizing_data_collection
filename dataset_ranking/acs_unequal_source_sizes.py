import argparse
import os
import pickle
import random
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, rankdata
from sklearn.ensemble import RandomForestRegressor

from dataset_ranking.experiment_equal_source_sizes import (
    DATA_DIR, RESULT_DIR, compute_baseline_scores, experiment_config,
    file_sha256, get_weight, index_hash, prepare_acs, procedure_sha256,
    summarize_results,
)


SOURCE_SIZE_BY_STATE = {
    "FL": 5000,
    "GA": 300,
    "IL": 5000,
    "MI": 4000,
    "NJ": 300,
    "NY": 300,
    "NC": 300,
    "OH": 3000,
    "PA": 3000,
    "TX": 5000,
}

DATA = None
SOURCE_SIZES = None
N_TARGET_MEAN = None


def run_trial(seed):
    target = DATA["target"]
    sources = DATA["sources"]
    target_duc = DATA["target_duc"]
    sources_duc = DATA["sources_duc"]
    response = DATA["response"]
    features = [column for column in target.columns if column != response]
    duc_features = [
        column for column in target_duc.columns if column != response
    ]

    py_rng = random.Random(seed)
    np_rng = np.random.RandomState(seed)
    test_index = py_rng.sample(target.index.tolist(), 1000)
    samples_test = target[target.index.isin(test_index)]
    target_train = target.drop(test_index)
    samples_t = target_train.sample(n=30, random_state=np_rng)
    samples_t_duc = target_duc.loc[samples_t.index]

    target_model_cov = DATA["target_model_cov"]
    target_cov = DATA["target_cov"]
    target_mean_index = None
    if N_TARGET_MEAN is not None:
        target_mean_pool = target_train.drop(samples_t.index)
        target_mean_index = py_rng.sample(
            target_mean_pool.index.tolist(), N_TARGET_MEAN
        )
        target_model_cov = (
            target.loc[target_mean_index, features].mean(axis=0).to_numpy()
        )
        target_cov = (
            target_duc.loc[target_mean_index, duc_features]
            .mean(axis=0).to_numpy()
        )

    xbar_t = samples_t[features].mean(axis=0).to_numpy()
    y_reg = target_model_cov - xbar_t
    xbar_t_duc = samples_t_duc[duc_features].mean(axis=0).to_numpy()
    y_reg_duc = target_cov - xbar_t_duc

    samples_sources = [
        source.sample(n=size, random_state=np_rng)
        for source, size in zip(sources, SOURCE_SIZES)
    ]
    weights = []
    duc_weights = []
    ducs = []
    kls = []
    scores = []
    source_hashes = []
    for i, samples_s in enumerate(samples_sources):
        x_source = samples_s[features].mean(axis=0).to_numpy() - xbar_t
        weight, _ = get_weight(y_reg, x_source)
        samples_s_duc = sources_duc[i].loc[samples_s.index]
        x_source_duc = (
            samples_s_duc[duc_features].mean(axis=0).to_numpy()
            - xbar_t_duc
        )
        duc_weight, duc = get_weight(y_reg_duc, x_source_duc)
        kl, score = compute_baseline_scores(
            samples_t[features], samples_s[features], seed
        )
        weights.append(weight)
        duc_weights.append(duc_weight)
        ducs.append(duc)
        kls.append(kl)
        scores.append(score)
        source_hashes.append(index_hash(samples_s.index))

    fit_target = RandomForestRegressor(random_state=np_rng, n_jobs=1).fit(
        samples_t[features], samples_t[response]
    )
    y_target = fit_target.predict(samples_test[features])

    pool_mse = []
    for samples_s in samples_sources:
        fit_pool = RandomForestRegressor(random_state=np_rng, n_jobs=1).fit(
            pd.concat([samples_s[features], samples_t[features]]),
            pd.concat([samples_s[response], samples_t[response]]),
        )
        y_pool = fit_pool.predict(samples_test[features])
        pool_mse.append(np.mean((samples_test[response] - y_pool) ** 2))

    weighted_mse = []
    for weight, samples_s in zip(weights, samples_sources):
        fit_source = RandomForestRegressor(random_state=np_rng, n_jobs=1).fit(
            samples_s[features], samples_s[response]
        )
        y_source = fit_source.predict(samples_test[features])
        prediction = weight * y_source + (1 - weight) * y_target
        weighted_mse.append(
            np.mean((samples_test[response] - prediction) ** 2)
        )

    return {
        "target_mse": np.mean((samples_test[response] - y_target) ** 2),
        "pool_mse_list": pool_mse,
        "weighted_mse_list": weighted_mse,
        "weights_list": weights,
        "duc_weights_list": duc_weights,
        "duc_list": ducs,
        "kl_list": kls,
        "score_x_list": scores,
        "seed": seed,
        "test_index_hash": index_hash(samples_test.index),
        "target_index_hash": index_hash(samples_t.index),
        "target_mean_index_hash": (
            None if target_mean_index is None else index_hash(target_mean_index)
        ),
        "source_index_hashes": source_hashes,
    }


def unequal_size_config(data, source_sizes, trials, seed_start,
                        n_target_mean=None):
    config = experiment_config(
        data, seed_start, trials, 1000, 30, None, n_target_mean
    )
    config["unequal_source_size_code_sha256"] = file_sha256(__file__)
    config["source_sizes"] = list(source_sizes)
    return config


def run_trials(data, source_sizes, trials=1000, seed_start=123, jobs=1,
               n_target_mean=None, checkpoint_path=None,
               checkpoint_every=4):
    if checkpoint_every < 1:
        raise ValueError("checkpoint_every must be positive.")

    global DATA, SOURCE_SIZES, N_TARGET_MEAN
    DATA = data
    SOURCE_SIZES = list(source_sizes)
    N_TARGET_MEAN = n_target_mean
    seeds = list(range(seed_start, seed_start + trials))
    config = unequal_size_config(
        data, source_sizes, trials, seed_start, n_target_mean
    )
    run_hash = procedure_sha256(config)

    results = []
    if checkpoint_path is not None:
        checkpoint_path = Path(checkpoint_path)
        if checkpoint_path.exists():
            with open(checkpoint_path, "rb") as file:
                checkpoint = pickle.load(file)
            if checkpoint["config"] != config:
                raise ValueError("Checkpoint settings do not match this run.")
            results = checkpoint["results"]
            if [result["seed"] for result in results] != seeds[:len(results)]:
                raise ValueError("Checkpoint seeds do not match this run.")
            if any(
                result.get("run_config_sha256")
                != procedure_sha256(checkpoint["config"])
                for result in results
            ):
                raise ValueError("Checkpoint provenance does not match this run.")

    pending = seeds[len(results):]
    pool = None if jobs == 1 else get_context("fork").Pool(jobs)
    try:
        for start in range(0, len(pending), checkpoint_every):
            batch_seeds = pending[start:start + checkpoint_every]
            batch = (
                [run_trial(seed) for seed in batch_seeds]
                if pool is None else pool.map(run_trial, batch_seeds)
            )
            for result in batch:
                result["run_config_sha256"] = run_hash
            results.extend(batch)
            if checkpoint_path is not None:
                checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
                temporary = checkpoint_path.with_suffix(".pkl.tmp")
                with open(temporary, "wb") as file:
                    pickle.dump({"config": config, "results": results}, file)
                os.replace(temporary, checkpoint_path)
            print(f"Checkpoint: {len(results)}/{trials} trials", flush=True)
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    if [result["seed"] for result in results] != seeds:
        raise ValueError("Completed trial seeds do not match this run.")
    return results, config


def summarize_unequal_results(results, labels):
    weighted_mse = np.asarray([
        result["weighted_mse_list"] for result in results
    ])
    ranks = np.asarray([
        rankdata(row, method="average") for row in weighted_mse
    ])
    return pd.DataFrame({
        "label": labels,
        "weighted_mse": weighted_mse.mean(axis=0),
        "avg_rank": ranks.mean(axis=0),
        "duc": np.asarray([
            result["duc_list"] for result in results
        ]).mean(axis=0),
        "neg_kl": -np.asarray([
            result["kl_list"] for result in results
        ]).mean(axis=0),
        "neg_score_x": -np.asarray([
            result["score_x_list"] for result in results
        ]).mean(axis=0),
    })


def save_results(data, source_sizes, results, config, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = summarize_results(results, data["labels"])
    summary.insert(1, "n_source", source_sizes)
    correlations = {
        column: pearsonr(summary[column], summary["avg_rank"])[0]
        for column in ["duc", "neg_kl", "neg_score_x"]
    }
    payload = {
        "config": config,
        "dataset": data["name"],
        "labels": data["labels"],
        "source_ids": data["source_ids"],
        "source_sizes": list(source_sizes),
        "trials": len(results),
        "seed_start": config["seed_start"],
        "seeds": [result["seed"] for result in results],
        "results": results,
        "summary": summary,
        "correlations": correlations,
    }
    with open(output_dir / "acs_unequal_source_sizes_results.pkl", "wb") as file:
        pickle.dump(payload, file)
    summary.to_csv(
        output_dir / "acs_unequal_source_sizes_summary.csv", index=False
    )
    return summary, correlations


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default=DATA_DIR)
    parser.add_argument(
        "--output-dir", default=RESULT_DIR
    )
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--seed-start", type=int, default=123)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--n-target-mean", type=int)
    parser.add_argument("--checkpoint-every", type=int, default=4)
    args = parser.parse_args()

    data = prepare_acs(args.data_dir, target_state=6)
    source_sizes = [SOURCE_SIZE_BY_STATE[label] for label in data["labels"]]
    output_dir = Path(args.output_dir)
    print(pd.DataFrame({
        "label": data["labels"],
        "n_source": source_sizes,
    }).to_string(index=False))
    results, config = run_trials(
        data,
        source_sizes,
        trials=args.trials,
        seed_start=args.seed_start,
        jobs=args.jobs,
        n_target_mean=args.n_target_mean,
        checkpoint_path=output_dir / "acs_unequal_source_sizes_checkpoint.pkl",
        checkpoint_every=args.checkpoint_every,
    )
    summary, correlations = save_results(
        data, source_sizes, results, config, output_dir
    )
    print(summary.to_string(index=False))
    print(correlations)


if __name__ == "__main__":
    main()
