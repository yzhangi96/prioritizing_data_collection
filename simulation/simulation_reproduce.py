import argparse
import hashlib
import os
import pickle
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn

from simulation_comparison_plot import plot_simulation
from simulation_score_x import compute_score_x


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_matrix(path, rows, columns=30):
    values = np.fromfile(path, dtype="<f8")
    if values.size != rows * columns:
        raise ValueError(f"Unexpected matrix size in {path}.")
    return values.reshape(rows, columns)


def run_trial(seed, rscript):
    generator = HERE / "generate_simulation_trial.R"
    with tempfile.TemporaryDirectory(prefix=f"simulation_{seed}_") as temporary:
        temporary = Path(temporary)
        subprocess.run(
            [rscript, str(generator), str(seed), str(temporary)],
            check=True, stdout=subprocess.DEVNULL,
        )
        manifest = pd.read_csv(temporary / "manifest.csv")
        reference = read_matrix(temporary / "reference.bin", 700)
        scores = []
        for row in manifest.itertuples(index=False):
            source = read_matrix(
                temporary / f"candidate_{int(row.candidate):02d}.bin",
                int(row.rows), int(row.columns),
            )
            scores.append(compute_score_x(
                reference, source,
                random_state=seed + int(row.candidate),
            ))
    return seed, scores


def save_checkpoint(path, config, results):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with open(temporary, "wb") as file:
        pickle.dump({"config": config, "results": results}, file)
    os.replace(temporary, path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--seed-start", type=int, default=1)
    parser.add_argument("--jobs", type=int, default=min(10, os.cpu_count() or 1))
    parser.add_argument("--checkpoint-every", type=int, default=10)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument(
        "--output-dir", type=Path, default=HERE / "results"
    )
    args = parser.parse_args()

    config = {
        "trials": args.trials,
        "seed_start": args.seed_start,
        "score_x_tolerance": 1e-4,
        "scikit_learn": sklearn.__version__,
        "generator_sha256": file_sha256(HERE / "generate_simulation_trial.R"),
        "score_x_sha256": file_sha256(HERE / "simulation_score_x.py"),
    }
    seeds = list(range(args.seed_start, args.seed_start + args.trials))
    checkpoint = args.output_dir / "simulation_score_x_checkpoint.pkl"
    results = {}
    if checkpoint.exists():
        with open(checkpoint, "rb") as file:
            saved = pickle.load(file)
        if saved["config"] != config:
            raise ValueError("Checkpoint configuration differs from this run.")
        results = saved["results"]

    pending = [seed for seed in seeds if seed not in results]
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(run_trial, seed, args.rscript): seed
            for seed in pending
        }
        for completed, future in enumerate(as_completed(futures), 1):
            seed, scores = future.result()
            results[seed] = scores
            if completed % args.checkpoint_every == 0 or completed == len(pending):
                save_checkpoint(checkpoint, config, results)
                print(f"Completed {len(results)}/{len(seeds)} trials", flush=True)

    ordered = np.asarray([results[seed] for seed in seeds])
    if ordered.shape != (args.trials, 15) or not np.isfinite(ordered).all():
        raise ValueError("Score X result matrix is invalid.")
    save_checkpoint(
        args.output_dir / "simulation_score_x_trials.pkl", config, results
    )

    kl = pd.read_csv(HERE / "simulation_kl_reference.csv")["kl"]
    pd.DataFrame({
        "kl": kl,
        "score_x": ordered.mean(axis=0),
    }).to_csv(HERE / "kl_score_x_avg.csv", index=False)

    subprocess.run(
        [args.rscript, str(HERE / "simulation_fig_comparison.R")],
        check=True, cwd=HERE,
    )
    data = pd.read_csv(HERE / "simulation_plot_data.csv")
    figure = plot_simulation(data, ROOT / "figures" / "sim.png")
    plt.close(figure)


if __name__ == "__main__":
    main()
