# Empirical experiment 

Use Python 3.12.1 and install requirements in `requirements.txt`.

Figures are written to `figures/`. The plotting code creates this directory
when it does not already exist.

Use `acs.ipynb` and `nhanes_chol.ipynb` to reproduce the manuscript figures
and tables.

Run the main ACS experiment:

```bash
python -m dataset_ranking.acs_unequal_source_sizes \
  --trials 1000 --seed-start 123 --jobs 8 --checkpoint-every 4
```

Run the main NHANES experiment:

```bash
python -m dataset_ranking.experiment_equal_source_sizes \
  --dataset nhanes --trials 1000 --seed-start 12000 \
  --jobs 8 --checkpoint-every 4
```

Score X uses the covariate-only implementation from
`the-chen-lab/data-addition-dilemma` at commit
`279777ff5ab8757b5da7430a788bf10b08014522`
