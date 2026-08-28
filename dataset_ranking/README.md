# Empirical experiment 

Use Python 3.12.1 and install requirements in `requirements.txt`.

Use `acs.ipynb` and `nhanes_chol.ipynb` to reproduce the manuscript figures
and tables.

Run an experiment with checkpointing from the repository root:

```bash
python dataset_ranking/experiment.py \
  --dataset acs --trials 1000 --seed-start 123 \
  --checkpoint-every 1 --output-dir results

python dataset_ranking/experiment.py \
  --dataset nhanes --trials 1000 --seed-start 12000 \
  --checkpoint-every 1 --output-dir results
```

Score X uses the covariate-only implementation from
`the-chen-lab/data-addition-dilemma` at commit
`279777ff5ab8757b5da7430a788bf10b08014522`