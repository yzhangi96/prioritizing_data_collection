# Simulation figure 

Run from the repository root:

```bash
python simulation/simulation_reproduce.py --trials 1000 --jobs 10
```

This generates every Score X trial from the original R/calinf seeds,
checkpoints under `simulation/results/`, writes `kl_score_x_avg.csv`, combines
it with the stored 1000-trial MSE/DUC results, and creates the manuscript figure
at `figures/sim.png`. Score X uses
`LogisticRegression(max_iter=10000, tol=1e-4)`.

The run with Python dependencies from `requirements.txt`, R
4.4.x, and `calinf` 0.1.0 at commit `5092dd8`. Install the optional R
dependencies used only for KDE, Score X, or parquet export with:

```r
install.packages(c("arrow", "ks", "remotes", "reticulate"))
remotes::install_github("rothenhaeusler/calinf@5092dd8")
```
