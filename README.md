# Employee Ownership and Entrenchment Strategies: replication code

This repository contains the R scripts that produce every numerical
result of the paper *Employee Ownership and Entrenchment Strategies*
(Gwenaël Piaser, IPAG Business School), submitted to *Games* (MDPI),
manuscript games-4527113.

## Contents

| File | Paper reference | What it computes |
|---|---|---|
| `R/S6_continuous_effort.R` | Script S1, Section 6, Table 2 | Separating equilibria with continuous effort: the good manager's contract, the bad manager's complete information optimum, and the set of contracts satisfying condition (30), for the sixteen specifications of the grid; robustness of Table 2 to a triangular return density |
| `R/A4_mixed_equilibrium_example.R` | Script S2, Appendix A.4 and Section 5 | The numerical example of Proposition A1, at full precision, with its robustness to a triangular return density; the separating equilibrium with a randomizing employee discussed after Corollary 1 |
| `R/A5_overgranting_range.R` | Script S3, Appendix A.5 | The magnitude of over-granting in the equilibria of Proposition A1 on a grid of 75 parameter points |
| `run_all.R` | | Runs the three scripts and writes their outputs to `results/` |
| `results/S6_results.csv` | | Output of Script S1 |
| `results/A4_example_output.csv` | | Output of Script S2, example of Appendix A.4 |
| `results/Corollary1_counterexample.csv` | | Output of Script S2, equilibrium discussed after Corollary 1 |
| `results/A5_grid_results.csv` | | Output of Script S3 |

## Requirements

Base R only; no package is needed. The results in `results/` were
produced with R 4.3.3.

## Usage

From the root of the repository:

```
Rscript run_all.R
```

Each script can also be run on its own from the `results/` folder, for
instance `Rscript ../R/A4_mixed_equilibrium_example.R`. Scripts S2 and
S3 run in a few seconds; Script S1 takes a few minutes.

## Numerical methods

The employee's incentive thresholds, beliefs and first-order
conditions are solved by root finding (`uniroot`, tolerances between
1e-12 and 1e-14). In Script S1, the managers' optima are located on a
grid of contracts on [0, Wd/Ws] and refined by local optimization
(`optimize`); the admissible set is computed on a grid of 5001
contracts and its endpoints are refined by root finding. Expectations
under the triangular density f(r) = (a - |r|) / a^2 on [-a, a],
a = 0.01, are computed with a composite Simpson rule. The return r is
degenerate at 0 in the reported results, which is a limiting
illustration; the scripts report the robustness of every result to the
triangular density.

## Checks reported

The output files contain, besides the results themselves, the residual
of the employee's first-order or indifference condition, the residuals
of the binding constraints at the endpoints of the admissible
intervals, the slack of the constraints at their midpoints, the number
of local maxima of each objective, the margins of all existence and
participation conditions, and a status distinguishing a nonempty
admissible set, an empty one, and a numerical failure.

## Citation

Piaser, G. Employee Ownership and Entrenchment Strategies. Submitted to
*Games* (MDPI).
