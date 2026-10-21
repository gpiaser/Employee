## =====================================================================
## Employee Ownership and Entrenchment Strategies
## Script S1. Section 6 (continuous effort): separating equilibria and
## Table 2
##
## Model of Section 6 with
##   u(x)   = x^(1 - gamma) / (1 - gamma) for gamma != 1, u(x) = log(x)
##            for gamma = 1 (only u' enters the computations),
##   mu(e)  = alpha * sqrt(e),   psi(e) = beta * e^2 / 2,
##   mu^B = 0, mu^G = delta, Wd = 10, Ws = 1, r0 = 0.02, Vm = 0.
## The return r is degenerate at 0 in Table 2, a limiting illustration;
## the script also recomputes the four rows of Table 2 with the
## triangular density f(r) = (a - |r|) / a^2 on [-a, a], a = 0.01.
##
## For each specification the script computes
##   cG     the contract of the type G manager (Proposition 3);
##   cFI_B  the complete information optimum of the type B manager;
##   the set of contracts cB satisfying condition (30):
##          fB(cB) >= MB  and  gG(cB) <= UG,
## and reports a status (nonempty, empty or numerical failure), the
## largest relative residual of the employee's first-order condition,
## the residuals of the binding constraints at the interval endpoints,
## the number of local maxima of each objective, and the slack of both
## constraints at the midpoint of the admissible interval. No
## refinement is applied, and the prior pG plays no role. Results are
## written to S6_results.csv.
## Numerical methods: the employee's first-order condition is solved by
## root finding on log effort (uniroot, tolerance 1e-13); the managers'
## optima are located on a grid of 2001 contracts on [0, Wd/Ws] and
## refined by local optimization (optimize, tolerance 1e-12); the
## admissible set is computed on a grid of 5001 contracts and its
## endpoints are refined by root finding (tolerance 1e-12).
## Expectations under the triangular density use a composite Simpson
## rule (51 nodes in the robustness check).
## =====================================================================

par_base <- list(Wd = 10, Ws = 1, r0 = 0.02, muB = 0)

quad_degenerate <- list(r = 0, w = 1)
quad_triangular <- function(a, n = 101) {
  r <- seq(-a, a, length.out = n)
  h <- r[2] - r[1]
  simpson <- c(1, rep(c(4, 2), length.out = n - 2), 1)
  w <- h / 3 * simpson * (a - abs(r)) / a^2
  list(r = r, w = w / sum(w))
}

u_prime <- function(x, gamma) x^(-gamma)

## Employee's first-order condition and optimal effort e*(c, q).
foc <- function(e, c, q, s, quad, p = par_base) {
  mu_k <- c(p$muB + s$delta, p$muB)                    # (G, B)
  w_k  <- c(q, 1 - q)
  Eu <- vapply(1:2, function(k) {
    x <- p$Ws * (1 + p$r0) + c * p$Ws * (1 + quad$r + s$alpha * sqrt(e) + mu_k[k])
    sum(quad$w * u_prime(x, s$gamma))
  }, numeric(1))
  c * p$Ws * s$alpha / (2 * sqrt(e)) * sum(w_k * Eu) - s$beta * e
}
estar <- function(c, q, s, quad, p = par_base) {
  if (c <= 0) return(0)                                # e*(0, q) = 0
  exp(uniroot(function(le) foc(exp(le), c, q, s, quad, p),
              c(-60, 15), tol = 1e-13)$root)
}

## Gross payoff (equation (25)): V_k(c, e) = (Wd - c Ws)(rho + mu(e) + mu^k).
V <- function(c, e, k, s, quad, p = par_base) {
  rho <- 1 + sum(quad$w * quad$r)
  mu_k <- if (k == "G") p$muB + s$delta else p$muB
  (p$Wd - c * p$Ws) * (rho + s$alpha * sqrt(e) + mu_k)
}

argmax_c <- function(f, p = par_base, n = 4001) {
  grid <- seq(0, p$Wd / p$Ws, length.out = n)
  v <- vapply(grid, f, numeric(1))
  i <- which.max(v)
  o <- optimize(f, c(grid[max(1, i - 1)], grid[min(n, i + 1)]),
                maximum = TRUE, tol = 1e-12)
  c(arg = o$maximum, value = o$objective)
}

separating_set <- function(s, quad, p = par_base, n = 5001, n_arg = 2001) {
  fG1 <- function(c) V(c, estar(c, 1, s, quad, p), "G", s, quad, p)
  fB1 <- function(c) V(c, estar(c, 1, s, quad, p), "B", s, quad, p)
  fB0 <- function(c) V(c, estar(c, 0, s, quad, p), "B", s, quad, p)
  gG0 <- function(c) V(c, estar(c, 0, s, quad, p), "G", s, quad, p)

  G  <- argmax_c(fG1, p, n_arg); UG <- G[["value"]]; cG <- G[["arg"]]
  MB <- argmax_c(fB1, p, n_arg)[["value"]]
  cFI_B <- argmax_c(fB0, p, n_arg)[["arg"]]

  slack <- function(c) min(fB0(c) - MB, UG - gG0(c))
  grid <- seq(0, p$Wd / p$Ws, length.out = n)
  ok <- vapply(grid, slack, numeric(1)) >= 0

  ## Diagnostics: residual of the employee's first-order condition on a
  ## sample of contracts, and number of local maxima of each objective.
  sample_c <- seq(0.05, p$Wd / p$Ws - 0.05, length.out = 10)
  foc_resid <- max(vapply(sample_c, function(c) {
    e <- estar(c, 0, s, quad, p)
    abs(foc(e, c, 0, s, quad, p)) / (s$beta * e)
  }, numeric(1)))
  sub <- grid[seq(1, n, by = 5)]
  n_peaks <- function(f) { v <- vapply(sub, f, numeric(1))
                           sum(diff(sign(diff(v))) < 0) }

  base_out <- data.frame(cG = cG, cFI_B = cFI_B, UG = UG, MB = MB,
                         foc_resid = foc_resid,
                         peaks_G = n_peaks(fG1), peaks_FIB = n_peaks(fB0),
                         peaks_MB = n_peaks(fB1))
  if (!any(ok)) {
    return(cbind(base_out, status = "empty", cB_low = NA, cB_high = NA,
                 n_intervals = 0, end_resid = NA, mid_slack_fB = NA,
                 mid_slack_gG = NA))
  }
  refine <- function(a, b) uniroot(slack, c(a, b), tol = 1e-12)$root
  idx <- which(ok); i_lo <- min(idx); i_hi <- max(idx)
  lo <- if (i_lo > 1) refine(grid[i_lo - 1], grid[i_lo]) else grid[i_lo]
  hi <- if (i_hi < n) refine(grid[i_hi], grid[i_hi + 1]) else grid[i_hi]
  mid <- (lo + hi) / 2
  cbind(base_out, status = "nonempty", cB_low = lo, cB_high = hi,
        n_intervals = sum(rle(ok)$values),
        end_resid = max(abs(slack(lo)), abs(slack(hi))),
        mid_slack_fB = fB0(mid) - MB, mid_slack_gG = UG - gG0(mid))
}

safe_set <- function(s, quad, ...) {
  tryCatch(separating_set(s, quad, ...), error = function(err) {
    data.frame(cG = NA, cFI_B = NA, UG = NA, MB = NA, foc_resid = NA,
               peaks_G = NA, peaks_FIB = NA, peaks_MB = NA,
               status = paste("numerical failure:", conditionMessage(err)),
               cB_low = NA, cB_high = NA, n_intervals = NA, end_resid = NA,
               mid_slack_fB = NA, mid_slack_gG = NA)
  })
}

## ---------------------------------------------------------------------
## Sweep over the sixteen specifications, degenerate return
## ---------------------------------------------------------------------
specs <- expand.grid(alpha = c(0.4, 0.6), beta = c(0.3, 0.5),
                     delta = c(0.2, 0.4), gamma = c(0.5, 1.0))
res <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i)
  cbind(specs[i, ], safe_set(as.list(specs[i, ]), quad_degenerate))))
res$over_low  <- res$cB_low  - res$cFI_B
res$over_high <- res$cB_high - res$cFI_B
## Over-granting in percentage points of ownership: s(c) = c Ws / Wd.
res$over_high_pp <- 100 * res$over_high * par_base$Ws / par_base$Wd
write.csv(res, "S6_results.csv", row.names = FALSE)

ne <- res$status == "nonempty"
cat("Status:", paste(names(table(res$status)), table(res$status),
                     collapse = ", "), "\n")
cat("Every admissible set is a single interval:",
    all(res$n_intervals[ne] == 1), "\n")
cat("Every admissible contract lies strictly above cFI_B:",
    all(res$cB_low[ne] > res$cFI_B[ne]), "\n")
cat("Unique maximizer of each objective (one local maximum):",
    all(res$peaks_G == 1 & res$peaks_FIB == 1 & res$peaks_MB == 1), "\n")
cat(sprintf("Largest relative FOC residual: %.1e; largest endpoint residual: %.1e\n",
            max(res$foc_resid), max(res$end_resid[ne])))
cat(sprintf("Smallest midpoint slack: %.3e\n",
            min(c(res$mid_slack_fB[ne], res$mid_slack_gG[ne]))))
cat(sprintf("Over-granting ranges from +%.3f to +%.3f in c (up to %.2f percentage points of ownership)\n",
            min(res$over_low[ne]), max(res$over_high[ne]),
            max(res$over_high_pp[ne])))

## ---------------------------------------------------------------------
## Table 2 and robustness to a nondegenerate return
## ---------------------------------------------------------------------
table_specs <- data.frame(alpha = c(0.4, 0.6, 0.6, 0.6),
                          beta  = c(0.5, 0.3, 0.5, 0.3),
                          delta = c(0.4, 0.4, 0.4, 0.2),
                          gamma = c(0.5, 0.5, 1.0, 1.0))
tab <- merge(table_specs, res, sort = FALSE)
cat("\nTable 2 (degenerate return)\n")
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  cat(sprintf("$%.1f$ & $%.1f$ & $%.1f$ & $%.1f$ & $%.3f$ & $%.3f$ & $[%.3f, %.3f]$ \\\\\n",
              r$alpha, r$beta, r$delta, r$gamma, r$cG, r$cFI_B,
              r$cB_low, r$cB_high))
}
## Coarser grids suffice here: the endpoints are refined by root finding.
tri <- do.call(rbind, lapply(seq_len(nrow(table_specs)), function(i)
  safe_set(as.list(table_specs[i, ]), quad_triangular(0.01, n = 51),
           n = 1001, n_arg = 201)))
d <- abs(as.matrix(tri[, c("cG", "cFI_B", "cB_low", "cB_high")]) -
         as.matrix(tab[, c("cG", "cFI_B", "cB_low", "cB_high")]))
cat(sprintf("Triangular return (a = 0.01): largest change in Table 2 entries = %.1e\n",
            max(d)))
