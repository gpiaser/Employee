## =====================================================================
## Employee Ownership and Entrenchment Strategies
## Script S3. Appendix A.5: magnitude of over-granting in the
## equilibria of Proposition A1
##
## Binary effort model with u(x) = sqrt(x), r degenerate at 0 (rho = 1),
## Wd = 10, Ws = 1, r0 = 0.02, mu^{L,B} = 0, psi^L = 0, and an effort
## increment Delta = mu^{H,k} - mu^{L,k} = 0.5 for both types, as in the
## example of Appendix A.4. The script sweeps the gap between types
## delta = mu^{L,G} - mu^{L,B}, the cost of effort psi^H and the prior pG.
##
## For each parameter point it computes the interval of dismissal
## thresholds Vm on which the equilibria of Proposition A1 exist under
## the maintained assumptions of Appendix A,
##   max{ V(c*,H,B), V(0,L,G), V(0,L,B) } < Vm
##     < min{ V(c*_B,H,B), V(c*_G,H,G) },
## and reports the relative excess (chat_B - c*_B) / c*_B of the largest
## contract a type B manager can offer over the complete information
## benchmark, both for Vm at the midpoint of that interval (the rule used
## in Appendix A.4) and in the limit where Vm approaches its lower bound.
## =====================================================================

base <- list(Wd = 10, Ws = 1, r0 = 0.02, rho = 1, muLB = 0, Delta = 0.5,
             psiL = 0)
u <- function(x) sqrt(x)

one_point <- function(delta, psiH, pG, b = base) {
  mu <- c(LB = b$muLB, HB = b$muLB + b$Delta,
          LG = b$muLB + delta, HG = b$muLB + delta + b$Delta)
  y <- b$Ws * (1 + b$r0)
  omega <- function(c, k) {
    u(y + c * b$Ws * (b$rho + mu[[paste0("H", k)]])) -
      u(y + c * b$Ws * (b$rho + mu[[paste0("L", k)]]))
  }
  dpsi <- psiH - b$psiL
  cmax <- b$Wd / b$Ws
  solve_c <- function(q) {
    f <- function(c) q * omega(c, "G") + (1 - q) * omega(c, "B") - dpsi
    if (f(cmax) <= 0) return(NA_real_)          # condition (16) fails
    uniroot(f, c(1e-12, cmax), tol = 1e-13)$root
  }
  cB <- solve_c(0); cs <- solve_c(pG); cG <- solve_c(1)
  out <- data.frame(delta = delta, psiH = psiH, pG = pG, exists = FALSE,
                    excess_mid = NA_real_, excess_sup = NA_real_)
  if (anyNA(c(cB, cs, cG))) return(out)
  A  <- function(c) b$Wd - c * b$Ws
  V0 <- function(c, j, k) A(c) * (b$rho + mu[[paste0(j, k)]])
  cbarG <- cmax * b$Delta / (b$rho + mu[["HG"]])
  if (!(cG < cbarG)) return(out)                  # Assumption 1 fails
  lo <- max(V0(cs, "H", "B"), V0(0, "L", "G"), V0(0, "L", "B"))
  hi <- min(V0(cB, "H", "B"), V0(cG, "H", "G"))
  if (!(lo < hi)) return(out)                     # no admissible Vm
  chat <- function(Vm) (b$Wd - Vm / (b$rho + mu[["HB"]])) / b$Ws
  out$exists <- TRUE
  out$excess_mid <- (chat((lo + hi) / 2) - cB) / cB
  out$excess_sup <- (chat(lo) - cB) / cB
  out
}

grid <- expand.grid(delta = c(0.1, 0.2, 0.3, 0.4, 0.5),
                    psiH  = c(0.05, 0.10, 0.15, 0.20, 0.25),
                    pG    = c(0.3, 0.5, 0.7))
res <- do.call(rbind, Map(one_point, grid$delta, grid$psiH, grid$pG))
write.csv(res, "A5_grid_results.csv", row.names = FALSE)

ok <- res[res$exists, ]
cat("Parameter points:", nrow(res), " with admissible thresholds:",
    nrow(ok), "\n")
cat(sprintf("Excess at the midpoint threshold: %.2f%% to %.2f%%\n",
            100 * min(ok$excess_mid), 100 * max(ok$excess_mid)))
cat(sprintf("Excess at the lowest admissible threshold: %.2f%% to %.2f%%\n",
            100 * min(ok$excess_sup), 100 * max(ok$excess_sup)))
cat("\nLargest excesses (midpoint rule)\n")
print(head(ok[order(-ok$excess_mid), ], 5), row.names = FALSE, digits = 3)
cat("\nSmallest excesses (midpoint rule)\n")
print(head(ok[order(ok$excess_mid), ], 5), row.names = FALSE, digits = 3)
cat("\nExample of Appendix A.4 (delta = 0.3, psiH = 0.15, pG = 0.5)\n")
print(one_point(0.3, 0.15, 0.5), row.names = FALSE, digits = 4)
