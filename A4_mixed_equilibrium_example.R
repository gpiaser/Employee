## =====================================================================
## Employee Ownership and Entrenchment Strategies
## Script S2. Appendix A.4: numerical example of Proposition A1 (n = 3),
## and the separating equilibrium with a randomizing employee discussed
## after Corollary 1 (Section 5).
##
## Binary effort model with u(x) = sqrt(x), Wd = 10, Ws = 1, r0 = 0.02,
## mu^{L,B} = 0, mu^{H,B} = 0.5, mu^{L,G} = 0.3, mu^{H,G} = 0.8,
## psi^L = 0, psi^H = 0.15, pG = 0.5, epsilon = 0.01. The dismissal
## threshold Vm is the midpoint of the interval ] V0(c*, H, B),
## V0(c*_B, H, B) [ on which the equilibria of Proposition A1 exist.
##
## The return r is degenerate at 0 in the paper's example, which is a
## limiting illustration since a degenerate return has no density. The
## script also solves the example with the triangular density
## f(r) = (a - |r|) / a^2 on [-a, a], a = 0.01, and reports the
## differences. Full-precision results are written to
## A4_example_output.csv, and those of the separating equilibrium with a
## randomizing employee to Corollary1_counterexample.csv.
## Numerical methods: thresholds and beliefs are solved by root finding
## (uniroot, tolerance 1e-14); expectations under the triangular density
## use a composite Simpson rule with 801 nodes.
## =====================================================================

base <- list(Wd = 10, Ws = 1, r0 = 0.02,
             mu = c(LB = 0, HB = 0.5, LG = 0.3, HG = 0.8),
             psiL = 0, psiH = 0.15, pG = 0.5, eps = 0.01)
u <- function(x) sqrt(x)

## Quadrature rules for the distribution of r: nodes and weights.
quad_degenerate <- list(r = 0, w = 1)
quad_triangular <- function(a, n = 801) {
  r <- seq(-a, a, length.out = n)
  h <- r[2] - r[1]
  simpson <- c(1, rep(c(4, 2), length.out = n - 2), 1)
  w <- h / 3 * simpson * (a - abs(r)) / a^2
  list(r = r, w = w / sum(w))
}

solve_example <- function(q, p = base) {
  y <- p$Ws * (1 + p$r0)
  rho <- 1 + sum(q$w * q$r)
  omega <- function(c, k) {
    muH <- p$mu[[paste0("H", k)]]; muL <- p$mu[[paste0("L", k)]]
    sum(q$w * (u(y + c * p$Ws * (1 + q$r + muH)) -
               u(y + c * p$Ws * (1 + q$r + muL))))
  }
  dpsi <- p$psiH - p$psiL
  Omega <- function(c, qq) qq * omega(c, "G") + (1 - qq) * omega(c, "B") - dpsi
  solve_c <- function(qq) uniroot(function(c) Omega(c, qq),
                                  c(1e-9, p$Wd / p$Ws), tol = 1e-14)$root
  cB <- solve_c(0); cs <- solve_c(p$pG); cG <- solve_c(1)
  A  <- function(c) p$Wd - c * p$Ws
  V0 <- function(c, j, k) A(c) * (rho + p$mu[[paste0(j, k)]])
  Vm <- (V0(cs, "H", "B") + V0(cB, "H", "B")) / 2
  chat <- (p$Wd - Vm / (rho + p$mu[["HB"]])) / p$Ws
  Phi <- A(cG)

  ## Condition (22), Proposition A1 and participation conditions
  checks <- c(
    cond22_G      = V0(cG, "H", "G") - Vm,          # > 0
    cond22_B      = Vm - V0(cs, "H", "B"),          # >= 0
    lowG_dismiss  = Vm - V0(0, "L", "G"),           # >= 0
    B_viable      = V0(cB, "H", "B") - Vm,          # > 0
    eps_G         = V0(cG, "H", "G") - p$eps,       # > 0
    eps_B         = Phi * (rho + p$mu[["HB"]]) - p$eps)  # >= 0
  stopifnot(all(checks[c("cond22_G", "B_viable", "eps_G")] > 0),
            all(checks[c("cond22_B", "lowG_dismiss", "eps_B")] >= 0))

  ci <- cB + c(0.15, 0.50, 0.85) * (chat - cB)
  qi <- vapply(ci, function(c) uniroot(function(qq) Omega(c, qq),
                                       c(0, 1), tol = 1e-14)$root, numeric(1))
  pii <- Phi / A(ci)
  lambda <- (1 - p$pG) / p$pG * qi / (1 - qi)
  PM_B <- rep(1 / 3, 3); PM_G <- lambda * PM_B
  stopifnot(sum(PM_G) < 1, all(pii < 1), all(qi < p$pG))

  ## Residual checks: employee indifference at each c_i, equal payoffs
  ## of the type B manager across his support, martingale property.
  resid_employee <- vapply(seq_along(ci), function(i) Omega(ci[i], qi[i]),
                           numeric(1))
  payB <- pii * V0(ci, "H", "B")
  post <- c(qi, 1)
  prob <- c(p$pG * PM_G + (1 - p$pG) * PM_B, p$pG * (1 - sum(PM_G)))
  list(rho = rho, cB = cB, cs = cs, cG = cG, Vm = Vm, chat = chat,
       Phi = Phi, ci = ci, qi = qi, pii = pii, PM_B = PM_B, PM_G = PM_G,
       UB = Phi * (rho + p$mu[["HB"]]), UG = Phi * (rho + p$mu[["HG"]]),
       checks = checks,
       max_resid_employee = max(abs(resid_employee)),
       spread_payB = diff(range(payB)),
       avg_posterior = sum(post * prob))
}

deg <- solve_example(quad_degenerate)
tri <- solve_example(quad_triangular(0.01))

cat(sprintf("Degenerate r: c*_B = %.9f, c* = %.9f, c*_G = %.9f\n",
            deg$cB, deg$cs, deg$cG))
cat(sprintf("              Vm = %.9f, chat_B = %.9f, Phi* = %.9f\n",
            deg$Vm, deg$chat, deg$Phi))
tab <- data.frame(contract = c("c1", "c2", "c3", "c*_G"),
                  c = c(deg$ci, deg$cG), pi = c(deg$pii, 1),
                  q = c(deg$qi, 1), PM_B = c(deg$PM_B, 0),
                  PM_G = c(deg$PM_G, 1 - sum(deg$PM_G)))
print(format(tab, digits = 9), row.names = FALSE)
cat(sprintf("UB = %.9f, UG = %.9f, average posterior = %.12f\n",
            deg$UB, deg$UG, deg$avg_posterior))
cat(sprintf("Max employee indifference residual = %.2e; spread of type B payoffs = %.2e\n",
            deg$max_resid_employee, deg$spread_payB))
cat("Condition margins (all must be >= 0, strict where required):\n")
print(signif(deg$checks, 6))

cat(sprintf("\nTriangular r (a = 0.01): c*_B = %.9f, c* = %.9f, c*_G = %.9f\n",
            tri$cB, tri$cs, tri$cG))
diffs <- c(cB = tri$cB - deg$cB, cs = tri$cs - deg$cs, cG = tri$cG - deg$cG,
           ci = tri$ci - deg$ci, pi = tri$pii - deg$pii, q = tri$qi - deg$qi,
           UB = tri$UB - deg$UB, UG = tri$UG - deg$UG)
cat(sprintf("Largest absolute difference with the degenerate case: %.2e\n",
            max(abs(diffs))))

out <- data.frame(
  quantity = c("c*_B", "c*", "c*_G", "Vm", "chat_B", "Phi*",
               paste0("c", 1:3), paste0("pi(c", 1:3, ")"),
               paste0("q(c", 1:3, ")"), paste0("PM_G(c", 1:3, ")"),
               "PM_G(c*_G)", "UB", "UG"),
  degenerate = c(deg$cB, deg$cs, deg$cG, deg$Vm, deg$chat, deg$Phi,
                 deg$ci, deg$pii, deg$qi, deg$PM_G, 1 - sum(deg$PM_G),
                 deg$UB, deg$UG),
  triangular = c(tri$cB, tri$cs, tri$cG, tri$Vm, tri$chat, tri$Phi,
                 tri$ci, tri$pii, tri$qi, tri$PM_G, 1 - sum(tri$PM_G),
                 tri$UB, tri$UG))
write.csv(out, "A4_example_output.csv", row.names = FALSE)

## =====================================================================
## Separating equilibrium with a randomizing employee in the regime
## where the threshold does not bind at c* (discussion after
## Corollary 1). Same primitives, Vm = 13.61, triangular density with
## a = 0.01 (and the degenerate case for comparison).
## Profile: type G offers c*_G, type B offers c*_B, the employee works
## with probability one at c*_G and with probability
## pi = A(c*_G) / A(c*_B) at c*_B, and attributes every other contract
## to a type G manager.
## =====================================================================
separating_mixing <- function(q, Vm = 13.61, p = base) {
  y <- p$Ws * (1 + p$r0)
  rho <- 1 + sum(q$w * q$r)
  omega <- function(c, k) {
    muH <- p$mu[[paste0("H", k)]]; muL <- p$mu[[paste0("L", k)]]
    sum(q$w * (u(y + c * p$Ws * (1 + q$r + muH)) -
               u(y + c * p$Ws * (1 + q$r + muL))))
  }
  Omega <- function(c, qq) qq * omega(c, "G") + (1 - qq) * omega(c, "B") -
    (p$psiH - p$psiL)
  solve_c <- function(qq) uniroot(function(c) Omega(c, qq),
                                  c(1e-9, p$Wd / p$Ws), tol = 1e-14)$root
  cB <- solve_c(0); cs <- solve_c(p$pG); cG <- solve_c(1)
  A  <- function(c) p$Wd - c * p$Ws
  V0 <- function(c, j, k) A(c) * (rho + p$mu[[paste0(j, k)]])
  chat <- (p$Wd - Vm / (rho + p$mu[["HB"]])) / p$Ws
  ## Conditions: c* < chat_B <= c*_G, V0(c*_G,H,G) > Vm >= V0(0,L,G)
  cond <- c(nonbinding_at_cstar = V0(cs, "H", "B") - Vm,       # > 0
            binding_at_cGstar   = Vm - V0(cG, "H", "B"),       # >= 0
            G_viable            = V0(cG, "H", "G") - Vm,       # > 0
            low_effort_zero     = Vm - V0(0, "L", "G"))        # >= 0
  stopifnot(cond[["nonbinding_at_cstar"]] > 0, cond[["binding_at_cGstar"]] >= 0,
            cond[["G_viable"]] > 0, cond[["low_effort_zero"]] >= 0)
  pi_B <- A(cG) / A(cB)
  UB <- pi_B * V0(cB, "H", "B"); UG <- V0(cG, "H", "G")
  ## Residuals: employee indifference at c*_B and c*_G, type G
  ## indifference between c*_G and c*_B.
  res <- c(employee_at_cB = Omega(cB, 0), employee_at_cG = Omega(cG, 1),
           G_indifference = pi_B * V0(cB, "H", "G") - UG)
  ## A deviation slightly above c* is profitable for both types under
  ## the prior, so the equilibrium does not survive Definition 4.
  cdev <- (cs + min(chat, cG)) / 2
  gains <- c(G = V0(cdev, "H", "G") - UG, B = V0(cdev, "H", "B") - UB)
  list(cB = cB, cs = cs, cG = cG, chat = chat, pi_B = pi_B, UB = UB,
       UG = UG, mimic_B = V0(cG, "H", "B"), cond = cond, res = res,
       cdev = cdev, gains = gains)
}

cat("\nSeparating equilibrium with a randomizing employee (Vm = 13.61)\n")
cex_rows <- list()
for (nm in c("triangular", "degenerate")) {
  q <- if (nm == "triangular") quad_triangular(0.01) else quad_degenerate
  z <- separating_mixing(q)
  cex_rows[[nm]] <- data.frame(
    distribution = nm, Vm = 13.61, cB_star = z$cB, c_star = z$cs,
    cG_star = z$cG, chat_B = z$chat, pi_cB = z$pi_B, UB = z$UB, UG = z$UG,
    B_mimic_gross = z$mimic_B,
    margin_nonbinding_at_cstar = z$cond[["nonbinding_at_cstar"]],
    margin_binding_at_cGstar = z$cond[["binding_at_cGstar"]],
    margin_G_viable = z$cond[["G_viable"]],
    margin_low_effort_zero = z$cond[["low_effort_zero"]],
    resid_employee_cB = z$res[["employee_at_cB"]],
    resid_employee_cG = z$res[["employee_at_cG"]],
    resid_G_indifference = z$res[["G_indifference"]],
    deviation_contract = z$cdev, deviation_gain_G = z$gains[["G"]],
    deviation_gain_B = z$gains[["B"]])
  cat(sprintf("%s: c*_B = %.9f, c* = %.9f, c*_G = %.9f, chat_B = %.6f\n",
              nm, z$cB, z$cs, z$cG, z$chat))
  cat(sprintf("   pi(c*_B) = %.9f, UB = %.9f, UG = %.9f, B mimicking c*_G: %.9f\n",
              z$pi_B, z$UB, z$UG, z$mimic_B))
  cat(sprintf("   largest residual = %.1e; gains from deviating to %.4f: G %+.4f, B %+.4f\n",
              max(abs(z$res)), z$cdev, z$gains[["G"]], z$gains[["B"]]))
}

write.csv(do.call(rbind, cex_rows), "Corollary1_counterexample.csv",
          row.names = FALSE)
