fit_maiz <- function(..., nIter = 800, chains = 1, seed = 1) {
  bayes_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = nIter, .burnin = 200,
             .chains = chains, .seed = seed, ...)
}

test_that("the compiled sampler draws the same numbers as the R reference sampler", {
  ns <- asNamespace("bayesammi")
  cd <- ns$.bam_prep("t", Maiz, quote(y), quote(entry), quote(site), quote(rep), environment(), allow_missing = TRUE)
  ls <- ns$.bam_ls(cd$Y, 2)
  s2 <- (cd$ssw + 4 * ls$rss) / 600
  pr <- bam_priors()
  pr$s2 <- s2
  st <- ns$.bam_start(ls, s2, 1)
  set.seed(5)
  r <- ns$.bam_chain(cd$Y, 4, cd$ssw, 2, pr, st, list(alpha = ls$alpha), 200, 1:200)
  set.seed(5)
  cc <- ns$.bam_run_chain(cd, 2, FALSE, FALSE, pr, st, ls$alpha, 200, rep(1L, 200), TRUE, FALSE)
  expect_equal(cc$scalar, r$scalar, tolerance = 1e-10)
  expect_equal(cc$alpha, r$alpha, tolerance = 1e-10)
  expect_equal(cc$gamma, r$gamma, tolerance = 1e-10)
})

test_that("WAIC and DIC agree with their definitions (and with loo)", {
  fit <- fit_maiz(nIter = 300)
  lay <- attr(fit, "meta")$layout
  all <- bayesammi:::.bam_all(fit)
  cell <- (Maiz$site - 1L) * 9L + Maiz$entry
  ll <- vapply(seq_len(nrow(all)), function(s) {
    eta <- vapply(1:20, function(j) bayesammi:::.bam_eta_env(all[s, , drop = FALSE], lay, j, 9), numeric(9))
    stats::dnorm(Maiz$y, eta[cell], sqrt(all[s, "sigma2"]), log = TRUE)
  }, numeric(nrow(Maiz)))
  ll <- t(ll)
  lppd <- sum(log(colMeans(exp(ll))))
  pw   <- sum(apply(ll, 2L, stats::var))
  expect_equal(fit$fit$lppd, lppd, tolerance = 1e-8)
  expect_equal(fit$fit$p_waic, pw, tolerance = 1e-8)
  expect_equal(fit$fit$Dbar, mean(-2 * rowSums(ll)), tolerance = 1e-8)
  skip_if_not_installed("loo")
  w <- suppressWarnings(loo::waic(ll))
  expect_equal(fit$fit$WAIC, w$estimates["waic", "Estimate"], tolerance = 1e-8)
})

test_that("environment-specific variances: convergence, sensible variances, better fit", {
  fh <- fit_maiz(.het = TRUE, nIter = 1500, chains = 2)
  s  <- fh$summary
  expect_equal(sum(s$block == "sigma2"), 20L)
  expect_lt(max(s$rhat[s$block %in% c("mu", "sigma2", "lambda", "tau", "delta")]), 1.1)
  # each variance is close to the within-site residual variance of the least-squares fit
  ls <- ge_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep)
  fitted <- ls$mu + outer(drop(ls$tao), drop(ls$delta), "+") + ls$alphas[, 1:2] %*% (ls$lambdas[1:2] * t(ls$gammas[, 1:2]))
  r  <- Maiz$y - fitted[cbind(Maiz$entry, Maiz$site)]
  v  <- tapply(r^2, Maiz$site, sum) / (36 - 9)
  expect_gt(cor(s[block == "sigma2", mean], v), 0.95)
  f0 <- fit_maiz(nIter = 1500, chains = 2)
  expect_lt(fh$fit$WAIC, f0$fit$WAIC)
  cp <- bam_compare(common = f0, het = fh)
  expect_equal(cp$table$model[1L], "het")
  expect_equal(cp$table$dWAIC_se[1L], 0)
  expect_equal(cp$table$dWAIC_se[2L], sqrt(nrow(Maiz) * var(f0$fit$waic_i - fh$fit$waic_i)))
  expect_equal(sum(fh$fit$waic_i), fh$fit$WAIC)
  expect_s3_class(plot(cp), "ggplot")
  # the posterior predictive check detects the unequal variances under a common variance
  expect_lt(bam_ppc(f0, .draws = 200, .seed = 1)$checks[statistic == "largest / smallest environment variance", p_value], 0.05)
  expect_gt(bam_ppc(fh, .draws = 200, .seed = 1)$checks[statistic == "largest / smallest environment variance", p_value], 0.05)
})

test_that("missing plots and cells: data augmentation recovers the truth", {
  set.seed(11)
  G <- 8; E <- 10; B <- 3
  a <- bayesammi:::.bam_orth1(matrix(rnorm(G * 2), G)); g <- bayesammi:::.bam_orth1(matrix(rnorm(E * 2), E))
  tau <- as.vector(scale(rnorm(G, 0, 5), scale = FALSE)); del <- as.vector(scale(rnorm(E, 0, 10), scale = FALSE))
  M <- 100 + outer(tau, del, "+") + a %*% (c(60, 30) * t(g))
  d <- expand.grid(rep = 1:B, gen = 1:G, env = 1:E)
  d$y <- M[cbind(d$gen, d$env)] + rnorm(nrow(d), 0, 4)
  drop <- sample(nrow(d), 30)
  dm <- d[-drop, ]
  dm <- dm[!(dm$gen == 2 & dm$env == 3) & !(dm$gen == 5 & dm$env == 7), ]
  fit <- bayes_ammi(.data = dm, .y = y, .gen = gen, .env = env, .rep = rep, .nIter = 2000, .burnin = 500,
                    .chains = 2, .seed = 2)
  expect_equal(attr(fit, "meta")$missing, 2L)
  s <- fit$summary
  expect_gte(mean(s[block == "tau", lower] <= tau & tau <= s[block == "tau", upper]), 0.75)
  expect_true(s[parameter == "sigma2", lower] < 16 & 16 < s[parameter == "sigma2", upper])
  pr <- bam_predict(fit)
  miss <- pr$cells[n == 0]
  expect_equal(nrow(miss), 2L)
  truth <- M[cbind(as.integer(as.character(miss$Gen)), as.integer(as.character(miss$Env)))]
  expect_true(all(abs(miss$mean - truth) < 4 * miss$sd))
  expect_s3_class(plot(pr), "ggplot")
  expect_output(print(fit), "empty cells")
})

test_that("GGE model: constraints, no genotype effects, winners", {
  fg <- fit_maiz(.model = "GGE", nIter = 600)
  expect_false(any(fg$summary$block == "tau"))
  expect_equal(ncol(fg$tao1), 0L)
  A <- fg$alphas1
  a <- A[A[, ".nIter"] == 10, c("alphas1", "alphas2")]
  expect_equal(colSums(a), c(0, 0), tolerance = 1e-8, ignore_attr = TRUE)
  g <- fg$gammas1
  g <- g[g[, ".nIter"] == 10, c("gammas1", "gammas2")]
  expect_equal(crossprod(g), diag(2), tolerance = 1e-8, ignore_attr = TRUE)
  expect_gt(abs(sum(colSums(g))), 1e-6)
  w <- bam_winners(fg)
  expect_equal(unname(colSums(w$probabilities)), rep(1, 20))
  expect_true(all(w$winners$Winner %in% as.character(1:9)))
  expect_s3_class(plot(w), "ggplot")
  bp <- biplots(model = fg, burnin = 0, thin = 1, pb = 0.05)
  expect_s3_class(bp, "bam_biplots")
})

test_that("probabilities of superior performance, ranks and stability are coherent", {
  fit <- fit_maiz(nIter = 800)
  sp <- bam_superiority(fit, .s = 0.25)
  ns <- attr(sp, "meta")$selected
  expect_equal(ns, 3L)
  expect_equal(sum(sp$marginal$prob_performance), ns)
  expect_equal(sum(sp$marginal$prob_stability), ns)
  expect_equal(unname(colSums(sp$conditional)), rep(ns, 20))
  P <- sp$pairwise_performance
  expect_equal(P[upper.tri(P)] + t(P)[upper.tri(P)], rep(1, sum(upper.tri(P))))
  expect_equal(sp$marginal$prob_joint, sp$marginal$prob_performance * sp$marginal$prob_stability)
  for (w in c("marginal", "conditional", "pairwise")) expect_s3_class(plot(sp, which = w), "ggplot")
  rk <- bam_ranks(fit)
  expect_equal(unname(rowSums(rk$probabilities)), rep(1, 9))
  expect_equal(unname(colSums(rk$probabilities)), rep(1, 9))
  expect_equal(rk$summary$Gen[1L], "6")
  expect_s3_class(bam_ranks(fit, .env = "12"), "bam_ranks")
  expect_error(bam_ranks(fit, .env = "nope"), "unknown environment")
  st <- bam_stability(fit)
  expect_true(all(st$stability$ASV >= st$stability$ASV_lower & st$stability$ASV <= st$stability$ASV_upper))
  expect_equal(sum(st$stability$prob_most_stable), 1)
  expect_s3_class(plot(st), "ggplot")
  for (x in list(sp, rk, st, bam_winners(fit), bam_predict(fit))) {
    txt <- bam_interpret(x)
    expect_s3_class(txt, "bam_interpretation")
    expect_true(nzchar(txt$headline))
  }
})

test_that("draws for posterior, bayesplot and coda", {
  fit <- fit_maiz(nIter = 300, chains = 2)
  expect_s3_class(coda::as.mcmc.list(fit), "mcmc.list")
  skip_if_not_installed("posterior")
  dr <- posterior::as_draws_df(fit)
  expect_equal(posterior::ndraws(dr), 600)
  expect_true("lambda[1]" %in% posterior::variables(dr))
  skip_if_not_installed("bayesplot")
  expect_s3_class(bayesplot::mcmc_trace(fit$draws, pars = "lambda[1]"), "ggplot")
})

test_that("with all terms the genotype scores move (exact pair moves) and match the R reference", {
  ns <- asNamespace("bayesammi")
  cd <- ns$.bam_prep("t", Maiz, quote(y), quote(entry), quote(site), quote(rep), environment(), allow_missing = TRUE)
  ls <- ns$.bam_ls(cd$Y, 8)
  pr <- bam_priors()
  pr$s2 <- 7e5
  st <- ns$.bam_start(ls, 7e5, 1)
  set.seed(5)
  r <- ns$.bam_chain(cd$Y, 4, cd$ssw, 8, pr, st, list(alpha = ls$alpha), 5, 1:5)
  set.seed(5)
  cc <- ns$.bam_run_chain(cd, 8, FALSE, FALSE, pr, st, ls$alpha, 5, rep(1L, 5), TRUE, FALSE)
  expect_equal(cc$alpha, r$alpha, tolerance = 1e-10)
  f <- fit_maiz(.k = 8, nIter = 400, chains = 2)
  expect_gt(min(f$summary[block == "alpha", sd]), 0.01)
  expect_lt(stats::median(f$summary[block == "alpha", rhat]), 1.05)
})
