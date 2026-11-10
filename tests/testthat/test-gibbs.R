test_that("von Mises-Fisher draws have the expected mean resultant length and direction", {
  set.seed(1)
  for (p in c(2L, 5L, 18L)) {
    mu <- rep(1, p) / sqrt(p)
    for (kap in c(0.5, 5, 50)) {
      X <- t(replicate(4000, bayesammi:::.bam_rvmf(kap * mu)))
      expect_equal(rowSums(X^2), rep(1, 4000))
      r <- drop(X %*% mu)
      A <- besselI(kap, p / 2, expon.scaled = TRUE) / besselI(kap, p / 2 - 1, expon.scaled = TRUE)
      expect_equal(mean(r), A, tolerance = 4 * sqrt(var(r) / 4000) / A + 1e-3, label = paste(p, kap))
      m <- colMeans(X)
      if (kap >= 5) expect_gt(sum(m * mu) / sqrt(sum(m^2)), 0.95)
    }
  }
  # one dimension: P(+1) = 1 / (1 + exp(-2 v))
  x <- replicate(20000, bayesammi:::.bam_rvmf(0.4))
  expect_equal(mean(x == 1), stats::plogis(0.8), tolerance = 0.02)
})

test_that("truncated normal draws stay in their interval and have the right mean", {
  set.seed(2)
  for (b in list(c(0, Inf), c(-1, 1), c(3, 4), c(8, Inf), c(-Inf, -9))) {
    x <- replicate(4000, bayesammi:::.bam_rtnorm(0, 1, b[1], b[2]))
    expect_true(all(x >= b[1] & x <= b[2]))
    m <- (dnorm(b[1]) - dnorm(b[2])) / (pnorm(b[2]) - pnorm(b[1]))
    if (is.finite(m)) expect_equal(mean(x), m, tolerance = 0.05 + 0.05 * abs(m))
  }
})

test_that("the sampler recovers the parameters of simulated data", {
  set.seed(3)
  G <- 8; E <- 10; B <- 3
  a <- bayesammi:::.bam_orth1(matrix(rnorm(G * 2), G)); g <- bayesammi:::.bam_orth1(matrix(rnorm(E * 2), E))
  lam <- c(60, 30); tau <- rnorm(G, 0, 5); tau <- tau - mean(tau); del <- rnorm(E, 0, 10); del <- del - mean(del)
  M <- 100 + outer(tau, del, "+") + a %*% (lam * t(g))
  d <- expand.grid(rep = 1:B, gen = 1:G, env = 1:E)
  d$y <- M[cbind(d$gen, d$env)] + rnorm(nrow(d), 0, 4)
  fit <- bayes_ammi(.data = d, .y = y, .gen = gen, .env = env, .rep = rep, .nIter = 3000, .k = 2, .burnin = 500,
                    .chains = 2, .seed = 4)
  s <- fit$summary
  expect_true(all(s[block %in% c("mu", "sigma2", "lambda", "tau", "delta"), rhat] < 1.05))
  covered <- function(b, truth) mean(s[block == b, lower] <= truth & truth <= s[block == b, upper])
  expect_gte(covered("tau", tau), 0.75)
  expect_gte(covered("delta", del), 0.75)
  expect_true(s[parameter == "sigma2", lower] < 16 & 16 < s[parameter == "sigma2", upper])
  expect_true(all(abs(s[block == "lambda", mean] - lam) < 4 * s[block == "lambda", sd]))
  # scores up to the sign of each axis
  sg <- sign(colSums(fit$ls$alpha * a))
  est <- matrix(s[block == "alpha", mean], G)
  expect_gt(abs(cor(est[, 1], a[, 1])), 0.95)
})

test_that("draws keep the constraints; seeds reproduce runs and leave the session's random numbers alone", {
  set.seed(10)
  before <- .Random.seed
  f1 <- bayes_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 60, .seed = 7)
  expect_identical(.Random.seed, before)
  f2 <- bayes_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 60, .seed = 7)
  expect_equal(f1$mu1, f2$mu1)
  expect_result(f1, "bayes_ammi", components = c("mu1", "tau1", "tao1", "delta1", "lambdas1", "alphas1", "gammas1",
                                                   "alphas0", "gammas0", "summary", "scores", "draws",
                                                   "diagnostics", "ls", "settings"))
  expect_equal(names(f1)[1:9], c("mu1", "tau1", "tao1", "delta1", "lambdas1", "alphas1", "gammas1", "alphas0",
                                 "gammas0"))
  A <- f1$alphas1
  for (it in c(1, 30, 60)) {
    a <- A[A[, ".nIter"] == it, c("alphas1", "alphas2")]
    expect_equal(crossprod(a), diag(2), tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(colSums(a), c(0, 0), tolerance = 1e-8, ignore_attr = TRUE)
  }
  expect_true(all(f1$lambdas1$lambdas1 >= f1$lambdas1$lambdas2))
  expect_equal(rowSums(f1$tao1), rep(0, 60), tolerance = 1e-8)
  expect_equal(nrow(f1$alphas1), 60 * 9)
  # thinning and burn-in: kept draws
  f3 <- bayes_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 100, .burnin = 20,
                   .thin = 5, .chains = 2, .seed = 1)
  expect_equal(nrow(f3$mu1), 40)
  expect_equal(coda::niter(f3$draws), 20)
  expect_equal(coda::nchain(f3$draws), 2)
  expect_false(anyNA(f3$summary$rhat))
  expect_output(print(f3), "Posterior means")
  for (w in c("trace", "density", "biplot", "effects")) expect_s3_class(plot(f3, which = w), "ggplot")
  expect_error(plot(f3, which = "trace", pars = "nope"), "unknown")
})

test_that("priors", {
  p <- bam_priors(nu = 3, s2 = 2)
  expect_s3_class(p, "bam_priors")
  expect_output(print(p), "Inv-chi")
  expect_error(bam_priors(mu_var = -1))
  # an informative prior on lambda pulls it towards the prior mean
  f <- bayes_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 400, .seed = 1,
                  .priors = bam_priors(lambda_mean = 3000, lambda_var = 1e3))
  expect_lt(abs(mean(f$lambdas1$lambdas1) - 3000), 600)
})
