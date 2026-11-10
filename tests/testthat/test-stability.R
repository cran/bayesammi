test_that("the compiled stability sampler draws the same numbers as the R reference", {
  ns <- asNamespace("bayesammi")
  cd <- ns$.bam_prep("t", Maiz, quote(y), quote(entry), quote(site), quote(rep), environment(), allow_missing = TRUE)
  for (code in 0:2) {
    ls <- ns$.bam_stab_ls(cd, code)
    pr <- ns$.bam_stab_fill_priors(bam_stab_priors(), cd, ls)
    st <- ns$.bam_stab_start(ls, 1L, code)
    set.seed(3)
    r <- ns$.bam_stab_chain(cd$S, cd$Q, cd$N * 1, code, pr, st, 40, rep(1L, 40))
    set.seed(3)
    cc <- ns$.bam_stab_cpp(cd$S, cd$Q, cd$N * 1, code, pr, st, 40L, rep(1L, 40), cd$data$y, cd$data$cell - 1L, FALSE)
    for (nm in c("mu", "sigma2", "g", "h", "b", "d")) expect_equal(drop(cc[[nm]]), drop(r[[nm]]), tolerance = 1e-9)
    expect_equal(drop(cc$sd2), drop(r$sd2), tolerance = 1e-9)
  }
})

sim_er <- function(seed, drop_cells = TRUE) {
  set.seed(seed)
  G <- 10; E <- 15; B <- 3
  g <- as.vector(scale(rnorm(G, 0, 300), scale = FALSE))
  h <- as.vector(scale(rnorm(E, 0, 1000), scale = FALSE))
  b <- as.vector(scale(runif(G, -0.5, 0.5), scale = FALSE))
  sd2 <- rep(c(2e4, 2e5), length.out = G)
  D <- matrix(rnorm(G * E, 0, sqrt(sd2)), G, E)
  M <- 5000 + outer(g, rep(1, E)) + outer(1 + b, h) + D
  d <- expand.grid(rep = 1:B, gen = 1:G, env = 1:E)
  d$y <- M[cbind(d$gen, d$env)] + rnorm(nrow(d), 0, 300)
  if (drop_cells) {
    d <- d[-sample(nrow(d), 25), ]
    d <- d[!(d$gen == 3 & d$env == 4), ]
  }
  list(data = d, slope = 1 + b, sd2 = sd2, g = g, h = h)
}

test_that("Eberhart-Russell: slopes and deviation variances are recovered, also with missing data", {
  s <- sim_er(21)
  fit <- bayes_stability(.data = s$data, .y = y, .gen = gen, .env = env, .rep = rep, .model = "ER", .nIter = 2500,
                         .burnin = 500, .chains = 2, .seed = 4)
  st <- fit$stability
  expect_gte(mean(st$slope_lower <= s$slope & s$slope <= st$slope_upper), 0.8)
  expect_gt(cor(st$slope, s$slope), 0.9)
  expect_gt(mean(st$dev_var[s$sd2 == 2e5]), 3 * mean(st$dev_var[s$sd2 == 2e4]))
  expect_lt(max(fit$summary[block %in% c("mu", "sigma2", "slope"), rhat]), 1.1)
  expect_equal(attr(fit, "meta")$missing, 1L)
  pr <- bam_predict(fit)
  expect_equal(nrow(pr$cells[n == 0]), 1L)
  expect_s3_class(fit$draws, "mcmc.list")
  expect_equal(dim(fit$deviations), c(5000L, 10L, 15L))
  expect_output(print(fit), "Sources:")
  for (w in c("slopes", "lines", "deviations", "trace", "density")) expect_s3_class(plot(fit, which = w), "ggplot")
  it <- bam_interpret(fit)
  expect_true(nzchar(it$headline))
  expect_true(length(it$sources) >= 3L)
})

test_that("Finlay-Wilkinson and Shukla models", {
  s <- sim_er(5, drop_cells = FALSE)
  fw <- bayes_stability(.data = s$data, .y = y, .gen = gen, .env = env, .rep = rep, .model = "FW", .nIter = 1500,
                        .burnin = 300, .chains = 1, .seed = 2)
  expect_equal(sum(fw$summary$block == "sigma2_d"), 1L)
  expect_gt(cor(fw$stability$slope, s$slope), 0.9)
  sh <- bayes_stability(.data = s$data, .y = y, .gen = gen, .env = env, .rep = rep, .model = "Shukla",
                        .nIter = 1500, .burnin = 300, .chains = 1, .seed = 2)
  expect_false("slope" %in% names(sh$stability))
  expect_false(any(sh$summary$block == "slope"))
  expect_error(plot(sh, which = "slopes"), "no slopes")
  expect_s3_class(plot(sh, which = "deviations"), "ggplot")
  # the safety-first variance of the Shukla model is the variance of the index plus the stability variance
  expect_true(all(sh$stability$safety_first < sh$stability$mean))
})

test_that("decisions, model comparison and checks accept stability fits", {
  fs <- bayes_stability(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 500, .burnin = 100,
                        .seed = 1)
  fa <- bayes_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 500, .burnin = 100,
                   .seed = 1)
  cp <- bam_compare(AMMI = fa, ER = fs)
  expect_equal(nrow(cp$table), 2L)
  expect_true(is.na(cp$table[model == "ER", terms]))
  expect_s3_class(bam_ppc(fs, .draws = 50, .seed = 1), "bam_ppc")
  expect_s3_class(bam_superiority(fs), "bam_superiority")
  expect_s3_class(bam_ranks(fs), "bam_ranks")
  expect_s3_class(bam_winners(fs), "bam_winners")
  expect_error(bam_stability(fs), "bilinear terms")
  expect_error(bam_ammi_indices(fs), "bilinear terms")
})

test_that("posterior stability indices agree with their formulas", {
  fa <- bayes_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 300, .burnin = 100,
                   .seed = 3)
  ind <- bam_indices(fa, .prob = 0.25, .alpha = 0.05)
  # the indices of the posterior mean draw by draw: check one draw by hand
  d <- bayesammi:::.bam_eta(fa)
  E1 <- d$eta[1L, , ]
  ge <- E1 - outer(rowMeans(E1), colMeans(E1), "+") + mean(E1)
  W <- rowSums(ge^2)
  M <- apply(E1, 2L, max)
  P <- rowSums((E1 - matrix(M, 9, 20, byrow = TRUE))^2) / 40
  Z <- 100 * E1 / matrix(colMeans(E1), 9, 20, byrow = TRUE)
  A <- rowMeans(Z) - qnorm(0.75) * apply(Z, 1L, sd)
  SF <- rowMeans(E1) - qnorm(0.95) * apply(E1, 1L, sd)
  Sh <- 9 / (7 * 19) * W - sum(W) / (8 * 7 * 19)
  expect_equal(unname(rowSums(d$interaction[1L, , ]^2)), unname(W), tolerance = 1e-10)
  expect_equal(sum(ind$indices[index == "wricke", mean] > 0), 9L)
  expect_true(all(c("wricke", "shukla", "lin_binns", "annicchiarico", "safety_first") %in% ind$indices$index))
  expect_equal(unname(colSums(matrix(ind$indices$prob_best, 9))), rep(1, 5), tolerance = 1e-10)
  # per-draw formulas against the function on a one-draw fit
  f1 <- fa
  f1$draws <- coda::mcmc.list(coda::mcmc(fa$draws[[1L]][c(1L, 1L), , drop = FALSE]))
  i1 <- bam_indices(f1, .prob = 0.25, .alpha = 0.05)$indices
  expect_equal(i1[index == "wricke", mean], unname(W), tolerance = 1e-8)
  expect_equal(i1[index == "shukla", mean], unname(Sh), tolerance = 1e-8)
  expect_equal(i1[index == "lin_binns", mean], unname(P), tolerance = 1e-8)
  expect_equal(i1[index == "annicchiarico", mean], unname(A), tolerance = 1e-8)
  expect_equal(i1[index == "safety_first", mean], unname(SF), tolerance = 1e-8)
  ai <- bam_ammi_indices(fa, .weight_y = 50)
  expect_equal(ai$indices[index == "ASV", mean], bam_stability(fa)$stability[order(as.integer(Gen)), ASV],
               tolerance = 1e-10)
  expect_true(all(ai$indices[index == "WAASY", mean] >= 0 & ai$indices[index == "WAASY", mean] <= 100))
  expect_output(print(ind), "Sources:")
  expect_s3_class(plot(ind), "ggplot")
  expect_s3_class(plot(ai), "ggplot")
  expect_true(nzchar(bam_interpret(ai)$headline))
})

test_that("every method cited by a result is in the reference registry", {
  refs <- bayesammi:::.bam_refs_list()
  src <- get(".bam_sources_of", asNamespace("bayesammi"))
  body_keys <- unique(regmatches(deparse(body(src)), gregexpr("\"[A-Za-z]+[0-9]{4}\"", deparse(body(src)))))
  keys <- gsub("\"", "", unlist(body_keys))
  expect_true(all(keys %in% names(refs)))
  b <- bam_cite()
  expect_equal(length(unclass(b)), length(refs))
  expect_output(print(bam_cite(NULL, style = "bibtex")), "@Article")
})
