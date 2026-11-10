test_that("the tables and effects agree with the released version", {
  f <- file.path(legacy_dir(), "helpers.rds")
  skip_if(!nzchar(legacy_dir()) || !file.exists(f), "no legacy snapshot (source tree only)")
  old <- readRDS(f)
  expect_equal(ge_mean(.data = Maiz, .y = y, .gen = entry, .env = site)$ge_means, old$ge_mean$ge_means,
               ignore_attr = TRUE)
  expect_equal(ge_mean(.data = Maiz, .y = y, .gen = entry, .env = site)$grand_mean, old$ge_mean$grand_mean)
  expect_equal(ge_var(.data = Maiz, .y = y, .gen = entry, .env = site)$ge_variances, old$ge_var$ge_variances,
               ignore_attr = TRUE)
  expect_equal(g_eff(.data = Maiz, .y = y, .gen = entry, .env = site), old$g_eff, ignore_attr = TRUE)
  expect_equal(e_eff(.data = Maiz, .y = y, .gen = entry, .env = site), old$e_eff, ignore_attr = TRUE)
  ge <- ge_eff(.data = Maiz, .y = y, .gen = entry, .env = site)
  expect_equal(ge$ge_effects, old$ge_eff$ge_effects, ignore_attr = TRUE)
  expect_equal(ge$ge_svd$d, old$ge_eff$ge_svd$d)
  expect_equal(matrix_k(5), old$matrix_k)
  expect_equal(orthnorm(cbind(1, 1:5, (1:5)^2)), old$orthnorm)
})

test_that("tables and effects from their definitions; columns by name, string or vector", {
  M <- tapply(Maiz$y, list(Maiz$entry, Maiz$site), mean)
  m <- ge_mean(.data = Maiz, .y = "y", .gen = "entry", .env = "site")
  expect_s3_class(m, "bam_ge_mean")
  expect_equal(m$ge_means, M, ignore_attr = TRUE)
  expect_equal(unname(drop(g_eff(.data = Maiz, .y = y, .gen = entry, .env = site))), unname(rowMeans(M) - mean(M)))
  expect_equal(unname(drop(e_eff(.data = Maiz, .y = Maiz$y, .gen = Maiz$entry, .env = Maiz$site))),
               unname(colMeans(M) - mean(M)))
  Z <- ge_eff(.data = Maiz, .y = y, .gen = entry, .env = site)$ge_effects
  expect_equal(unname(rowSums(Z)), rep(0, 9))
  expect_equal(unname(colSums(Z)), rep(0, 20))
  v <- ge_var(.data = Maiz, .y = y, .gen = entry, .env = site)$ge_variances
  expect_equal(unname(v), unname(tapply(Maiz$y, list(Maiz$entry, Maiz$site), var)))
})

test_that("ge_ammi: closed-form error variance equals the REML fit of ge_model()", {
  a <- ge_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep)
  expect_s3_class(a, "bam_ge_ammi")
  expect_equal(a$k, 8L)
  expect_equal(round(a$lambdas[1:2]), c(5923, 3070))
  expect_equal(a$tau, 1 / a$sigma2)
  skip_if_not_installed("lme4")
  fm <- ge_model(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep)
  expect_equal(a$sigma2, stats::sigma(fm)^2, tolerance = 1e-6)
  # the boundary case: replicate effects removed, so the REML replicate variance is zero
  d <- Maiz
  d$y <- d$y - ave(d$y, d$site, d$rep) + ave(d$y, d$site)
  a2 <- ge_ammi(.data = d, .y = y, .gen = entry, .env = site, .rep = rep)
  expect_equal(a2$sigma2, stats::sigma(ge_model(.data = d, .y = y, .gen = entry, .env = site, .rep = rep))^2,
               tolerance = 1e-4)
})

test_that("matrix_k and orthnorm", {
  for (n in c(2, 3, 9, 20)) {
    K <- matrix_k(n)
    expect_equal(dim(K), c(n - 1, n))
    expect_equal(K %*% t(K), diag(n - 1))
    expect_equal(drop(K %*% rep(1, n)), rep(0, n - 1))
    expect_equal(t(K) %*% K, diag(n) - 1 / n)
  }
  expect_error(matrix_k(1))
  V <- orthnorm(cbind(1, 1:6), basis = TRUE, norm = TRUE)
  expect_equal(crossprod(V), diag(6))
  V2 <- orthnorm(cbind(1, 1:6), basis = FALSE, norm = FALSE)
  expect_equal(sum(V2[, 1] * V2[, 2]), 0)
  expect_null(orthnorm(NULL))
  expect_warning(orthnorm(matrix(1:6, 2)))
})

test_that("input errors are reported by name", {
  expect_error(ge_mean(.data = Maiz, .y = yield, .gen = entry, .env = site), "\\.y")
  d <- Maiz[!(Maiz$entry == 1 & Maiz$site == 1), ]
  expect_error(ge_mean(.data = d, .y = y, .gen = entry, .env = site), "genotype 1 in environment 1")
  d <- Maiz[-1, ]
  # unequal replicates are completed by data augmentation since version 1.0.0
  expect_s3_class(bayes_ammi(.data = d, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 10), "bam_bayes_ammi")
  expect_error(bayes_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 10, .het = TRUE,
                          .model = "nope"))
  expect_error(bayes_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 10, .k = 9))
})
