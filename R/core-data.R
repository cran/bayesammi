#-------------------------------------------------------------------------------------------------------------------
# One core for the data of every function: the cells of the genotype-by-environment table (sums, sums of squares,
# counts and means of the observed plots), the within-cell sum of squares and the least-squares AMMI or GGE fit used
# as the starting point of the sampler (with missing cells filled by the EM algorithm).
#-------------------------------------------------------------------------------------------------------------------

# factor with sorted levels (as in the earlier versions) when it is not already a factor
.bam_factor <- function(x) if (is.factor(x)) droplevels(x) else factor(x)

# Columns of an analysis, given as bare names, strings or vectors; `rep_expr` NULL when the function has no .rep
.bam_prep <- function(fn, data, y_expr, gen_expr, env_expr, rep_expr, env_caller, need_rep = TRUE,
                      allow_missing = FALSE) {
  data <- .bam_check_data(data, fn)
  y    <- .bam_column(data, y_expr, ".y", fn, env_caller)
  gen  <- .bam_column(data, gen_expr, ".gen", fn, env_caller)
  env  <- .bam_column(data, env_expr, ".env", fn, env_caller)
  rep  <- if (need_rep) .bam_column(data, rep_expr, ".rep", fn, env_caller)
  cd   <- .bam_ge_data(y$value, gen$value, env$value, if (need_rep) rep$value, fn, allow_missing)
  cd$names <- c(y = y$name, gen = gen$name, env = env$name, rep = if (need_rep) rep$name else NA_character_)
  cd
}

# Genotype-by-environment data. Returns list(Y = G x E means (NA for an empty cell), S, Q, N (sums, sums of squares
# and counts of the plots), G, E (labels), R (replicates per cell when equal, else NA), B (largest count), ssw
# (within-cell sum of squares), n (plots), balanced, missing (number of empty cells), data (data.table of the plots,
# with the cell number))
.bam_ge_data <- function(y, gen, env, rep = NULL, fn, allow_missing = FALSE) {
  if (!is.numeric(y)) .bam_stop(fn, "`.y` must be numeric.")
  ok <- !is.na(y) & !is.na(gen) & !is.na(env) & (if (is.null(rep)) TRUE else !is.na(rep))
  if (!all(ok)) .bam_warn(fn, sum(!ok), " row(s) with missing values are left out.")
  y   <- as.numeric(y[ok])
  gen <- .bam_factor(gen[ok])
  env <- .bam_factor(env[ok])
  G   <- levels(gen)
  E   <- levels(env)
  if (length(G) < 3L) .bam_stop(fn, "at least 3 genotypes are needed.")
  if (length(E) < 3L) .bam_stop(fn, "at least 3 environments are needed.")
  cell <- as.integer(gen) + (as.integer(env) - 1L) * length(G)
  nb   <- length(G) * length(E)
  N    <- matrix(tabulate(cell, nbins = nb), length(G), length(E), dimnames = list(G, E))
  if (any(N == 0L) && !allow_missing) {
    miss <- which(N == 0L, arr.ind = TRUE)
    .bam_stop(fn, sum(N == 0L), " genotype-environment combination(s) have no data (for example genotype ",
              G[miss[1L, 1L]], " in environment ", E[miss[1L, 2L]], "); every genotype must be grown in every ",
              "environment.")
  }
  if (any(rowSums(N) == 0L) || any(colSums(N) == 0L)) .bam_stop(fn, "every genotype and environment needs data.")
  S <- matrix(0, length(G), length(E), dimnames = list(G, E))
  Q <- S
  used <- sort(unique(cell))
  S[used] <- collapse::fsum(y, cell)
  Q[used] <- collapse::fsum(y^2, cell)
  Y <- S / N
  Y[N == 0L] <- NA_real_
  ssw <- sum((y - Y[cell])^2)
  counts <- unique(as.vector(N))
  out <- list(Y = Y, S = S, Q = Q, N = N, G = G, E = E,
              R = if (length(counts) == 1L) N[1L, 1L] else NA_integer_, B = max(N), ssw = ssw, n = length(y),
              balanced = length(counts) == 1L, missing = sum(N == 0L),
              data = data.table::data.table(Gen = gen, Env = env, y = y, cell = cell))
  if (!is.null(rep)) out$data[, Rep := .bam_factor(rep[ok])]
  out
}

# the replicate-within-environment sum of squares (randomized complete blocks in every environment)
.bam_ss_rep <- function(cd) {
  d <- cd$data
  if (is.null(d$Rep)) return(NA_real_)
  er  <- interaction(d$Env, d$Rep, drop = TRUE)
  mer <- collapse::fmean(d$y, er)
  ner <- tabulate(as.integer(er), nbins = nlevels(er))
  me  <- colMeans(cd$Y)
  e_er <- as.integer(d$Env[match(levels(er), as.character(er))])
  sum(ner * (mer - me[e_er])^2)
}

# Least-squares AMMI (gge = FALSE) or GGE (gge = TRUE) fit with K multiplicative terms: grand mean, genotype effects
# (zero for GGE) and environment effects, singular values and vectors, residual sum of squares of the means. Empty
# cells are filled by the EM algorithm (fitted values of the model, iterated).
.bam_ls <- function(Y, K, gge = FALSE, iter = 200L, tol = 1e-9) {
  miss <- is.na(Y)
  fit1 <- function(Yc) {
    mu  <- mean(Yc)
    del <- colMeans(Yc) - mu
    tau <- if (gge) rep(0, nrow(Yc)) else rowMeans(Yc) - mu
    Z   <- Yc - outer(tau, del, `+`) - mu
    s   <- svd(Z)
    kk  <- seq_len(K)
    fit <- s$u[, kk, drop = FALSE] %*% (s$d[kk] * t(s$v[, kk, drop = FALSE]))
    list(mu = mu, tau = tau, delta = del, Z = Z, d = s$d, lambda = s$d[kk], alpha = s$u[, kk, drop = FALSE],
         gamma = s$v[, kk, drop = FALSE], fitted = mu + outer(tau, del, `+`) + fit, rss = sum((Z - fit)^2))
  }
  Yc <- Y
  if (any(miss)) {
    add <- outer(rowMeans(Y, na.rm = TRUE), colMeans(Y, na.rm = TRUE), `+`) - mean(Y, na.rm = TRUE)
    Yc[miss] <- add[miss]
    for (i in seq_len(iter)) {
      f <- fit1(Yc)
      d <- max(abs(Yc[miss] - f$fitted[miss]))
      Yc[miss] <- f$fitted[miss]
      if (d < tol * (1 + max(abs(Yc)))) break
    }
  }
  f <- fit1(Yc)
  # residual sum of squares of the observed cells only
  f$rss <- sum(((Yc - f$fitted)^2)[!miss])
  f$filled <- Yc
  f
}

utils::globalVariables(c("Rep", ":="))
