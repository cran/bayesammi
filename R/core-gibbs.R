#-------------------------------------------------------------------------------------------------------------------
# The Gibbs sampler of Crossa et al. (2011) for the AMMI model with K bilinear terms, written for the table of cell
# means Y (G x E) of a balanced trial with R replicates per cell. With the sum-to-zero contrasts of the main effects
# and the orthonormality of the singular vectors, the conditional distributions of the paper (eq. 5-10) reduce to
#   mu        ~ N((n ybar + c mu0) / (n + c), s2 / (n + c)),                             c = s2 / v_mu
#   t (G - 1) ~ N((R E K_G ybar_i. + c K_G m_tau) / (R E + c), s2 / (R E + c) I),       tau = K_G' t
#   d (E - 1) ~ N((R G K_E ybar_.j + c K_E m_del) / (R G + c), s2 / (R G + c) I),       delta = K_E' d
#   lambda_k  ~ N((R a_k' Y g_k + c l0) / (R + c), s2 / (R + c)) truncated to (lambda_k+1, lambda_k-1)
#   a_k       ~ von Mises-Fisher(R lambda_k / s2 * Y g_k) on the unit sphere orthogonal to 1 and a_k*, k* != k
#   g_k       ~ von Mises-Fisher(R lambda_k / s2 * Y' a_k) on the unit sphere orthogonal to 1 and g_k*
#   s2        ~ (nu s0^2 + SSW + R sum (Y - fit)^2) / chi^2(nu + n)
#-------------------------------------------------------------------------------------------------------------------

#' Prior distributions of the Bayesian AMMI model
#'
#' @description
#' \code{bam_priors()} sets the hyperparameters of the prior distributions of \code{\link{bayes_ammi}()}. The
#' defaults are the vague but proper priors of Crossa et al. (2011): normal priors with mean zero and variance
#' \eqn{10^8}{10^8} for the overall mean, the main effects and the singular values, and a scaled inverse
#' chi-square prior with one degree of freedom for the error variance.
#'
#' @details
#' With \eqn{G}{G} genotypes, \eqn{E}{E} environments and \eqn{K}{K} bilinear terms,
#' \deqn{\mu \sim N(\mu_0, \sigma^2_\mu), \quad \boldsymbol{\tau} \sim N(\mathbf{m}_\tau, \sigma^2_\tau
#'   \mathbf{I}_G), \quad \boldsymbol{\delta} \sim N(\mathbf{m}_\delta, \sigma^2_\delta \mathbf{I}_E),}{mu ~ N(mu0,
#'   v_mu); tau ~ N(m_tau, v_tau I); delta ~ N(m_delta, v_delta I);}
#' each main effect restricted to sum zero,
#' \deqn{\lambda_k \sim N^{+}(\lambda_0, \sigma^2_\lambda), \quad \lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_K > 0,
#'   \qquad \sigma^2 \sim \mathrm{Inv}\textrm{-}\chi^2(\nu, S^2_0),}{lambda_k ~ N+(lambda0, v_lambda), lambda_1 >= ...
#'   >= lambda_K > 0; sigma^2 ~ Inv-chi^2(nu, s0^2),}
#' and spherical uniform priors for the singular vectors on the subspaces orthogonal to \eqn{\mathbf{1}}{1} and to
#' the other vectors. \eqn{S^2_0}{s0^2} = \code{NULL} (default) takes the residual mean square of the
#' least-squares AMMI fit, so that the prior is equally vague whatever the unit of the trait; with \eqn{\nu = 1}{nu
#' = 1} it carries the information of a single observation. Crossa et al. (2011) report \eqn{S_0 = 10^8}{s0 = 10^8};
#' their published posterior of \eqn{\sigma^2}{sigma^2} is reproduced with the vague default (see the vignette on
#' validation).
#'
#' @param mu_mean,mu_var prior mean and variance of the overall mean \eqn{\mu}{mu}.
#' @param tau_mean,tau_var prior mean (a scalar or a vector of length \eqn{G}{G}) and variance of the genotype
#'   effects \eqn{\tau_i}{tau_i}.
#' @param delta_mean,delta_var prior mean (scalar or length \eqn{E}{E}) and variance of the environment effects
#'   \eqn{\delta_j}{delta_j}.
#' @param lambda_mean,lambda_var prior mean and variance of the singular values \eqn{\lambda_k}{lambda_k}.
#' @param nu,s2 degrees of freedom \eqn{\nu}{nu} and scale \eqn{S^2_0}{s0^2} of the prior of the error variance
#'   (\code{s2 = NULL}: the least-squares residual mean square).
#' @return A list of class \code{"bam_priors"} with these elements.
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#'
#' Hobert, J. P. and Casella, G. (1996). The effect of improper priors on Gibbs sampling in hierarchical linear mixed
#' models. \emph{Journal of the American Statistical Association}, \strong{91}(436), 1461--1473.
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bayes_ammi}}
#' @examples
#' bam_priors(
#'     mu_mean     = 0
#'   , mu_var      = 1e8
#'   , tau_mean    = 0
#'   , tau_var     = 1e8
#'   , delta_mean  = 0
#'   , delta_var   = 1e8
#'   , lambda_mean = 0
#'   , lambda_var  = 1e8
#'   , nu          = 1
#'   , s2          = NULL
#'   )
#' @export
bam_priors <- function(mu_mean = 0, mu_var = 1e8, tau_mean = 0, tau_var = 1e8, delta_mean = 0, delta_var = 1e8,
                       lambda_mean = 0, lambda_var = 1e8, nu = 1, s2 = NULL) {
  fn <- "bam_priors"
  for (nm in c("mu_var", "tau_var", "delta_var", "lambda_var")) .bam_check_number(get(nm), nm, fn, min = 1e-12)
  .bam_check_number(nu, "nu", fn, min = 1e-12)
  if (!is.null(s2)) .bam_check_number(s2, "s2", fn, min = 1e-12)
  for (nm in c("mu_mean", "lambda_mean")) .bam_check_number(get(nm), nm, fn)
  for (nm in c("tau_mean", "delta_mean")) .bam_check_number(get(nm), nm, fn, scalar = FALSE)
  structure(list(mu_mean = mu_mean, mu_var = mu_var, tau_mean = tau_mean, tau_var = tau_var,
                 delta_mean = delta_mean, delta_var = delta_var, lambda_mean = lambda_mean,
                 lambda_var = lambda_var, nu = nu, s2 = s2), class = "bam_priors")
}

#' @rdname bam_priors
#' @param x a \code{"bam_priors"} object.
#' @param ... unused.
#' @export
print.bam_priors <- function(x, ...) {
  f <- function(v) if (is.null(v)) "least-squares residual mean square" else
    paste(format(v, scientific = TRUE, digits = 3L), collapse = ", ")
  cat("Priors of the Bayesian AMMI model\n")
  cat(sprintf("  mu            ~ N(%s, %s)\n", f(x$mu_mean), f(x$mu_var)))
  cat(sprintf("  tau_i         ~ N(%s, %s), sum zero\n", f(x$tau_mean), f(x$tau_var)))
  cat(sprintf("  delta_j       ~ N(%s, %s), sum zero\n", f(x$delta_mean), f(x$delta_var)))
  cat(sprintf("  lambda_k      ~ N+(%s, %s), ordered\n", f(x$lambda_mean), f(x$lambda_var)))
  cat(sprintf("  sigma^2       ~ Inv-chi^2(nu = %s, s2 = %s)\n", f(x$nu), f(x$s2)))
  cat("  alpha_k, gamma_k uniform on the constrained unit spheres\n")
  invisible(x)
}

# one draw from N(m, s^2) truncated to (lo, hi), by the inverse distribution function on the tail side that keeps
# precision; an exponential tail approximation when the interval is beyond double precision
.bam_rtnorm <- function(m, s, lo, hi) {
  a <- (lo - m) / s
  b <- (hi - m) / s
  if (a > 0) {
    pa <- stats::pnorm(a, lower.tail = FALSE)
    pb <- stats::pnorm(b, lower.tail = FALSE)
    if (pa - pb <= 0) return(lo + s * stats::rexp(1L, a))
    return(m + s * stats::qnorm(stats::runif(1L, pb, pa), lower.tail = FALSE))
  }
  pa <- stats::pnorm(a)
  pb <- stats::pnorm(b)
  if (pb - pa <= 0) return(hi - s * stats::rexp(1L, -b))
  m + s * stats::qnorm(stats::runif(1L, pa, pb))
}

# one draw from the von Mises-Fisher distribution on the unit sphere of R^p with parameter v = kappa * direction
# (Wood 1994; p = 1: the two points +1 and -1)
.bam_rvmf <- function(v) {
  p   <- length(v)
  kap <- sqrt(sum(v^2))
  if (p == 1L) return(if (stats::runif(1L) < stats::plogis(2 * v)) 1 else -1)
  if (kap < 1e-12) {
    x <- stats::rnorm(p)
    return(x / sqrt(sum(x^2)))
  }
  mu <- v / kap
  m1 <- p - 1
  b  <- m1 / (2 * kap + sqrt(4 * kap^2 + m1^2))
  x0 <- (1 - b) / (1 + b)
  cc <- kap * x0 + m1 * log(1 - x0^2)
  repeat {
    z <- stats::rbeta(1L, m1 / 2, m1 / 2)
    w <- (1 - (1 + b) * z) / (1 - (1 - b) * z)
    if (kap * w + m1 * log(1 - x0 * w) - cc >= log(stats::runif(1L))) break
  }
  u <- stats::rnorm(p)
  u <- u - sum(u * mu) * mu
  u <- u / sqrt(sum(u^2))
  w * mu + sqrt(max(0, 1 - w^2)) * u
}

# orthonormal basis of the orthogonal complement of the columns of A: Gram-Schmidt (two passes) of the columns of A,
# then of the unit vectors in order (the algorithm of null_basis() in src/gibbs.cpp, so that the R reference sampler
# and the compiled one draw the same numbers)
.bam_null <- function(A) {
  p <- nrow(A)
  Q <- matrix(0, p, 0L)
  cand <- c(lapply(seq_len(ncol(A)), function(c) A[, c]), lapply(seq_len(p), function(i) replace(numeric(p), i, 1)))
  tol  <- c(rep(1e-10, ncol(A)), rep(1e-6, p))
  r <- NA_integer_
  for (h in seq_along(cand)) {
    if (h == ncol(A) + 1L) r <- ncol(Q)
    if (h > ncol(A) && ncol(Q) == p) break
    v <- cand[[h]]
    for (pass in 1:2) if (ncol(Q)) v <- v - Q %*% crossprod(Q, v)
    nv <- sqrt(sum(v^2))
    if (nv > tol[h]) Q <- cbind(Q, v / nv)
  }
  if (is.na(r)) r <- ncol(Q)
  Q[, -seq_len(r), drop = FALSE]
}

# orthonormal columns orthogonal to 1, from a starting matrix (Gram-Schmidt through QR)
.bam_orth1 <- function(A) {
  Q <- qr.Q(qr(cbind(1, A)))
  Q[, -1L, drop = FALSE] * rep(sign(colSums(Q[, -1L, drop = FALSE] * A)), each = nrow(A))
}

# One chain. Y: G x E means; R: replicates; ssw: within-cell sum of squares; K: bilinear terms; pr: priors with s2
# filled in; start: list(sigma2, lambda, alpha, gamma); ref: list(alpha) for the signs; keep: indices of the
# iterations kept. Returns the kept draws: `scalar` (mu, sigma2, tau, delta, lambda), `alpha`, `gamma` arrays.
.bam_chain <- function(Y, R, ssw, K, pr, start, ref, n_total, keep, align = TRUE) {
  G  <- nrow(Y)
  E  <- ncol(Y)
  n  <- R * G * E
  KG <- matrix_k(G)
  KE <- matrix_k(E)
  ybar <- mean(Y)
  yr   <- drop(KG %*% rowMeans(Y))
  yc   <- drop(KE %*% colMeans(Y))
  mtau <- drop(KG %*% rep_len(pr$tau_mean, G))
  mdel <- drop(KE %*% rep_len(pr$delta_mean, E))
  s2  <- start$sigma2
  lam <- start$lambda
  A   <- start$alpha
  Gm  <- start$gamma
  nk  <- length(keep)
  sc  <- matrix(NA_real_, nk, 2L + G + E + K)
  aa  <- array(NA_real_, c(nk, G, K))
  gg  <- array(NA_real_, c(nk, E, K))
  slot <- integer(n_total)
  slot[keep] <- seq_len(nk)
  for (it in seq_len(n_total)) {
    # overall mean and main effects (eq. 5)
    c_mu <- s2 / pr$mu_var
    mu   <- stats::rnorm(1L, (n * ybar + c_mu * pr$mu_mean) / (n + c_mu), sqrt(s2 / (n + c_mu)))
    c_t  <- s2 / pr$tau_var
    tau  <- drop(crossprod(KG, stats::rnorm(G - 1L, (R * E * yr + c_t * mtau) / (R * E + c_t),
                                            sqrt(s2 / (R * E + c_t)))))
    c_d  <- s2 / pr$delta_var
    del  <- drop(crossprod(KE, stats::rnorm(E - 1L, (R * G * yc + c_d * mdel) / (R * G + c_d),
                                            sqrt(s2 / (R * G + c_d)))))
    # bilinear terms (eq. 7-9), one after the other
    c_l <- s2 / pr$lambda_var
    for (k in seq_len(K)) {
      yg <- drop(Y %*% Gm[, k])
      lo <- if (k < K) lam[k + 1L] else 0
      hi <- if (k > 1L) lam[k - 1L] else Inf
      lam[k] <- .bam_rtnorm((R * sum(A[, k] * yg) + c_l * pr$lambda_mean) / (R + c_l), sqrt(s2 / (R + c_l)), lo, hi)
      H <- .bam_null(cbind(1, A[, -k, drop = FALSE]))
      A[, k] <- H %*% .bam_rvmf(drop(crossprod(H, yg)) * (R * lam[k] / s2))
      H <- .bam_null(cbind(1, Gm[, -k, drop = FALSE]))
      Gm[, k] <- H %*% .bam_rvmf(drop(crossprod(H, crossprod(Y, A[, k]))) * (R * lam[k] / s2))
      # the pair (alpha_k, gamma_k) and its negative give the same model: keep the sign of the reference
      if (align && sum(A[, k] * ref$alpha[, k]) < 0) {
        A[, k]  <- -A[, k]
        Gm[, k] <- -Gm[, k]
      }
    }
    # pair moves for score matrices whose columns fill the constrained space
    if (K > 1L) {
      if (G - 1L == K) A <- .bam_pair_sweep(A, (Y %*% Gm %*% diag(lam * R, K)) / s2, TRUE)
      if (E - 1L == K) Gm <- .bam_pair_sweep(Gm, (crossprod(Y, A) %*% diag(lam * R, K)) / s2, TRUE)
      if (align) for (k in seq_len(K)) if (sum(A[, k] * ref$alpha[, k]) < 0) {
        A[, k]  <- -A[, k]
        Gm[, k] <- -Gm[, k]
      }
    }
    # error variance (eq. 10)
    res <- Y - mu - outer(tau, del, `+`) - A %*% (lam * t(Gm))
    s2  <- (pr$nu * pr$s2 + ssw + R * sum(res^2)) / stats::rchisq(1L, pr$nu + n)
    j <- slot[it]
    if (j > 0L) {
      sc[j, ]   <- c(mu, s2, tau, del, lam)
      aa[j, , ] <- A
      gg[j, , ] <- Gm
    }
  }
  list(scalar = sc, alpha = aa, gamma = gg)
}

# starting values of chain `ch`: the least-squares fit for the first chain, dispersed values for the others
.bam_start <- function(ls, s2, ch, het = FALSE, E = 1L) {
  s2v <- if (het) rep_len(s2, E) else s2[1L]
  base <- list(mu = ls$mu, tau = ls$tau, delta = ls$delta)
  if (ch == 1L) return(c(base, list(sigma2 = s2v, lambda = ls$lambda, alpha = ls$alpha, gamma = ls$gamma)))
  K   <- length(ls$lambda)
  lam <- sort(ls$lambda * stats::runif(K, 0.5, 1.5), decreasing = TRUE)
  sd_a <- 0.5 / sqrt(nrow(ls$alpha))
  sd_g <- 0.5 / sqrt(nrow(ls$gamma))
  gam <- ls$gamma + stats::rnorm(length(ls$gamma), 0, sd_g)
  c(base, list(sigma2 = s2v * stats::runif(1L, 0.5, 2), lambda = lam,
               alpha = .bam_orth1(ls$alpha + stats::rnorm(length(ls$alpha), 0, sd_a)),
               gamma = if (isTRUE(ls$gge)) qr.Q(qr(gam)) * rep(sign(colSums(qr.Q(qr(gam)) * gam)), each = nrow(gam))
                       else .bam_orth1(gam)))
}

# one draw from the von Mises distribution (Best and Fisher 1979), the random numbers in the order of the compiled
# sampler
.bam_rvonmises <- function(m, kap) {
  if (kap < 1e-8) return(pi * (2 * stats::runif(1L) - 1))
  t <- 1 + sqrt(1 + 4 * kap^2)
  rho <- (t - sqrt(2 * t)) / (2 * kap)
  r <- (1 + rho^2) / (2 * rho)
  repeat {
    u1 <- stats::runif(1L)
    u2 <- stats::runif(1L)
    z <- cos(pi * u1)
    f <- (1 + r * z) / (r + z)
    cc <- kap * (r - f)
    if (cc * (2 - cc) - u2 > 0 || log(cc / u2) + 1 - cc >= 0) break
  }
  u3 <- stats::runif(1L)
  (if (u3 > 0.5) 1 else -1) * acos(f) + m
}

# two columns at a time (Hoff 2009) for score matrices whose columns fill the constrained space: given the other
# columns the pair is a rotation (and reflection) in a plane, drawn exactly (see src/gibbs.cpp)
.bam_pair_sweep <- function(X, C, one) {
  K <- ncol(X)
  for (k in seq_len(K - 1L)) for (l in (k + 1L):K) {
    O <- cbind(if (one) 1, X[, -c(k, l), drop = FALSE])
    N <- .bam_null(O)
    if (ncol(N) != 2L) next
    a <- drop(crossprod(N, C[, k]))
    b <- drop(crossprod(N, C[, l]))
    ap <- a[1] + b[2]; bp <- a[2] - b[1]
    am <- a[1] - b[2]; bm <- a[2] + b[1]
    kp <- sqrt(ap^2 + bp^2); km <- sqrt(am^2 + bm^2)
    lp <- log(besselI(kp, 0, TRUE)) + kp
    lm <- log(besselI(km, 0, TRUE)) + km
    e <- if (stats::runif(1L) < 1 / (1 + exp(lm - lp))) 1 else -1
    th <- if (e > 0) .bam_rvonmises(atan2(bp, ap), kp) else .bam_rvonmises(atan2(bm, am), km)
    X[, k] <- N[, 1L] * cos(th) + N[, 2L] * sin(th)
    X[, l] <- e * (-N[, 1L] * sin(th) + N[, 2L] * cos(th))
  }
  X
}

# one chain of the compiled sampler; cd: the data core, pr: priors (s2 filled in), keep: 0/1 per iteration
.bam_run_chain <- function(cd, K, gge, het, pr, start, ref_alpha, n_total, keep, align, waic) {
  G <- length(cd$G)
  E <- length(cd$E)
  prc <- list(mu_mean = pr$mu_mean, mu_var = pr$mu_var, tau_mean = rep_len(pr$tau_mean, G), tau_var = pr$tau_var,
              delta_mean = rep_len(pr$delta_mean, E), delta_var = pr$delta_var, lambda_mean = pr$lambda_mean,
              lambda_var = pr$lambda_var, nu = pr$nu, s2 = as.numeric(pr$s2))
  .bam_gibbs_cpp(cd$S, cd$Q, matrix(as.numeric(cd$N), G, E), as.integer(cd$B), as.integer(K), gge, het, prc, start,
                 ref_alpha, align, as.integer(n_total), as.integer(keep), cd$data$y,
                 as.integer(cd$data$cell - 1L), waic)
}
