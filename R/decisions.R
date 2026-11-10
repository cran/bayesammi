#-------------------------------------------------------------------------------------------------------------------
# Decisions from the posterior: predictions of every cell, which genotype wins where, probabilities of superior
# performance and stability (Dias et al. 2022), rank distributions and Bayesian stability values. All of them are
# functions of the draws of the cell means eta_ij = mu + tau_i + delta_j + sum_k lambda_k alpha_ik gamma_jk.
#-------------------------------------------------------------------------------------------------------------------

# ranks within every row of X (smallest = 1), as t(apply(X, 1, rank, ties.method = ties)) but by comparisons of
# whole columns, which is much faster for many draws and few genotypes
.bam_row_rank <- function(X, ties = c("first", "average")) {
  ties <- match.arg(ties)
  G <- ncol(X)
  R <- matrix(1, nrow(X), G)
  for (i in seq_len(G)) for (k in seq_len(G)[-i]) {
    R[, i] <- R[, i] + (X[, k] < X[, i]) + (if (ties == "first") (k < i) else 0.5) * (X[, k] == X[, i])
  }
  R
}

# variance over the environments of every draw and genotype of an array draws x G x E
.bam_row_var <- function(A) {
  dm <- dim(A)
  X  <- matrix(A, dm[1L] * dm[2L], dm[3L])
  matrix(rowSums((X - rowMeans(X))^2) / (dm[3L] - 1L), dm[1L], dm[2L])
}

# draws of the cell means: array (draws x G x E), and the draws of the marginal genotype values and interaction
.bam_eta <- function(fit) {
  m   <- attr(fit, "meta")
  G   <- m$genotypes
  E   <- m$environments
  if (inherits(fit, "bam_bayes_stability")) {
    lay <- fit$layout
    all <- do.call(rbind, lapply(fit$draws, unclass))
    eta <- .bam_stab_eta(fit)
    dimnames(eta) <- NULL
  } else {
    lay <- m$layout
    all <- .bam_all(fit)
    eta <- array(0, c(nrow(all), G, E))
    for (j in seq_len(E)) eta[, , j] <- .bam_eta_env(all, lay, j, G)
  }
  S    <- nrow(all)
  marg <- matrix(rowMeans(matrix(eta, S * G, E)), S, G)
  envm <- vapply(seq_len(E), function(j) rowMeans(eta[, , j, drop = FALSE]), numeric(S))
  envm <- matrix(envm, S, E)
  ge   <- eta - as.vector(marg) - as.vector(envm[, rep(seq_len(E), each = G), drop = FALSE]) + rowMeans(marg)
  list(eta = eta, marginal = marg, interaction = ge, G = m$gen_labels, E = m$env_labels, all = all, lay = lay,
       level = m$level)
}

.bam_need_fit <- function(x, fn) {
  if (!inherits(x, c("bam_bayes_ammi", "bam_bayes_stability")))
    .bam_stop(fn, "`fit` must be a result of bayes_ammi() or bayes_stability().")
  invisible(TRUE)
}

# functions of the bilinear terms need an AMMI or GGE fit
.bam_need_ammi <- function(x, fn) {
  if (!inherits(x, "bam_bayes_ammi"))
    .bam_stop(fn, "`fit` must be a result of bayes_ammi(): this analysis uses the singular values and scores of ",
              "the bilinear terms. For bayes_stability() fits see the component stability and bam_indices().")
  invisible(TRUE)
}

#' Predictions of every genotype-environment cell
#'
#' @description
#' \code{bam_predict()} gives the posterior distribution of the expected response of every genotype in every
#' environment, including the cells that were not observed, and of the genotype means over the environments.
#'
#' @details
#' For draw \eqn{s}{s},
#' \deqn{\eta^{(s)}_{ij} = \mu^{(s)} + \tau^{(s)}_i + \delta^{(s)}_j + \sum_{k=1}^{K}\lambda^{(s)}_k
#'   \alpha^{(s)}_{ik}\gamma^{(s)}_{jk}}{eta_ij = mu + tau_i + delta_j + sum_k lambda_k alpha_ik gamma_jk}
#' (without \eqn{\tau_i}{tau_i} in the GGE model); the genotype mean over the environments is \eqn{\bar{\eta}_{i\cdot}
#' = E^{-1}\sum_j\eta_{ij}}{eta_i. = sum_j eta_ij / E}. A missing cell is predicted from the main effects and from
#' the interaction pattern its genotype and environment show in the other cells; the data augmentation of
#' \code{\link{bayes_ammi}()} has already propagated its uncertainty into all parameters.
#'
#' @param fit a result of \code{\link{bayes_ammi}()} or \code{\link{bayes_stability}()} (\code{bam_stability()}
#'   needs \code{bayes_ammi()}).
#' @param .level probability of the HPD intervals (default: that of the fit).
#' @return An object of class \code{c("bam_predict", "bam_result", "list")}:
#' \describe{
#'   \item{\code{cells}}{\code{data.table}: \code{Gen}, \code{Env}, \code{n} (plots observed), \code{observed}
#'     (their mean, \code{NA} when missing), \code{mean}, \code{sd}, \code{lower}, \code{upper}.}
#'   \item{\code{genotypes}}{\code{data.table}: posterior mean, sd and HPD interval of every genotype mean.}
#' }
#'   Methods: \code{print()}, \code{plot()} (\code{which = c("cells", "genotypes")}), \code{\link{bam_interpret}()}.
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bam_winners}}, \code{\link{bam_superiority}}
#' @family decisions
#' @examples
#' data(Maiz)
#' d <- Maiz[!(Maiz$entry == 8 & Maiz$site == 8), ]
#' fit <-
#'   bayes_ammi(
#'       .data   = d
#'     , .y      = y
#'     , .gen    = entry
#'     , .env    = site
#'     , .rep    = rep
#'     , .nIter  = 1000
#'     , .k      = 2
#'     , .burnin = 200
#'     , .thin   = 1
#'     , .chains = 1
#'     , .priors = bam_priors(
#'           mu_mean     = 0
#'         , mu_var      = 1e+8
#'         , tau_mean    = 0
#'         , tau_var     = 1e+8
#'         , delta_mean  = 0
#'         , delta_var   = 1e+8
#'         , lambda_mean = 0
#'         , lambda_var  = 1e+8
#'         , nu          = 1
#'         , s2          = NULL
#'         )
#'     , .level  = 0.95
#'     , .seed   = 1
#'     , .align  = TRUE
#'     , .model  = "AMMI"
#'     , .het    = FALSE
#'     , .waic   = TRUE
#'     )
#' pr <-
#'   bam_predict(
#'       fit    = fit
#'     , .level = 0.95
#'     )
#' pr$cells[n == 0]
#' @export
bam_predict <- function(fit, .level = NULL) {
  fn <- "bam_predict"
  .bam_need_fit(fit, fn)
  lev <- if (is.null(.level)) attr(fit, "meta")$level else .bam_check_number(.level, ".level", fn, min = 0.5, max = 0.999)
  d   <- .bam_eta(fit)
  m   <- attr(fit, "meta")
  G <- length(d$G)
  E <- length(d$E)
  M   <- matrix(d$eta, dim(d$eta)[1L])
  h   <- .bam_hpd(M, lev)
  cells <- data.table::data.table(Gen = factor(rep(d$G, E), levels = d$G), Env = factor(rep(d$E, each = G), levels = d$E),
                                  n = as.vector(m$cells$N), observed = as.vector(m$cells$Y), mean = colMeans(M),
                                  sd = apply(M, 2L, stats::sd), lower = h[, 1L], upper = h[, 2L])
  hg <- .bam_hpd(d$marginal, lev)
  gen <- data.table::data.table(Gen = factor(d$G, levels = d$G), mean = colMeans(d$marginal),
                                sd = apply(d$marginal, 2L, stats::sd), lower = hg[, 1L], upper = hg[, 2L])
  .bam_result(list(cells = cells, genotypes = gen), "predict",
              label = sprintf("Posterior predictions of %d cells (%d not observed), %g%% HPD", G * E,
                              sum(cells$n == 0L), 100 * lev),
              call = match.call(), meta = list(level = lev, response = m$response))
}

#' Which genotype wins where: posterior probabilities
#'
#' @description
#' \code{bam_winners()} gives, for every environment, the posterior probability that each genotype has the best
#' expected response, the most probable winner, and the groups of environments with the same most probable winner
#' (candidate mega-environments). This is the Bayesian counterpart of the which-won-where view of the GGE biplot,
#' with a probability attached to every statement.
#'
#' @details
#' \deqn{P_{ij} = \Pr\left(\eta_{ij} = \max_{i'}\eta_{i'j} \mid \mathbf{y}\right) \approx \frac{1}{S}\sum_{s=1}^{S}
#'   I\left(\eta^{(s)}_{ij} = \max_{i'}\eta^{(s)}_{i'j}\right),}{P_ij = Pr(eta_ij = max_i' eta_i'j | y), the share of
#'   draws in which genotype i is best in environment j,}
#' with the minimum instead of the maximum when \code{higher_better = FALSE}.
#'
#' @inheritParams bam_predict
#' @param higher_better \code{TRUE} (default) when larger values of the trait are better.
#' @return An object of class \code{c("bam_winners", "bam_result", "list")}:
#' \describe{
#'   \item{\code{probabilities}}{matrix \eqn{G \times E}{G x E} of the probabilities \eqn{P_{ij}}{P_ij}.}
#'   \item{\code{winners}}{\code{data.table}: \code{Env}, \code{Winner} (most probable), \code{Probability},
#'     \code{Second} and \code{Probability2} (the runner-up).}
#'   \item{\code{mega}}{\code{data.table}: groups of environments with the same most probable winner.}
#' }
#'   Methods: \code{print()}, \code{plot()}, \code{\link{bam_interpret}()}.
#' @references
#' Yan, W., Hunt, L. A., Sheng, Q. and Szlavnics, Z. (2000). Cultivar evaluation and mega-environment investigation
#' based on the GGE biplot. \emph{Crop Science}, \strong{40}(3), 597--605. \doi{10.2135/cropsci2000.403597x}
#'
#' Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A general Bayesian estimation method of linear-bilinear
#' models applied to plant breeding trials with genotype x environment interaction. \emph{Journal of Agricultural,
#' Biological, and Environmental Statistics}, \strong{17}, 15--37. \doi{10.1007/s13253-011-0063-9}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bam_superiority}}, \code{\link{biplots}}
#' @family decisions
#' @examples
#' data(Maiz)
#' fit <-
#'   bayes_ammi(
#'       .data   = Maiz
#'     , .y      = y
#'     , .gen    = entry
#'     , .env    = site
#'     , .rep    = rep
#'     , .nIter  = 1000
#'     , .k      = 2
#'     , .burnin = 200
#'     , .thin   = 1
#'     , .chains = 1
#'     , .priors = bam_priors(
#'           mu_mean     = 0
#'         , mu_var      = 1e+8
#'         , tau_mean    = 0
#'         , tau_var     = 1e+8
#'         , delta_mean  = 0
#'         , delta_var   = 1e+8
#'         , lambda_mean = 0
#'         , lambda_var  = 1e+8
#'         , nu          = 1
#'         , s2          = NULL
#'         )
#'     , .level  = 0.95
#'     , .seed   = 1
#'     , .align  = TRUE
#'     , .model  = "GGE"
#'     , .het    = FALSE
#'     , .waic   = TRUE
#'     )
#' w <-
#'   bam_winners(
#'       fit           = fit
#'     , higher_better = TRUE
#'     )
#' w$winners
#' plot(x = w)
#' @export
bam_winners <- function(fit, higher_better = TRUE) {
  fn <- "bam_winners"
  .bam_need_fit(fit, fn)
  higher_better <- .bam_check_flag(higher_better, "higher_better", fn)
  d <- .bam_eta(fit)
  G <- length(d$G)
  E <- length(d$E)
  P <- vapply(seq_len(E), function(j) {
    best <- if (higher_better) max.col(d$eta[, , j], ties.method = "first") else
      max.col(-d$eta[, , j], ties.method = "first")
    tabulate(best, nbins = G) / length(best)
  }, numeric(G))
  dimnames(P) <- list(d$G, d$E)
  o  <- apply(P, 2L, order, decreasing = TRUE)
  win <- data.table::data.table(Env = d$E, Winner = d$G[o[1L, ]], Probability = P[cbind(o[1L, ], seq_len(E))],
                                Second = d$G[o[2L, ]], Probability2 = P[cbind(o[2L, ], seq_len(E))])
  mega <- win[, list(Environments = paste(Env, collapse = ", "), n = .N, mean_probability = mean(Probability)),
              by = Winner][order(-n)]
  .bam_result(list(probabilities = P, winners = win, mega = mega), "winners",
              label = sprintf("Which genotype wins where: posterior probabilities (%s is better)",
                              if (higher_better) "higher" else "lower"),
              call = match.call(), meta = list(higher_better = higher_better))
}

#' Probabilities of superior performance and stability
#'
#' @description
#' \code{bam_superiority()} computes the probabilities of Dias et al. (2022) for recommending cultivars: the
#' probability that a genotype is among the best selected fraction across environments and within each environment,
#' the probability that it is among the most stable, their joint probability, and pairwise probabilities that one
#' genotype outperforms or is more stable than another.
#'
#' @details
#' With \eqn{\Omega^{(s)}}{Omega(s)} the \eqn{\lceil sG \rceil}{ceiling(s G)} genotypes with the best marginal values
#' \eqn{\bar{\eta}^{(s)}_{i\cdot}}{eta_i.} in draw \eqn{s}{s} (selection intensity \eqn{s}{s} = \code{.s}),
#' \deqn{P^{\mathrm{perf}}_i = \Pr(i \in \Omega \mid \mathbf{y}), \qquad P^{\mathrm{perf}}_{ij} = \Pr(i \in
#'   \Omega_j \mid \mathbf{y}),}{P_perf_i = Pr(i in Omega | y), P_perf_ij = Pr(i in Omega_j | y),}
#' where \eqn{\Omega_j}{Omega_j} selects on \eqn{\eta_{ij}}{eta_ij} within environment \eqn{j}{j}. Stability is the
#' variance over environments of the interaction effects \eqn{ge_{ij} = \eta_{ij} - \bar{\eta}_{i\cdot} -
#' \bar{\eta}_{\cdot j} + \bar{\eta}_{\cdot\cdot}}{ge_ij = eta_ij - eta_i. - eta_.j + eta..} (Shukla 1972):
#' \eqn{P^{\mathrm{stab}}_i}{P_stab_i} is the probability that genotype \eqn{i}{i} is among the
#' \eqn{\lceil sG \rceil}{ceiling(s G)} genotypes with the smallest variance, and the joint probability is
#' \eqn{P^{\mathrm{perf}}_i P^{\mathrm{stab}}_i}{P_perf_i P_stab_i}, as in Dias et al. (2022). Pairwise,
#' \eqn{\Pr(\bar{\eta}_{i\cdot} > \bar{\eta}_{i'\cdot} \mid \mathbf{y})}{Pr(eta_i. > eta_i'. | y)} and
#' \eqn{\Pr(\mathrm{var}_i < \mathrm{var}_{i'} \mid \mathbf{y})}{Pr(var_i < var_i' | y)}.
#'
#' @inheritParams bam_winners
#' @param .s selection intensity: the proportion of genotypes selected (default 0.2).
#' @return An object of class \code{c("bam_superiority", "bam_result", "list")}:
#' \describe{
#'   \item{\code{marginal}}{\code{data.table} sorted by the joint probability: \code{Gen}, \code{mean} (posterior
#'     mean of the marginal value), \code{prob_performance}, \code{prob_stability}, \code{prob_joint}.}
#'   \item{\code{conditional}}{matrix \eqn{G \times E}{G x E}: probability of superior performance in every
#'     environment.}
#'   \item{\code{pairwise_performance}, \code{pairwise_stability}}{matrices \eqn{G \times G}{G x G}: entry
#'     \eqn{(i, i')}{(i, i')} is the probability that row genotype \eqn{i}{i} is better (more stable) than
#'     \eqn{i'}{i'}.}
#' }
#'   Methods: \code{print()}, \code{plot()} (\code{which = c("marginal", "conditional", "pairwise")}),
#'   \code{\link{bam_interpret}()}.
#' @references
#' Dias, K. O. G., dos Santos, J. P. R., Krause, M. D., Piepho, H.-P., Guimaraes, L. J. M., Pastina, M. M. and
#' Garcia, A. A. F. (2022). Leveraging probability concepts for cultivar recommendation in multi-environment trials.
#' \emph{Theoretical and Applied Genetics}. \doi{10.1007/s00122-022-04041-y}
#'
#' Shukla, G. K. (1972). Some statistical aspects of partitioning genotype-environmental components of variability.
#' \emph{Heredity}, \strong{29}(2), 237--245. \doi{10.1038/hdy.1972.87}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bam_winners}}, \code{\link{bam_ranks}}, \code{\link{bam_stability}}
#' @family decisions
#' @examples
#' data(Maiz)
#' fit <-
#'   bayes_ammi(
#'       .data   = Maiz
#'     , .y      = y
#'     , .gen    = entry
#'     , .env    = site
#'     , .rep    = rep
#'     , .nIter  = 1000
#'     , .k      = 2
#'     , .burnin = 200
#'     , .thin   = 1
#'     , .chains = 1
#'     , .priors = bam_priors(
#'           mu_mean     = 0
#'         , mu_var      = 1e+8
#'         , tau_mean    = 0
#'         , tau_var     = 1e+8
#'         , delta_mean  = 0
#'         , delta_var   = 1e+8
#'         , lambda_mean = 0
#'         , lambda_var  = 1e+8
#'         , nu          = 1
#'         , s2          = NULL
#'         )
#'     , .level  = 0.95
#'     , .seed   = 1
#'     , .align  = TRUE
#'     , .model  = "AMMI"
#'     , .het    = FALSE
#'     , .waic   = TRUE
#'     )
#' sp <-
#'   bam_superiority(
#'       fit           = fit
#'     , .s            = 0.2
#'     , higher_better = TRUE
#'     )
#' sp$marginal
#' plot(x = sp, which = "marginal")
#' @export
bam_superiority <- function(fit, .s = 0.2, higher_better = TRUE) {
  fn <- "bam_superiority"
  .bam_need_fit(fit, fn)
  .s <- .bam_check_number(.s, ".s", fn, min = 1e-6, max = 1)
  higher_better <- .bam_check_flag(higher_better, "higher_better", fn)
  d  <- .bam_eta(fit)
  G  <- length(d$G)
  E  <- length(d$E)
  ns <- max(1L, ceiling(.s * G))
  sgn <- if (higher_better) 1 else -1
  top <- function(M) {
    r <- .bam_row_rank(-sgn * M, "first")
    colMeans(r <= ns)
  }
  pp <- top(d$marginal)
  cond <- vapply(seq_len(E), function(j) top(d$eta[, , j]), numeric(G))
  dimnames(cond) <- list(d$G, d$E)
  V  <- .bam_row_var(d$interaction)
  rs <- .bam_row_rank(V, "first")
  ps <- colMeans(rs <= ns)
  pw <- function(M, better) {
    out <- vapply(seq_len(G), function(i) colMeans(better(M[, i], M)), numeric(G))
    out <- t(out)
    diag(out) <- NA_real_
    dimnames(out) <- list(d$G, d$G)
    out
  }
  ppf <- pw(d$marginal, function(a, M) sgn * a > sgn * M)
  pst <- pw(V, function(a, M) a < M)
  marg <- data.table::data.table(Gen = d$G, mean = colMeans(d$marginal), prob_performance = pp,
                                 prob_stability = ps, prob_joint = pp * ps)
  data.table::setorder(marg, -prob_joint, -prob_performance)
  .bam_result(list(marginal = marg, conditional = cond, pairwise_performance = ppf, pairwise_stability = pst),
              "superiority",
              label = sprintf("Probabilities of superior performance and stability (selection of %d of %d genotypes)",
                              ns, G),
              call = match.call(), meta = list(s = .s, selected = ns, higher_better = higher_better))
}

#' Posterior distribution of the ranks of the genotypes
#'
#' @description
#' \code{bam_ranks()} gives the posterior probability that each genotype takes each rank, by its mean over the
#' environments (the marginal value) or within one environment, with the expected rank and its HPD interval.
#'
#' @details
#' In draw \eqn{s}{s} the genotypes are ranked by \eqn{\bar{\eta}^{(s)}_{i\cdot}}{eta_i.} (1 = best); the
#' probability of rank \eqn{r}{r} is the share of draws in which genotype \eqn{i}{i} takes it. Ranks summarise the
#' uncertainty of a selection decision better than intervals of the means, since they account for all the
#' genotypes jointly.
#'
#' @inheritParams bam_winners
#' @param .env \code{NULL} (default: marginal values) or the label of an environment.
#' @return An object of class \code{c("bam_ranks", "bam_result", "list")}: \code{probabilities} (matrix genotypes x
#'   ranks) and \code{summary} (\code{data.table}: \code{Gen}, \code{mean_rank}, \code{median_rank}, \code{lower},
#'   \code{upper} of the HPD interval of the rank, \code{prob_first}).
#'   Methods: \code{print()}, \code{plot()}, \code{\link{bam_interpret}()}.
#' @references
#' Dias, K. O. G., dos Santos, J. P. R., Krause, M. D., Piepho, H.-P., Guimaraes, L. J. M., Pastina, M. M. and
#' Garcia, A. A. F. (2022). Leveraging probability concepts for cultivar recommendation in multi-environment trials.
#' \emph{Theoretical and Applied Genetics}. \doi{10.1007/s00122-022-04041-y}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bam_superiority}}
#' @family decisions
#' @examples
#' data(Maiz)
#' fit <-
#'   bayes_ammi(
#'       .data   = Maiz
#'     , .y      = y
#'     , .gen    = entry
#'     , .env    = site
#'     , .rep    = rep
#'     , .nIter  = 1000
#'     , .k      = 2
#'     , .burnin = 200
#'     , .thin   = 1
#'     , .chains = 1
#'     , .priors = bam_priors(
#'           mu_mean     = 0
#'         , mu_var      = 1e+8
#'         , tau_mean    = 0
#'         , tau_var     = 1e+8
#'         , delta_mean  = 0
#'         , delta_var   = 1e+8
#'         , lambda_mean = 0
#'         , lambda_var  = 1e+8
#'         , nu          = 1
#'         , s2          = NULL
#'         )
#'     , .level  = 0.95
#'     , .seed   = 1
#'     , .align  = TRUE
#'     , .model  = "AMMI"
#'     , .het    = FALSE
#'     , .waic   = TRUE
#'     )
#' rk <-
#'   bam_ranks(
#'       fit           = fit
#'     , higher_better = TRUE
#'     , .env          = NULL
#'     )
#' rk$summary
#' plot(x = rk)
#' @export
bam_ranks <- function(fit, higher_better = TRUE, .env = NULL) {
  fn <- "bam_ranks"
  .bam_need_fit(fit, fn)
  higher_better <- .bam_check_flag(higher_better, "higher_better", fn)
  d <- .bam_eta(fit)
  G <- length(d$G)
  M <- if (is.null(.env)) d$marginal else {
    j <- match(as.character(.env), d$E)
    if (is.na(j)) .bam_stop(fn, "unknown environment \"", .env, "\".")
    d$eta[, , j]
  }
  sgn <- if (higher_better) -1 else 1
  R <- .bam_row_rank(sgn * M, "first")
  P <- vapply(seq_len(G), function(i) tabulate(R[, i], nbins = G) / nrow(R), numeric(G))
  P <- t(P)
  dimnames(P) <- list(d$G, seq_len(G))
  h <- .bam_hpd(R, d$level)
  sm <- data.table::data.table(Gen = d$G, mean_rank = colMeans(R), median_rank = apply(R, 2L, stats::median),
                               lower = h[, 1L], upper = h[, 2L], prob_first = P[, 1L])
  data.table::setorder(sm, mean_rank)
  .bam_result(list(probabilities = P, summary = sm), "ranks",
              label = sprintf("Posterior ranks of the genotypes (%s)",
                              if (is.null(.env)) "means over the environments" else paste("environment", .env)),
              call = match.call(), meta = list(env = .env, higher_better = higher_better, level = d$level))
}

#' Bayesian stability of the genotypes
#'
#' @description
#' \code{bam_stability()} gives the posterior distributions of two stability measures of every genotype: the AMMI
#' stability value (ASV) of Purchase et al. (2000) computed draw by draw, and the variance of the interaction effects
#' over the environments (Shukla's stability variance of the model), with posterior means, HPD intervals, the
#' probability of being the most stable genotype and the posterior expected ranks, combined with the mean into a
#' Bayesian yield stability index.
#'
#' @details
#' With the biplot scores \eqn{u_{ik} = \lambda_k^{1/2}\alpha_{ik}}{u_ik = lambda_k^(1/2) alpha_ik} and the sums of
#' squares of the terms proportional to \eqn{\lambda_k^2}{lambda_k^2},
#' \deqn{\mathrm{ASV}_i = \sqrt{\left(\frac{\bar{\lambda}_1^2}{\bar{\lambda}_2^2}u_{i1}\right)^2 + u_{i2}^2},
#'   \qquad V_i = \frac{1}{E - 1}\sum_j ge_{ij}^2,}{ASV_i = sqrt((lambda_1^2 / lambda_2^2 u_i1)^2 + u_i2^2); V_i =
#'   sum_j ge_ij^2 / (E - 1),}
#' (\eqn{\mathrm{ASV}_i = |u_{i1}|}{ASV_i = |u_i1|} with one term), evaluated in every draw. The weight of the first
#' axis uses the posterior means \eqn{\bar{\lambda}_k}{lambda_k} of the singular values: the ratio of the draws
#' themselves has no finite mean when the posterior of \eqn{\lambda_2}{lambda_2} reaches zero. The Bayesian yield
#' stability index is the sum of the posterior expected ranks of the mean and of the ASV.
#'
#' @inheritParams bam_winners
#' @return An object of class \code{c("bam_stability", "bam_result", "list")} with \code{stability}
#'   (\code{data.table} sorted by the index: \code{Gen}, \code{mean}, \code{ASV}, \code{ASV_lower},
#'   \code{ASV_upper}, \code{variance}, \code{prob_most_stable} (smallest ASV), \code{rank_mean}, \code{rank_ASV},
#'   \code{YSI}).
#'   Methods: \code{print()}, \code{plot()}, \code{\link{bam_interpret}()}.
#' @references
#' Purchase, J. L., Hatting, H. and van Deventer, C. S. (2000). Genotype x environment interaction of winter wheat in
#' South Africa: II. Stability analysis of yield performance. \emph{South African Journal of Plant and Soil},
#' \strong{17}(3), 101--107.
#'
#' Shukla, G. K. (1972). Some statistical aspects of partitioning genotype-environmental components of variability.
#' \emph{Heredity}, \strong{29}(2), 237--245. \doi{10.1038/hdy.1972.87}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bam_superiority}}
#' @family decisions
#' @examples
#' data(Maiz)
#' fit <-
#'   bayes_ammi(
#'       .data   = Maiz
#'     , .y      = y
#'     , .gen    = entry
#'     , .env    = site
#'     , .rep    = rep
#'     , .nIter  = 1000
#'     , .k      = 2
#'     , .burnin = 200
#'     , .thin   = 1
#'     , .chains = 1
#'     , .priors = bam_priors(
#'           mu_mean     = 0
#'         , mu_var      = 1e+8
#'         , tau_mean    = 0
#'         , tau_var     = 1e+8
#'         , delta_mean  = 0
#'         , delta_var   = 1e+8
#'         , lambda_mean = 0
#'         , lambda_var  = 1e+8
#'         , nu          = 1
#'         , s2          = NULL
#'         )
#'     , .level  = 0.95
#'     , .seed   = 1
#'     , .align  = TRUE
#'     , .model  = "AMMI"
#'     , .het    = FALSE
#'     , .waic   = TRUE
#'     )
#' st <-
#'   bam_stability(
#'       fit           = fit
#'     , higher_better = TRUE
#'     )
#' st$stability
#' @export
bam_stability <- function(fit, higher_better = TRUE) {
  fn <- "bam_stability"
  .bam_need_ammi(fit, fn)
  higher_better <- .bam_check_flag(higher_better, "higher_better", fn)
  d   <- .bam_eta(fit)
  lam <- .bam_block(d$all, d$lay, "lambda")
  A   <- .bam_block(d$all, d$lay, "alpha")
  u1  <- A[, , 1L] * sqrt(lam[, 1L])
  # the weight of the first axis from the posterior means of the singular values: a ratio evaluated draw by draw
  # explodes when lambda_2 comes close to zero
  wt  <- if (ncol(lam) > 1L) mean(lam[, 1L])^2 / mean(lam[, 2L])^2 else 1
  asv <- if (ncol(lam) > 1L) sqrt((wt * u1)^2 + (A[, , 2L] * sqrt(lam[, 2L]))^2) else abs(u1)
  V   <- .bam_row_var(d$interaction)
  h   <- .bam_hpd(asv, d$level)
  best <- tabulate(max.col(-asv, ties.method = "first"), nbins = length(d$G)) / nrow(asv)
  rm  <- colMeans(.bam_row_rank(if (higher_better) -d$marginal else d$marginal, "average"))
  ra  <- colMeans(.bam_row_rank(asv, "average"))
  out <- data.table::data.table(Gen = d$G, mean = colMeans(d$marginal), ASV = colMeans(asv), ASV_lower = h[, 1L],
                                ASV_upper = h[, 2L], variance = colMeans(V), prob_most_stable = best,
                                rank_mean = rm, rank_ASV = ra, YSI = rm + ra)
  data.table::setorder(out, YSI)
  .bam_result(list(stability = out), "stability", label = "Bayesian stability: AMMI stability value and interaction variance",
              call = match.call(), meta = list(higher_better = higher_better, level = d$level))
}

#-------------------------------------------------------------------------------------------------------------------
# print and plot methods of the decision results
#-------------------------------------------------------------------------------------------------------------------
.bam_print_predict <- function(x, digits = 4L, ...) {
  cat(attr(x, "label"), "\n\nGenotype means\n")
  print(x$genotypes, digits = digits)
  cat("\nCells (first rows; all in x$cells)\n")
  print(utils::head(x$cells, 10L), digits = digits)
  invisible(x)
}

.bam_plot_predict <- function(x, which = c("cells", "genotypes"), ...) {
  if (match.arg(which) == "genotypes") {
    d <- data.table::copy(x$genotypes)[order(mean)]
    d[, Gen := factor(Gen, levels = Gen)]
    return(ggplot2::ggplot(d, ggplot2::aes(x = mean, y = Gen)) +
             ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper), height = 0.3, colour = .bam_palette[1L]) +
             ggplot2::geom_point(colour = .bam_palette[1L], size = 2) + .bam_x_cont() +
             ggplot2::labs(x = "posterior mean over the environments with HPD interval", y = NULL,
                           title = "Genotype means") + bam_theme())
  }
  d <- data.table::copy(x$cells)
  ggplot2::ggplot(d, ggplot2::aes(x = Env, y = Gen, fill = mean)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_point(data = d[n == 0L], shape = 4, size = 3, colour = "black") +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = .bam_palette[1L], name = "posterior\nmean",
                                 breaks = scales::pretty_breaks(5)) +
    ggplot2::labs(x = "environment", y = "genotype", title = "Predicted cell means",
                  subtitle = "x: cell not observed (predicted from the model)") +
    bam_theme() + ggplot2::theme(legend.position = "right")
}

.bam_print_winners <- function(x, digits = 3L, ...) {
  cat(attr(x, "label"), "\n\n")
  print(x$winners, digits = digits)
  cat("\nGroups of environments with the same most probable winner\n")
  print(x$mega, digits = digits)
  invisible(x)
}

.bam_plot_winners <- function(x, ...) {
  P <- x$probabilities
  d <- data.table::data.table(Gen = factor(rep(rownames(P), ncol(P)), levels = rownames(P)),
                              Env = factor(rep(colnames(P), each = nrow(P)), levels = colnames(P)), p = as.vector(P))
  ggplot2::ggplot(d, ggplot2::aes(x = Env, y = Gen, fill = p)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(data = d[p >= 0.1], ggplot2::aes(label = sprintf("%.2f", p)), size = 2.6, family = "serif") +
    ggplot2::scale_fill_gradient(low = "white", high = .bam_palette[3L], limits = c(0, 1),
                                 breaks = scales::pretty_breaks(5), name = "probability\nof winning") +
    ggplot2::labs(x = "environment", y = "genotype", title = "Which genotype wins where") +
    bam_theme() + ggplot2::theme(legend.position = "right")
}

.bam_print_superiority <- function(x, digits = 3L, ...) {
  cat(attr(x, "label"), "\n\n")
  print(x$marginal, digits = digits)
  invisible(x)
}

.bam_plot_superiority <- function(x, which = c("marginal", "conditional", "pairwise"), ...) {
  which <- match.arg(which)
  if (which == "marginal") {
    d <- data.table::melt(x$marginal[, list(Gen, performance = prob_performance, stability = prob_stability,
                                            joint = prob_joint)], id.vars = "Gen", variable.name = "Probability",
                          value.name = "p")
    d[, Gen := factor(Gen, levels = rev(x$marginal$Gen))]
    return(ggplot2::ggplot(d, ggplot2::aes(x = p, y = Gen, fill = Probability)) +
             ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.75) +
             ggplot2::scale_fill_manual(values = .bam_palette[c(1L, 3L, 2L)], name = NULL) +
             ggplot2::scale_x_continuous(limits = c(0, 1), breaks = scales::pretty_breaks(5),
                                         expand = ggplot2::expansion(mult = c(0, 0.02))) +
             ggplot2::labs(x = "probability", y = NULL,
                           title = sprintf("Probabilities of superior performance and stability (top %d)",
                                           attr(x, "meta")$selected)) +
             bam_theme())
  }
  P <- if (which == "conditional") x$conditional else x$pairwise_performance
  d <- data.table::data.table(r = factor(rep(rownames(P), ncol(P)), levels = rownames(P)),
                              c = factor(rep(colnames(P), each = nrow(P)), levels = colnames(P)), p = as.vector(P))
  ggplot2::ggplot(d, ggplot2::aes(x = c, y = r, fill = p)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_gradient(low = "white", high = .bam_palette[1L], limits = c(0, 1), na.value = "grey90",
                                 breaks = scales::pretty_breaks(5), name = "probability") +
    ggplot2::labs(x = if (which == "conditional") "environment" else "compared with genotype", y = "genotype",
                  title = if (which == "conditional") "Probability of superior performance in each environment" else
                    "Probability that the row genotype outperforms the column genotype") +
    bam_theme() + ggplot2::theme(legend.position = "right")
}

.bam_print_ranks <- function(x, digits = 3L, ...) {
  cat(attr(x, "label"), "\n\n")
  print(x$summary, digits = digits)
  invisible(x)
}

.bam_plot_ranks <- function(x, ...) {
  P <- x$probabilities[x$summary$Gen, , drop = FALSE]
  d <- data.table::data.table(Gen = factor(rep(rownames(P), ncol(P)), levels = rev(rownames(P))),
                              rank = rep(seq_len(ncol(P)), each = nrow(P)), p = as.vector(P))
  ggplot2::ggplot(d, ggplot2::aes(x = rank, y = Gen, fill = p)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_gradient(low = "white", high = .bam_palette[4L], limits = c(0, 1),
                                 breaks = scales::pretty_breaks(5), name = "probability") +
    ggplot2::scale_x_continuous(breaks = .bam_breaks(TRUE), expand = c(0, 0)) +
    ggplot2::labs(x = "rank (1 = best)", y = NULL, title = "Posterior distribution of the ranks") +
    bam_theme() + ggplot2::theme(legend.position = "right")
}

.bam_print_stability <- function(x, digits = 3L, ...) {
  cat(attr(x, "label"), "\n\n")
  print(x$stability, digits = digits)
  invisible(x)
}

.bam_plot_stability <- function(x, ...) {
  d <- data.table::copy(x$stability)
  ggplot2::ggplot(d, ggplot2::aes(x = mean, y = ASV)) +
    ggplot2::geom_vline(xintercept = mean(d$mean), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ASV_lower, ymax = ASV_upper), width = 0, colour = "grey55") +
    ggplot2::geom_point(colour = .bam_palette[1L], size = 2) +
    ggrepel::geom_text_repel(ggplot2::aes(label = Gen), size = 3, family = "serif") +
    .bam_x_cont() + .bam_y_cont() +
    ggplot2::labs(x = "posterior mean over the environments", y = "AMMI stability value (HPD interval)",
                  title = "Mean against stability", subtitle = "lower right: high and stable") +
    bam_theme()
}

utils::globalVariables(c("Gen", "Env", "n", "observed", "Winner", "Probability", "prob_joint", "prob_performance",
                         "prob_stability", "p", "r", "mean_rank", "YSI", "ASV", "ASV_lower", "ASV_upper", "rank"))
