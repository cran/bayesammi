#' Compare Bayesian AMMI and GGE models
#'
#' @description
#' \code{bam_compare()} compares fits of \code{\link{bayes_ammi}()} that differ in the number of multiplicative
#' terms, the model (AMMI or GGE), the error variances or the priors, by the deviance information criterion (DIC)
#' and the widely applicable information criterion (WAIC). Smaller values predict new data better.
#'
#' @details
#' With \eqn{D(\boldsymbol{\theta}) = -2\log p(\mathbf{y} \mid \boldsymbol{\theta})}{D(theta) = -2 log p(y |
#' theta)} over the observed plots,
#' \deqn{\mathrm{DIC} = D(\bar{\boldsymbol{\theta}}) + 2p_D, \qquad p_D = \overline{D(\boldsymbol{\theta})} -
#'   D(\bar{\boldsymbol{\theta}})}{DIC = D(theta_bar) + 2 pD, pD = mean D(theta) - D(theta_bar)}
#' (Spiegelhalter et al. 2002), the deviance evaluated at the posterior means of the cell means and variances; and
#' \deqn{\mathrm{WAIC} = -2\left(\sum_o \log\frac{1}{S}\sum_s p(y_o \mid \boldsymbol{\theta}^{(s)}) -
#'   \sum_o \mathrm{Var}_s\log p(y_o \mid \boldsymbol{\theta}^{(s)})\right)}{WAIC = -2 (lppd - p_waic)}
#' (Watanabe 2010; Gelman et al. 2014), accumulated over all draws during sampling. All fits must use the same data.
#'
#' The standard error of a difference in WAIC comes from the plot-by-plot differences of the contributions
#' \eqn{w_o}{w_o} to WAIC: \eqn{\mathrm{se} = \sqrt{n\,\mathrm{Var}_o(w_o^{(m)} - w_o^{(\mathrm{best})})}}{se = sqrt(n
#' var(w_m - w_best))} (Vehtari et al. 2017). Because the models are fitted to the same plots, it is smaller
#' than the standard error of either WAIC. A difference of less than about two of its standard errors does not
#' separate two models.
#'
#' @param ... fits of \code{\link{bayes_ammi}()} or \code{\link{bayes_stability}()} to the same data, named or not
#'   (names label the rows).
#' @return An object of class \code{c("bam_compare", "bam_result", "list")} with \code{table}
#'   (\code{data.table} sorted by WAIC: \code{model}, \code{terms}, \code{variances}, \code{DIC}, \code{pD},
#'   \code{WAIC}, \code{p_waic}, \code{se} (standard error of WAIC), \code{lambda_last_lower} (lower HPD limit of the
#'   last singular value), \code{lambda_last_mean}, \code{dWAIC} (difference from the best) and \code{dWAIC_se} (its
#'   standard error)).
#'   Methods: \code{print()}, \code{plot()}, \code{\link{bam_interpret}()}.
#' @references
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A. (2002). Bayesian measures of model
#' complexity and fit. \emph{Journal of the Royal Statistical Society B}, \strong{64}(4), 583--639.
#' \doi{10.1111/1467-9868.00353}
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion
#' in singular learning theory. \emph{Journal of Machine Learning Research}, \strong{11}, 3571--3594.
#'
#' Gelman, A., Hwang, J. and Vehtari, A. (2014). Understanding predictive information criteria for Bayesian models.
#' \emph{Statistics and Computing}, \strong{24}(6), 997--1016. \doi{10.1007/s11222-013-9416-2}
#'
#' Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out
#' cross-validation and WAIC. \emph{Statistics and Computing}, \strong{27}(5), 1413--1432.
#' \doi{10.1007/s11222-016-9696-4}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bam_ppc}}
#' @family model choice
#' @examples
#' data(Maiz)
#' fits <-
#'   lapply(
#'       X   = 1:2
#'     , FUN = function(K) {
#'         bayes_ammi(
#'             .data   = Maiz
#'           , .y      = y
#'           , .gen    = entry
#'           , .env    = site
#'           , .rep    = rep
#'           , .nIter  = 600
#'           , .k      = K
#'           , .burnin = 100
#'           , .thin   = 1
#'           , .chains = 1
#'           , .priors = bam_priors(
#'                 mu_mean     = 0
#'               , mu_var      = 1e+8
#'               , tau_mean    = 0
#'               , tau_var     = 1e+8
#'               , delta_mean  = 0
#'               , delta_var   = 1e+8
#'               , lambda_mean = 0
#'               , lambda_var  = 1e+8
#'               , nu          = 1
#'               , s2          = NULL
#'               )
#'           , .level  = 0.95
#'           , .seed   = 1
#'           , .align  = TRUE
#'           , .model  = "AMMI"
#'           , .het    = FALSE
#'           , .waic   = TRUE
#'           )
#'       }
#'     )
#' bam_compare(
#'     AMMI1 = fits[[1]]
#'   , AMMI2 = fits[[2]]
#'   )
#' @export
bam_compare <- function(...) {
  fn <- "bam_compare"
  fits <- list(...)
  if (!length(fits)) .bam_stop(fn, "give at least one fit.")
  for (f in fits) if (!inherits(f, c("bam_bayes_ammi", "bam_bayes_stability")))
    .bam_stop(fn, "every fit must be a result of bayes_ammi() or bayes_stability().")
  nm <- names(fits)
  if (is.null(nm)) nm <- rep("", length(fits))
  nobs <- vapply(fits, function(f) sum(attr(f, "meta")$cells$N), numeric(1L))
  if (length(unique(nobs)) > 1L) .bam_stop(fn, "the fits use different data (", paste(nobs, collapse = ", "), " plots).")
  tab <- data.table::rbindlist(lapply(seq_along(fits), function(i) {
    f <- fits[[i]]
    s <- f$settings
    stab <- inherits(f, "bam_bayes_stability")
    lam <- if (stab) data.table::data.table(lower = NA_real_, mean = NA_real_) else f$summary[block == "lambda"]
    data.table::data.table(
        model = if (nzchar(nm[i])) nm[i] else if (stab) s$model else
          sprintf("%s%d%s", s$model, s$k, if (s$het) "-het" else "")
      , terms = if (stab) NA_integer_ else s$k
      , variances = if (stab) c(ER = "deviation per genotype", FW = "one deviation variance",
                                Shukla = "interaction per genotype")[[s$model]] else
          if (s$het) "per environment" else "common"
      , DIC = f$fit$DIC, pD = f$fit$pD
      , WAIC = if (is.null(f$fit$WAIC)) NA_real_ else f$fit$WAIC
      , p_waic = if (is.null(f$fit$p_waic)) NA_real_ else f$fit$p_waic
      , se = if (is.null(f$fit$waic_se)) NA_real_ else f$fit$waic_se
      , lambda_last_lower = lam$lower[nrow(lam)], lambda_last_mean = lam$mean[nrow(lam)])
  }))
  key <- if (all(is.na(tab$WAIC))) tab$DIC else tab$WAIC
  tab[, dWAIC := key - min(key)]
  # standard error of the difference from the best model, from the plot-by-plot contributions to WAIC
  wi <- lapply(fits, function(f) f$fit$waic_i)
  best <- which.min(key)
  tab[, dWAIC_se := if (any(vapply(wi, is.null, NA))) NA_real_ else
    vapply(wi, function(w) sqrt(length(w) * stats::var(w - wi[[best]])), numeric(1L))]
  tab <- tab[order(key)]
  .bam_result(list(table = tab), "compare", label = "Comparison of Bayesian models (smaller is better)",
              call = match.call())
}

#' Posterior predictive checks
#'
#' @description
#' \code{bam_ppc()} compares statistics of the observed trial with their distribution in trials replicated from the
#' posterior predictive distribution (Gelman et al. 1996). A posterior predictive p-value near 0 or 1 shows a feature
#' of the data that the model does not reproduce, for instance unequal error variances or a missing interaction term.
#'
#' @details
#' For a subset of \code{.draws} posterior draws, a trial \eqn{\mathbf{y}^{\mathrm{rep}}}{y_rep} with the observed
#' plots is simulated from \eqn{N(\eta_{ij}, \sigma^2_j)}{N(eta_ij, sigma_j^2)}, and for each statistic \eqn{T}{T}
#' \eqn{p = \Pr(T(\mathbf{y}^{\mathrm{rep}}) \ge T(\mathbf{y}) \mid \mathbf{y})}{p = Pr(T(y_rep) >= T(y) | y)}. The
#' statistics: variance of the cell means, largest and smallest cell mean, variance of the genotype means, first
#' singular value of the interaction effects of the cell means, pooled within-cell variance, and the ratio of the
#' largest to the smallest within-environment variance.
#'
#' @param fit a result of \code{\link{bayes_ammi}()} or \code{\link{bayes_stability}()}.
#' @param .draws number of posterior draws used (default 1000).
#' @param .seed seed for the replicated trials (the random number state of the session is restored).
#' @return An object of class \code{c("bam_ppc", "bam_result", "list")} with \code{checks} (\code{data.table}:
#'   \code{statistic}, \code{observed}, \code{mean_rep}, \code{p_value}) and \code{replicated} (\code{data.table} of
#'   the replicated statistics).
#'   Methods: \code{print()}, \code{plot()}, \code{\link{bam_interpret}()}.
#' @references
#' Gelman, A., Meng, X.-L. and Stern, H. (1996). Posterior predictive assessment of model fitness via realized
#' discrepancies. \emph{Statistica Sinica}, \strong{6}(4), 733--807.
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bam_compare}}
#' @family model choice
#' @examples
#' data(Maiz)
#' fit <-
#'   bayes_ammi(
#'       .data   = Maiz
#'     , .y      = y
#'     , .gen    = entry
#'     , .env    = site
#'     , .rep    = rep
#'     , .nIter  = 600
#'     , .k      = 2
#'     , .burnin = 100
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
#' pc <-
#'   bam_ppc(
#'       fit    = fit
#'     , .draws = 300
#'     , .seed  = 1
#'     )
#' pc$checks
#' plot(x = pc)
#' @export
bam_ppc <- function(fit, .draws = 1000L, .seed = NULL) {
  fn <- "bam_ppc"
  .bam_need_fit(fit, fn)
  .draws <- .bam_check_number(.draws, ".draws", fn, min = 10, integer = TRUE)
  m   <- attr(fit, "meta")
  de  <- .bam_eta(fit)
  all <- de$all
  lay <- de$lay
  N   <- m$cells$N
  G <- nrow(N)
  E <- ncol(N)
  if (!is.null(.seed)) {
    had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else
      rm(".Random.seed", envir = globalenv()), add = TRUE)
    set.seed(.seed)
  }
  idx <- unique(round(seq(1, nrow(all), length.out = min(.draws, nrow(all)))))
  # the observed plots: cell number of every plot, from the counts
  cell <- rep(seq_len(G * E), as.vector(N))
  env  <- (cell - 1L) %/% G + 1L
  stat <- function(yv) {
    Ym <- matrix(NA_real_, G, E)
    Ym[sort(unique(cell))] <- collapse::fmean(yv, cell)
    full <- !is.na(Ym)
    Yf <- Ym
    if (any(!full)) Yf[!full] <- (outer(rowMeans(Ym, na.rm = TRUE), colMeans(Ym, na.rm = TRUE), `+`) -
                                    mean(Ym, na.rm = TRUE))[!full]
    Z  <- Yf - outer(rowMeans(Yf), colMeans(Yf), `+`) + mean(Yf)
    within <- yv - Ym[cell]
    wv <- vapply(seq_len(E), function(j) {
      w <- env == j
      sum(within[w]^2) / max(1, sum(w) - sum(N[, j] > 0))
    }, numeric(1L))
    c(`variance of cell means` = stats::var(Ym[full]), `largest cell mean` = max(Ym[full]),
      `smallest cell mean` = min(Ym[full]), `variance of genotype means` = stats::var(rowMeans(Yf)),
      `first singular value` = svd(Z, nu = 0, nv = 0)$d[1L],
      `within-cell variance` = sum(within^2) / max(1, length(yv) - sum(N > 0)),
      `largest / smallest environment variance` = max(wv) / max(min(wv), 1e-12))
  }
  obs <- .bam_obs_stats(fit, stat, cell)
  rep_stats <- t(vapply(idx, function(s) {
    eta <- de$eta[s, , ]
    sg  <- all[s, lay$parameter[lay$block == "sigma2"]]
    sd  <- sqrt(if (length(sg) == 1L) rep(sg, E) else sg)
    yv  <- eta[cell] + stats::rnorm(length(cell)) * sd[env]
    stat(yv)
  }, numeric(length(obs))))
  checks <- data.table::data.table(statistic = names(obs), observed = unname(obs), mean_rep = colMeans(rep_stats),
                                   p_value = colMeans(sweep(rep_stats, 2L, obs, `>=`)))
  .bam_result(list(checks = checks, replicated = data.table::as.data.table(rep_stats)), "ppc",
              label = sprintf("Posterior predictive checks (%d replicated trials)", length(idx)), call = match.call())
}

# the observed statistics: the plot values are kept in the fit's meta data (sufficient statistics are not enough for
# the within-environment variances)
.bam_obs_stats <- function(fit, stat, cell) {
  y <- attr(fit, "meta")$y_by_cell
  if (is.null(y)) .bam_stop("bam_ppc", "the fit has no plot data; refit with bayesammi 1.0.0.")
  stat(y)
}

#-------------------------------------------------------------------------------------------------------------------
# draws for the packages posterior and bayesplot, and coda
#-------------------------------------------------------------------------------------------------------------------

#' The draws of a fit for posterior, bayesplot and coda
#'
#' @description
#' The draws of \code{\link{bayes_ammi}()} are a \code{coda::mcmc.list} (component \code{draws}). These methods hand
#' them to the packages \pkg{posterior} and \pkg{coda}; every plot of \pkg{bayesplot} (\code{mcmc_trace()},
#' \code{mcmc_areas()}, \code{mcmc_rank_overlay()}, ...) then works on them, as does \code{posterior::summarise_draws()}
#' (rank-normalised \eqn{\hat{R}}{R-hat}, bulk and tail effective sample sizes).
#'
#' @param x a result of \code{\link{bayes_ammi}()}.
#' @param ... unused.
#' @return \code{as_draws_df()}: a \code{posterior::draws_df}; \code{as_draws()}: a \code{posterior::draws_list};
#'   \code{as.mcmc.list()}: the \code{coda::mcmc.list}.
#' @references
#' Gabry, J., Simpson, D., Vehtari, A., Betancourt, M. and Gelman, A. (2019). Visualization in Bayesian workflow.
#' \emph{Journal of the Royal Statistical Society A}, \strong{182}(2), 389--402. \doi{10.1111/rssa.12378}
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C. (2021). Rank-normalization, folding, and
#' localization: an improved R-hat for assessing convergence of MCMC. \emph{Bayesian Analysis}, \strong{16}(2),
#' 667--718. \doi{10.1214/20-BA1221}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @examples
#' data(Maiz)
#' fit <-
#'   bayes_ammi(
#'       .data   = Maiz
#'     , .y      = y
#'     , .gen    = entry
#'     , .env    = site
#'     , .rep    = rep
#'     , .nIter  = 500
#'     , .k      = 2
#'     , .burnin = 100
#'     , .thin   = 1
#'     , .chains = 2
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
#' if (requireNamespace("posterior", quietly = TRUE)) {
#'   dr <- posterior::as_draws_df(x = fit)
#'   posterior::summarise_draws(
#'       .x = posterior::subset_draws(
#'           x        = dr
#'         , variable = c("mu", "lambda[1]", "lambda[2]")
#'         )
#'     )
#' }
#' if (requireNamespace("bayesplot", quietly = TRUE)) {
#'   bayesplot::mcmc_trace(x = fit$draws, pars = c("lambda[1]", "lambda[2]"))
#' }
#' @name bam-draws
NULL

#' @rdname bam-draws
#' @exportS3Method posterior::as_draws_df
as_draws_df.bam_bayes_ammi <- function(x, ...) posterior::as_draws_df(x$draws)

#' @rdname bam-draws
#' @exportS3Method posterior::as_draws
as_draws.bam_bayes_ammi <- function(x, ...) posterior::as_draws(x$draws)

#' @rdname bam-draws
#' @exportS3Method coda::as.mcmc.list
as.mcmc.list.bam_bayes_ammi <- function(x, ...) x$draws

.bam_print_compare <- function(x, digits = 1L, ...) {
  cat(attr(x, "label"), "\n\n")
  t <- data.table::copy(x$table)
  num <- c("DIC", "pD", "WAIC", "p_waic", "se", "dWAIC", "dWAIC_se", "lambda_last_lower", "lambda_last_mean")
  t[, (num) := lapply(.SD, round, digits = digits), .SDcols = num]
  print(t)
  invisible(x)
}

.bam_plot_compare <- function(x, ...) {
  d <- data.table::copy(x$table)
  d[, model := factor(model, levels = rev(model))]
  d[, dDIC := DIC - min(DIC)]
  d[, half := 2 * ifelse(is.na(dWAIC_se), se, dWAIC_se)]
  ggplot2::ggplot(d, ggplot2::aes(x = dWAIC, y = model)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = dWAIC - half, xmax = dWAIC + half), height = 0.25,
                            colour = .bam_palette[1L]) +
    ggplot2::geom_point(size = 2.2, colour = .bam_palette[1L]) +
    ggplot2::geom_point(ggplot2::aes(x = dDIC), shape = 4, size = 2.2, colour = .bam_palette[2L]) +
    .bam_x_cont() +
    ggplot2::labs(x = "Difference from the best model: WAIC (dot, with 2 standard errors) and DIC (cross)", y = NULL,
                  title = "Model comparison (smaller is better)") +
    bam_theme()
}

.bam_print_ppc <- function(x, digits = 4L, ...) {
  cat(attr(x, "label"), "\n\n")
  print(x$checks, digits = digits)
  invisible(x)
}

.bam_plot_ppc <- function(x, ...) {
  r <- data.table::melt(x$replicated[, lapply(.SD, as.numeric)], measure.vars = names(x$replicated),
                        variable.name = "statistic", value.name = "value")
  o <- data.table::copy(x$checks)
  o[, statistic := factor(statistic, levels = names(x$replicated))]
  o[, lab := sprintf("p = %.2f", p_value)]
  ggplot2::ggplot(r, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 30L, fill = "grey75", colour = "white", linewidth = 0.1) +
    ggplot2::geom_vline(data = o, ggplot2::aes(xintercept = observed), colour = .bam_palette[2L], linewidth = 0.7) +
    ggplot2::geom_text(data = o, ggplot2::aes(x = Inf, y = Inf, label = lab), hjust = 1.1, vjust = 1.5, size = 3,
                       family = "serif") +
    ggplot2::facet_wrap(ggplot2::vars(statistic), scales = "free", ncol = 2L) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(4)) +
    .bam_y_cont(bars = TRUE) +
    ggplot2::labs(x = NULL, y = "replicated trials", title = "Posterior predictive checks",
                  subtitle = "histogram: replicated trials; line: observed trial") +
    bam_theme()
}

utils::globalVariables(c("dWAIC", "dWAIC_se", "dDIC", "half", "model", "WAIC", "se", "DIC", "statistic", "observed", "p_value", "lab"))
