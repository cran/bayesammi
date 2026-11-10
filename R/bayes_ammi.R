#' Bayesian estimation of the AMMI and GGE models
#'
#' @description
#' \code{bayes_ammi()} estimates the additive main effects and multiplicative interaction (AMMI) model of a
#' multi-environment trial by the Gibbs sampler of Crossa et al. (2011), or the GGE (sites regression, SREG) model
#' of Perez-Elizalde et al. (2012), for balanced or unbalanced trials, with one error variance or one for every
#' environment. The sampler is compiled (C++). Every parameter (overall mean, genotype and
#' environment effects, singular values, genotype and environment scores of the bilinear terms, error variance) has a
#' posterior distribution, so that the biplot carries statistical inference: highest posterior density (HPD)
#' regions around the scores show which genotypes and environments interact significantly and which form
#' homogeneous groups (\code{\link{biplots}}).
#'
#' @details
#' \strong{Model.} The response of genotype \eqn{i}{i} in replicate \eqn{l}{l} of environment \eqn{j}{j} is
#' \deqn{y_{ijl} = \mu + \tau_i + \delta_j + \sum_{k=1}^{K} \lambda_k \alpha_{ik} \gamma_{jk} + e_{ijl}, \qquad
#'   e_{ijl} \sim N(0, \sigma^2),}{y_ijl = mu + tau_i + delta_j + sum_k lambda_k alpha_ik gamma_jk + e_ijl, e_ijl ~
#'   N(0, sigma^2),}
#' with \eqn{\sum_i \tau_i = \sum_j \delta_j = 0}{sum tau_i = sum delta_j = 0}, \eqn{\lambda_1 \ge \cdots \ge
#' \lambda_K > 0}{lambda_1 >= ... >= lambda_K > 0}, and genotype and environment score vectors
#' \eqn{\boldsymbol{\alpha}_k}{alpha_k}, \eqn{\boldsymbol{\gamma}_k}{gamma_k} of unit length, orthogonal to
#' \eqn{\mathbf{1}}{1} and to each other (\eqn{K}{K} = \code{.k}). The priors are those of \code{\link{bam_priors}()}.
#'
#' \strong{Conditional distributions.} With \eqn{\mathbf{Y}}{Y} the \eqn{G \times E}{G x E} table of cell means,
#' \eqn{B}{B} replicates, \eqn{n = BGE}{n = BGE} observations, the orthonormal contrasts
#' \eqn{\mathbf{K}_G}{K_G} of \code{\link{matrix_k}()} (\eqn{\boldsymbol{\tau} = \mathbf{K}_G^\top
#' \mathbf{t}}{tau = K_G' t}) and \eqn{c_x = \sigma^2 / \sigma^2_x}{c_x = sigma^2 / v_x}, each iteration draws
#' \deqn{\mu \sim N\left(\frac{n\bar{y} + c_\mu \mu_0}{n + c_\mu}, \frac{\sigma^2}{n + c_\mu}\right), \qquad
#'   \mathbf{t} \sim N\left(\frac{BE\,\mathbf{K}_G\bar{\mathbf{y}}_{G} + c_\tau \mathbf{K}_G\mathbf{m}_\tau}{BE +
#'   c_\tau}, \frac{\sigma^2}{BE + c_\tau}\mathbf{I}\right),}{mu ~ N((n ybar + c_mu mu0)/(n + c_mu), sigma^2/(n +
#'   c_mu)); t ~ N((BE K_G ybar_G + c_tau K_G m_tau)/(BE + c_tau), sigma^2/(BE + c_tau) I),}
#' the same for \eqn{\boldsymbol{\delta}}{delta} with \eqn{BG}{BG}, and for \eqn{k = 1, \ldots, K}{k = 1, ..., K}
#' \deqn{\lambda_k \sim N\left(\frac{B\,\boldsymbol{\alpha}_k^\top\mathbf{Y}\boldsymbol{\gamma}_k + c_\lambda
#'   \lambda_0}{B + c_\lambda}, \frac{\sigma^2}{B + c_\lambda}\right) \textrm{ on } (\lambda_{k+1},
#'   \lambda_{k-1}),}{lambda_k ~ N((B alpha_k' Y gamma_k + c_lambda lambda0)/(B + c_lambda), sigma^2/(B + c_lambda))
#'   on (lambda_k+1, lambda_k-1),}
#' \deqn{p(\boldsymbol{\alpha}_k \mid \cdot) \propto \exp\left(\frac{B\lambda_k}{\sigma^2}
#'   \boldsymbol{\alpha}_k^\top\mathbf{Y}\boldsymbol{\gamma}_k\right), \qquad p(\boldsymbol{\gamma}_k \mid \cdot)
#'   \propto \exp\left(\frac{B\lambda_k}{\sigma^2}\boldsymbol{\gamma}_k^\top\mathbf{Y}^\top
#'   \boldsymbol{\alpha}_k\right),}{p(alpha_k | .) ~ exp(B lambda_k / sigma^2 alpha_k' Y gamma_k), p(gamma_k | .) ~
#'   exp(B lambda_k / sigma^2 gamma_k' Y' alpha_k),}
#' von Mises-Fisher distributions on the unit sphere of the subspace orthogonal to \eqn{\mathbf{1}}{1} and to the
#' other score vectors (Crossa et al. 2011, eq. 8-9), and finally
#' \deqn{\sigma^2 \sim \frac{\nu S^2_0 + \mathrm{SSW} + B\sum_{ij}(\bar{y}_{ij} - \hat{y}_{ij})^2}{\chi^2_{\nu +
#'   n}},}{sigma^2 ~ (nu s0^2 + SSW + B sum (y_ij - fit_ij)^2) / chi^2(nu + n),}
#' with SSW the within-cell sum of squares. The chains start at the least-squares AMMI fit (the first chain) or at
#' dispersed values (the others).
#'
#' \strong{A variance for every environment} (\code{.het = TRUE}): \eqn{e_{ijl} \sim N(0, \sigma^2_j)}{e_ijl ~ N(0,
#' sigma_j^2)} with independent priors \eqn{\mathrm{Inv}\textrm{-}\chi^2(\nu, S^2_{0j})}{Inv-chi^2(nu, s0_j^2)}.
#' With weights \eqn{w_j = 1/\sigma^2_j}{w_j = 1 / sigma_j^2} the genotype effects keep a closed normal
#' conditional (precision \eqn{B\sum_j w_j}{B sum w_j}); the overall mean and the environment effects are drawn
#' jointly; \eqn{\lambda_k}{lambda_k} has precision \eqn{B\sum_j w_j\gamma_{jk}^2}{B sum w_j gamma_jk^2}; the
#' genotype scores stay von Mises-Fisher; the environment scores have the Bingham-von Mises-Fisher conditional
#' \deqn{p(\boldsymbol{\gamma}_k \mid \cdot) \propto \exp\left(B\lambda_k\sum_j w_j\gamma_{jk}r_{\cdot jk} -
#'   \frac{B\lambda_k^2}{2}\sum_j w_j\gamma_{jk}^2\right),}{p(gamma_k | .) ~ exp(B lambda_k sum_j w_j gamma_jk r_jk -
#'   B lambda_k^2 / 2 sum_j w_j gamma_jk^2),}
#' (\eqn{r_{\cdot jk} = \sum_i\alpha_{ik}r_{ijk}}{r_jk = sum_i alpha_ik r_ijk} with \eqn{r_{ijk}}{r_ijk} the cell
#' residual without term \eqn{k}{k}), drawn by a Gibbs sweep in the eigenbasis of its quadratic term (Hoff 2009);
#' every \eqn{\sigma^2_j}{sigma_j^2} is a scaled inverse chi-square with \eqn{\nu + BG}{nu + BG} degrees of freedom.
#'
#' \strong{GGE model} (\code{.model = "GGE"}): \eqn{y_{ijl} = \mu + \delta_j + \sum_k\lambda_k\alpha_{ik}\gamma_{jk}
#' + e_{ijl}}{y_ijl = mu + delta_j + sum_k lambda_k alpha_ik gamma_jk + e_ijl}: the genotype main effect is part of the
#' multiplicative terms, the genotype scores are orthogonal to \eqn{\mathbf{1}}{1}, the environment scores only to
#' each other (Perez-Elizalde et al. 2012). The conditional distributions are those above without
#' \eqn{\boldsymbol{\tau}}{tau} and without the constraint \eqn{\boldsymbol{\gamma}_k \perp \mathbf{1}}{gamma_k
#' orthogonal to 1}.
#'
#' \strong{Unbalanced trials.} Missing plots and whole missing genotype-environment cells are completed by data
#' augmentation: with \eqn{B}{B} the largest number of replicates of a cell, every iteration first draws the missing
#' plots from \eqn{N(\eta_{ij}, \sigma^2_j)}{N(eta_ij, sigma_j^2)}, then the conditional distributions of the balanced
#' design apply exactly (the plots are assumed missing at random). Every genotype and every environment needs some
#' data.
#'
#' \strong{Model choice.} With \code{.waic = TRUE} the sampler accumulates the pointwise log-likelihood of the observed
#' plots, giving WAIC; the deviance of every kept draw gives DIC (component \code{fit}; \code{\link{bam_compare}()}).
#'
#' \strong{Signs.} The pair \eqn{(\boldsymbol{\alpha}_k, \boldsymbol{\gamma}_k)}{(alpha_k, gamma_k)} and its
#' negative give the same model. With \code{.align = TRUE} (default) every draw takes the sign that agrees with the
#' least-squares vectors, so that posterior means and HPD regions refer to one orientation of each axis; with
#' \code{.align = FALSE} a weak axis can switch sign during the run and its posterior means shrink towards zero.
#'
#' \strong{Changes from the earlier versions.} The singular values were computed once and never sampled; the error
#' variance was not returned (\code{tau1} held the constant precision of the least-squares fit); the residuals for
#' \eqn{\sigma^2}{sigma^2} were computed against the wrong environments; the within-cell sum of squares was a sum of
#' variances; the covariance matrices of the main effects were square roots of variances; the scores were not drawn
#' from their constrained conditional distributions; all \eqn{\min(G, E) - 1}{min(G, E) - 1} axes were fitted.
#' \code{tau1} now holds the draws of the error variance \eqn{\sigma^2}{sigma^2}.
#'
#' @param .data a data frame with one row per plot.
#' @param .y,.gen,.env,.rep the columns of the response, genotypes, environments and replicates: bare names, strings
#'   or vectors.
#' @param .nIter number of iterations kept per chain before thinning (after the burn-in).
#' @param .k number of bilinear terms \eqn{K}{K} (default 2, as in Crossa et al. 2011).
#' @param .burnin number of iterations discarded at the start of every chain (default 0).
#' @param .thin keep every \code{.thin}-th iteration (default 1).
#' @param .chains number of chains (default 1; two or more give the Gelman-Rubin diagnostic).
#' @param .priors the priors, from \code{\link{bam_priors}()}.
#' @param .level probability of the HPD intervals (default 0.95).
#' @param .seed integer seed for a reproducible run, or \code{NULL}; the random number state of the session is
#'   restored afterwards.
#' @param .align keep the sign of every bilinear term (see Signs; default \code{TRUE}).
#' @param .model \code{"AMMI"} (default) or \code{"GGE"} (sites regression).
#' @param .het \code{FALSE} (default): one error variance; \code{TRUE}: one for every environment.
#' @param .waic compute WAIC during sampling (default \code{TRUE}; costs one pass over the plots per kept draw).
#' @return An object of class \code{c("bam_bayes_ammi", "bam_result", "list")}; the components of the earlier
#'   versions come first, the draws of all chains one after the other:
#' \describe{
#'   \item{\code{mu1}}{\code{data.table}, column \code{mu}: draws of the overall mean.}
#'   \item{\code{tau1}}{\code{data.table}, column \code{tau}: draws of the error variance \eqn{\sigma^2}{sigma^2}
#'     (columns \code{tau1}, ..., \code{tauE} with a variance for every environment).}
#'   \item{\code{tao1}}{\code{data.table}, columns \code{tao1}, ..., \code{taoG}: draws of the genotype effects (no
#'     columns for the GGE model).}
#'   \item{\code{delta1}}{\code{data.table}, columns \code{delta1}, ...: draws of the environment effects.}
#'   \item{\code{lambdas1}}{\code{data.table}, columns \code{lambdas1}, ...: draws of the singular values.}
#'   \item{\code{alphas1}, \code{gammas1}}{matrices with columns \code{.nIter} (draw), \code{Gen} or \code{Env}
#'     (number) and \code{alphas1}, ... or \code{gammas1}, ...: draws of the genotype and environment scores.}
#'   \item{\code{alphas0}, \code{gammas0}}{\code{data.table}s: least-squares scores (the starting values).}
#'   \item{\code{summary}}{\code{data.table}: posterior mean, standard deviation, Monte Carlo standard error, HPD
#'     interval, effective sample size and potential scale reduction factor of every parameter (Tables 1 and B1
#'     of Crossa et al. 2011).}
#'   \item{\code{scores}}{\code{data.table}: posterior means of the biplot coordinates
#'     \eqn{\lambda_k^{1/2}\alpha_{ik}}{lambda_k^(1/2) alpha_ik} and \eqn{\lambda_k^{1/2}\gamma_{jk}}{lambda_k^(1/2)
#'     gamma_jk}, with the mean of every genotype and environment.}
#'   \item{\code{draws}}{\code{coda::mcmc.list} of all parameters, one element per chain.}
#'   \item{\code{diagnostics}}{list: \code{gelman} (potential scale reduction factors, two or more chains),
#'     \code{raftery} (Raftery-Lewis dependence factors), \code{heidel} (Heidelberger-Welch stationarity and
#'     half-width tests), for the overall mean, the error variance, the main effects and the singular values.}
#'   \item{\code{ls}}{list: the least-squares AMMI estimates (\code{mu}, \code{tau}, \code{delta}, \code{lambda},
#'     \code{alpha}, \code{gamma}, \code{sigma2}).}
#'   \item{\code{fit}}{list: \code{deviance} of every kept draw, \code{Dbar}, \code{Dhat}, \code{pD}, \code{DIC},
#'     and with \code{.waic = TRUE} \code{lppd}, \code{p_waic}, \code{WAIC}, \code{waic_se} and \code{waic_i}
#'     (the contribution of every observed plot to WAIC, in the order of the data).}
#'   \item{\code{settings}}{list: iterations, burn-in, thinning, chains, number of terms, model, variances, priors
#'     (with the scale of the error variance filled in), level, seed, alignment, time used.}
#' }
#'   Methods: \code{print()}, \code{summary()}, \code{plot()} (\code{which = c("trace", "density", "biplot",
#'   "effects")}), \code{\link{biplots}()}, \code{\link{bam_interpret}()}, the decision tools
#'   (\code{\link{bam_predict}()}, \code{\link{bam_winners}()}, \code{\link{bam_superiority}()},
#'   \code{\link{bam_ranks}()}, \code{\link{bam_stability}()}), \code{\link{bam_compare}()}, \code{\link{bam_ppc}()} and
#'   \code{posterior::as_draws_df()} (\code{\link{bam-draws}}).
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#'
#' Viele, K. and Srinivasan, C. (2000). Parsimonious estimation of multiplicative interaction in analysis of variance
#' using Kullback-Leibler information. \emph{Journal of Statistical Planning and Inference}, \strong{84}(1--2),
#' 201--219.
#'
#' Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A general Bayesian estimation method of linear-bilinear
#' models applied to plant breeding trials with genotype x environment interaction. \emph{Journal of Agricultural,
#' Biological, and Environmental Statistics}, \strong{17}, 15--37. \doi{10.1007/s13253-011-0063-9}
#'
#' Hoff, P. D. (2009). Simulation of the matrix Bingham-von Mises-Fisher distribution, with applications to
#' multivariate and relational data. \emph{Journal of Computational and Graphical Statistics}, \strong{18}(2),
#' 438--456. \doi{10.1198/jcgs.2009.07177}
#'
#' Wood, A. T. A. (1994). Simulation of the von Mises Fisher distribution. \emph{Communications in Statistics -
#' Simulation and Computation}, \strong{23}(1), 157--164.
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com}), Jose Crossa (\email{j.crossa@@cgiar.org}), Sergio
#'   Perez-Elizalde (\email{sergiop@@colpos.mx}), Diego Jarquin (\email{diego.jarquin@@gmail.com}), Jose Miguel
#'   Cotes, Kert Viele, Genzhou Liu, Paul L. Cornelius
#' @seealso \code{\link{biplots}}, \code{\link{bam_priors}}, \code{\link{bam_interpret}}, \code{\link{ge_ammi}}
#' @examples
#' data(Maiz)
#' fm1 <-
#'   bayes_ammi(
#'       .data   = Maiz
#'     , .y      = y
#'     , .gen    = entry
#'     , .env    = site
#'     , .rep    = rep
#'     , .nIter  = 200
#'     , .k      = 2
#'     , .burnin = 50
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
#' fm1
#' head(fm1$summary)
#' fm1$fit$WAIC
#' plot(x = fm1, which = "trace")
#' plot(x = fm1, which = "biplot")
#' @export
bayes_ammi <- function(.data, .y, .gen, .env, .rep, .nIter = 1000L, .k = 2L, .burnin = 0L, .thin = 1L,
                       .chains = 1L, .priors = bam_priors(), .level = 0.95, .seed = NULL, .align = TRUE,
                       .model = c("AMMI", "GGE"), .het = FALSE, .waic = TRUE) {
  UseMethod("bayes_ammi")
}

#' @rdname bayes_ammi
#' @export
bayes_ammi.default <- function(.data, .y, .gen, .env, .rep, .nIter = 1000L, .k = 2L, .burnin = 0L, .thin = 1L,
                               .chains = 1L, .priors = bam_priors(), .level = 0.95, .seed = NULL, .align = TRUE,
                               .model = c("AMMI", "GGE"), .het = FALSE, .waic = TRUE) {
  fn  <- "bayes_ammi"
  t0  <- proc.time()[["elapsed"]]
  cd  <- .bam_prep(fn, .data, substitute(.y), substitute(.gen), substitute(.env), substitute(.rep), parent.frame(),
                   allow_missing = TRUE)
  G <- length(cd$G)
  E <- length(cd$E)
  .model  <- .bam_check_choice(.model[1L], c("AMMI", "GGE"), ".model", fn)
  gge     <- .model == "GGE"
  .nIter  <- .bam_check_number(.nIter, ".nIter", fn, min = 2, integer = TRUE)
  .k      <- .bam_check_number(.k, ".k", fn, min = 1, max = if (gge) min(G - 1L, E) else min(G, E) - 1L,
                               integer = TRUE)
  .burnin <- .bam_check_number(.burnin, ".burnin", fn, min = 0, integer = TRUE)
  .thin   <- .bam_check_number(.thin, ".thin", fn, min = 1, max = .nIter, integer = TRUE)
  .chains <- .bam_check_number(.chains, ".chains", fn, min = 1, integer = TRUE)
  .level  <- .bam_check_number(.level, ".level", fn, min = 0.5, max = 0.999)
  .align  <- .bam_check_flag(.align, ".align", fn)
  .het    <- .bam_check_flag(.het, ".het", fn)
  .waic   <- .bam_check_flag(.waic, ".waic", fn)
  .bam_check_class(.priors, "bam_priors", ".priors", fn)
  if (!is.null(.seed)) .seed <- .bam_check_number(.seed, ".seed", fn, integer = TRUE)
  if (.het && cd$B < 2L) .bam_stop(fn, "environment-specific variances (.het = TRUE) need replicates.")
  K  <- .k
  ls <- .bam_ls(cd$Y, K, gge)
  ls$gge <- gge
  # residual mean square of the least-squares fit: start of the chains and default scale of the prior
  df_res <- cd$n - 1 - (E - 1) - (if (gge) 0 else G - 1) - K * (G + E - 1 - 2 * K + gge) - cd$missing
  if (df_res < 1) .bam_stop(fn, "too few observations for ", K, " multiplicative terms.")
  resid  <- cd$ssw + sum((cd$N * (cd$Y - ls$fitted)^2)[cd$N > 0])
  s2_ls  <- resid / df_res
  s2_env <- if (.het) {
    vapply(seq_len(E), function(j) {
      w <- cd$data$Env == cd$E[j]
      r <- cd$data$y[w] - ls$fitted[cd$data$cell[w]]
      max(sum(r^2) / max(1, sum(w) - G), 1e-8 * s2_ls)
    }, numeric(1L))
  } else s2_ls
  pr <- .priors
  if (is.null(pr$s2)) pr$s2 <- s2_env
  if (.het) pr$s2 <- rep_len(pr$s2, E)
  # reproducible runs leave the random number state of the session as it was
  if (!is.null(.seed)) {
    had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else
      rm(".Random.seed", envir = globalenv()), add = TRUE)
    set.seed(.seed)
  }
  n_total <- .burnin + .nIter
  keep    <- seq.int(.burnin + .thin, n_total, by = .thin)
  flag    <- integer(n_total)
  flag[keep] <- 1L
  chains  <- lapply(seq_len(.chains), function(ch) {
    .bam_run_chain(cd, K, gge, .het, pr, .bam_start(ls, s2_env, ch, .het, E), ls$alpha, n_total, flag, .align,
                   .waic)
  })
  lay <- .bam_layout(cd$G, cd$E, K, gge, .het)
  out <- .bam_assemble(chains, cd, ls, K, keep, .level, lay)
  out$ls <- list(mu = ls$mu, tau = stats::setNames(ls$tau, cd$G), delta = stats::setNames(ls$delta, cd$E),
                 lambda = ls$lambda, alpha = ls$alpha, gamma = ls$gamma, sigma2 = s2_ls)
  dimnames(out$ls$alpha) <- list(cd$G, paste0("axis", seq_len(K)))
  dimnames(out$ls$gamma) <- list(cd$E, paste0("axis", seq_len(K)))
  out$fit <- .bam_fit_stats(chains, out, cd, lay)
  out$settings <- list(nIter = .nIter, burnin = .burnin, thin = .thin, chains = .chains, kept = length(keep),
                       k = K, model = .model, het = .het, priors = pr, level = .level, seed = .seed, align = .align,
                       seconds = proc.time()[["elapsed"]] - t0)
  .bam_result(out, "bayes_ammi",
              label = sprintf("Bayesian %s with %d multiplicative term(s)%s: %d genotypes, %d environments, %s",
                              .model, K, if (.het) " and a variance per environment" else "", G, E,
                              if (cd$balanced) sprintf("%d replicates", cd$R) else
                                sprintf("%d plots (%d empty cells)", cd$n, cd$missing)),
              call = match.call(),
              meta = list(response = cd$names[["y"]], genotypes = G, environments = E, replicates = cd$B,
                          gen_labels = cd$G, env_labels = cd$E, k = K, level = .level, chains = .chains,
                          draws = length(keep) * .chains, model = .model, het = .het, layout = lay,
                          cells = list(Y = cd$Y, N = cd$N), balanced = cd$balanced, missing = cd$missing,
                          y_by_cell = cd$data$y[order(cd$data$cell)]))
}

# the components of the result from the chains
.bam_assemble <- function(chains, cd, ls, K, keep, level, lay) {
  G  <- length(cd$G)
  E  <- length(cd$E)
  nk <- length(keep)
  mats <- lapply(chains, function(ch) {
    m <- cbind(ch$scalar, matrix(ch$alpha, nk), matrix(ch$gamma, nk))
    colnames(m) <- lay$parameter
    coda::mcmc(m, start = keep[1L], thin = if (nk > 1L) keep[2L] - keep[1L] else 1L)
  })
  draws <- coda::mcmc.list(mats)
  all   <- do.call(rbind, lapply(mats, unclass))
  N     <- nrow(all)
  long  <- function(Arr, n) {
    out <- matrix(0, N * n, K)
    for (k in seq_len(K)) out[, k] <- as.vector(t(Arr[, , k, drop = TRUE]))
    cbind(.nIter = rep(seq_len(N), each = n), rep(seq_len(n), N), out)
  }
  alphas1 <- long(.bam_block(all, lay, "alpha"), G)
  gammas1 <- long(.bam_block(all, lay, "gamma"), E)
  colnames(alphas1) <- c(".nIter", "Gen", paste0("alphas", seq_len(K)))
  colnames(gammas1) <- c(".nIter", "Env", paste0("gammas", seq_len(K)))
  dt <- function(M, nm) data.table::setnames(data.table::as.data.table(M), nm)
  sig <- .bam_block(all, lay, "sigma2")
  tau <- if (any(lay$block == "tau")) .bam_block(all, lay, "tau") else matrix(numeric(0), N, 0L)
  out <- list(
      mu1      = dt(all[, "mu", drop = FALSE], "mu")
    , tau1     = dt(sig, if (ncol(sig) == 1L) "tau" else paste0("tau", seq_len(ncol(sig))))
    , tao1     = if (ncol(tau)) dt(tau, paste0("tao", seq_len(G))) else data.table::data.table()
    , delta1   = dt(.bam_block(all, lay, "delta"), paste0("delta", seq_len(E)))
    , lambdas1 = dt(.bam_block(all, lay, "lambda"), paste0("lambdas", seq_len(K)))
    , alphas1  = alphas1
    , gammas1  = gammas1
    , alphas0  = dt(ls$alpha, paste0("alphas", seq_len(K)))
    , gammas0  = dt(ls$gamma, paste0("gammas", seq_len(K)))
    )
  out$summary     <- .bam_summary_table(draws, all, lay, level)
  out$scores      <- .bam_scores(all, lay, cd, K)
  out$draws       <- draws
  out$diagnostics <- .bam_diagnostics(draws, sum(!lay$block %in% c("alpha", "gamma")))
  out
}

# posterior summaries of every parameter: mean, sd, Monte Carlo standard error, HPD interval, effective sample size,
# potential scale reduction factor
.bam_summary_table <- function(draws, all, lay, level) {
  hpd <- .bam_hpd(all, level)
  ess <- coda::effectiveSize(draws)
  rhat <- if (length(draws) > 1L) {
    tryCatch(coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1L],
             error = function(e) rep(NA_real_, ncol(all)))
  } else rep(NA_real_, ncol(all))
  sdv <- apply(all, 2L, stats::sd)
  cbind(data.table::copy(lay),
        data.table::data.table(mean = colMeans(all), sd = sdv, mc_se = sdv / sqrt(pmax(ess, 1)), lower = hpd[, 1L],
                               upper = hpd[, 2L], ess = unname(ess), rhat = unname(rhat)))
}

# deviance information criterion and WAIC from the accumulators of the chains
.bam_fit_stats <- function(chains, out, cd, lay) {
  dev  <- unlist(lapply(chains, `[[`, "deviance"))
  all  <- .bam_all(out)
  pm   <- colMeans(all)
  # deviance at the posterior means of mu, tau, delta, sigma2 and of the cell means eta
  G <- length(cd$G)
  E <- length(cd$E)
  eta <- vapply(seq_len(E), function(j) colMeans(.bam_eta_env(all, lay, j, G)), numeric(G))
  s2  <- pm[lay$parameter[lay$block == "sigma2"]]
  s2  <- if (length(s2) == 1L) rep(s2, E) else s2
  N <- cd$N
  dhat <- sum(N * log(2 * pi * rep(s2, each = G)) + (cd$Q - 2 * eta * cd$S + N * eta^2) / rep(s2, each = G))
  pd   <- mean(dev) - dhat
  c(list(deviance = dev, Dbar = mean(dev), Dhat = dhat, pD = pd, DIC = dhat + 2 * pd), .bam_waic(chains))
}

# WAIC from the accumulators of the chains (log-sum-exp, running mean and sum of squares of the pointwise
# log-likelihood of every observed plot); an empty list when WAIC was not computed
.bam_waic <- function(chains) {
  res  <- list()
  lse  <- lapply(chains, `[[`, "lse")
  if (all(is.finite(lse[[1L]]))) {
    S  <- vapply(chains, function(ch) length(ch$deviance), numeric(1L))
    lt <- do.call(cbind, lse)
    mx <- apply(lt, 1L, max)
    lppd_i <- mx + log(rowSums(exp(lt - mx))) - log(sum(S))
    # pooled mean and variance of the pointwise log-likelihood over the chains
    mt <- do.call(cbind, lapply(chains, `[[`, "lmean"))
    m2 <- do.call(cbind, lapply(chains, `[[`, "lm2"))
    gm <- drop(mt %*% S) / sum(S)
    ss <- rowSums(m2) + drop(((mt - gm)^2) %*% S)
    p_i <- ss / (sum(S) - 1)
    res$lppd   <- sum(lppd_i)
    res$p_waic <- sum(p_i)
    res$WAIC   <- -2 * (res$lppd - res$p_waic)
    res$waic_i  <- -2 * (lppd_i - p_i)
    res$waic_se <- sqrt(length(p_i) * stats::var(res$waic_i))
  }
  res
}

# shortest intervals containing the proportion `level` of the draws (the HPD interval of a unimodal posterior;
# Chen and Shao 1999, the algorithm of Crossa et al. 2011)
.bam_hpd <- function(M, level) {
  M <- as.matrix(M)
  n <- nrow(M)
  s <- max(1L, min(n - 1L, ceiling(level * n) - 1L))
  t(apply(M, 2L, function(v) {
    v <- sort(v)
    w <- v[(s + 1L):n] - v[seq_len(n - s)]
    j <- which.min(w)
    c(v[j], v[j + s])
  }))
}

# biplot coordinates: posterior means of lambda_k^(1/2) alpha_ik and lambda_k^(1/2) gamma_jk
.bam_scores <- function(all, lay, cd, K) {
  lam <- .bam_block(all, lay, "lambda")
  A   <- .bam_block(all, lay, "alpha")
  Gm  <- .bam_block(all, lay, "gamma")
  ga <- vapply(seq_len(K), function(k) colMeans(A[, , k, drop = TRUE] * sqrt(lam[, k])), numeric(dim(A)[2L]))
  ea <- vapply(seq_len(K), function(k) colMeans(Gm[, , k, drop = TRUE] * sqrt(lam[, k])), numeric(dim(Gm)[2L]))
  ga <- matrix(ga, ncol = K)
  ea <- matrix(ea, ncol = K)
  colnames(ga) <- colnames(ea) <- paste0("PC", seq_len(K))
  data.table::data.table(type = rep(c("Genotype", "Environment"), c(nrow(ga), nrow(ea))), label = c(cd$G, cd$E),
                         mean = c(rowMeans(cd$Y, na.rm = TRUE), colMeans(cd$Y, na.rm = TRUE)), rbind(ga, ea))
}

# convergence diagnostics of the scalar parameters (Crossa et al. 2011: Raftery-Lewis, Heidelberger-Welch,
# Gelman-Rubin)
.bam_diagnostics <- function(draws, n_scalar) {
  sub <- coda::as.mcmc.list(lapply(draws, function(m) m[, seq_len(n_scalar), drop = FALSE]))
  safe <- function(expr) tryCatch(expr, error = function(e) NULL)
  list(
      gelman  = if (length(sub) > 1L) safe(coda::gelman.diag(sub, autoburnin = FALSE, multivariate = FALSE))
    , raftery = safe(lapply(sub, coda::raftery.diag))
    , heidel  = safe(lapply(sub, coda::heidel.diag))
    )
}

.bam_print_bayes_ammi <- function(x, digits = 4L, ...) {
  s  <- x$summary
  st <- x$settings
  cat(attr(x, "label"), "\n")
  cat(sprintf("%d chain(s) x %d kept draws (iterations %d after a burn-in of %d, thinning %d); %.1f s\n\n",
              st$chains, st$kept, st$nIter, st$burnin, st$thin, st$seconds))
  cat(sprintf("Posterior means, standard deviations and %g%% HPD intervals\n", 100 * st$level))
  tab <- s[block %in% c("mu", "sigma2", "lambda", "tau", "delta"),
           list(parameter, mean, sd, lower, upper, ess = round(ess), rhat)]
  if (all(is.na(tab$rhat))) tab[, rhat := NULL]
  print(tab, digits = digits, nrows = 200L)
  cat("\nScores of the bilinear terms: x$summary (blocks alpha, gamma); biplot coordinates: x$scores\n")
  invisible(x)
}

.bam_plot_bayes_ammi <- function(x, which = c("trace", "density", "biplot", "effects"), pars = NULL, ...) {
  which <- match.arg(which)
  if (which == "biplot") return(.bam_plot_mean_biplot(x))
  if (which == "effects") return(.bam_plot_effects(x))
  if (is.null(pars)) pars <- c("mu", "sigma2", sprintf("lambda[%d]", seq_len(x$settings$k)))
  .bam_plot_draws(x, which, pars)
}

# traces or histograms of the draws of chosen parameters (a result with $draws and $summary)
.bam_plot_draws <- function(x, which, pars) {
  bad <- setdiff(pars, x$summary$parameter)
  if (length(bad)) .bam_stop("plot", "unknown parameter(s): ", paste(bad, collapse = ", "),
                             "; see x$summary$parameter.")
  d <- data.table::rbindlist(lapply(seq_along(x$draws), function(ch) {
    m  <- x$draws[[ch]][, pars, drop = FALSE]
    # traces show at most 1000 evenly spaced draws per chain (enough for the eye, small files)
    ix <- if (which == "trace" && nrow(m) > 1000L) unique(round(seq(1, nrow(m), length.out = 1000L))) else
      seq_len(nrow(m))
    m  <- m[ix, , drop = FALSE]
    data.table::data.table(chain = factor(ch), draw = rep(ix, length(pars)),
                           parameter = factor(rep(pars, each = nrow(m)), levels = pars), value = as.vector(m))
  }))
  if (which == "trace") {
    return(ggplot2::ggplot(data = d, mapping = ggplot2::aes(x = draw, y = value, colour = chain)) +
             ggplot2::geom_line(linewidth = 0.25) +
             ggplot2::facet_wrap(ggplot2::vars(parameter), scales = "free_y", ncol = 1L) +
             ggplot2::scale_colour_manual(values = .bam_palette, name = "chain") +
             .bam_x_cont() +
             ggplot2::labs(x = "kept draw", y = NULL, title = "Traces of the Markov chains") +
             bam_theme())
  }
  ggplot2::ggplot(data = d, mapping = ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 40L, fill = .bam_palette[1L], colour = "white", linewidth = 0.1) +
    ggplot2::facet_wrap(ggplot2::vars(parameter), scales = "free", ncol = 2L) +
    .bam_x_cont() + .bam_y_cont(bars = TRUE) +
    ggplot2::labs(x = NULL, y = "draws", title = "Marginal posterior distributions") +
    bam_theme()
}

# the biplot of the posterior means of the coordinates (Crossa et al. 2011, Fig. 2)
.bam_plot_mean_biplot <- function(x) {
  s <- x$scores
  if (x$settings$k < 2L) .bam_stop("plot", "the biplot needs at least two bilinear terms (.k >= 2).")
  lim <- max(abs(c(s$PC1, s$PC2))) * 1.08
  sm  <- x$summary[block == "lambda"]
  e   <- s[type == "Environment"]
  ggplot2::ggplot(data = s, mapping = ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_segment(data = e, mapping = ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2),
                          colour = .bam_palette[2L], linewidth = 0.3,
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"))) +
    ggplot2::geom_point(mapping = ggplot2::aes(colour = type, shape = type), size = 2) +
    ggrepel::geom_text_repel(mapping = ggplot2::aes(label = label, colour = type), size = 3, family = "serif",
                             max.overlaps = 50L, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(Genotype = .bam_palette[1L], Environment = .bam_palette[2L]), name = NULL) +
    ggplot2::scale_shape_manual(values = c(Genotype = 16L, Environment = 17L), name = NULL) +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::labs(x = expression(lambda[1]^{1/2} * alpha[i1] ~ ", " ~ lambda[1]^{1/2} * gamma[j1]),
                  y = expression(lambda[2]^{1/2} * alpha[i2] ~ ", " ~ lambda[2]^{1/2} * gamma[j2]),
                  title = "Bayesian AMMI biplot of the posterior means",
                  subtitle = sprintf("lambda_1 = %s, lambda_2 = %s (posterior means)",
                                     format(signif(sm$mean[1L], 4L)), format(signif(sm$mean[2L], 4L)))) +
    bam_theme()
}

# posterior means and HPD intervals of the main effects
.bam_plot_effects <- function(x) {
  s <- x$summary[block %in% c("tau", "delta")]
  s[, block := factor(block, levels = c("tau", "delta"), labels = c("genotype effects", "environment effects"))]
  s[, label := factor(label, levels = unique(label))]
  ggplot2::ggplot(data = s, mapping = ggplot2::aes(x = mean, y = label)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_errorbarh(mapping = ggplot2::aes(xmin = lower, xmax = upper), height = 0.3,
                            colour = .bam_palette[1L]) +
    ggplot2::geom_point(size = 1.6, colour = .bam_palette[1L]) +
    ggplot2::facet_wrap(ggplot2::vars(block), scales = "free", ncol = 2L) +
    .bam_x_cont() +
    ggplot2::labs(x = sprintf("posterior mean and %g%% HPD interval", 100 * x$settings$level), y = NULL,
                  title = "Main effects") +
    bam_theme()
}

utils::globalVariables(c("sd", "block", "parameter", "lower", "upper", "ess", "rhat", "draw", "value", "chain", "type",
                         "PC1", "PC2", "label", ".nIter"))
