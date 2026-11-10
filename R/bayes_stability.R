#' Priors of the Bayesian stability models
#'
#' @description
#' \code{bam_stab_priors()} sets the prior distributions of \code{\link{bayes_stability}()}: normal priors for the
#' overall mean and the contrasts of the genotype effects, environment effects (the environmental index) and slope
#' deviations, and scaled inverse chi-squared priors for the error variance and the deviation variances.
#'
#' @details
#' \deqn{\mu \sim N(0, \sigma^2_\mu), \quad \mathbf{K}_G\mathbf{g} \sim N(\mathbf{0}, \sigma^2_g\mathbf{I}), \quad
#'   \mathbf{K}_E\mathbf{h} \sim N(\mathbf{0}, \sigma^2_h\mathbf{I}), \quad \mathbf{K}_G\mathbf{b} \sim
#'   N(\mathbf{0}, \sigma^2_b\mathbf{I}),}{mu ~ N(0, v_mu), K_G g ~ N(0, v_g I), K_E h ~ N(0, v_h I), K_G b ~ N(0,
#'   v_b I),}
#' \deqn{\sigma^2 \sim \mathrm{Inv}\textrm{-}\chi^2(\nu, S^2_0), \qquad \sigma^2_{d,i} \sim
#'   \mathrm{Inv}\textrm{-}\chi^2(\nu_d, S^2_d),}{sigma^2 ~ Inv-chi^2(nu, S0^2), sigma_d,i^2 ~ Inv-chi^2(nu_d,
#'   S_d^2),}
#' where \eqn{\mathbf{K}_n}{K_n} are the orthonormal contrasts of \code{\link{matrix_k}()}, so that the effects and
#' the slope deviations sum to zero. Variances left \code{NULL} are set from the data: \eqn{10^4}{10^4} times the
#' variance of the response for the effects, and the least-squares estimates for \eqn{S^2_0}{S0^2} and
#' \eqn{S^2_d}{S_d^2}, so that the priors carry the information of \eqn{\nu}{nu} and \eqn{\nu_d}{nu_d}
#' observations whatever the unit of the trait.
#'
#' @param mu_var,g_var,h_var prior variances of the overall mean and of the contrasts of the genotype and environment
#'   effects (\code{NULL}: \eqn{10^4}{10^4} times the variance of the response).
#' @param b_var prior variance of the contrasts of the slope deviations \eqn{b_i}{b_i} (default 100; the slopes are
#'   \eqn{1 + b_i}{1 + b_i}, without unit).
#' @param nu,s2 degrees of freedom and scale of the prior of the error variance (\code{s2 = NULL}: the within-cell
#'   mean square).
#' @param nu_d,s2_d degrees of freedom and scale of the prior of the deviation variances (\code{s2_d = NULL}: the
#'   pooled least-squares deviation variance).
#' @return An object of class \code{c("bam_stab_priors", "list")}.
#' @references
#' Lian, L. and de los Campos, G. (2016). FW: an R package for Finlay-Wilkinson regression that incorporates
#' genomic/pedigree information and covariance structures between environments. \emph{G3: Genes, Genomes,
#' Genetics}, \strong{6}(3), 589--597. \doi{10.1534/g3.115.026328}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bayes_stability}}
#' @family stability
#' @examples
#' bam_stab_priors(
#'     mu_var = NULL
#'   , g_var  = NULL
#'   , h_var  = NULL
#'   , b_var  = 100
#'   , nu     = 1
#'   , s2     = NULL
#'   , nu_d   = 1
#'   , s2_d   = NULL
#'   )
#' @export
bam_stab_priors <- function(mu_var = NULL, g_var = NULL, h_var = NULL, b_var = 100, nu = 1, s2 = NULL, nu_d = 1,
                            s2_d = NULL) {
  fn <- "bam_stab_priors"
  chk <- function(v, nm) if (is.null(v)) NULL else .bam_check_number(v, nm, fn, min = 1e-12)
  structure(list(mu_var = chk(mu_var, "mu_var"), g_var = chk(g_var, "g_var"), h_var = chk(h_var, "h_var"),
                 b_var = chk(b_var, "b_var"), nu = chk(nu, "nu"), s2 = chk(s2, "s2"), nu_d = chk(nu_d, "nu_d"),
                 s2_d = chk(s2_d, "s2_d")),
            class = c("bam_stab_priors", "list"))
}

#' @export
print.bam_stab_priors <- function(x, ...) {
  f <- function(v, txt) if (is.null(v)) txt else format(v)
  cat("Priors of the Bayesian stability models\n")
  cat("  mu       ~ N(0, ", f(x$mu_var, "1e4 x var(y)"), ")\n", sep = "")
  cat("  K_G g    ~ N(0, ", f(x$g_var, "1e4 x var(y)"), " I)\n", sep = "")
  cat("  K_E h    ~ N(0, ", f(x$h_var, "1e4 x var(y)"), " I)\n", sep = "")
  cat("  K_G b    ~ N(0, ", f(x$b_var, ""), " I)   (slopes 1 + b_i)\n", sep = "")
  cat("  sigma^2  ~ Inv-chi^2(", f(x$nu, ""), ", ", f(x$s2, "within-cell mean square"), ")\n", sep = "")
  cat("  sigma_d^2 ~ Inv-chi^2(", f(x$nu_d, ""), ", ", f(x$s2_d, "least-squares deviation variance"), ")\n", sep = "")
  invisible(x)
}

#' Bayesian stability analysis: Finlay-Wilkinson, Eberhart-Russell and Shukla models
#'
#' @description
#' \code{bayes_stability()} fits by Gibbs sampling the joint regression models of Finlay and Wilkinson (1963) and
#' Eberhart and Russell (1966), and the heterogeneous interaction variance model of Shukla (1972), to the plots of a
#' multi-environment trial. Every stability parameter (slope, deviation variance, stability variance, safety-first
#' index) then has a posterior distribution, so genotypes can be compared by probabilities rather than by
#' significance tests. The data may be unbalanced: plots or whole genotype-environment cells may be missing.
#'
#' @details
#' For genotype \eqn{i}{i}, environment \eqn{j}{j} and plot \eqn{l}{l},
#' \deqn{y_{ijl} = \mu + g_i + (1 + b_i)h_j + d_{ij} + e_{ijl}, \qquad d_{ij} \sim N(0, \sigma^2_{d,i}), \quad
#'   e_{ijl} \sim N(0, \sigma^2),}{y_ijl = mu + g_i + (1 + b_i) h_j + d_ij + e_ijl, d_ij ~ N(0, sigma_d,i^2), e_ijl ~
#'   N(0, sigma^2),}
#' with \eqn{\sum_i g_i = \sum_j h_j = \sum_i b_i = 0}{sum g = sum h = sum b = 0}. The environment effects
#' \eqn{h_j}{h_j} are the environmental index, estimated jointly with the other parameters as in the Bayesian
#' Finlay-Wilkinson regression of Lian and de los Campos (2016) instead of being fixed at the environment means; the
#' slope of genotype \eqn{i}{i} is \eqn{1 + b_i}{1 + b_i}, with average 1.
#' \describe{
#'   \item{\code{"ER"}}{Eberhart-Russell: slopes and a deviation variance \eqn{\sigma^2_{d,i}}{sigma_d,i^2} for every
#'     genotype (the Bayesian version of Oliveira et al. 2018). A stable genotype has a slope near 1 and a small
#'     deviation variance.}
#'   \item{\code{"FW"}}{Finlay-Wilkinson: slopes and one deviation variance for all genotypes.}
#'   \item{\code{"Shukla"}}{no slopes (\eqn{b_i = 0}{b_i = 0}) and an interaction variance
#'     \eqn{\sigma^2_{d,i}}{sigma_d,i^2} for every genotype: Shukla's stability variance, as in the Bayesian
#'     heterogeneous variance model of Edwards and Jannink (2006).}
#' }
#' The sampler draws the overall mean, the genotype effects, the environmental index and the slopes from their
#' normal conditional distributions with the deviations \eqn{d_{ij}}{d_ij} integrated out (each observed cell mean
#' then has variance \eqn{\sigma^2_{d,i} + \sigma^2/n_{ij}}{sigma_d,i^2 + sigma^2 / n_ij}); then the deviations;
#' then the variances from their scaled inverse chi-squared conditionals. This partially collapsed Gibbs sampler
#' mixes well even when the deviations are large. Empty cells contribute no data; their deviations are drawn from the
#' prior, and \code{\link{bam_predict}()} predicts them.
#'
#' \code{stability} reports, for every genotype, the posterior mean and HPD interval of its mean
#' \eqn{\mu + g_i}{mu + g_i}, slope, deviation variance and the share of its variance across environments explained by
#' the regression,
#' \eqn{R^2_i = (1 + b_i)^2 s^2_h / ((1 + b_i)^2 s^2_h + \sigma^2_{d,i})}{R2_i = (1 + b_i)^2 s_h^2 / ((1 + b_i)^2
#' s_h^2 + sigma_d,i^2)} with \eqn{s^2_h}{s_h^2} the variance of the environmental index; the posterior probability
#' that the slope exceeds 1; and the safety-first index of Eskridge (1990),
#' \eqn{\mu + g_i - z_{1-\alpha}\sqrt{V_i}}{mu + g_i - z_(1-alpha) sqrt(V_i)}, with \eqn{V_i = b_i^2 s^2_h(1 - 1/E) +
#' \sigma^2_{d,i}}{V_i = b_i^2 s_h^2 (1 - 1/E) + sigma_d,i^2} (Eberhart-Russell), \eqn{b_i^2 s^2_h(1 - 1/E)}{b_i^2
#' s_h^2 (1 - 1/E)} (Finlay-Wilkinson) or \eqn{s^2_h + \sigma^2_{d,i}}{s_h^2 + sigma_d,i^2} (Shukla). When smaller
#' values are better (\code{higher_better = FALSE} in \code{\link{bam_interpret}()} and the decision functions), the
#' index is read the other way: \eqn{\mu + g_i + z_{1-\alpha}\sqrt{V_i}}{mu + g_i + z sqrt(V_i)}, smaller better.
#'
#' The result works with \code{\link{bam_predict}()}, \code{\link{bam_winners}()}, \code{\link{bam_superiority}()},
#' \code{\link{bam_ranks}()}, \code{\link{bam_indices}()}, \code{\link{bam_compare}()} and
#' \code{\link{bam_ppc}()}, which use the draws of the cell means \eqn{\mu + g_i + (1 + b_i)h_j + d_{ij}}{mu + g_i +
#' (1 + b_i) h_j + d_ij}.
#'
#' @param .data data frame with one row per plot.
#' @param .y,.gen,.env,.rep columns of the response, genotype, environment and replicate (bare names or strings).
#' @param .model \code{"ER"} (default), \code{"FW"} or \code{"Shukla"}; see Details.
#' @param .nIter number of iterations kept per chain before thinning (default 1000).
#' @param .burnin number of iterations discarded at the start of every chain (default 0).
#' @param .thin keep every \code{.thin}-th iteration (default 1).
#' @param .chains number of chains (default 1); the first starts at the least-squares estimates, the others at
#'   dispersed values.
#' @param .priors priors from \code{\link{bam_stab_priors}()}.
#' @param .level probability of the highest posterior density intervals (default 0.95).
#' @param .alpha probability of the safety-first index (default 0.05, as in Eskridge 1990).
#' @param .seed seed of the random numbers (the random number state of the session is restored afterwards).
#' @param .waic compute WAIC during sampling (default \code{TRUE}).
#' @return An object of class \code{c("bam_bayes_stability", "bam_result", "list")}:
#' \describe{
#'   \item{\code{stability}}{\code{data.table}, one row per genotype: \code{Gen}; \code{mean}, \code{mean_lower},
#'     \code{mean_upper}; \code{slope}, \code{slope_lower}, \code{slope_upper}, \code{prob_slope_above_1};
#'     \code{dev_var}, \code{dev_var_lower}, \code{dev_var_upper}; \code{R2}; \code{safety_first},
#'     \code{safety_first_lower}, \code{safety_first_upper}, \code{prob_best_safety}.}
#'   \item{\code{environments}}{\code{data.table}: \code{Env}, posterior mean and HPD interval of the environmental
#'     index \eqn{h_j}{h_j}, and the observed environment mean.}
#'   \item{\code{summary}}{\code{data.table} of every parameter (\code{mu}, \code{sigma2}, \code{g[i]},
#'     \code{h[j]}, \code{slope[i]}, \code{sigma2_d[i]}): posterior mean, sd, Monte Carlo error, HPD interval,
#'     effective sample size, potential scale reduction factor.}
#'   \item{\code{draws}}{\code{coda::mcmc.list} of these parameters.}
#'   \item{\code{deviations}}{array (draws x genotypes x environments) of the deviations \eqn{d_{ij}}{d_ij}.}
#'   \item{\code{diagnostics}}{Gelman-Rubin, Raftery-Lewis and Heidelberger-Welch diagnostics of the overall mean,
#'     error variance, slopes and deviation variances.}
#'   \item{\code{fit}}{DIC and WAIC, as in \code{\link{bayes_ammi}()}.}
#'   \item{\code{ls}}{least-squares estimates (the start of the first chain).}
#'   \item{\code{settings}}{iterations, burn-in, thinning, chains, model, priors, level, alpha, seed, time used.}
#' }
#' Methods: \code{print()}, \code{plot(x, which = c("slopes", "lines", "deviations", "trace", "density"))},
#' \code{\link{bam_interpret}()}, \code{\link{bam_cite}()}.
#' @references
#' Finlay, K. W. and Wilkinson, G. N. (1963). The analysis of adaptation in a plant-breeding programme.
#' \emph{Australian Journal of Agricultural Research}, \strong{14}(6), 742--754. \doi{10.1071/AR9630742}
#'
#' Eberhart, S. A. and Russell, W. A. (1966). Stability parameters for comparing varieties. \emph{Crop Science},
#' \strong{6}(1), 36--40. \doi{10.2135/cropsci1966.0011183X000600010011x}
#'
#' Shukla, G. K. (1972). Some statistical aspects of partitioning genotype-environmental components of variability.
#' \emph{Heredity}, \strong{29}(2), 237--245. \doi{10.1038/hdy.1972.87}
#'
#' Eskridge, K. M. (1990). Selection of stable cultivars using a safety-first rule. \emph{Crop Science},
#' \strong{30}, 369--374.
#'
#' Edwards, J. W. and Jannink, J.-L. (2006). Bayesian modeling of heterogeneous error and genotype x environment
#' interaction variances. \emph{Crop Science}, \strong{46}(2), 820. \doi{10.2135/cropsci2005.0164}
#'
#' Lian, L. and de los Campos, G. (2016). FW: an R package for Finlay-Wilkinson regression that incorporates
#' genomic/pedigree information and covariance structures between environments. \emph{G3: Genes, Genomes,
#' Genetics}, \strong{6}(3), 589--597. \doi{10.1534/g3.115.026328}
#'
#' Oliveira, T. R. A., Carvalho, H. W. L., Nascimento, M., Costa, E. F. N., Amaral Junior, A. T., Gravina, G. A. and
#' Carvalho Filho, J. L. S. (2018). The Eberhart and Russell's Bayesian method used as an instrument to select maize
#' hybrids. \emph{Euphytica}, \strong{214}, 64. \doi{10.1007/s10681-018-2146-y}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bam_stab_priors}}, \code{\link{bam_indices}}, \code{\link{bayes_ammi}}
#' @family stability
#' @examples
#' data(Maiz)
#' fs <-
#'   bayes_stability(
#'       .data   = Maiz
#'     , .y      = y
#'     , .gen    = entry
#'     , .env    = site
#'     , .rep    = rep
#'     , .model  = "ER"
#'     , .nIter  = 1000
#'     , .burnin = 200
#'     , .thin   = 1
#'     , .chains = 1
#'     , .priors = bam_stab_priors(
#'           mu_var = NULL
#'         , g_var  = NULL
#'         , h_var  = NULL
#'         , b_var  = 100
#'         , nu     = 1
#'         , s2     = NULL
#'         , nu_d   = 1
#'         , s2_d   = NULL
#'         )
#'     , .level  = 0.95
#'     , .alpha  = 0.05
#'     , .seed   = 1
#'     , .waic   = TRUE
#'     )
#' fs
#' plot(x = fs, which = "slopes")
#' @export
bayes_stability <- function(.data, .y, .gen, .env, .rep, .model = c("ER", "FW", "Shukla"), .nIter = 1000L,
                            .burnin = 0L, .thin = 1L, .chains = 1L, .priors = bam_stab_priors(), .level = 0.95,
                            .alpha = 0.05, .seed = NULL, .waic = TRUE) {
  fn <- "bayes_stability"
  t0 <- proc.time()[["elapsed"]]
  cd <- .bam_prep(fn, .data, substitute(.y), substitute(.gen), substitute(.env), substitute(.rep), parent.frame(),
                  allow_missing = TRUE)
  G <- length(cd$G)
  E <- length(cd$E)
  .model  <- .bam_check_choice(.model[1L], c("ER", "FW", "Shukla"), ".model", fn)
  .nIter  <- .bam_check_number(.nIter, ".nIter", fn, min = 2, integer = TRUE)
  .burnin <- .bam_check_number(.burnin, ".burnin", fn, min = 0, integer = TRUE)
  .thin   <- .bam_check_number(.thin, ".thin", fn, min = 1, max = .nIter, integer = TRUE)
  .chains <- .bam_check_number(.chains, ".chains", fn, min = 1, integer = TRUE)
  .level  <- .bam_check_number(.level, ".level", fn, min = 0.5, max = 0.999)
  .alpha  <- .bam_check_number(.alpha, ".alpha", fn, min = 1e-4, max = 0.5)
  .waic   <- .bam_check_flag(.waic, ".waic", fn)
  .bam_check_class(.priors, "bam_stab_priors", ".priors", fn)
  if (!is.null(.seed)) .seed <- .bam_check_number(.seed, ".seed", fn, integer = TRUE)
  if (cd$B < 2L) .bam_stop(fn, "the deviations and the error variance are separated by the replicates: every ",
                           "environment needs plots with replicates.")
  if (E < 4L) .bam_stop(fn, "at least 4 environments are needed for a regression on the environmental index.")
  code <- match(.model, c("ER", "FW", "Shukla")) - 1L
  ls <- .bam_stab_ls(cd, code)
  pr <- .bam_stab_fill_priors(.priors, cd, ls)
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
  ord <- order(cd$data$cell)
  chains <- lapply(seq_len(.chains), function(ch) {
    .bam_stab_cpp(cd$S, cd$Q, cd$N * 1, code, pr, .bam_stab_start(ls, ch, code), n_total, flag,
                  cd$data$y, cd$data$cell - 1L, .waic)
  })
  out <- .bam_stab_assemble(chains, cd, ls, code, keep, .level, .alpha)
  out$fit <- .bam_stab_fit_stats(chains, out, cd)
  out$settings <- list(nIter = .nIter, burnin = .burnin, thin = .thin, chains = .chains, kept = length(keep),
                       model = .model, priors = pr, level = .level, alpha = .alpha, seed = .seed,
                       seconds = proc.time()[["elapsed"]] - t0)
  nm <- c(ER = "Eberhart-Russell", FW = "Finlay-Wilkinson", Shukla = "Shukla")[[.model]]
  .bam_result(out, "bayes_stability",
              label = sprintf("Bayesian %s stability model: %d genotypes, %d environments, %s", nm, G, E,
                              if (cd$balanced) sprintf("%d replicates", cd$R) else
                                sprintf("%d plots (%d empty cells)", cd$n, cd$missing)),
              call = match.call(),
              meta = list(response = cd$names[["y"]], genotypes = G, environments = E, replicates = cd$B,
                          gen_labels = cd$G, env_labels = cd$E, level = .level, chains = .chains,
                          draws = length(keep) * .chains, model = .model, cells = list(Y = cd$Y, N = cd$N),
                          balanced = cd$balanced, missing = cd$missing, y_by_cell = cd$data$y[ord],
                          het = FALSE))
}

# least-squares start: additive fit (empty cells filled by EM), slopes by regression of the interaction on the
# environment effects, deviation variances from the residuals of that regression
.bam_stab_ls <- function(cd, code) {
  f <- .bam_ls(cd$Y, 0L)
  G <- length(cd$G)
  E <- length(cd$E)
  Z <- f$filled - f$fitted
  h <- f$delta
  b <- if (code == 2L) rep(0, G) else drop(Z %*% h) / sum(h^2)
  b <- b - mean(b)
  Dv <- Z - outer(b, h)
  obs <- cd$N > 0
  df_w <- cd$n - sum(obs)
  s2 <- if (df_w > 0) cd$ssw / df_w else stats::var(as.vector(Dv[obs])) / 2
  mean_inv_n <- mean(1 / cd$N[obs])
  sd2 <- vapply(seq_len(G), function(i) {
    o <- obs[i, ]
    max(sum(Dv[i, o]^2) / max(1, sum(o) - 1 - (code != 2L)) - s2 * mean_inv_n, 0.1 * s2 * mean_inv_n)
  }, numeric(1L))
  list(mu = f$mu, g = f$tau, h = h, b = b, sigma2 = s2, sd2 = sd2, Z = Z)
}

.bam_stab_fill_priors <- function(pr, cd, ls) {
  vy <- stats::var(cd$data$y)
  big <- 1e4 * vy
  list(mu_var = if (is.null(pr$mu_var)) big else pr$mu_var,
       g_var  = if (is.null(pr$g_var)) big else pr$g_var,
       h_var  = if (is.null(pr$h_var)) big else pr$h_var,
       b_var  = pr$b_var, nu = pr$nu,
       s2     = if (is.null(pr$s2)) ls$sigma2 else pr$s2,
       nu_d   = pr$nu_d,
       s2_d   = if (is.null(pr$s2_d)) mean(ls$sd2) else pr$s2_d)
}

# starting values of chain `ch`: the least-squares estimates for the first chain, dispersed values for the others
.bam_stab_start <- function(ls, ch, code) {
  G <- length(ls$g)
  E <- length(ls$h)
  common <- code == 1L
  sd2 <- if (common) mean(ls$sd2) else ls$sd2
  if (ch == 1L) return(list(mu = ls$mu, g = ls$g, h = ls$h, b = ls$b, sigma2 = ls$sigma2, sd2 = sd2))
  ctr <- function(v) v - mean(v)
  se_g <- sqrt(mean(ls$sd2) / E + ls$sigma2 / E)
  list(mu = ls$mu + 3 * stats::rnorm(1L) * se_g / sqrt(G),
       g = ls$g + ctr(stats::rnorm(G, 0, 2 * se_g)),
       h = ls$h + ctr(stats::rnorm(E, 0, 2 * sqrt(mean(ls$sd2) / G))),
       b = if (code == 2L) rep(0, G) else ls$b + ctr(stats::rnorm(G, 0, 0.1)),
       sigma2 = ls$sigma2 * exp(stats::rnorm(1L, 0, 0.3)),
       sd2 = sd2 * exp(stats::rnorm(length(sd2), 0, 0.3)))
}

# The R reference of the compiled sampler: the same conditional distributions with the random numbers drawn in the
# same order, written directly from the formulas of the help page (used by the tests)
.bam_stab_chain <- function(S, Q, N, code, pr, start, n_total, keep) {
  G <- nrow(S)
  E <- ncol(S)
  KG <- matrix_k(G)
  KE <- matrix_k(E)
  slopes <- code != 2L
  common <- code == 1L
  obs <- N > 0
  Yb <- ifelse(obs, S / pmax(N, 1), 0)
  ssw <- sum((Q - S^2 / pmax(N, 1))[obs])
  ntot <- sum(N)
  mu <- start$mu
  g <- start$g
  h <- start$h
  b <- if (slopes) start$b else rep(0, G)
  s2 <- start$sigma2
  sd2 <- start$sd2
  D <- matrix(0, G, E)
  draw_ct <- function(K, w, s, v) {
    P <- K %*% (w * t(K))
    diag(P) <- diag(P) + 1 / v
    R <- chol(P)
    m <- backsolve(R, forwardsolve(t(R), K %*% s))
    z <- stats::rnorm(nrow(K))
    drop(t(K) %*% (m + backsolve(R, z)))
  }
  out <- list()
  for (it in seq_len(n_total)) {
    sdv <- if (common) rep(sd2, G) else sd2
    W <- ifelse(obs, 1 / (sdv + s2 / pmax(N, 1)), 0)
    cc <- 1 + b
    r <- Yb - g - outer(cc, h)
    prec <- sum(W) + 1 / pr$mu_var
    mu <- sum(W * r) / prec + stats::rnorm(1L) / sqrt(prec)
    r <- Yb - mu - outer(cc, h)
    g <- draw_ct(KG, rowSums(W), rowSums(W * r), pr$g_var)
    r <- Yb - mu - g
    h <- draw_ct(KE, colSums(W * cc^2), colSums(W * r * cc), pr$h_var)
    if (slopes) {
      r <- Yb - mu - outer(g, h, `+`)
      b <- draw_ct(KG, drop(W %*% h^2), drop((W * r) %*% h), pr$b_var)
      cc <- 1 + b
    }
    eta <- mu + g + outer(cc, h)
    for (j in seq_len(E)) for (i in seq_len(G)) {
      if (obs[i, j]) {
        a <- N[i, j] / s2
        p <- a + 1 / sdv[i]
        D[i, j] <- a * (Yb[i, j] - eta[i, j]) / p + stats::rnorm(1L) / sqrt(p)
      } else {
        D[i, j] <- stats::rnorm(1L) * sqrt(sdv[i])
      }
    }
    m <- eta + D
    ss <- ssw + sum((N * (Yb - m)^2)[obs])
    s2 <- (pr$nu * pr$s2 + ss) / stats::rchisq(1L, pr$nu + ntot)
    if (common) {
      sd2 <- (pr$nu_d * pr$s2_d + sum(D^2)) / stats::rchisq(1L, pr$nu_d + G * E)
    } else {
      for (i in seq_len(G)) sd2[i] <- (pr$nu_d * pr$s2_d + sum(D[i, ]^2)) / stats::rchisq(1L, pr$nu_d + E)
    }
    if (keep[it]) out[[length(out) + 1L]] <- list(mu = mu, sigma2 = s2, g = g, h = h, b = b, sd2 = sd2,
                                                  d = as.vector(D))
  }
  pick <- function(nm) do.call(rbind, lapply(out, `[[`, nm))
  list(mu = pick("mu")[, 1L], sigma2 = pick("sigma2")[, 1L], g = pick("g"), h = pick("h"), b = pick("b"),
       sd2 = pick("sd2"), d = pick("d"))
}

# the components of the result from the chains
.bam_stab_assemble <- function(chains, cd, ls, code, keep, level, alpha) {
  G <- length(cd$G)
  E <- length(cd$E)
  nk <- length(keep)
  common <- code == 1L
  nsd <- if (common) 1L else G
  par <- c("mu", "sigma2", sprintf("g[%d]", seq_len(G)), sprintf("h[%d]", seq_len(E)),
           if (code != 2L) sprintf("slope[%d]", seq_len(G)),
           if (common) "sigma2_d" else sprintf("sigma2_d[%d]", seq_len(G)))
  lay <- data.table::data.table(
      parameter = par
    , block     = c("mu", "sigma2", rep("g", G), rep("h", E), if (code != 2L) rep("slope", G), rep("sigma2_d", nsd))
    , label     = c(NA, NA, cd$G, cd$E, if (code != 2L) cd$G, if (common) NA else cd$G)
    , axis      = NA_integer_
    )
  mats <- lapply(chains, function(ch) {
    m <- cbind(ch$mu, ch$sigma2, ch$g, ch$h, if (code != 2L) 1 + ch$b, ch$sd2)
    colnames(m) <- par
    coda::mcmc(m, start = keep[1L], thin = if (nk > 1L) keep[2L] - keep[1L] else 1L)
  })
  draws <- coda::mcmc.list(mats)
  all <- do.call(rbind, lapply(mats, unclass))
  dev <- do.call(rbind, lapply(chains, `[[`, "d"))
  S <- nrow(all)
  Darr <- array(dev, c(S, G, E), dimnames = list(NULL, cd$G, cd$E))
  blk <- function(b) all[, lay$parameter[lay$block == b], drop = FALSE]
  mu <- all[, "mu"]
  gm <- mu + blk("g")
  hh <- blk("h")
  sl <- if (code != 2L) blk("slope") else matrix(1, S, G)
  sd2 <- blk("sigma2_d")
  if (common) sd2 <- sd2[, rep(1L, G), drop = FALSE]
  sh2 <- rowSums(hh^2) / (E - 1)
  bb <- sl - 1
  V <- switch(code + 1L,
              bb^2 * sh2 * (1 - 1 / E) + sd2,
              bb^2 * sh2 * (1 - 1 / E),
              sh2 + sd2)
  z <- stats::qnorm(1 - alpha)
  sf <- gm - z * sqrt(V)
  r2 <- sl^2 * sh2 / (sl^2 * sh2 + sd2)
  hp <- function(M) .bam_hpd(M, level)
  hg <- hp(gm)
  hs <- hp(sl)
  hd <- hp(sd2)
  hf <- hp(sf)
  stab <- data.table::data.table(
      Gen = cd$G, mean = colMeans(gm), mean_lower = hg[, 1L], mean_upper = hg[, 2L]
    , slope = colMeans(sl), slope_lower = hs[, 1L], slope_upper = hs[, 2L], prob_slope_above_1 = colMeans(sl > 1)
    , dev_var = colMeans(sd2), dev_var_lower = hd[, 1L], dev_var_upper = hd[, 2L], R2 = colMeans(r2)
    , safety_first = colMeans(sf), safety_first_lower = hf[, 1L], safety_first_upper = hf[, 2L]
    , prob_best_safety = tabulate(max.col(sf, ties.method = "first"), nbins = G) / S
    )
  if (code == 2L) stab[, c("slope", "slope_lower", "slope_upper", "prob_slope_above_1", "R2") := NULL]
  hh_h <- hp(hh)
  envm <- colSums(cd$S) / colSums(cd$N)
  envs <- data.table::data.table(Env = cd$E, index = colMeans(hh), lower = hh_h[, 1L], upper = hh_h[, 2L],
                                 env_mean = envm)
  list(stability = stab, environments = envs, summary = .bam_summary_table(draws, all, lay, level), draws = draws,
       deviations = Darr,
       diagnostics = .bam_diagnostics(coda::mcmc.list(lapply(mats, function(m)
         coda::mcmc(m[, lay$block %in% c("mu", "sigma2", "slope", "sigma2_d"), drop = FALSE], start = stats::start(m),
                    thin = coda::thin(m)))), sum(lay$block %in% c("mu", "sigma2", "slope", "sigma2_d"))),
       ls = list(mu = ls$mu, g = stats::setNames(ls$g, cd$G), h = stats::setNames(ls$h, cd$E),
                 slope = stats::setNames(1 + ls$b, cd$G), dev_var = stats::setNames(ls$sd2, cd$G),
                 sigma2 = ls$sigma2),
       layout = lay)
}

# short axis labels for large numbers: 500K, 1M
.bam_short <- function(x) scales::label_number(scale_cut = scales::cut_short_scale())(x)

# DIC (deviance at the posterior means of the cell means and the error variance) and WAIC
.bam_stab_fit_stats <- function(chains, out, cd) {
  dev <- unlist(lapply(chains, `[[`, "deviance"))
  eta <- apply(.bam_eta(.bam_stab_as_fit(out, cd))$eta, c(2L, 3L), mean)
  s2 <- mean(out$summary[parameter == "sigma2", mean])
  N <- cd$N
  dhat <- sum((N * log(2 * pi * s2) + (cd$Q - 2 * eta * cd$S + N * eta^2) / s2)[N > 0])
  pd <- mean(dev) - dhat
  c(list(deviance = dev, Dbar = mean(dev), Dhat = dhat, pD = pd, DIC = dhat + 2 * pd), .bam_waic(chains))
}

# a minimal fit object for .bam_eta() while the result is being assembled
.bam_stab_as_fit <- function(out, cd) {
  structure(out, class = c("bam_bayes_stability", "bam_result", "list"),
            meta = list(genotypes = length(cd$G), environments = length(cd$E), gen_labels = cd$G, env_labels = cd$E,
                        level = 0.95))
}

# draws of the cell means of a stability fit: mu + g_i + slope_i h_j + d_ij
.bam_stab_eta <- function(fit) {
  all <- do.call(rbind, lapply(fit$draws, unclass))
  lay <- fit$layout
  blk <- function(b) all[, lay$parameter[lay$block == b], drop = FALSE]
  G <- length(lay$parameter[lay$block == "g"])
  E <- length(lay$parameter[lay$block == "h"])
  gm <- all[, "mu"] + blk("g")
  sl <- if (any(lay$block == "slope")) blk("slope") else matrix(1, nrow(all), G)
  hh <- blk("h")
  eta <- fit$deviations
  for (j in seq_len(E)) eta[, , j] <- eta[, , j] + gm + sl * hh[, j]
  eta
}

#-------------------------------------------------------------------------------------------------------------------
# print and plot
#-------------------------------------------------------------------------------------------------------------------
.bam_print_bayes_stability <- function(x, digits = 4L, ...) {
  st <- x$settings
  cat(attr(x, "label"), "\n")
  cat(sprintf("%d chain(s) x %d kept draws (iterations %d after a burn-in of %d, thinning %d); %.1f s\n\n",
              st$chains, st$kept, st$nIter, st$burnin, st$thin, st$seconds))
  cat(sprintf("Genotypes: posterior means and %g%% HPD intervals; safety-first index at alpha = %g\n",
              100 * st$level, st$alpha))
  s <- data.table::copy(x$stability)
  num <- names(s)[vapply(s, is.numeric, NA)]
  s[, (num) := lapply(.SD, signif, digits = digits), .SDcols = num]
  print(s, nrows = 500L)
  sm <- x$summary[block %in% c("mu", "sigma2")]
  cat(sprintf("\nOverall mean %s, error variance %s; largest potential scale reduction factor %s\n",
              format(signif(sm$mean[1L], digits)), format(signif(sm$mean[2L], digits)),
              if (all(is.na(x$summary$rhat))) "(one chain)" else format(round(max(x$summary$rhat, na.rm = TRUE), 3L))))
  cat("All parameters: x$summary; environmental index: x$environments\n")
  invisible(x)
}

.bam_plot_bayes_stability <- function(x, which = c("slopes", "lines", "deviations", "trace", "density"), pars = NULL,
                                      ...) {
  which <- match.arg(which)
  s <- x$stability
  if (which %in% c("trace", "density")) {
    if (is.null(pars)) pars <- intersect(c("mu", "sigma2", "slope[1]", "sigma2_d[1]", "sigma2_d"), x$summary$parameter)
    return(.bam_plot_draws(x, which, pars))
  }
  if (which == "slopes") {
    if (!"slope" %in% names(s)) .bam_stop("plot", "the Shukla model has no slopes; use which = \"deviations\".")
    gm <- mean(s$mean)
    return(ggplot2::ggplot(data = s, mapping = ggplot2::aes(x = mean, y = slope)) +
             ggplot2::geom_hline(yintercept = 1, linewidth = 0.3, colour = "grey50", linetype = 2) +
             ggplot2::geom_vline(xintercept = gm, linewidth = 0.3, colour = "grey50", linetype = 2) +
             ggplot2::geom_errorbar(mapping = ggplot2::aes(ymin = slope_lower, ymax = slope_upper), width = 0,
                                    colour = .bam_palette[1L], alpha = 0.6) +
             ggplot2::geom_errorbarh(mapping = ggplot2::aes(xmin = mean_lower, xmax = mean_upper), height = 0,
                                     colour = .bam_palette[1L], alpha = 0.6) +
             ggplot2::geom_point(mapping = ggplot2::aes(size = dev_var), colour = .bam_palette[1L]) +
             ggrepel::geom_text_repel(mapping = ggplot2::aes(label = Gen), size = 3, family = "serif",
                                      max.overlaps = 50L) +
             ggplot2::scale_size_continuous(name = "deviation\nvariance", range = c(1.5, 5),
                                            breaks = scales::pretty_breaks(n = 3), labels = .bam_short) +
             .bam_x_cont() + .bam_y_cont() +
             ggplot2::labs(x = "genotype mean", y = "slope on the environmental index",
                           title = "Adaptability and stability (posterior means, HPD intervals)",
                           subtitle = "Slope 1: average response; point size: deviation variance") +
             bam_theme())
  }
  if (which == "deviations") {
    d <- data.table::copy(s)[order(dev_var)]
    d[, Gen := factor(Gen, levels = Gen)]
    return(ggplot2::ggplot(data = d, mapping = ggplot2::aes(x = dev_var, y = Gen)) +
             ggplot2::geom_errorbarh(mapping = ggplot2::aes(xmin = dev_var_lower, xmax = dev_var_upper), height = 0.25,
                                     colour = .bam_palette[1L]) +
             ggplot2::geom_point(size = 2, colour = .bam_palette[1L]) +
             .bam_x_cont(labels = .bam_short) +
             ggplot2::labs(x = if (x$settings$model == "Shukla") "stability variance" else "deviation variance",
                           y = NULL, title = "Genotype-specific variances (posterior means, HPD intervals)") +
             bam_theme())
  }
  # regression lines on the environmental index with the observed cell means
  e <- x$environments
  Y <- attr(x, "meta")$cells$Y
  mu <- x$summary[parameter == "mu", mean]
  sl <- if ("slope" %in% names(s)) s$slope else rep(1, nrow(s))
  pts <- data.table::data.table(Gen = rep(s$Gen, nrow(e)), index = rep(e$index, each = nrow(s)), y = as.vector(Y))
  lines <- data.table::data.table(Gen = s$Gen, intercept = s$mean, slope = sl)
  ggplot2::ggplot(data = pts, mapping = ggplot2::aes(x = index, y = y, colour = Gen)) +
    ggplot2::geom_abline(data = lines, mapping = ggplot2::aes(intercept = intercept, slope = slope, colour = Gen),
                         linewidth = 0.5) +
    ggplot2::geom_point(size = 1.2, alpha = 0.7, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = grDevices::colorRampPalette(.bam_palette)(nrow(s)), name = "genotype") +
    .bam_x_cont() + .bam_y_cont() +
    ggplot2::labs(x = "environmental index (posterior mean)", y = "cell mean",
                  title = "Regression on the environmental index",
                  subtitle = sprintf("Lines: posterior mean of genotype mean + slope x index (overall mean %s)",
                                     format(signif(mu, 4L)))) +
    bam_theme()
}

utils::globalVariables(c("slope", "slope_lower", "slope_upper", "mean_lower", "mean_upper", "dev_var",
                         "dev_var_lower", "dev_var_upper", "Gen", "index", "intercept", "parameter", "safety_first",
                         "prob_best_safety"))
