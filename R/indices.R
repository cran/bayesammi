#-------------------------------------------------------------------------------------------------------------------
# Classical stability indices made Bayesian: every index is a function of the cell means eta_ij, so its posterior
# distribution is the index computed on every draw of the cell means (of an AMMI, GGE or stability fit). AMMI-based
# indices are functions of the singular values and scores of an AMMI or GGE fit.
#-------------------------------------------------------------------------------------------------------------------

#' Posterior distributions of classical stability indices
#'
#' @description
#' \code{bam_indices()} computes five classical stability and superiority indices on every posterior draw of the cell
#' means of a fit, so that each index has a posterior mean, an HPD interval and, for every genotype, a probability of
#' being the best. The indices are those of Wricke (1962), Shukla (1972), Lin and Binns (1988), Annicchiarico (1992)
#' and Eskridge (1990).
#'
#' @details
#' With \eqn{\eta_{ij}}{eta_ij} the cell mean of genotype \eqn{i}{i} in environment \eqn{j}{j} (one posterior draw),
#' \eqn{\bar{\eta}_{i\cdot}}{eta_i.}, \eqn{\bar{\eta}_{\cdot j}}{eta_.j} and \eqn{\bar{\eta}_{\cdot\cdot}}{eta_..}
#' the genotype, environment and overall means, and the interaction
#' \eqn{\epsilon_{ij} = \eta_{ij} - \bar{\eta}_{i\cdot} - \bar{\eta}_{\cdot j} + \bar{\eta}_{\cdot\cdot}}{ge_ij =
#' eta_ij - eta_i. - eta_.j + eta_..}:
#' \describe{
#'   \item{\code{wricke}}{ecovalence \eqn{W_i = \sum_j \epsilon_{ij}^2}{W_i = sum_j ge_ij^2}; smaller is more stable.}
#'   \item{\code{shukla}}{stability variance \eqn{\sigma^2_i = \frac{G}{(G-2)(E-1)}W_i -
#'     \frac{\sum_k W_k}{(G-1)(G-2)(E-1)}}{sigma_i^2 = G W_i / ((G - 2)(E - 1)) - sum(W) / ((G - 1)(G - 2)(E - 1))};
#'     smaller is more stable.}
#'   \item{\code{lin_binns}}{superiority \eqn{P_i = \sum_j(\eta_{ij} - M_j)^2/(2E)}{P_i = sum_j (eta_ij - M_j)^2 /
#'     (2E)}, \eqn{M_j}{M_j} the best cell mean of environment \eqn{j}{j}; smaller is better.}
#'   \item{\code{annicchiarico}}{reliability \eqn{I_i = \bar{Z}_i - z_{1-p}s_{Z_i}}{I_i = mean(Z_i) - z_(1-p)
#'     sd(Z_i)} with \eqn{Z_{ij} = 100\,\eta_{ij}/\bar{\eta}_{\cdot j}}{Z_ij = 100 eta_ij / eta_.j}; larger is better
#'     (with \code{higher_better = FALSE}: \eqn{\bar{Z}_i + z_{1-p}s_{Z_i}}{mean + z sd}, smaller better).}
#'   \item{\code{safety_first}}{Eskridge's index with the variance across environments,
#'     \eqn{\bar{\eta}_{i\cdot} - z_{1-\alpha}s_i}{eta_i. - z_(1-alpha) s_i}; larger is better (with
#'     \code{higher_better = FALSE}: \eqn{\bar{\eta}_{i\cdot} + z_{1-\alpha}s_i}{eta_i. + z s_i}, smaller better).
#'     \code{\link{bayes_stability}()} gives the model-based versions of the regression and Shukla models.}
#' }
#' Because the indices are computed from the posterior of the cell means, they measure the stability of the expected
#' performance, free of the plot error; the HPD intervals show how well the trial determines them.
#'
#' @param fit a result of \code{\link{bayes_ammi}()} or \code{\link{bayes_stability}()}.
#' @param higher_better \code{TRUE} (default) when larger values of the trait are better.
#' @param .prob error probability \eqn{p}{p} of Annicchiarico's index (default 0.25).
#' @param .alpha probability \eqn{\alpha}{alpha} of the safety-first index (default 0.05).
#' @return An object of class \code{c("bam_indices", "bam_result", "list")} with
#' \describe{
#'   \item{\code{indices}}{\code{data.table} in long form: \code{index}, \code{Gen}, \code{mean}, \code{lower},
#'     \code{upper} (HPD interval), \code{prob_best} (posterior probability that the genotype is the best on this
#'     index) and \code{rank_mean} (expected rank, 1 = best).}
#'   \item{\code{means}}{\code{data.table}: one row per genotype with the posterior means of the indices and the mean
#'     of their expected ranks (\code{rank_overall}), sorted by it.}
#' }
#' Methods: \code{print()}, \code{plot()}, \code{\link{bam_interpret}()}, \code{\link{bam_cite}()}.
#' @references
#' Wricke, G. (1962). Uber eine Methode zur Erfassung der okologischen Streubreite in Feldversuchen.
#' \emph{Zeitschrift fur Pflanzenzuchtung}, \strong{47}, 92--96.
#'
#' Shukla, G. K. (1972). Some statistical aspects of partitioning genotype-environmental components of variability.
#' \emph{Heredity}, \strong{29}(2), 237--245. \doi{10.1038/hdy.1972.87}
#'
#' Lin, C. S. and Binns, M. R. (1988). A superiority measure of cultivar performance for cultivar x location data.
#' \emph{Canadian Journal of Plant Science}, \strong{68}(1), 193--198. \doi{10.4141/cjps88-018}
#'
#' Annicchiarico, P. (1992). Cultivar adaptation and recommendation from alfalfa trials in Northern Italy.
#' \emph{Journal of Genetics and Breeding}, \strong{46}, 269--278.
#'
#' Eskridge, K. M. (1990). Selection of stable cultivars using a safety-first rule. \emph{Crop Science},
#' \strong{30}, 369--374.
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bam_ammi_indices}}, \code{\link{bayes_stability}}, \code{\link{bam_superiority}}
#' @family stability
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
#' ind <-
#'   bam_indices(
#'       fit           = fit
#'     , higher_better = TRUE
#'     , .prob         = 0.25
#'     , .alpha        = 0.05
#'     )
#' ind$means
#' plot(x = ind)
#' @export
bam_indices <- function(fit, higher_better = TRUE, .prob = 0.25, .alpha = 0.05) {
  fn <- "bam_indices"
  .bam_need_fit(fit, fn)
  higher_better <- .bam_check_flag(higher_better, "higher_better", fn)
  .prob  <- .bam_check_number(.prob, ".prob", fn, min = 1e-4, max = 0.5)
  .alpha <- .bam_check_number(.alpha, ".alpha", fn, min = 1e-4, max = 0.5)
  d <- .bam_eta(fit)
  eta <- d$eta
  S <- dim(eta)[1L]
  G <- dim(eta)[2L]
  E <- dim(eta)[3L]
  if (G < 3L) .bam_stop(fn, "Shukla's stability variance needs at least 3 genotypes.")
  sgn <- if (higher_better) 1 else -1
  W <- matrix(rowSums(matrix(d$interaction^2, S * G, E)), S, G)
  shk <- G / ((G - 2) * (E - 1)) * W - rowSums(W) / ((G - 1) * (G - 2) * (E - 1))
  P <- matrix(0, S, G)
  for (j in seq_len(E)) {
    ej <- eta[, , j, drop = FALSE]
    dim(ej) <- c(S, G)
    best <- if (higher_better) apply(ej, 1L, max) else apply(ej, 1L, min)
    P <- P + (ej - best)^2
  }
  P <- P / (2 * E)
  envm <- vapply(seq_len(E), function(j) rowMeans(eta[, , j, drop = FALSE]), numeric(S))
  envm <- matrix(envm, S, E)
  Z <- eta
  for (j in seq_len(E)) Z[, , j] <- 100 * eta[, , j] / envm[, j]
  zm <- matrix(rowMeans(matrix(Z, S * G, E)), S, G)
  zs <- sqrt(.bam_row_var(Z))
  ann <- zm - sgn * stats::qnorm(1 - .prob) * zs
  sf  <- d$marginal - sgn * stats::qnorm(1 - .alpha) * sqrt(.bam_row_var(eta))
  idx <- list(wricke = W, shukla = shk, lin_binns = P, annicchiarico = ann, safety_first = sf)
  larger <- c(wricke = FALSE, shukla = FALSE, lin_binns = FALSE, annicchiarico = higher_better,
              safety_first = higher_better)
  .bam_index_result(idx, larger, d, "indices",
                    sprintf("Posterior stability indices (%s values of the trait are better)",
                            if (higher_better) "larger" else "smaller"),
                    match.call(), list(higher_better = higher_better, prob = .prob, alpha = .alpha))
}

#' Posterior distributions of AMMI-based stability indices
#'
#' @description
#' \code{bam_ammi_indices()} computes the AMMI stability value and the other stability indices built from the
#' singular values and genotype scores of an AMMI or GGE fit on every posterior draw: ASV (Purchase et al. 2000),
#' SIPC (Sneller et al. 1997), EV (Zobel 1994), Za (Zali et al. 2012), WAAS and the mean-stability index WAASY
#' (Olivoto et al. 2019).
#'
#' @details
#' With \eqn{K}{K} multiplicative terms, singular values \eqn{\lambda_k}{lambda_k}, genotype scores
#' \eqn{\alpha_{ik}}{alpha_ik}, coordinates \eqn{\mathrm{PC}_{ik} = \lambda_k^{1/2}\alpha_{ik}}{PC_ik =
#' lambda_k^(1/2) alpha_ik} and the share of term \eqn{k}{k} in the fitted interaction
#' \eqn{\theta_k = \lambda_k^2/\sum_m\lambda_m^2}{theta_k = lambda_k^2 / sum lambda^2}, on every draw:
#' \describe{
#'   \item{\code{ASV}}{\eqn{\sqrt{(\lambda_1^2/\lambda_2^2\,\mathrm{PC}_{i1})^2 + \mathrm{PC}_{i2}^2}}{sqrt((lambda_1^2
#'     / lambda_2^2 PC_i1)^2 + PC_i2^2)}, with the ratio of the posterior means of the singular values (as in
#'     \code{\link{bam_stability}()}; a ratio of draws explodes when \eqn{\lambda_2}{lambda_2} comes near zero);
#'     needs \eqn{K \ge 2}{K >= 2}.}
#'   \item{\code{SIPC}}{\eqn{\sum_k|\mathrm{PC}_{ik}|}{sum_k |PC_ik|}.}
#'   \item{\code{EV}}{\eqn{\sum_k\alpha_{ik}^2/K}{sum_k alpha_ik^2 / K}.}
#'   \item{\code{Za}}{\eqn{\sum_k\theta_k|\alpha_{ik}|}{sum_k theta_k |alpha_ik|}.}
#'   \item{\code{WAAS}}{\eqn{\sum_k\theta_k|\mathrm{PC}_{ik}|}{sum_k theta_k |PC_ik|} (the weights sum to one).}
#'   \item{\code{WAASY}}{\eqn{(\theta_Y R^Y_i + (100 - \theta_Y)R^W_i)/100}{(w R^Y_i + (100 - w) R^W_i) / 100}, where
#'     \eqn{R^Y_i}{R^Y_i} rescales the genotype means to 0-100 (100 for the best) and \eqn{R^W_i}{R^W_i} rescales WAAS
#'     to 0-100 (100 for the most stable); larger is better.}
#' }
#' All but WAASY are smaller for more stable genotypes. The shares \eqn{\theta_k}{theta_k} are those of the
#' \eqn{K}{K} fitted terms; Zali et al. (2012) and Olivoto et al. (2019) use the shares of the significant terms in
#' the whole interaction, which differ by a constant factor for Za and not at all for WAAS.
#'
#' @param fit a result of \code{\link{bayes_ammi}()}.
#' @param higher_better \code{TRUE} (default) when larger values of the trait are better (for WAASY).
#' @param .weight_y weight \eqn{\theta_Y}{theta_Y} (0-100) of the mean in WAASY (default 50, equal weights).
#' @return An object of class \code{c("bam_ammi_indices", "bam_result", "list")} with \code{indices} and
#'   \code{means}, as in \code{\link{bam_indices}()}.
#'   Methods: \code{print()}, \code{plot()}, \code{\link{bam_interpret}()}, \code{\link{bam_cite}()}.
#' @references
#' Purchase, J. L., Hatting, H. and van Deventer, C. S. (2000). Genotype x environment interaction of winter wheat
#' (Triticum aestivum L.) in South Africa: II. Stability analysis of yield performance. \emph{South African Journal
#' of Plant and Soil}, \strong{17}(3), 101--107.
#'
#' Sneller, C. H., Kilgore-Norquest, L. and Dombek, D. (1997). Repeatability of yield stability statistics in
#' soybean. \emph{Crop Science}, \strong{37}(2), 383--390.
#'
#' Zobel, R. W. (1994). Stress resistance and root systems. In \emph{Proceedings of the Workshop on Adaptation of
#' Plants to Soil Stress}, INTSORMIL Publication 94-2, 80--99. University of Nebraska-Lincoln.
#'
#' Zali, H., Farshadfar, E., Sabaghpour, S. H. and Karimizadeh, R. (2012). Evaluation of genotype x environment
#' interaction in chickpea using measures of stability from AMMI model. \emph{Annals of Biological Research},
#' \strong{3}(7), 3126--3136.
#'
#' Olivoto, T., Lucio, A. D. C., da Silva, J. A. G., Marchioro, V. S., de Souza, V. Q. and Jost, E. (2019). Mean
#' performance and stability in multi-environment trials I: combining features of AMMI and BLUP techniques.
#' \emph{Agronomy Journal}, \strong{111}(6), 2949--2960. \doi{10.2134/agronj2019.03.0220}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bam_indices}}, \code{\link{bam_stability}}
#' @family stability
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
#' ai <-
#'   bam_ammi_indices(
#'       fit           = fit
#'     , higher_better = TRUE
#'     , .weight_y     = 50
#'     )
#' ai$means
#' plot(x = ai)
#' @export
bam_ammi_indices <- function(fit, higher_better = TRUE, .weight_y = 50) {
  fn <- "bam_ammi_indices"
  .bam_need_ammi(fit, fn)
  higher_better <- .bam_check_flag(higher_better, "higher_better", fn)
  .weight_y <- .bam_check_number(.weight_y, ".weight_y", fn, min = 0, max = 100)
  d <- .bam_eta(fit)
  lam <- .bam_block(d$all, d$lay, "lambda")
  A <- .bam_block(d$all, d$lay, "alpha")
  S <- nrow(lam)
  K <- ncol(lam)
  G <- dim(A)[2L]
  th <- lam^2 / rowSums(lam^2)
  PC <- A
  for (k in seq_len(K)) PC[, , k] <- A[, , k] * sqrt(lam[, k])
  sumk <- function(f) {
    out <- matrix(0, S, G)
    for (k in seq_len(K)) out <- out + f(k)
    out
  }
  sipc <- sumk(function(k) abs(PC[, , k]))
  ev   <- sumk(function(k) A[, , k]^2) / K
  za   <- sumk(function(k) th[, k] * abs(A[, , k]))
  waas <- sumk(function(k) th[, k] * abs(PC[, , k]))
  resc <- function(M, larger) {
    lo <- apply(M, 1L, min)
    hi <- apply(M, 1L, max)
    r <- (M - lo) / pmax(hi - lo, .Machine$double.eps)
    100 * (if (larger) r else 1 - r)
  }
  waasy <- (.weight_y * resc(d$marginal, higher_better) + (100 - .weight_y) * resc(waas, FALSE)) / 100
  idx <- list(SIPC = sipc, EV = ev, Za = za, WAAS = waas, WAASY = waasy)
  larger <- c(SIPC = FALSE, EV = FALSE, Za = FALSE, WAAS = FALSE, WAASY = TRUE)
  if (K >= 2L) {
    wt <- mean(lam[, 1L])^2 / mean(lam[, 2L])^2
    asv <- sqrt((wt * PC[, , 1L])^2 + PC[, , 2L]^2)
    idx <- c(list(ASV = asv), idx)
    larger <- c(ASV = FALSE, larger)
  }
  .bam_index_result(idx, larger, d, "ammi_indices",
                    sprintf("Posterior AMMI-based stability indices (%d term(s); WAASY with weight %g on the mean)",
                            K, .weight_y),
                    match.call(), list(higher_better = higher_better, weight_y = .weight_y, k = K))
}

# the result of an index function: posterior mean, HPD interval, probability of being best and expected rank of
# every index (a list of draws x genotypes matrices) for every genotype
.bam_index_result <- function(idx, larger, d, cls, label, call, meta) {
  G <- length(d$G)
  long <- data.table::rbindlist(lapply(names(idx), function(nm) {
    M <- idx[[nm]]
    h <- .bam_hpd(M, d$level)
    X <- if (larger[[nm]]) -M else M
    best <- tabulate(max.col(-X, ties.method = "first"), nbins = G) / nrow(M)
    data.table::data.table(index = nm, Gen = d$G, mean = colMeans(M), lower = h[, 1L], upper = h[, 2L],
                           prob_best = best, rank_mean = colMeans(.bam_row_rank(X, "average")))
  }))
  long[, index := factor(index, levels = names(idx))]
  wide <- data.table::dcast(long, Gen ~ index, value.var = "mean")
  rk <- long[, list(rank_overall = mean(rank_mean)), by = Gen]
  wide <- merge(wide, rk, by = "Gen", sort = FALSE)
  data.table::setorder(wide, rank_overall)
  .bam_result(list(indices = long, means = wide), cls, label = label, call = call,
              meta = c(meta, list(level = d$level, larger = larger)))
}

.bam_print_indices <- function(x, digits = 4L, ...) {
  cat(attr(x, "label"), "\n\n")
  cat("Posterior means; rank_overall is the mean of the expected ranks over the indices (1 = best)\n")
  w <- data.table::copy(x$means)
  num <- names(w)[vapply(w, is.numeric, NA)]
  w[, (num) := lapply(.SD, signif, digits = digits), .SDcols = num]
  print(w, nrows = 500L)
  b <- x$indices[, .SD[which.max(prob_best)], by = index]
  cat("\nMost probably the best on each index:\n")
  print(b[, list(index, Gen, prob_best = round(prob_best, 3L), mean = signif(mean, digits), lower = signif(lower, digits),
                 upper = signif(upper, digits))], nrows = 50L)
  cat("\nAll HPD intervals and probabilities: x$indices\n")
  invisible(x)
}
.bam_print_ammi_indices <- .bam_print_indices

.bam_plot_indices <- function(x, ...) {
  d <- data.table::copy(x$indices)
  ord <- x$means$Gen
  d[, Gen := factor(Gen, levels = rev(ord))]
  ggplot2::ggplot(data = d, mapping = ggplot2::aes(x = mean, y = Gen)) +
    ggplot2::geom_errorbarh(mapping = ggplot2::aes(xmin = lower, xmax = upper), height = 0.3,
                            colour = .bam_palette[1L]) +
    ggplot2::geom_point(mapping = ggplot2::aes(size = prob_best), colour = .bam_palette[1L]) +
    ggplot2::facet_wrap(ggplot2::vars(index), scales = "free_x") +
    ggplot2::scale_size_continuous(name = "probability\nof being best", range = c(1, 4),
                                   breaks = scales::pretty_breaks(n = 3)) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 4), labels = .bam_short) +
    ggplot2::labs(x = "posterior mean and HPD interval", y = "genotype (ordered by the mean expected rank)",
                  title = attr(x, "label")) +
    bam_theme()
}
.bam_plot_ammi_indices <- .bam_plot_indices

utils::globalVariables(c("index", "prob_best", "rank_mean", "rank_overall", "lower", "upper"))
