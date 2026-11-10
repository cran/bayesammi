#' Least-squares AMMI model
#'
#' @description
#' \code{ge_ammi()} computes the least-squares estimates of the AMMI model with all interaction axes: grand mean,
#' genotype and environment effects, singular values and vectors of the interaction, and the error variance. They
#' are the starting values of \code{\link{bayes_ammi}()} and the classical estimates it is compared with.
#'
#' @details
#' With the cell means \eqn{\bar{y}_{ij}}{y_ij} of \eqn{B}{B} replicates, \eqn{\hat{\mu} = \bar{y}}{mu = ybar},
#' \eqn{\hat{\tau}_i = \bar{y}_{i\cdot} - \bar{y}}{tau_i = y_i. - ybar}, \eqn{\hat{\delta}_j = \bar{y}_{\cdot j} -
#' \bar{y}}{delta_j = y_.j - ybar} and the singular value decomposition of the interaction effects
#' (\code{\link{ge_eff}()}), \eqn{K = \min(G, E) - 1}{K = min(G, E) - 1} axes. The error variance is the residual
#' variance of the model with random replicates within environments (\code{\link{ge_model}()}), in closed form for
#' balanced data:
#' \deqn{\hat{\sigma}^2 = \frac{\mathrm{SSW} - \mathrm{SS}_{B(E)}}{E(B - 1)(G - 1)}}{sigma^2 = (SSW - SS_B(E)) /
#'   (E (B - 1) (G - 1))}
#' when the replicate mean square is not smaller, otherwise (the replicate variance component estimated as zero)
#' \eqn{\mathrm{SSW}/(GE(B - 1))}{SSW / (GE (B - 1))}; SSW is the within-cell and \eqn{\mathrm{SS}_{B(E)}}{SS_B(E)}
#' the replicate-within-environment sum of squares.
#'
#' \strong{Change from the earlier versions.} The error variance came from \code{lme4::lmer()}; the closed form gives
#' the same value for balanced data without needing \pkg{lme4}.
#'
#' @inheritParams ge_model
#' @return A list (class \code{c("bam_ge_ammi", "bam_result", "list")}) with the components of the earlier versions:
#' \describe{
#'   \item{\code{g}, \code{e}, \code{Rep}, \code{k}}{numbers of genotypes, environments, replicates and axes.}
#'   \item{\code{ge_means}}{the result of \code{\link{ge_mean}()}.}
#'   \item{\code{mu}, \code{sigma2}, \code{tau}}{grand mean, error variance and its reciprocal (precision).}
#'   \item{\code{tao}, \code{delta}}{genotype and environment effects (\code{\link{g_eff}()}, \code{\link{e_eff}()}).}
#'   \item{\code{lambdas}, \code{alphas}, \code{gammas}}{singular values and genotype and environment singular
#'     vectors of the interaction.}
#' }
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#'
#' Gauch, H. G. (1992). \emph{Statistical Analysis of Regional Yield Trials: AMMI Analysis of Factorial Designs}.
#' Elsevier, Amsterdam.
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bayes_ammi}}, \code{\link{ge_eff}}
#' @examples
#' data(Maiz)
#' fit <-
#'   ge_ammi(
#'       .data = Maiz
#'     , .y    = y
#'     , .gen  = entry
#'     , .env  = site
#'     , .rep  = rep
#'     )
#' fit$lambdas
#' fit$sigma2
#' @export
ge_ammi <- function(.data, .y, .gen, .env, .rep) UseMethod("ge_ammi")

#' @rdname ge_ammi
#' @export
ge_ammi.default <- function(.data, .y, .gen, .env, .rep) {
  fn <- "ge_ammi"
  cd <- .bam_prep(fn, .data, substitute(.y), substitute(.gen), substitute(.env), substitute(.rep), parent.frame())
  if (!cd$balanced) .bam_stop(fn, "every genotype must have the same number of replicates in every environment.")
  Y  <- cd$Y
  G  <- nrow(Y)
  E  <- ncol(Y)
  B  <- cd$R
  if (B < 2L) .bam_stop(fn, "at least 2 replicates are needed for the error variance.")
  K  <- min(G, E) - 1L
  ls <- .bam_ls(Y, K)
  ssr  <- .bam_ss_rep(cd)
  df_e <- E * (B - 1) * (G - 1)
  mse  <- (cd$ssw - ssr) / df_e
  s2   <- if (ssr / (E * (B - 1)) >= mse) mse else cd$ssw / (G * E * (B - 1))
  .bam_result(list(
      g        = G
    , e        = E
    , Rep      = B
    , k        = K
    , ge_means = list(ge_means = Y, grand_mean = mean(Y))
    , mu       = mean(Y)
    , sigma2   = s2
    , tau      = 1 / s2
    , tao      = matrix(ls$tau, ncol = 1L, dimnames = list(cd$G, "GEffs"))
    , delta    = matrix(ls$delta, ncol = 1L, dimnames = list(cd$E, "EEffs"))
    , lambdas  = ls$d[seq_len(K)]
    , alphas   = ls$alpha
    , gammas   = ls$gamma
    ), "ge_ammi", label = sprintf("Least-squares AMMI: %d genotypes, %d environments, %d replicates", G, E, B),
    call = match.call(), meta = list(response = cd$names[["y"]]))
}
