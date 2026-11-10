#' Genotype effects
#'
#' @description
#' \code{g_eff()} computes the genotype main effects, the deviations of the genotype means from the grand mean.
#'
#' @details
#' \deqn{\hat{\tau}_i = \bar{y}_{i\cdot} - \bar{y}, \qquad \bar{y}_{i\cdot} = \frac{1}{E}\sum_j \bar{y}_{ij},}{tau_i
#'   = y_i. - ybar, y_i. = sum_j y_ij / E,}
#' the least-squares estimates of the AMMI model (equal to the deviations of the plot means when the data are
#' balanced).
#'
#' @inheritParams ge_mean
#' @return A \eqn{G \times 1}{G x 1} matrix with column \code{GEffs} (row names: genotypes).
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{e_eff}}, \code{\link{ge_eff}}
#' @examples
#' data(Maiz)
#' g_eff(
#'     .data = Maiz
#'   , .y    = y
#'   , .gen  = entry
#'   , .env  = site
#'   )
#' @export
g_eff <- function(.data, .y, .gen, .env) UseMethod("g_eff")

#' @rdname g_eff
#' @export
g_eff.default <- function(.data, .y, .gen, .env) {
  cd <- .bam_prep("g_eff", .data, substitute(.y), substitute(.gen), substitute(.env), NULL, parent.frame(),
                  need_rep = FALSE)
  matrix(rowMeans(cd$Y) - mean(cd$Y), ncol = 1L, dimnames = list(cd$G, "GEffs"))
}
