#' Environment effects
#'
#' @description
#' \code{e_eff()} computes the environment main effects, the deviations of the environment means from the grand mean.
#'
#' @details
#' \deqn{\hat{\delta}_j = \bar{y}_{\cdot j} - \bar{y}, \qquad \bar{y}_{\cdot j} = \frac{1}{G}\sum_i
#'   \bar{y}_{ij}.}{delta_j = y_.j - ybar, y_.j = sum_i y_ij / G.}
#'
#' @inheritParams ge_mean
#' @return An \eqn{E \times 1}{E x 1} matrix with column \code{EEffs} (row names: environments).
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{g_eff}}, \code{\link{ge_eff}}
#' @examples
#' data(Maiz)
#' e_eff(
#'     .data = Maiz
#'   , .y    = y
#'   , .gen  = entry
#'   , .env  = site
#'   )
#' @export
e_eff <- function(.data, .y, .gen, .env) UseMethod("e_eff")

#' @rdname e_eff
#' @export
e_eff.default <- function(.data, .y, .gen, .env) {
  cd <- .bam_prep("e_eff", .data, substitute(.y), substitute(.gen), substitute(.env), NULL, parent.frame(),
                  need_rep = FALSE)
  matrix(colMeans(cd$Y) - mean(cd$Y), ncol = 1L, dimnames = list(cd$E, "EEffs"))
}
