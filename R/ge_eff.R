#' Genotype-by-environment interaction effects
#'
#' @description
#' \code{ge_eff()} computes the interaction effects (the cell means centred on genotype and environment means) and
#' their singular value decomposition, the least-squares estimate of the multiplicative part of the AMMI model.
#'
#' @details
#' \deqn{\mathbf{Z} = (z_{ij}), \quad z_{ij} = \bar{y}_{ij} - \bar{y}_{i\cdot} - \bar{y}_{\cdot j} + \bar{y}, \qquad
#'   \mathbf{Z} = \sum_k \lambda_k\boldsymbol{\alpha}_k\boldsymbol{\gamma}_k^\top,}{z_ij = y_ij - y_i. - y_.j + ybar;
#'   Z = sum_k lambda_k alpha_k gamma_k',}
#' with singular values \eqn{\lambda_1 \ge \lambda_2 \ge \cdots}{lambda_1 >= lambda_2 >= ...} and orthonormal
#' singular vectors orthogonal to \eqn{\mathbf{1}}{1}.
#'
#' @inheritParams ge_mean
#' @return A list (class \code{c("bam_ge_eff", "bam_result", "list")}) with \code{ge_effects} (the \eqn{G \times
#'   E}{G x E} matrix \eqn{\mathbf{Z}}{Z}) and \code{ge_svd} (\code{svd()} of it: \code{d}, \code{u}, \code{v}).
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#'
#' Gauch, H. G. (1992). \emph{Statistical Analysis of Regional Yield Trials: AMMI Analysis of Factorial Designs}.
#' Elsevier, Amsterdam.
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{ge_ammi}}
#' @examples
#' data(Maiz)
#' fit <-
#'   ge_eff(
#'       .data = Maiz
#'     , .y    = y
#'     , .gen  = entry
#'     , .env  = site
#'     )
#' fit$ge_svd$d
#' @export
ge_eff <- function(.data, .y, .gen, .env) UseMethod("ge_eff")

#' @rdname ge_eff
#' @export
ge_eff.default <- function(.data, .y, .gen, .env) {
  cd <- .bam_prep("ge_eff", .data, substitute(.y), substitute(.gen), substitute(.env), NULL, parent.frame(),
                  need_rep = FALSE)
  Y  <- cd$Y
  Z  <- Y - outer(rowMeans(Y), colMeans(Y), `+`) + mean(Y)
  .bam_result(list(ge_effects = Z, ge_svd = svd(Z)), "ge_eff", label = "Interaction effects", call = match.call(),
              meta = list(response = cd$names[["y"]]))
}
