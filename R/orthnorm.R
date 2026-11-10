#' Gram-Schmidt orthonormalisation
#'
#' @description
#' \code{orthnorm()} orthogonalises the columns of a matrix by the Gram-Schmidt process, optionally completes them to
#' a basis of the whole space and scales them to unit length.
#'
#' @details
#' For columns \eqn{\mathbf{u}_1, \ldots, \mathbf{u}_m}{u_1, ..., u_m},
#' \deqn{\mathbf{v}_1 = \mathbf{u}_1, \qquad \mathbf{v}_i = \mathbf{u}_i - \sum_{h < i}\frac{\mathbf{u}_i^\top
#'   \mathbf{v}_h}{\mathbf{v}_h^\top\mathbf{v}_h}\mathbf{v}_h,}{v_1 = u_1; v_i = u_i - sum_{h < i} (u_i' v_h / v_h'
#'   v_h) v_h,}
#' each divided by its length when \code{norm = TRUE}. With \code{basis = TRUE} and fewer columns than rows, unit
#' vectors are added (those least explained by the given columns) to complete a basis.
#'
#' @param u a matrix (or vector) whose columns are orthogonalised.
#' @param basis complete the columns to a basis of the whole space (default \code{TRUE}).
#' @param norm scale the columns to unit length (default \code{TRUE}).
#' @return The matrix of orthogonal (orthonormal with \code{norm = TRUE}) columns; \code{NULL} for \code{u = NULL}.
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{matrix_k}}
#' @examples
#' V <-
#'   orthnorm(
#'       u     = cbind(1, 1:4)
#'     , basis = TRUE
#'     , norm  = TRUE
#'     )
#' round(x = crossprod(V), digits = 10)
#' @export
orthnorm <- function(u = NULL, basis = TRUE, norm = TRUE) UseMethod("orthnorm")

#' @rdname orthnorm
#' @export
orthnorm.default <- function(u = NULL, basis = TRUE, norm = TRUE) {
  if (is.null(u)) return(NULL)
  if (!is.matrix(u)) u <- as.matrix(u)
  p <- nrow(u)
  n <- ncol(u)
  if (n > p) {
    .bam_warn("orthnorm", "more columns than rows: only the first ", p, " are used.")
    u <- u[, seq_len(p), drop = FALSE]
    n <- p
  }
  if (basis && p > n) {
    base <- diag(p)
    coef_proj <- crossprod(u, base) / diag(crossprod(u))
    base2 <- base - u %*% matrix(coef_proj, nrow = n, ncol = p)
    norm_base2 <- diag(crossprod(base2))
    base <- base[, order(norm_base2) > n, drop = FALSE]
    u <- cbind(u, base)
    n <- p
  }
  v <- u
  if (n > 1L) {
    for (i in 2:n) {
      prev <- v[, seq_len(i - 1L), drop = FALSE]
      v[, i] <- u[, i] - prev %*% (drop(crossprod(prev, u[, i])) / colSums(prev^2))
    }
  }
  if (norm) v <- sweep(v, 2L, sqrt(colSums(v^2)), `/`)
  v
}
