#' @name    matrix_k
#' @aliases matrix_k
#' @title   k Matrix
#' @description Gives k matrix
#'
#' @param n  Number of columns
#'
#' @return Matrix
#'
#' @author
#' \enumerate{
#'     \item Muhammad Yaseen (\email{myaseen208@gmail.com})
#'    }
#'
#' @references
#'  Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J.M., Viele, K., Liu, G., and Cornelius, P.L. (2011)
#'  Bayesian Estimation of the Additive Main Effects and Multiplicative Interaction Model
#'  \emph{Crop Science}, \strong{51}, 1458–1469.
#'  (\href{https://acsess.onlinelibrary.wiley.com/doi/abs/10.2135/cropsci2010.06.0343}{doi: 10.2135/cropsci2010.06.0343})
#'
#' @export


matrix_k <- function(n){
  UseMethod("matrix_k")
}


#' @export
#' @rdname matrix_k


matrix_k.default <- function(n){
  m <- matrix(0, nrow = (n-1), ncol = n)
  long <- sqrt((n - (0:(n-2)) - 1) * (n - (0:(n-2))))
  diag(m) <- (n - (0:(n-2)) - 1)
  m[upper.tri(m)] <- -1
  m / long[row(m)]
}
