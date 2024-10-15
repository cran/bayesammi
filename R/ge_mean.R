#' @name    ge_mean
#' @aliases ge_mean
#' @title   Genotype by Environment Interaction Means
#' @description Calcuates Genotype by Environment Interaction Means
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#'
#' @return Genotype by Environment Interaction Means
#'
#' @author
#' \enumerate{
#'     \item Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'    }
#'
#' @references
#'  Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J.M., Viele, K., Liu, G., and Cornelius, P.L. (2011)
#'  Bayesian Estimation of the Additive Main Effects and Multiplicative Interaction Model
#'  \emph{Crop Science}, \strong{51}, 1458–1469.
#'  (\href{https://acsess.onlinelibrary.wiley.com/doi/abs/10.2135/cropsci2010.06.0343}{doi: 10.2135/cropsci2010.06.0343})
#'
#' @import rlang
#' @import tidyr
#' @importFrom  dplyr select group_by count
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#'
#' data(Maiz)
#' ge_mean(
#'     .data  = Maiz
#'    , .y    = y
#'    , .gen  = entry
#'    , .env  = site
#'    )
#'
#'


if(getRversion() >= "2.15.1"){
  utils::globalVariables(
    c(
      "Gen"
      , "Env"
      , "Y"
      , "GEMean"
      , "."
    )
  )
}

ge_mean <- function(.data, .y, .gen, .env) {
  UseMethod("ge_mean")
}


#' @export
#' @rdname ge_mean

ge_mean.default <-
  function(.data, .y, .gen, .env){

    .y    <- deparse(substitute(.y))
    .gen  <- deparse(substitute(.gen))
    .env  <- deparse(substitute(.env))

    df1 <- tibble::as_tibble(data.frame(
        Env = factor(.data[[.env]])
      , Gen = factor(.data[[.gen]])
      , Y   = .data[[.y]]
    ))

    ge_means <-
      df1 %>%
      dplyr::group_by(Gen, Env) %>%
      dplyr::summarise(GEMean = mean(Y)) %>%
      tidyr::spread(key = Env, value = GEMean) %>%
    #  magrittr::set_rownames(.$Gen) %>%
      dplyr::ungroup() %>%
      dplyr::select(- Gen) %>%
      as.matrix()

    grand_mean <- mean(ge_means)

    return(list(
        ge_means   = ge_means
      , grand_mean = grand_mean
      ))
}
