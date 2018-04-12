#' @name    ge_var
#' @aliases ge_var
#' @title   Genotype by Environment Interaction Variances
#' @description Calcuates Genotype by Environment Interaction Variances
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#'
#' @return Genotype by Environment Interaction Variances
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
#'  (\href{https://dl.sciencesocieties.org/publications/cs/abstracts/51/4/1458?access=0&view=pdf}{doi: 10.2135/cropsci2010.06.0343})
#'
#' @import rlang
#' @import tidyr
#' @importFrom  dplyr select group_by count
#' @importFrom magrittr %>%
#' @importFrom stats var
#'
#' @export
#'
#' @examples
#'
#' data(Maiz)
#' ge_var(
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
      , "GEVar"
      , "."
    )
  )
}

ge_var <- function(.data, .y, .gen, .env) {
  UseMethod("ge_var")
}


#' @export
#' @rdname ge_var

ge_var.default <-
  function(.data, .y, .gen, .env){

    .y    <- deparse(substitute(.y))
    .gen  <- deparse(substitute(.gen))
    .env  <- deparse(substitute(.env))

    df1 <- tibble::as_tibble(data.frame(
      Env = factor(.data[[.env]])
      , Gen = factor(.data[[.gen]])
      , Y   = .data[[.y]]
    ))


    ge_variances <-
      df1 %>%
      dplyr::group_by(Gen, Env) %>%
      dplyr::summarise(GEVar = var(Y)) %>%
      tidyr::spread(key = Env, value = GEVar) %>%
     # magrittr::set_rownames(.$Gen) %>%
      dplyr::ungroup() %>%
      dplyr::select(- Gen) %>%
      as.matrix()

    return(list(
      ge_variances   = ge_variances
      ))
}
