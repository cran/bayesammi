#' @name    e_eff
#' @aliases e_eff
#' @title   Environment Effects
#' @description Calcuates Environment Effects
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#'
#' @return Environment Effects
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
#'
#' @export
#'
#' @examples
#'
#' data(Maiz)
#' e_eff(
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
      "Env"
      , "Y"
      , "EMean"
      , "."
      , "EEffs"
    )
  )
}


e_eff <- function(.data, .y, .gen, .env) {
  UseMethod("e_eff")
}


#' @export
#' @rdname e_eff

e_eff.default <-
  function(.data, .y, .gen, .env){


    .y    <- deparse(substitute(.y))
    .gen  <- deparse(substitute(.gen))
    .env  <- deparse(substitute(.env))

    df1 <- tibble::as_tibble(data.frame(
      Env = factor(.data[[.env]])
      , Gen = factor(.data[[.gen]])
      , Y   = .data[[.y]]
    ))


    e_effects <-
      df1 %>%
      dplyr::group_by(Env) %>%
      dplyr::summarise(EMean = mean(Y)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(EEffs = EMean - mean(EMean)) %>%
    #  magrittr::set_rownames(.$Env) %>%
      dplyr::ungroup() %>%
      dplyr::select(EEffs) %>%
      as.matrix()

    return(e_effects)
}
