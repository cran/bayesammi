#' @name    ge_eff
#' @aliases ge_eff
#' @title   Genotype by Environment Interaction Effects
#' @description Calcuates Genotype by Environment Interaction Effects
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#'
#' @return Genotype by Environment Interaction Effects
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
#' ge_eff(
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
      , "GMean"
      , "EMean"
      , "Mean"
      , "GEEffs"
      , "."
    )
  )
}


ge_eff <- function(.data, .y, .gen, .env) {
  UseMethod("ge_eff")
}


#' @export
#' @rdname ge_eff

ge_eff.default <-
  function(.data, .y, .gen, .env){

    .y    <- deparse(substitute(.y))
    .gen  <- deparse(substitute(.gen))
    .env  <- deparse(substitute(.env))

    df1 <- tibble::as_tibble(data.frame(
      Env = factor(.data[[.env]])
      , Gen = factor(.data[[.gen]])
      , Y   = .data[[.y]]
    ))


    ge_effects <-
      df1 %>%
      dplyr::group_by(Gen, Env) %>%
      dplyr::mutate(GEMean = mean(Y)) %>%
      dplyr::group_by(Gen) %>%
      dplyr::mutate(GMean = mean(Y)) %>%
      dplyr::group_by(Env) %>%
      dplyr::mutate(EMean = mean(Y)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Mean = mean(Y)) %>%
      dplyr::group_by(Gen, Env) %>%
      dplyr::summarize(GEEffs = mean(GEMean - GMean - EMean + Mean)) %>%
      tidyr::spread(key = Env, value = GEEffs) %>%
  #    magrittr::set_rownames(.$Gen) %>%
      dplyr::ungroup() %>%
      dplyr::select(- Gen) %>%
      as.matrix()

    ge_svd <- svd(ge_effects)

    return(list(
        ge_effects = ge_effects
      , ge_svd     = ge_svd
      ))
}
