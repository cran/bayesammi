#' @name    g_eff
#' @aliases g_eff
#' @title   Genotype Effects
#' @description Calcuates Genotype Effects
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#'
#' @return Genotype Effects
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
#' g_eff(
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
      , "Y"
      , "GMean"
      , "."
      , "GEffs"
    )
  )
}
g_eff <- function(.data, .y, .gen, .env) {
  UseMethod("g_eff")
}


#' @export
#' @rdname g_eff

g_eff.default <-
  function(.data, .y, .gen, .env){

    .y    <- deparse(substitute(.y))
    .gen  <- deparse(substitute(.gen))
    .env  <- deparse(substitute(.env))

    df1 <- tibble::as_tibble(data.frame(
      Env = factor(.data[[.env]])
      , Gen = factor(.data[[.gen]])
      , Y   = .data[[.y]]
    ))

    g_effects <-
      df1 %>%
        dplyr::group_by(Gen) %>%
        dplyr::summarise(GMean = mean(Y)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(GEffs = GMean - mean(GMean)) %>%
      #  magrittr::set_rownames(.$Gen) %>%
        dplyr::ungroup() %>%
        dplyr::select(GEffs) %>%
        as.matrix()

    return(g_effects)
}

