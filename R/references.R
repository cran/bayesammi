#-------------------------------------------------------------------------------------------------------------------
# The sources of every method: one registry gives the short citation printed under each result, the references of
# the interpretations and the application, bam_cite(), inst/CITATION and inst/REFERENCES.bib (written from it by
# tools/write-citation.R).
#-------------------------------------------------------------------------------------------------------------------

.bam_ref <- function(key, short, text, ...) {
  list(key = key, short = short, text = text, bib = utils::bibentry(key = key, ...))
}

.bam_p <- function(given, family) utils::person(given = given, family = family)

.bam_refs_list <- function() {
  p <- .bam_p
  refs <- list(
    .bam_ref("Crossa2011", "Crossa et al. (2011)",
             paste("Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L.",
                   "(2011). Bayesian estimation of the additive main effects and multiplicative interaction model.",
                   "Crop Science, 51(4), 1458-1469. doi:10.2135/cropsci2010.06.0343"),
             bibtype = "Article",
             title = "Bayesian Estimation of the Additive Main Effects and Multiplicative Interaction Model",
             author = c(p("Jose", "Crossa"), p("Sergio", "Perez-Elizalde"), p("Diego", "Jarquin"),
                        p(c("Jose", "Miguel"), "Cotes"), p("Kert", "Viele"), p("Genzhou", "Liu"),
                        p(c("Paul", "L."), "Cornelius")),
             journal = "Crop Science", year = "2011", volume = "51", number = "4", pages = "1458--1469",
             doi = "10.2135/cropsci2010.06.0343"),
    .bam_ref("PerezElizalde2012", "Perez-Elizalde et al. (2012)",
             paste("Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A general Bayesian estimation method of",
                   "linear-bilinear models applied to plant breeding trials with genotype x environment interaction.",
                   "Journal of Agricultural, Biological, and Environmental Statistics, 17(1), 15-37.",
                   "doi:10.1007/s13253-011-0063-9"),
             bibtype = "Article",
             title = paste("A General {Bayesian} Estimation Method of Linear-Bilinear Models Applied to Plant Breeding",
                           "Trials with Genotype x Environment Interaction"),
             author = c(p("Sergio", "Perez-Elizalde"), p("Diego", "Jarquin"), p("Jose", "Crossa")),
             journal = "Journal of Agricultural, Biological, and Environmental Statistics", year = "2012",
             volume = "17", number = "1", pages = "15--37", doi = "10.1007/s13253-011-0063-9"),
    .bam_ref("Viele2000", "Viele and Srinivasan (2000)",
             paste("Viele, K. and Srinivasan, C. (2000). Parsimonious estimation of multiplicative interaction in",
                   "analysis of variance using Kullback-Leibler information. Journal of Statistical Planning and",
                   "Inference, 84(1-2), 201-219."),
             bibtype = "Article",
             title = paste("Parsimonious Estimation of Multiplicative Interaction in Analysis of Variance Using",
                           "{Kullback-Leibler} Information"),
             author = c(p("Kert", "Viele"), p("C.", "Srinivasan")),
             journal = "Journal of Statistical Planning and Inference", year = "2000", volume = "84",
             number = "1-2", pages = "201--219"),
    .bam_ref("Hoff2009", "Hoff (2009)",
             paste("Hoff, P. D. (2009). Simulation of the matrix Bingham-von Mises-Fisher distribution, with",
                   "applications to multivariate and relational data. Journal of Computational and Graphical",
                   "Statistics, 18(2), 438-456. doi:10.1198/jcgs.2009.07177"),
             bibtype = "Article",
             title = paste("Simulation of the Matrix {Bingham-von Mises-Fisher} Distribution, with Applications to",
                           "Multivariate and Relational Data"),
             author = p(c("Peter", "D."), "Hoff"),
             journal = "Journal of Computational and Graphical Statistics", year = "2009", volume = "18",
             number = "2", pages = "438--456", doi = "10.1198/jcgs.2009.07177"),
    .bam_ref("Yan2000", "Yan et al. (2000)",
             paste("Yan, W., Hunt, L. A., Sheng, Q. and Szlavnics, Z. (2000). Cultivar evaluation and",
                   "mega-environment investigation based on the GGE biplot. Crop Science, 40(3), 597-605.",
                   "doi:10.2135/cropsci2000.403597x"),
             bibtype = "Article",
             title = "Cultivar Evaluation and Mega-Environment Investigation Based on the {GGE} Biplot",
             author = c(p("W.", "Yan"), p(c("L.", "A."), "Hunt"), p("Q.", "Sheng"), p("Z.", "Szlavnics")),
             journal = "Crop Science", year = "2000", volume = "40", number = "3", pages = "597--605",
             doi = "10.2135/cropsci2000.403597x"),
    .bam_ref("Dias2022", "Dias et al. (2022)",
             paste("Dias, K. O. G., dos Santos, J. P. R., Krause, M. D., Piepho, H.-P., Guimaraes, L. J. M., Pastina,",
                   "M. M. and Garcia, A. A. F. (2022). Leveraging probability concepts for cultivar recommendation in",
                   "multi-environment trials. Theoretical and Applied Genetics.",
                   "doi:10.1007/s00122-022-04041-y"),
             bibtype = "Article",
             title = "Leveraging Probability Concepts for Cultivar Recommendation in Multi-Environment Trials",
             author = c(p(c("K.", "O.", "G."), "Dias"), p(c("J.", "P.", "R."), "dos Santos"), p(c("M.", "D."), "Krause"),
                        p("H.-P.", "Piepho"), p(c("L.", "J.", "M."), "Guimaraes"), p(c("M.", "M."), "Pastina"),
                        p(c("A.", "A.", "F."), "Garcia")),
             journal = "Theoretical and Applied Genetics", year = "2022",
             doi = "10.1007/s00122-022-04041-y"),
    .bam_ref("FinlayWilkinson1963", "Finlay and Wilkinson (1963)",
             paste("Finlay, K. W. and Wilkinson, G. N. (1963). The analysis of adaptation in a plant-breeding",
                   "programme. Australian Journal of Agricultural Research, 14(6), 742-754. doi:10.1071/AR9630742"),
             bibtype = "Article", title = "The Analysis of Adaptation in a Plant-Breeding Programme",
             author = c(p(c("K.", "W."), "Finlay"), p(c("G.", "N."), "Wilkinson")),
             journal = "Australian Journal of Agricultural Research", year = "1963", volume = "14", number = "6",
             pages = "742--754", doi = "10.1071/AR9630742"),
    .bam_ref("EberhartRussell1966", "Eberhart and Russell (1966)",
             paste("Eberhart, S. A. and Russell, W. A. (1966). Stability parameters for comparing varieties. Crop",
                   "Science, 6(1), 36-40. doi:10.2135/cropsci1966.0011183X000600010011x"),
             bibtype = "Article", title = "Stability Parameters for Comparing Varieties",
             author = c(p(c("S.", "A."), "Eberhart"), p(c("W.", "A."), "Russell")),
             journal = "Crop Science", year = "1966", volume = "6", number = "1", pages = "36--40",
             doi = "10.2135/cropsci1966.0011183X000600010011x"),
    .bam_ref("Shukla1972", "Shukla (1972)",
             paste("Shukla, G. K. (1972). Some statistical aspects of partitioning genotype-environmental components",
                   "of variability. Heredity, 29(2), 237-245. doi:10.1038/hdy.1972.87"),
             bibtype = "Article",
             title = "Some Statistical Aspects of Partitioning Genotype-Environmental Components of Variability",
             author = p(c("G.", "K."), "Shukla"), journal = "Heredity", year = "1972", volume = "29",
             number = "2", pages = "237--245", doi = "10.1038/hdy.1972.87"),
    .bam_ref("Wricke1962", "Wricke (1962)",
             paste("Wricke, G. (1962). Uber eine Methode zur Erfassung der okologischen Streubreite in Feldversuchen.",
                   "Zeitschrift fur Pflanzenzuchtung, 47, 92-96."),
             bibtype = "Article",
             title = "{\\\"U}ber eine Methode zur Erfassung der {\\\"o}kologischen Streubreite in Feldversuchen",
             author = p("G.", "Wricke"), journal = "Zeitschrift f{\\\"u}r Pflanzenz{\\\"u}chtung", year = "1962",
             volume = "47", pages = "92--96"),
    .bam_ref("LinBinns1988", "Lin and Binns (1988)",
             paste("Lin, C. S. and Binns, M. R. (1988). A superiority measure of cultivar performance for cultivar x",
                   "location data. Canadian Journal of Plant Science, 68(1), 193-198. doi:10.4141/cjps88-018"),
             bibtype = "Article",
             title = "A Superiority Measure of Cultivar Performance for Cultivar x Location Data",
             author = c(p(c("C.", "S."), "Lin"), p(c("M.", "R."), "Binns")),
             journal = "Canadian Journal of Plant Science", year = "1988", volume = "68", number = "1",
             pages = "193--198", doi = "10.4141/cjps88-018"),
    .bam_ref("Annicchiarico1992", "Annicchiarico (1992)",
             paste("Annicchiarico, P. (1992). Cultivar adaptation and recommendation from alfalfa trials in Northern",
                   "Italy. Journal of Genetics and Breeding, 46, 269-278."),
             bibtype = "Article",
             title = "Cultivar Adaptation and Recommendation from Alfalfa Trials in {Northern Italy}",
             author = p("P.", "Annicchiarico"), journal = "Journal of Genetics and Breeding", year = "1992",
             volume = "46", pages = "269--278"),
    .bam_ref("Eskridge1990", "Eskridge (1990)",
             paste("Eskridge, K. M. (1990). Selection of stable cultivars using a safety-first rule. Crop Science,",
                   "30, 369-374."),
             bibtype = "Article", title = "Selection of Stable Cultivars Using a Safety-First Rule",
             author = p(c("K.", "M."), "Eskridge"), journal = "Crop Science", year = "1990", volume = "30",
             pages = "369--374"),
    .bam_ref("EdwardsJannink2006", "Edwards and Jannink (2006)",
             paste("Edwards, J. W. and Jannink, J.-L. (2006). Bayesian modeling of heterogeneous error and genotype x",
                   "environment interaction variances. Crop Science, 46(2), 820. doi:10.2135/cropsci2005.0164"),
             bibtype = "Article",
             title = "{Bayesian} Modeling of Heterogeneous Error and Genotype x Environment Interaction Variances",
             author = c(p(c("J.", "W."), "Edwards"), p("J.-L.", "Jannink")), journal = "Crop Science",
             year = "2006", volume = "46", number = "2", pages = "820", doi = "10.2135/cropsci2005.0164"),
    .bam_ref("LianCampos2016", "Lian and de los Campos (2016)",
             paste("Lian, L. and de los Campos, G. (2016). FW: an R package for Finlay-Wilkinson regression that",
                   "incorporates genomic/pedigree information and covariance structures between environments. G3:",
                   "Genes, Genomes, Genetics, 6(3), 589-597. doi:10.1534/g3.115.026328"),
             bibtype = "Article",
             title = paste("{FW}: An {R} Package for {Finlay-Wilkinson} Regression that Incorporates",
                           "Genomic/Pedigree Information and Covariance Structures Between Environments"),
             author = c(p("L.", "Lian"), p("G.", "de los Campos")),
             journal = "G3: Genes, Genomes, Genetics", year = "2016", volume = "6", number = "3",
             pages = "589--597", doi = "10.1534/g3.115.026328"),
    .bam_ref("Oliveira2018", "Oliveira et al. (2018)",
             paste("Oliveira, T. R. A., Carvalho, H. W. L., Nascimento, M., Costa, E. F. N., Amaral Junior, A. T.,",
                   "Gravina, G. A. and Carvalho Filho, J. L. S. (2018). The Eberhart and Russell's Bayesian method used",
                   "as an instrument to select maize hybrids. Euphytica, 214, 64. doi:10.1007/s10681-018-2146-y"),
             bibtype = "Article",
             title = "The {Eberhart} and {Russell}'s {Bayesian} Method Used as an Instrument to Select Maize Hybrids",
             author = c(p(c("T.", "R.", "A."), "de Oliveira"), p(c("H.", "W.", "L."), "de Carvalho"),
                        p("M.", "Nascimento"), p(c("E.", "F.", "N."), "Costa"), p(c("A.", "T."), "do Amaral Junior"),
                        p(c("G.", "A."), "Gravina"), p(c("J.", "L.", "S."), "de Carvalho Filho")),
             journal = "Euphytica", year = "2018", volume = "214", pages = "64", doi = "10.1007/s10681-018-2146-y"),
    .bam_ref("Purchase2000", "Purchase et al. (2000)",
             paste("Purchase, J. L., Hatting, H. and van Deventer, C. S. (2000). Genotype x environment interaction",
                   "of winter wheat (Triticum aestivum L.) in South Africa: II. Stability analysis of yield",
                   "performance. South African Journal of Plant and Soil, 17(3), 101-107."),
             bibtype = "Article",
             title = paste("Genotype x Environment Interaction of Winter Wheat ({Triticum aestivum} {L}.) in South",
                           "Africa: {II}. Stability Analysis of Yield Performance"),
             author = c(p(c("J.", "L."), "Purchase"), p("H.", "Hatting"), p(c("C.", "S."), "van Deventer")),
             journal = "South African Journal of Plant and Soil", year = "2000", volume = "17", number = "3",
             pages = "101--107"),
    .bam_ref("Sneller1997", "Sneller et al. (1997)",
             paste("Sneller, C. H., Kilgore-Norquest, L. and Dombek, D. (1997). Repeatability of yield stability",
                   "statistics in soybean. Crop Science, 37(2), 383-390."),
             bibtype = "Article", title = "Repeatability of Yield Stability Statistics in Soybean",
             author = c(p(c("C.", "H."), "Sneller"), p("L.", "Kilgore-Norquest"), p("D.", "Dombek")),
             journal = "Crop Science", year = "1997", volume = "37", number = "2", pages = "383--390"),
    .bam_ref("Zobel1994", "Zobel (1994)",
             paste("Zobel, R. W. (1994). Stress resistance and root systems. In Proceedings of the Workshop on",
                   "Adaptation of Plants to Soil Stress, INTSORMIL Publication 94-2, 80-99. Institute of Agriculture",
                   "and Natural Resources, University of Nebraska-Lincoln."),
             bibtype = "InProceedings", title = "Stress Resistance and Root Systems",
             author = p(c("R.", "W."), "Zobel"),
             booktitle = "Proceedings of the Workshop on Adaptation of Plants to Soil Stress",
             series = "INTSORMIL Publication 94-2", year = "1994", pages = "80--99",
             publisher = "Institute of Agriculture and Natural Resources, University of Nebraska-Lincoln"),
    .bam_ref("Zali2012", "Zali et al. (2012)",
             paste("Zali, H., Farshadfar, E., Sabaghpour, S. H. and Karimizadeh, R. (2012). Evaluation of genotype x",
                   "environment interaction in chickpea using measures of stability from AMMI model. Annals of",
                   "Biological Research, 3(7), 3126-3136."),
             bibtype = "Article",
             title = paste("Evaluation of Genotype x Environment Interaction in Chickpea Using Measures of Stability",
                           "from {AMMI} Model"),
             author = c(p("H.", "Zali"), p("E.", "Farshadfar"), p(c("S.", "H."), "Sabaghpour"), p("R.", "Karimizadeh")),
             journal = "Annals of Biological Research", year = "2012", volume = "3", number = "7",
             pages = "3126--3136"),
    .bam_ref("Olivoto2019", "Olivoto et al. (2019)",
             paste("Olivoto, T., Lucio, A. D. C., da Silva, J. A. G., Marchioro, V. S., de Souza, V. Q. and Jost, E.",
                   "(2019). Mean performance and stability in multi-environment trials I: combining features of AMMI",
                   "and BLUP techniques. Agronomy Journal, 111(6), 2949-2960. doi:10.2134/agronj2019.03.0220"),
             bibtype = "Article",
             title = paste("Mean Performance and Stability in Multi-Environment Trials {I}: Combining Features of",
                           "{AMMI} and {BLUP} Techniques"),
             author = c(p("T.", "Olivoto"), p(c("A.", "D.", "C."), "Lucio"), p(c("J.", "A.", "G."), "da Silva"),
                        p(c("V.", "S."), "Marchioro"), p(c("V.", "Q."), "de Souza"), p("E.", "Jost")),
             journal = "Agronomy Journal", year = "2019", volume = "111", number = "6", pages = "2949--2960",
             doi = "10.2134/agronj2019.03.0220"),
    .bam_ref("Farshadfar2008", "Farshadfar (2008)",
             paste("Farshadfar, E. (2008). Incorporation of AMMI stability value and grain yield in a single",
                   "non-parametric index (GSI) in bread wheat. Pakistan Journal of Biological Sciences, 11(14),",
                   "1791-1796."),
             bibtype = "Article",
             title = paste("Incorporation of {AMMI} Stability Value and Grain Yield in a Single Non-Parametric Index",
                           "({GSI}) in Bread Wheat"),
             author = p("E.", "Farshadfar"), journal = "Pakistan Journal of Biological Sciences",
             year = "2008", volume = "11", number = "14", pages = "1791--1796"),
    .bam_ref("Spiegelhalter2002", "Spiegelhalter et al. (2002)",
             paste("Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A. (2002). Bayesian measures",
                   "of model complexity and fit. Journal of the Royal Statistical Society B, 64(4), 583-639.",
                   "doi:10.1111/1467-9868.00353"),
             bibtype = "Article", title = "{Bayesian} Measures of Model Complexity and Fit",
             author = c(p(c("David", "J."), "Spiegelhalter"), p(c("Nicola", "G."), "Best"),
                        p(c("Bradley", "P."), "Carlin"), p("Angelika", "van der Linde")),
             journal = "Journal of the Royal Statistical Society B", year = "2002", volume = "64", number = "4",
             pages = "583--639", doi = "10.1111/1467-9868.00353"),
    .bam_ref("Watanabe2010", "Watanabe (2010)",
             paste("Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable",
                   "information criterion in singular learning theory. Journal of Machine Learning Research, 11,",
                   "3571-3594."),
             bibtype = "Article",
             title = paste("Asymptotic Equivalence of {Bayes} Cross Validation and Widely Applicable Information",
                           "Criterion in Singular Learning Theory"),
             author = p("Sumio", "Watanabe"), journal = "Journal of Machine Learning Research", year = "2010",
             volume = "11", pages = "3571--3594"),
    .bam_ref("Vehtari2017", "Vehtari et al. (2017)",
             paste("Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model evaluation using",
                   "leave-one-out cross-validation and WAIC. Statistics and Computing, 27(5), 1413-1432.",
                   "doi:10.1007/s11222-016-9696-4"),
             bibtype = "Article",
             title = "Practical {Bayesian} Model Evaluation Using Leave-One-Out Cross-Validation and {WAIC}",
             author = c(p("Aki", "Vehtari"), p("Andrew", "Gelman"), p("Jonah", "Gabry")),
             journal = "Statistics and Computing", year = "2017", volume = "27", number = "5",
             pages = "1413--1432", doi = "10.1007/s11222-016-9696-4"),
    .bam_ref("Gelman1996", "Gelman et al. (1996)",
             paste("Gelman, A., Meng, X.-L. and Stern, H. (1996). Posterior predictive assessment of model fitness",
                   "via realized discrepancies. Statistica Sinica, 6(4), 733-807."),
             bibtype = "Article", title = "Posterior Predictive Assessment of Model Fitness via Realized Discrepancies",
             author = c(p("Andrew", "Gelman"), p("Xiao-Li", "Meng"), p("Hal", "Stern")),
             journal = "Statistica Sinica", year = "1996", volume = "6", number = "4", pages = "733--807"),
    .bam_ref("Gabry2019", "Gabry et al. (2019)",
             paste("Gabry, J., Simpson, D., Vehtari, A., Betancourt, M. and Gelman, A. (2019). Visualization in",
                   "Bayesian workflow. Journal of the Royal Statistical Society A, 182(2), 389-402.",
                   "doi:10.1111/rssa.12378"),
             bibtype = "Article", title = "Visualization in {Bayesian} Workflow",
             author = c(p("Jonah", "Gabry"), p("Daniel", "Simpson"), p("Aki", "Vehtari"), p("Michael", "Betancourt"),
                        p("Andrew", "Gelman")),
             journal = "Journal of the Royal Statistical Society A", year = "2019", volume = "182", number = "2",
             pages = "389--402", doi = "10.1111/rssa.12378")
    )
  stats::setNames(refs, vapply(refs, `[[`, "", "key"))
}

# the sources of a result, by its class (and model)
.bam_sources_of <- function(x) {
  cls <- .bam_cls(x)
  s <- x$settings
  switch(cls,
    bayes_ammi = c("Crossa2011", if (identical(s$model, "GGE")) c("PerezElizalde2012", "Yan2000"),
                   if (isTRUE(s$het)) c("PerezElizalde2012", "Hoff2009"), "Spiegelhalter2002", "Watanabe2010"),
    bayes_stability = c(switch(s$model,
                               ER = c("EberhartRussell1966", "Oliveira2018", "LianCampos2016"),
                               FW = c("FinlayWilkinson1963", "LianCampos2016"),
                               Shukla = c("Shukla1972", "EdwardsJannink2006")),
                        "Eskridge1990", "Spiegelhalter2002", "Watanabe2010"),
    biplots = "Crossa2011",
    ge_ammi = "Crossa2011",
    predict = "Crossa2011",
    winners = "Yan2000",
    superiority = "Dias2022",
    ranks = "Dias2022",
    stability = c("Purchase2000", "Shukla1972", "Farshadfar2008"),
    indices = c("Wricke1962", "Shukla1972", "LinBinns1988", "Annicchiarico1992", "Eskridge1990"),
    ammi_indices = c("Purchase2000", "Sneller1997", "Zobel1994", "Zali2012", "Olivoto2019"),
    compare = c("Spiegelhalter2002", "Watanabe2010", "Vehtari2017"),
    ppc = "Gelman1996",
    character(0))
}

# "Sources: Crossa et al. (2011); ..." for the end of a print
.bam_sources_line <- function(x) {
  k <- unique(.bam_sources_of(x))
  if (!length(k)) return(invisible(NULL))
  r <- .bam_refs_list()[k]
  txt <- strwrap(paste0("Sources: ", paste(vapply(r, `[[`, "", "short"), collapse = "; "),
                        ". Full references: bam_cite(x)."), width = getOption("width", 80L))
  cat("\n", paste(txt, collapse = "\n"), "\n", sep = "")
  invisible(NULL)
}

#' References of a result or of the package
#'
#' @description
#' \code{bam_cite()} gives the references of the methods behind a result of the package, or the references of
#' every method in the package, as text or BibTeX. \code{citation("bayesammi")} cites the package and its main
#' methods.
#'
#' @param x a result of a function of the package, or \code{NULL} (default) for every reference of the package.
#' @param style \code{"text"} (default) or \code{"bibtex"}.
#' @return An object of class \code{bibentry} (printed in the chosen style).
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @examples
#' bam_cite(
#'     x     = NULL
#'   , style = "text"
#'   )
#' data(Maiz)
#' g <-
#'   ge_ammi(
#'       .data = Maiz
#'     , .y    = y
#'     , .gen  = entry
#'     , .env  = site
#'     , .rep  = rep
#'     )
#' bam_cite(
#'     x     = g
#'   , style = "bibtex"
#'   )
#' @export
bam_cite <- function(x = NULL, style = c("text", "bibtex")) {
  style <- match.arg(style)
  refs <- .bam_refs_list()
  keys <- if (is.null(x)) names(refs) else {
    if (!inherits(x, "bam_result")) .bam_stop("bam_cite", "`x` must be a result of a function of the package or NULL.")
    unique(.bam_sources_of(x))
  }
  out <- do.call(c, lapply(refs[keys], `[[`, "bib"))
  if (is.null(out)) out <- utils::bibentry("Misc", title = "none", author = "none", year = "none")[0L]
  structure(out, bam_style = style, class = c("bam_cite", class(out)))
}

#' @export
print.bam_cite <- function(x, ...) {
  style <- attr(x, "bam_style")
  y <- x
  class(y) <- setdiff(class(y), "bam_cite")
  if (!length(y)) {
    cat("No references.\n")
    return(invisible(x))
  }
  if (identical(style, "bibtex")) print(y, style = "Bibtex") else {
    refs <- .bam_refs_list()
    keys <- vapply(unclass(y), function(e) attr(e, "key"), "")
    cat(paste(vapply(refs[keys], `[[`, "", "text"), collapse = "\n\n"), "\n")
  }
  invisible(x)
}
