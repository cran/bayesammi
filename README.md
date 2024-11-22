
## `bayesammi`: Bayesian Estimation of the Additive Main Effects and Multiplicative Interaction Model

###### Version : [0.3.0](https://myaseen208.com/bayesammi/); Copyright (C) 2018-2024: License: [GPL-2\|GPL-3](https://www.r-project.org/Licenses/)

##### *Muhammad Yaseen<sup>1,2</sup>, Jose Crossa<sup>3,4,5,6</sup>, Sergio Perez-Elizalde<sup>7</sup>, Diego Jarquin<sup>8</sup>, Jose Miguel Cotes<sup>9</sup>, Kert Viele<sup>10</sup>, Genzhou Liu<sup>11</sup>, Paul L. Cornelius<sup>12</sup>, Julian Garcia Abadillo Velasco<sup>12</sup>*

1.  School of Mathematical & Statistical Sciences, Clemson University,
    Clemson, South Carolina, USA
2.  Department of Mathematics & Statistics, University of Agriculture
    Faisalabad, Pakistan
3.  Department of Statistics and Operations Research, and Distinguish
    Scientist Fellowship Program, King Saud University, Riyah, Saudi
    Arabia
4.  AgCenter, Louisiana State University, Baton Rouge, Louisiana, USA
5.  Colegio de Postgraduados, Montecillos, Mexico
6.  International Maize and Wheat Improvement Center (CIMMYT),
    Mexico-Veracruz, Mexico
7.  Colegio de Postgraduados, Montecillo,Estado de México 56230, México
8.  Agronomy Department, University of Florida, Gainesville, FL, United
    States
9.  Dep. de Ciencias Agrono ́ micas, Facultad de CienciasAgropecuarias,
    Univ. Nacional de Colombia, Calle 59A No 63–20 B11101-07, Medellı́n,
    Colombia
10. Department of Statistics, University of Kentucky, Lexington, KY,
    40546-03121, USA
11. Auxilium Pharmaceuticals, Inc., PA, USA
12. Department of Plant and Soil Sciences and Department of Statistics,
    University of Kentucky, Lexington, KY, 40546-03121, USA

------------------------------------------------------------------------

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version-last-release/bayesammi)](https://cran.r-project.org/package=bayesammi)
[![rstudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/bayesammi?color=green)](https://CRAN.R-project.org/package=bayesammi)
<!-- [![packageversion](https://img.shields.io/badge/Package%20version-0.2.3.3-orange.svg)](https://github.com/myaseen208/bayesammi) -->

[![develVersion](https://img.shields.io/badge/devel%20version-0.3.0-orange.svg)](https://github.com/myaseen208/bayesammi)

<!-- [![GitHub Download Count](https://github-basic-badges.herokuapp.com/downloads/myaseen208/bayesammi/total.svg)] -->

[![Project Status:
WIP](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Last-changedate](https://img.shields.io/badge/last%20change-2024--11--22-yellowgreen.svg)](https://github.com/myaseen208/bayesammi)
\*\*\*

## Description

Performs Bayesian estimation of the additive main effects and
multiplicative interaction (AMMI) model. The method is explained in
Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J.M., Viele, K.,
Liu, G. and Cornelius, P.L. (2011)
([doi:10.2135/cropsci2010.06.0343](https://doi.org/10.2135/cropsci2010.06.0343)).

## Installation

The package can be installed from CRAN as follows:

``` r
install.packages("bayesammi", dependencies = TRUE)
```

The development version can be installed from github as follows:

``` r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("myaseen208/bayesammi")
```

## What’s new

To know whats new in this version type:

``` r
news(package = "bayesammi")
```

## Links

[CRAN page](https://cran.r-project.org/package=bayesammi)

[Github page](https://github.com/myaseen208/bayesammi)

[Documentation website](https://myaseen208.com/bayesammi/)

## Citing `bayesammi`

To cite the methods in the package use:

``` r
citation("bayesammi")
Please, support this project by citing it in your publications!

  Yaseen M, Crossa J, Perez-Elizalde S, Jarquin D, Cotes JM, Viele K,
  Liu G, Cornelius PL (2018). _bayesammi: Bayesian Estimation of the
  Additive Main Effects and Multiplicative Interaction Model_.

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {bayesammi: Bayesian Estimation of the Additive Main Effects and Multiplicative Interaction Model},
    author = {Muhammad Yaseen and Jose Crossa and Sergio Perez-Elizalde and Diego Jarquin and Jose Miguel Cotes and Kert Viele and Genzhou Liu and Paul L. Cornelius},
    year = {2018},
    journal = {The Comprehensive R Archive Network (CRAN)},
  }
```
