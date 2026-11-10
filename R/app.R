#' Bayesian AMMI and GGE analysis in a web browser
#'
#' @description
#' \code{bam_app()} starts the application of \pkg{bayesammi}: load a trial (the example \code{\link{Maiz}} or a CSV
#' file with one row per plot), fit the Bayesian AMMI or GGE model, check convergence, draw the biplot with HPD
#' regions, and turn the posterior into decisions (which genotype wins where, probabilities of superior performance
#' and stability, ranks, Bayesian stability, predictions of missing cells), and compare models. Every page shows an
#' interpretation computed from the results (\code{\link{bam_interpret}()}), the output as R prints it, formatted
#' tables, plots and the R code that repeats the analysis, each with a Copy button.
#'
#' @details
#' The application needs the suggested packages \pkg{shiny}, \pkg{bslib} and \pkg{DT}. Plots are sent to the browser
#' as WebP images when \pkg{magick} is installed (smaller than PNG; PNG otherwise) and can be downloaded as PNG, PDF
#' or WebP. Text is set in IBM Plex Sans and output in IBM Plex Mono when the browser can reach Google Fonts, in the
#' system fonts otherwise.
#'
#' @param launch.browser passed to \code{shiny::runApp()}: \code{TRUE} (default in an interactive session) opens the
#'   application in the web browser.
#' @param ... further arguments to \code{shiny::runApp()}, such as \code{port}.
#' @return \code{bam_app()} is called for its side effect and returns the value of \code{shiny::runApp()} when the
#'   application is closed; \code{bam_app_object()} returns the application object (for \code{shiny::runApp()},
#'   \code{shiny::testServer()} or deployment on a server).
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bayes_ammi}}, \code{\link{bam_interpret}}
#' @examples
#' if (interactive()) {
#'   bam_app(
#'       launch.browser = TRUE
#'     )
#' }
#' @export
bam_app <- function(launch.browser = interactive(), ...) {
  .bam_app_need()
  shiny::runApp(bam_app_object(), launch.browser = launch.browser, ...)
}

#' @rdname bam_app
#' @export
bam_app_object <- function() {
  .bam_app_need()
  shiny::shinyApp(ui = .bam_app_ui(), server = .bam_app_server)
}

.bam_app_need <- function() {
  miss <- c("shiny", "bslib", "DT")[!vapply(c("shiny", "bslib", "DT"), requireNamespace, NA, quietly = TRUE)]
  if (length(miss)) .bam_stop("bam_app", "the application needs the package(s) ", paste(miss, collapse = ", "),
                              ": install.packages(c(", paste0("\"", miss, "\"", collapse = ", "), ")).")
}

#-------------------------------------------------------------------------------------------------------------------
# look: theme, fonts, styles, Copy buttons
#-------------------------------------------------------------------------------------------------------------------
.bam_app_font <- function(family, kind = c("sans", "mono")) {
  fb <- switch(match.arg(kind),
               sans = c("system-ui", "-apple-system", "Segoe UI", "Helvetica Neue", "Arial", "sans-serif"),
               mono = c("SFMono-Regular", "Menlo", "Consolas", "Liberation Mono", "DejaVu Sans Mono", "monospace"))
  do.call(bslib::font_collection, c(list(bslib::font_google(family, local = FALSE)), as.list(fb)))
}

.bam_app_theme <- function() {
  bslib::bs_theme(version = 5, bootswatch = "flatly", primary = "#1f4e9e",
                  base_font = .bam_app_font("IBM Plex Sans", "sans"),
                  heading_font = .bam_app_font("IBM Plex Sans", "sans"),
                  code_font = .bam_app_font("IBM Plex Mono", "mono"))
}

.bam_app_head <- function() {
  shiny::tags$head(
      shiny::tags$style(shiny::HTML(paste(
          "pre.shiny-text-output{font-family:var(--bs-font-monospace);font-size:.85em;line-height:1.5;white-space:pre;"
        , "  overflow-x:auto;background:#f6f8fa;border:1px solid #dde3ea;border-radius:6px;padding:10px 12px;}"
        , "table.dataTable thead th{background:#1f4e9e;color:#fff;font-weight:600;}"
        , ".bam-bar{display:flex;gap:.4rem;justify-content:flex-end;margin-bottom:.4rem}"
        , ".dataTables_wrapper .dt-buttons{float:left;margin:0 1rem .5rem 0}"
        , ".dataTables_wrapper .dataTables_length{float:left;margin-bottom:.5rem}"
        , ".dataTables_wrapper .dataTables_filter{float:right;margin-bottom:.5rem}"
        , ".dataTables_wrapper .dataTables_length label,.dataTables_wrapper .dataTables_filter label{display:inline-flex;"
        , "  align-items:center;gap:.4rem;font-weight:400;margin:0}"
        , ".dataTables_wrapper select,.dataTables_wrapper input[type=search]{display:inline-block;width:auto}"
        , ".dataTables_wrapper .dataTables_info{float:left;padding-top:.6rem;font-size:.9em}"
        , ".dataTables_wrapper .dataTables_paginate{float:right;padding-top:.4rem}"
        , ".dataTables_wrapper table.dataTable{clear:both}"
        , ".dataTables_wrapper::after{content:'';display:table;clear:both}"
        , ".bam-plot img{max-width:100%;height:auto}")))
    , shiny::tags$script(shiny::HTML(paste(
          "$(document).on('click', '.bam-copy', function() {"
        , "  var el = document.getElementById($(this).data('target')); if (!el) return;"
        , "  var txt = el.innerText || el.textContent, b = $(this);"
        , "  var done = function() { b.find('.bam-copy-label').text('Copied'); setTimeout(function() { b.find('.bam-copy-label').text('Copy'); }, 1500); };"
        , "  if (navigator.clipboard && window.isSecureContext) { navigator.clipboard.writeText(txt).then(done); return; }"
        , "  var ta = document.createElement('textarea'); ta.value = txt; document.body.appendChild(ta); ta.select();"
        , "  document.execCommand('copy'); document.body.removeChild(ta); done();"
        , "});", sep = "\n")))
    )
}

.bam_copy_bar <- function(target) {
  shiny::div(class = "bam-bar",
             shiny::tags$button(type = "button", class = "btn btn-outline-secondary btn-sm bam-copy",
                                `data-target` = target, shiny::icon("copy"),
                                shiny::span(class = "bam-copy-label", "Copy")))
}

.bam_btn <- "btn btn-sm btn-outline-secondary me-1"

# a table for the application: Copy, CSV and Excel buttons, all rows on request, numbers rounded
.bam_dt <- function(d, title = "table") {
  d <- as.data.frame(d, check.names = FALSE)
  num <- which(vapply(d, is.numeric, NA))
  t <- DT::datatable(d, rownames = FALSE, extensions = "Buttons", style = "auto",
                     class = "compact stripe hover cell-border",
                     options = list(pageLength = 15L, scrollX = TRUE, dom = "Blfrtip",
                                    lengthMenu = list(c(15L, 50L, -1L), c("15", "50", "All")),
                                    buttons = list(list(extend = "copy", title = NULL, className = .bam_btn),
                                                   list(extend = "csv", title = title, className = .bam_btn),
                                                   list(extend = "excel", title = title, className = .bam_btn))))
  pcol  <- num[grepl("^p_|^prob|p_value|^Probability", names(d)[num])]
  whole <- setdiff(num[vapply(d[num], function(v) all(is.na(v) | v == round(v)), NA)], pcol)
  dec   <- setdiff(num, c(whole, pcol))
  if (length(pcol)) t <- DT::formatRound(t, pcol, digits = 3L)
  if (length(dec)) t <- DT::formatSignif(t, dec, digits = 5L)
  if (length(whole)) t <- DT::formatRound(t, whole, digits = 0L)
  t
}

# a ggplot as an image for renderImage(): WebP through magick when installed, PNG otherwise
.bam_plot_file <- function(p, width = 900, height = 600, res = 110) {
  png <- tempfile(fileext = ".png")
  if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png(png, width = width, height = height, res = res) else
    grDevices::png(png, width = width, height = height, res = res)
  print(p)
  grDevices::dev.off()
  if (.bam_webp_ok()) {
    webp <- sub("\\.png$", ".webp", png)
    ok <- tryCatch({
      magick::image_write(magick::image_read(png), webp, format = "webp", quality = 100,
                          defines = c(`webp:lossless` = "true"))
      TRUE
    }, error = function(e) FALSE)
    if (ok) return(list(src = webp, contentType = "image/webp", width = "100%", alt = "plot"))
  }
  list(src = png, contentType = "image/png", width = "100%", alt = "plot")
}

# a card with tabs: Interpretation, Output, Table, Plot (each with a Copy button where it is text)
.bam_result_tabs <- function(id, plot = TRUE, table = TRUE, code = FALSE, extra = NULL) {
  ns <- shiny::NS(id)
  tabs <- list(
      bslib::nav_panel("Interpretation", .bam_copy_bar(ns("int_text")), shiny::uiOutput(ns("int")),
                       shiny::tags$pre(id = ns("int_text"), class = "shiny-text-output", style = "display:none"))
    , if (plot) bslib::nav_panel("Plot", extra, shiny::div(class = "bam-plot", shiny::imageOutput(ns("plot"), height = "auto")),
                                 shiny::div(shiny::downloadButton(ns("png"), "PNG", class = "btn-sm"),
                                            shiny::downloadButton(ns("pdf"), "PDF", class = "btn-sm"),
                                            shiny::downloadButton(ns("webp"), "WebP", class = "btn-sm")))
    , if (table) bslib::nav_panel("Table", DT::DTOutput(ns("table")))
    , bslib::nav_panel("Output", .bam_copy_bar(ns("out")), shiny::verbatimTextOutput(ns("out")))
    , if (code) bslib::nav_panel("R code", .bam_copy_bar(ns("code")), shiny::verbatimTextOutput(ns("code")))
    )
  do.call(bslib::navset_card_tab, c(list(id = ns("tabs")), Filter(Negate(is.null), tabs)))
}

.bam_int_html <- function(it) {
  li <- function(s) if (length(s)) shiny::tags$ul(lapply(s, shiny::tags$li)) else NULL
  shiny::tagList(
      shiny::div(class = "alert alert-primary", shiny::strong(it$headline))
    , if (length(it$findings)) shiny::tagList(shiny::h6("Findings"), li(it$findings))
    , if (length(it$recommendations)) shiny::tagList(shiny::h6("Recommendations"), li(it$recommendations))
    , if (length(it$cautions)) shiny::tagList(shiny::h6("Cautions"), li(it$cautions))
    , shiny::tags$details(shiny::tags$summary("How this is read"), shiny::p(it$guide))
    , if (length(it$sources)) shiny::tags$details(class = "bam-sources", shiny::tags$summary("Sources"),
                                                  shiny::tags$ul(lapply(it$sources, shiny::tags$li)))
    )
}

# the server side of a result card: res() a bam_result, table_fun(res) a data frame, plot_fun(res) a ggplot
.bam_result_server <- function(id, res, table_fun = NULL, plot_fun = NULL, code_fun = NULL, hb = shiny::reactive(TRUE)) {
  force(res)
  force(table_fun)
  force(plot_fun)
  force(code_fun)
  shiny::moduleServer(id, function(input, output, session) {
    it <- shiny::reactive({
      shiny::req(res())
      bam_interpret(res(), higher_better = hb())
    })
    output$int <- shiny::renderUI(.bam_int_html(it()))
    output$int_text <- shiny::renderText(paste(format(it(), style = "text", width = 100L), collapse = "\n"))
    output$out <- shiny::renderPrint({
      shiny::req(res())
      op <- options(width = 120L)
      on.exit(options(op))
      print(res())
    })
    if (!is.null(table_fun)) output$table <- DT::renderDT({
      shiny::req(res())
      .bam_dt(table_fun(res()), id)
    })
    if (!is.null(plot_fun)) {
      pl <- shiny::reactive({
        shiny::req(res())
        plot_fun(res(), input)
      })
      output$plot <- shiny::renderImage(.bam_plot_file(pl()), deleteFile = TRUE)
      output$png <- shiny::downloadHandler(paste0(id, ".png"), function(file)
        ggplot2::ggsave(file, pl(), width = 10, height = 7, dpi = 200))
      output$pdf <- shiny::downloadHandler(paste0(id, ".pdf"), function(file)
        ggplot2::ggsave(file, pl(), width = 10, height = 7))
      output$webp <- shiny::downloadHandler(paste0(id, ".webp"), function(file) bam_save(pl(), file, width = 10, height = 7))
    }
    if (!is.null(code_fun)) output$code <- shiny::renderText(code_fun())
  })
}

#-------------------------------------------------------------------------------------------------------------------
# pages
#-------------------------------------------------------------------------------------------------------------------
.bam_app_ui <- function() {
  side <- function(...) bslib::sidebar(width = 300, ...)
  bslib::page_navbar(
      title = "bayesammi", id = "nav", theme = .bam_app_theme(), fillable = FALSE, header = .bam_app_head()
    , bslib::nav_panel("Home", .bam_app_home())
    , bslib::nav_panel("Data", bslib::layout_sidebar(
          sidebar = side(
              shiny::radioButtons("source", "Data", c(`Example: Maiz (9 genotypes, 20 sites)` = "maiz", `CSV file` = "csv"))
            , shiny::conditionalPanel("input.source == 'csv'", shiny::fileInput("file", "CSV file (one row per plot)", accept = ".csv"))
            , shiny::selectInput("col_y", "Response", NULL), shiny::selectInput("col_gen", "Genotype", NULL)
            , shiny::selectInput("col_env", "Environment", NULL), shiny::selectInput("col_rep", "Replicate", NULL)
            , shiny::checkboxInput("hb", "Higher values are better", TRUE))
        , bslib::card(bslib::card_header("The trial"), .bam_copy_bar("data_summary"),
                      shiny::verbatimTextOutput("data_summary"), DT::DTOutput("data_table"))))
    , bslib::nav_panel("Model", bslib::layout_sidebar(
          sidebar = side(
              shiny::selectInput("model", "Model", c(AMMI = "AMMI", `GGE (SREG)` = "GGE"))
            , shiny::numericInput("k", "Multiplicative terms", 2, min = 1, max = 6, step = 1)
            , shiny::checkboxInput("het", "A variance for every environment", FALSE)
            , shiny::numericInput("iter", "Iterations kept per chain", 3000, min = 100, step = 500)
            , shiny::numericInput("burnin", "Burn-in", 1000, min = 0, step = 500)
            , shiny::numericInput("thin", "Thinning", 1, min = 1, step = 1)
            , shiny::numericInput("chains", "Chains", 2, min = 1, max = 6, step = 1)
            , shiny::numericInput("seed", "Seed", 1, step = 1)
            , shiny::numericInput("level", "HPD probability", 0.95, min = 0.5, max = 0.999, step = 0.01)
            , shiny::tags$details(shiny::tags$summary("Priors"),
                                  shiny::numericInput("lambda_var", "Prior variance of the singular values", 1e8),
                                  shiny::numericInput("nu", "Degrees of freedom of the error variance prior", 1, min = 0.01))
            , shiny::actionButton("run", "Fit the model", class = "btn-primary w-100", icon = shiny::icon("play")))
        , .bam_result_tabs("fit", code = TRUE, extra = shiny::selectInput("fit_plot", NULL,
                           c(`traces` = "trace", `posterior histograms` = "density", `biplot of posterior means` = "biplot",
                             `main effects` = "effects")))))
    , bslib::nav_panel("Biplot", bslib::layout_sidebar(
          sidebar = side(
              shiny::checkboxGroupInput("pb", "HPD regions", c(`95%` = 0.05, `99%` = 0.01), selected = c(0.05, 0.01))
            , shiny::radioButtons("which", "Show", c(both = "both", genotypes = "genotype", environments = "environment"))
            , shiny::checkboxInput("stable", "Regions containing the origin", TRUE)
            , shiny::sliderInput("bp_thin", "Share of the draws used", 0.05, 1, 0.25, step = 0.05))
        , .bam_result_tabs("biplot")))
    , bslib::nav_menu("Decisions",
          bslib::nav_panel("Which wins where", .bam_result_tabs("winners"))
        , bslib::nav_panel("Superior performance", bslib::layout_sidebar(
              sidebar = side(shiny::sliderInput("s", "Selection intensity", 0.05, 0.6, 0.2, step = 0.05),
                             shiny::selectInput("sp_plot", "Plot", c(probabilities = "marginal", `by environment` = "conditional",
                                                                     pairwise = "pairwise")))
            , .bam_result_tabs("superiority")))
        , bslib::nav_panel("Ranks", .bam_result_tabs("ranks"))
        , bslib::nav_panel("Stability", .bam_result_tabs("stability"))
        , bslib::nav_panel("Predictions", .bam_result_tabs("predict", extra = shiny::selectInput("pr_plot", NULL,
                           c(cells = "cells", `genotype means` = "genotypes")))))
    , bslib::nav_menu("Stability",
          bslib::nav_panel("Regression models", bslib::layout_sidebar(
              sidebar = side(
                  shiny::selectInput("st_model", "Model", c(`Eberhart-Russell` = "ER", `Finlay-Wilkinson` = "FW",
                                                            `Shukla (stability variance)` = "Shukla"))
                , shiny::numericInput("st_alpha", "Probability of the safety-first index", 0.05, min = 0.001,
                                      max = 0.5, step = 0.01)
                , shiny::actionButton("st_run", "Fit the model", class = "btn-primary w-100", icon = shiny::icon("play"))
                , shiny::helpText("Uses the data of the Data page and the chain settings of the Model page."))
            , .bam_result_tabs("stab", code = TRUE, extra = shiny::selectInput("st_plot", NULL,
                               c(`slopes and means` = "slopes", `regression lines` = "lines",
                                 `deviation variances` = "deviations", traces = "trace",
                                 `posterior histograms` = "density")))))
        , bslib::nav_panel("Stability indices", bslib::layout_sidebar(
              sidebar = side(
                  shiny::radioButtons("ind_src", "Cell means of the fit of", c(`the Model page (AMMI or GGE)` = "ammi",
                                                                              `the regression models page` = "stab"))
                , shiny::numericInput("ind_prob", "Error probability of Annicchiarico's index", 0.25, min = 0.01,
                                      max = 0.5, step = 0.05)
                , shiny::numericInput("ind_alpha", "Probability of the safety-first index", 0.05, min = 0.001,
                                      max = 0.5, step = 0.01))
            , .bam_result_tabs("indices")))
        , bslib::nav_panel("AMMI indices", bslib::layout_sidebar(
              sidebar = side(shiny::sliderInput("wy", "Weight of the mean in WAASY", 0, 100, 50, step = 5))
            , .bam_result_tabs("ammi_indices"))))
    , bslib::nav_menu("Model choice",
          bslib::nav_panel("Compare models", bslib::layout_sidebar(
              sidebar = side(shiny::checkboxGroupInput("cmp_k", "Terms", 1:4, selected = 1:3, inline = TRUE),
                             shiny::checkboxInput("cmp_het", "Also with a variance per environment", TRUE),
                             shiny::actionButton("compare", "Compare", class = "btn-primary w-100"),
                             shiny::helpText("Uses the settings of the Model page; each model is fitted once."))
            , .bam_result_tabs("compare")))
        , bslib::nav_panel("Predictive checks", .bam_result_tabs("ppc", table = TRUE)))
    )
}

.bam_app_home <- function() {
  card <- function(title, text) bslib::card(bslib::card_header(title), shiny::p(text))
  shiny::tagList(
      shiny::h3("Bayesian AMMI and GGE analysis of multi-environment trials")
    , shiny::p("The Gibbs sampler of Crossa et al. (2011), extended to unbalanced trials, a variance for every",
               "environment and the GGE (SREG) model. Every parameter has a posterior distribution, so every statement",
               "about genotypes and environments comes with a probability.")
    , bslib::layout_column_wrap(width = 1 / 3,
          card("1. Data", "The example Maiz or your CSV file, one row per plot; missing plots and cells are allowed.")
        , card("2. Model", "AMMI or GGE, number of terms, variances, chains; traces and convergence diagnostics.")
        , card("3. Biplot", "Highest posterior density regions: stable genotypes and groups of environments.")
        , card("4. Decisions", "Probabilities of winning, of superior performance and stability, ranks, predictions.")
        , card("5. Stability", "Bayesian Eberhart-Russell, Finlay-Wilkinson and Shukla models; classical and AMMI-based stability indices with HPD intervals.")
        , card("6. Model choice", "WAIC and DIC of several models; posterior predictive checks.")
        , card("Every page", "Interpretation, output as R prints it, tables and plots, all with Copy or download."))
    )
}

#-------------------------------------------------------------------------------------------------------------------
# server
#-------------------------------------------------------------------------------------------------------------------
.bam_app_server <- function(input, output, session) {
  raw <- shiny::reactive({
    if (identical(input$source, "csv")) {
      shiny::req(input$file)
      data.table::fread(input$file$datapath, data.table = FALSE)
    } else {
      e <- new.env()
      utils::data("Maiz", package = "bayesammi", envir = e)
      e$Maiz
    }
  })
  shiny::observeEvent(raw(), {
    nm <- names(raw())
    pick <- function(pat, default) { h <- grep(pat, nm, ignore.case = TRUE, value = TRUE); if (length(h)) h[1L] else default }
    shiny::updateSelectInput(session, "col_y", choices = nm, selected = pick("^y$|yield", nm[length(nm)]))
    shiny::updateSelectInput(session, "col_gen", choices = nm, selected = pick("gen|entry|variety", nm[1L]))
    shiny::updateSelectInput(session, "col_env", choices = nm, selected = pick("env|site|loc", nm[2L]))
    shiny::updateSelectInput(session, "col_rep", choices = nm, selected = pick("rep|block", nm[3L]))
  })
  output$data_summary <- shiny::renderPrint({
    d <- raw()
    shiny::req(input$col_y %in% names(d), input$col_gen %in% names(d), input$col_env %in% names(d))
    N <- table(d[[input$col_gen]], d[[input$col_env]])
    cat(sprintf("%d plots; %d genotypes; %d environments; replicates per cell from %d to %d; %d empty cell(s)\n\n",
                nrow(d), nrow(N), ncol(N), min(N), max(N), sum(N == 0)))
    print(summary(d[[input$col_y]]))
  })
  output$data_table <- DT::renderDT(.bam_dt(raw(), "data"))

  fit <- shiny::reactiveVal(NULL)
  args <- shiny::reactive(list(
      .data = raw(), .y = input$col_y, .gen = input$col_gen, .env = input$col_env, .rep = input$col_rep
    , .nIter = as.integer(input$iter), .k = as.integer(input$k), .burnin = as.integer(input$burnin)
    , .thin = as.integer(input$thin), .chains = as.integer(input$chains)
    , .priors = bam_priors(lambda_var = input$lambda_var, nu = input$nu), .level = input$level
    , .seed = as.integer(input$seed), .align = TRUE, .model = input$model, .het = isTRUE(input$het), .waic = TRUE))
  shiny::observeEvent(input$run, {
    r <- tryCatch(shiny::withProgress(message = "Sampling", value = 0.5, do.call(bayes_ammi, args())),
                  error = function(e) e)
    if (inherits(r, "error")) {
      shiny::showNotification(conditionMessage(r), type = "error", duration = 10)
      return()
    }
    fit(r)
  })
  hb <- shiny::reactive(isTRUE(input$hb))
  code <- shiny::reactive(.bam_app_code(args(), input$source))
  # stability models
  stab <- shiny::reactiveVal(NULL)
  stab_args <- shiny::reactive(list(
      .data = raw(), .y = input$col_y, .gen = input$col_gen, .env = input$col_env, .rep = input$col_rep
    , .model = input$st_model, .nIter = as.integer(input$iter), .burnin = as.integer(input$burnin)
    , .thin = as.integer(input$thin), .chains = as.integer(input$chains), .priors = bam_stab_priors()
    , .level = input$level, .alpha = input$st_alpha, .seed = as.integer(input$seed), .waic = TRUE))
  shiny::observeEvent(input$st_run, {
    r <- tryCatch(shiny::withProgress(message = "Sampling", value = 0.5, do.call(bayes_stability, stab_args())),
                  error = function(e) e)
    if (inherits(r, "error")) {
      shiny::showNotification(conditionMessage(r), type = "error", duration = 10)
      return()
    }
    stab(r)
  })
  .bam_result_server("stab", stab, function(x) x$stability, function(x, i) plot(x, which = input$st_plot),
                     shiny::reactive(.bam_app_code_stab(stab_args(), input$source)), hb)
  .bam_result_server("indices", shiny::reactive({
    src <- if (identical(input$ind_src, "stab")) stab() else fit()
    shiny::req(src)
    bam_indices(src, higher_better = hb(), .prob = input$ind_prob, .alpha = input$ind_alpha)
  }), function(x) x$indices, function(x, i) plot(x), hb = hb)
  .bam_result_server("ammi_indices", shiny::reactive({
    shiny::req(fit())
    bam_ammi_indices(fit(), higher_better = hb(), .weight_y = input$wy)
  }), function(x) x$indices, function(x, i) plot(x), hb = hb)
  .bam_result_server("fit", fit, function(x) x$summary, function(x, i) plot(x, which = input$fit_plot), code, hb)
  bp <- shiny::reactive({
    shiny::req(fit(), length(input$pb) > 0L)
    if (fit()$settings$k < 2L) return(NULL)
    biplots(fit(), burnin = 0, thin = input$bp_thin, pb = as.numeric(input$pb), plot_stable = isTRUE(input$stable),
            plot_unstable = TRUE, .which = input$which)
  })
  .bam_result_server("biplot", bp, function(x) x$biplot_data, function(x, i) x$biplot, hb = hb)
  .bam_result_server("winners", shiny::reactive({ shiny::req(fit()); bam_winners(fit(), higher_better = hb()) }),
                     function(x) x$winners, function(x, i) plot(x), hb = hb)
  .bam_result_server("superiority", shiny::reactive({ shiny::req(fit()); bam_superiority(fit(), .s = input$s, higher_better = hb()) }),
                     function(x) x$marginal, function(x, i) plot(x, which = input$sp_plot), hb = hb)
  .bam_result_server("ranks", shiny::reactive({ shiny::req(fit()); bam_ranks(fit(), higher_better = hb()) }),
                     function(x) x$summary, function(x, i) plot(x), hb = hb)
  .bam_result_server("stability", shiny::reactive({ shiny::req(fit()); bam_stability(fit(), higher_better = hb()) }),
                     function(x) x$stability, function(x, i) plot(x), hb = hb)
  .bam_result_server("predict", shiny::reactive({ shiny::req(fit()); bam_predict(fit()) }),
                     function(x) x$cells, function(x, i) plot(x, which = input$pr_plot), hb = hb)
  .bam_result_server("ppc", shiny::reactive({ shiny::req(fit()); bam_ppc(fit(), .draws = 500, .seed = 1) }),
                     function(x) x$checks, function(x, i) plot(x), hb = hb)
  cmp <- shiny::reactiveVal(NULL)
  shiny::observeEvent(input$compare, {
    a <- args()
    ks <- as.integer(input$cmp_k)
    het <- if (isTRUE(input$cmp_het)) c(FALSE, TRUE) else FALSE
    grid <- expand.grid(k = ks, het = het)
    fits <- list()
    shiny::withProgress(message = "Fitting models", value = 0, {
      for (r in seq_len(nrow(grid))) {
        a$.k <- grid$k[r]
        a$.het <- grid$het[r]
        nm <- sprintf("%s%d%s", a$.model, a$.k, if (a$.het) "-het" else "")
        fits[[nm]] <- tryCatch(do.call(bayes_ammi, a), error = function(e) NULL)
        shiny::incProgress(1 / nrow(grid))
      }
    })
    fits <- Filter(Negate(is.null), fits)
    if (length(fits)) cmp(do.call(bam_compare, fits))
  })
  .bam_result_server("compare", cmp, function(x) x$table, function(x, i) plot(x), hb = hb)
}

# the R code of the current fit, every argument on its own line
# the R code of the stability fit of the application
.bam_app_code_stab <- function(a, source) {
  vals <- c(.data = if (identical(source, "csv")) "trial" else "Maiz", .y = a$.y, .gen = a$.gen, .env = a$.env,
            .rep = a$.rep, .model = sprintf("\"%s\"", a$.model), .nIter = a$.nIter, .burnin = a$.burnin,
            .thin = a$.thin, .chains = a$.chains, .priors = "bam_stab_priors()", .level = a$.level,
            .alpha = a$.alpha, .seed = a$.seed, .waic = "TRUE")
  w <- max(nchar(names(vals)))
  body <- paste0("    ", c(" ", rep(",", length(vals) - 1L)), " ", formatC(names(vals), width = -w), " = ", vals,
                 collapse = "\n")
  paste0("library(bayesammi)\n", if (identical(source, "csv")) "trial <- data.table::fread(\"trial.csv\")\n" else
           "data(Maiz)\n", "fs <-\n  bayes_stability(\n", body, "\n    )\nfs\nbam_interpret(\n    x             = fs\n",
         "  , higher_better = TRUE\n  , top           = 5\n  )\nbam_cite(\n    x     = fs\n  , style = \"text\"\n  )\n")
}

.bam_app_code <- function(a, source) {
  vals <- c(.data = if (identical(source, "csv")) "trial" else "Maiz", .y = a$.y, .gen = a$.gen, .env = a$.env,
            .rep = a$.rep, .nIter = a$.nIter, .k = a$.k, .burnin = a$.burnin, .thin = a$.thin, .chains = a$.chains,
            .priors = sprintf("bam_priors(lambda_var = %s, nu = %s)", format(a$.priors$lambda_var), format(a$.priors$nu)),
            .level = a$.level, .seed = a$.seed, .align = "TRUE", .model = sprintf("\"%s\"", a$.model),
            .het = as.character(a$.het), .waic = "TRUE")
  w <- max(nchar(names(vals)))
  body <- paste0("    ", c(" ", rep(",", length(vals) - 1L)), " ", formatC(names(vals), width = -w), " = ", vals,
                 collapse = "\n")
  paste0("library(bayesammi)\n", if (identical(source, "csv")) "trial <- data.table::fread(\"trial.csv\")\n" else
           "data(Maiz)\n", "fit <-\n  bayes_ammi(\n", body, "\n    )\nfit\nbam_interpret(\n    x             = fit\n",
         "  , higher_better = TRUE\n  , top           = 5\n  )\n")
}
