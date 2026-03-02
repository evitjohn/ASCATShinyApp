library(shiny)
library(bslib)
library(shinycssloaders)
library(ASCAT)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load Bioconductor packages for ASCAT function
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  stop("Please install GenomicRanges: BiocManager::install('GenomicRanges')")
}
if (!requireNamespace("IRanges", quietly = TRUE)) {
  stop("Please install IRanges: BiocManager::install('IRanges')")
}

# Source the new ASCAT function
source("ascat_runAscat_ForShiny.R")

# ─── UI ───────────────────────────────────────────────────────────────────────

ui <- page_fillable(
  title = "ASCAT Analysis",

  theme = bs_theme(
    bootswatch   = "flatly",
    primary      = "#2C6E49",
    base_font    = font_google("IBM Plex Sans"),
    heading_font = font_google("IBM Plex Mono")
  ),

  tags$head(tags$style(HTML("
    .card-header { font-weight: 600; font-size: 0.85rem; letter-spacing: 0.04em;
                   text-transform: uppercase; color: #555; }
    .sidebar { background: #f8f9fa !important; }
    #runBtn  { width: 100%; background-color: #2C6E49; border-color: #2C6E49;
               color: white; font-weight: 600; }
    #runBtn:hover { background-color: #1e4d33; border-color: #1e4d33; }
    .status-box     { padding: 8px 12px; border-radius: 6px; font-size: 0.85rem; margin-top: 8px; }
    .status-idle    { background: #e9ecef; color: #6c757d; }
    .status-running { background: #fff3cd; color: #856404; }
    .status-done    { background: #d1e7dd; color: #0a3622; }
    .status-error   { background: #f8d7da; color: #842029; }
  "))),

  layout_sidebar(
    sidebar = sidebar(
      width = 280,
      card(
        card_header("Input"),
        fileInput("ascatFile", "Upload ASCAT .Rdata file",
                  accept   = c(".Rdata", ".RData", ".rda"),
                  multiple = FALSE,
                  placeholder = "ascat.bc.Rdata"),
        hr(style = "margin: 8px 0;"),
        checkboxInput("manualRhoPsi",
                      "Enable manual rho / psi selection",
                      value = FALSE),
        conditionalPanel(
          condition = "input.manualRhoPsi == true",
          sliderInput("rho", "Rho (tumour purity, %)",
                      min = 1, max = 100, value = 50, step = 1),
          sliderInput("psi", "Psi (tumour ploidy)",
                      min = 1.0, max = 10.0, value = 2.0, step = 0.1)
        ),
        hr(style = "margin: 8px 0;"),
        numericInput("gamma", "Gamma",
                     value = 1, min = 0.1, max = 2, step = 0.05),
        hr(style = "margin: 8px 0;"),
        actionButton("runBtn",
                     label = tagList(shiny::icon("play"), "Run Analysis"),
                     class = "btn-primary"),
        uiOutput("statusUI")
      )
    ),

    layout_columns(
      # Left: three plots stacked
      layout_columns(
        card(card_header("Segmented Data"),
             plotOutput("segmentedPlot", height = "300px") |>
               withSpinner(type = 4, color = "#2C6E49")),
        card(card_header("Raw / Non-rounded Profile"),
             plotOutput("nonroundedPlot", height = "250px") |>
               withSpinner(type = 4, color = "#2C6E49")),
        card(card_header("ASCAT Copy-number Profile"),
             plotOutput("ascatProfilePlot", height = "250px") |>
               withSpinner(type = 4, color = "#2C6E49")),
        col_widths = c(12, 12, 12)
      ),
      # Right: sunrise
      card(card_header("Sunrise / Goodness-of-fit"),
           plotOutput("sunrisePlot", height = "400px") |>
             withSpinner(type = 4, color = "#2C6E49")),
      col_widths = c(8, 4)
    )
  )
)


# ─── SERVER ───────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  options(shiny.maxRequestSize = 100 * 1024^2)

  # Reactive values: store recorded plot objects
  plots <- reactiveValues(
    segmented    = NULL,   # path to PNG (segmented plot not changed)
    nonrounded   = NULL,   # recorded plot object
    ascatProfile = NULL,   # recorded plot object
    sunrise      = NULL,   # recorded plot object
    status       = "idle",
    statusMsg    = "Waiting for input."
  )

  output$statusUI <- renderUI({
    cls <- switch(plots$status,
                  idle    = "status-idle",
                  running = "status-running",
                  done    = "status-done",
                  error   = "status-error",
                  "status-idle")
    div(class = paste("status-box", cls), plots$statusMsg)
  })

  # Segmented plot: still using PNG approach
  output$segmentedPlot <- renderPlot({
    req(plots$segmented)
    img <- png::readPNG(plots$segmented)
    grid::grid.raster(img)
  })

  # New approach: use replayPlot for sunrise, nonrounded, and ASCAT profile
  output$sunrisePlot <- renderPlot({
    req(plots$sunrise)
    replayPlot(plots$sunrise)
  })

  output$nonroundedPlot <- renderPlot({
    req(plots$nonrounded)
    replayPlot(plots$nonrounded)
  })

  output$ascatProfilePlot <- renderPlot({
    req(plots$ascatProfile)
    replayPlot(plots$ascatProfile)
  })

  # ── Main run ──────────────────────────────────────────────────────────────
  observeEvent(input$runBtn, {
    req(input$ascatFile)

    plots$segmented    <- NULL
    plots$nonrounded   <- NULL
    plots$ascatProfile <- NULL
    plots$sunrise      <- NULL
    plots$status       <- "running"
    plots$statusMsg    <- "Loading ASCAT object…"

    # Load file
    env <- new.env(parent = emptyenv())
    tryCatch(
      load(input$ascatFile$datapath, envir = env),
      error = function(e) {
        plots$status    <<- "error"
        plots$statusMsg <<- paste("Error loading file:", conditionMessage(e))
      }
    )
    if (plots$status == "error") return()

    obj_name <- intersect(ls(env), c("ascat.bc", "ASCATobj", "ascat_bc"))
    if (length(obj_name) == 0) obj_name <- ls(env)[1]
    ascat.bc <- get(obj_name[1], envir = env)

    tmp_dir <- file.path(tempdir(), paste0("ascat_run_", as.integer(Sys.time())))
    dir.create(tmp_dir, showWarnings = FALSE)

    # ── Segmented data plot ─────────────────────────────────────────────────
    plots$statusMsg <- "Plotting segmented data…"
    tryCatch(
      ascat.plotSegmentedData(ascat.bc, img.dir = tmp_dir, img.prefix = "shiny_"),
      error = function(e) message("segmentedData error: ", conditionMessage(e))
    )
    seg_file <- file.path(tmp_dir, paste0("shiny_", ascat.bc$samples[1], ".ASPCF.png"))
    if (file.exists(seg_file)) plots$segmented <- seg_file

    # ── Determine rho/psi ───────────────────────────────────────────────────
    if (input$manualRhoPsi) {
      rho_final <- input$rho / 100
      psi_final <- input$psi
      plots$statusMsg <- "Running ASCAT with manual rho/psi…"
    } else {
      rho_final <- NA
      psi_final <- NA
      plots$statusMsg <- "Running ASCAT optimisation…"
    }

    # ── Run ASCAT with new function ─────────────────────────────────────────
    ascat.output <- tryCatch(
      ascat.runAscat(
        ASCATobj       = ascat.bc,
        gamma          = input$gamma,
        rho_manual     = rho_final,
        psi_manual     = psi_final,
        pdfPlot        = FALSE,
        write_segments = FALSE,
        img.dir        = tmp_dir,
        img.prefix     = "shiny_"
      ),
      error = function(e) {
        plots$status    <<- "error"
        plots$statusMsg <<- paste("ASCAT error:", conditionMessage(e))
        NULL
      }
    )
    if (is.null(ascat.output)) return()

    # Check if ASCAT found a solution
    if (is.null(ascat.output$purity) || is.na(ascat.output$purity[1])) {
      plots$status    <- "error"
      plots$statusMsg <- "ASCAT could not find a solution. Enable manual rho/psi and try again."
      return()
    }

    # ── Extract plots from ascat.output ─────────────────────────────────────
    plots$statusMsg <- "Extracting plots…"
    
    # The plots are stored in ascat.output$plots[[1]] for the first sample
    if (!is.null(ascat.output$plots) && length(ascat.output$plots) > 0) {
      sample_plots <- ascat.output$plots[[1]]
      
      plots$sunrise      <- sample_plots$sunrise_plot
      plots$nonrounded   <- sample_plots$nonrounded_plot
      plots$ascatProfile <- sample_plots$ascat_profile_plot
    }

    # ── Extract metrics ─────────────────────────────────────────────────────
    rho_final   <- ascat.output$purity[1]
    psi_final   <- ascat.output$psi[1]
    ploidy      <- ascat.output$ploidy[1]
    gof         <- ascat.output$goodnessOfFit[1]

    plots$status    <- "done"
    plots$statusMsg <- paste0(
      "Done.  Purity: ",  sprintf("%.0f%%", rho_final * 100), "  |  ",
      "Ploidy: ",         sprintf("%.2f",   ploidy),           "  |  ",
      "GoF: ",            sprintf("%.1f%%", gof)
    )
  })
}

# ─── RUN ──────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
