library(shiny)
library(Seurat)
library(qs)
library(ggplot2)
source("/Users/chen.13589/Documents/Tigerrr07/scripts/R_scripts/Dotplot.R")
source("/Users/chen.13589/Documents/Tigerrr07/scripts/R_scripts/Heatmap.R")

DEFAULT_SEURAT_PATH <- "/Users/chen.13589/Documents/Tigerrr07/scripts/scRNA-seq/write/seurat.qs"

ui <- fluidPage(
  titlePanel("Seurat Viewer"),
  sidebarLayout(
    sidebarPanel(
      textInput("seurat_path", "Seurat .qs path", value = DEFAULT_SEURAT_PATH),
      actionButton("load_btn", "Load object"),
      tags$hr(),
      conditionalPanel(
        "input.main_tab == 'UMAP/TSNE'",
        uiOutput("assay_ui"),
        uiOutput("reduction_ui"),
        selectInput("color_by", "Color by", choices = c("metadata", "gene")),
        uiOutput("feature_ui"),
        uiOutput("split_ui"),
        numericInput("split_nrow", "Split rows", value = 1, min = 1, max = 12, step = 1),
        numericInput("split_ncol", "Split columns (0 = auto)", value = 0, min = 0, max = 12, step = 1),
        checkboxInput("label", "Label clusters (metadata only)", value = FALSE),
        sliderInput("pt_size", "Point size", min = 0.1, max = 2.0, value = 0.4, step = 0.1)
      ),
      conditionalPanel(
        "input.main_tab == 'Dotplot'",
        uiOutput("assay_ui"),
        uiOutput("dotplot_genes_ui"),
        uiOutput("dotplot_group_ui")
      )
    ),
    mainPanel(
      verbatimTextOutput("status"),
      tabsetPanel(
        id = "main_tab",
        selected = "UMAP/TSNE",
        tabPanel(
          "UMAP/TSNE",
          tabsetPanel(
            id = "umap_subtab",
            selected = "UMAP/TSNE",
            tabPanel(
              "UMAP/TSNE",
              fluidRow(
                column(6, uiOutput("subset_ui")),
                column(2, numericInput("plot_width", "UMAP width (px)", value = 500, min = 50, max = 2000, step = 50)),
                column(2, numericInput("plot_height", "UMAP height (px)", value = 350, min = 50, max = 2000, step = 50))
              ),
              downloadButton("download_umap", "Save UMAP/TSNE (PDF)"),
              plotOutput("umap_plot")
            ),
            tabPanel(
              "Violin",
              fluidRow(
                column(2, numericInput("violin_width", "Violin width (px)", value = 800, min = 50, max = 2000, step = 50)),
                column(2, numericInput("violin_height", "Violin height (px)", value = 250, min = 50, max = 800, step = 25))
              ),
              downloadButton("download_violin", "Save Violin (PDF)"),
              plotOutput("violin_plot")
            )
          )
        ),
        tabPanel(
          "Dotplot",
          tabsetPanel(
            id = "dotplot_subtab",
            selected = "Dotplot",
            tabPanel(
              "Dotplot",
              fluidRow(
                column(2, numericInput("dotplot_width", "Dotplot width (px)", value = 900, min = 50, max = 2000, step = 50)),
                column(2, numericInput("dotplot_height", "Dotplot height (px)", value = 500, min = 50, max = 2000, step = 50))
              ),
              downloadButton("download_dotplot", "Save Dotplot (PDF)"),
              plotOutput("dotplot_plot")
            ),
            tabPanel(
              "Heatmap",
              fluidRow(
                column(2, numericInput("heatmap_width", "Heatmap width (px)", value = 900, min = 50, max = 2000, step = 50)),
                column(2, numericInput("heatmap_height", "Heatmap height (px)", value = 500, min = 50, max = 2000, step = 50))
              ),
              downloadButton("download_heatmap", "Save Heatmap (PDF)"),
              plotOutput("heatmap_plot")
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  seu <- reactiveVal(NULL)
  gene_choices <- reactiveVal(character())
  dotplot_gene_choices <- reactiveVal(character())

  px_to_in <- function(px) {
    px / 72
  }

  get_assay_features <- function(obj, assay) {
    tryCatch({
      rownames(GetAssayData(obj, assay = assay, layer = "data"))
    }, error = function(e) {
      rownames(GetAssayData(obj, assay = assay, slot = "data"))
    })
  }

  observeEvent(input$load_btn, {
    path <- input$seurat_path
    if (!file.exists(path)) {
      showNotification(paste("File not found:", path), type = "error")
      seu(NULL)
      return()
    }
    obj <- tryCatch(qread(path), error = function(e) e)
    if (inherits(obj, "error")) {
      showNotification(paste("Failed to read:", obj$message), type = "error")
      seu(NULL)
      return()
    }
    if (!inherits(obj, "Seurat")) {
      showNotification("Loaded object is not a Seurat object", type = "error")
      seu(NULL)
      return()
    }
    seu(obj)
  })

  observeEvent(list(seu(), input$assay, input$color_by), {
    obj <- seu()
    if (is.null(obj) || is.null(input$assay) || input$color_by != "gene") {
      gene_choices(character())
      return()
    }
    feats <- get_assay_features(obj, input$assay)
    gene_choices(feats)
    updateSelectizeInput(session, "gene", choices = feats, selected = feats[1], server = TRUE)
  })

  output$status <- renderText({
    obj <- seu()
    if (is.null(obj)) {
      return("No Seurat object loaded.")
    }
    paste0("Loaded Seurat object with ", ncol(obj), " cells and ", nrow(obj), " features.")
  })

  output$assay_ui <- renderUI({
    obj <- seu()
    if (is.null(obj)) return(NULL)
    selectInput("assay", "Assay", choices = Assays(obj), selected = DefaultAssay(obj))
  })

  output$reduction_ui <- renderUI({
    obj <- seu()
    if (is.null(obj)) return(NULL)
    reds <- Reductions(obj)
    if (length(reds) == 0) {
      return(tags$div("No reductions found in object."))
    }
    preferred <- if ("umap" %in% reds) "umap" else reds[[1]]
    selectInput("reduction", "Reduction", choices = reds, selected = preferred)
  })

  output$feature_ui <- renderUI({
    obj <- seu()
    req(obj)
    if (input$color_by == "metadata") {
      meta_cols <- colnames(obj@meta.data)
      selectInput("meta", "Metadata", choices = meta_cols)
    } else {
      selectizeInput("gene", "Gene", choices = character(0), options = list(maxOptions = 2000))
    }
  })

  output$split_ui <- renderUI({
    obj <- seu()
    req(obj)
    meta_cols <- colnames(obj@meta.data)
    choices <- c("none", meta_cols)
    selectInput("split_by", "Split by (metadata)", choices = choices, selected = "none")
  })

  output$subset_ui <- renderUI({
    obj <- seu()
    req(obj)
    meta_cols <- colnames(obj@meta.data)
    tagList(
      selectInput("subset_col", "Subset by (metadata)", choices = c("none", meta_cols), selected = "none"),
      uiOutput("subset_val_ui")
    )
  })

  output$subset_val_ui <- renderUI({
    obj <- seu()
    req(obj)
    if (is.null(input$subset_col) || input$subset_col == "none") return(NULL)
    vals <- sort(unique(as.character(obj@meta.data[[input$subset_col]])))
    selectInput("subset_val", "Subset value", choices = vals, selected = vals[1])
  })

  output$dotplot_genes_ui <- renderUI({
    obj <- seu()
    req(obj)
    req(input$assay)
    selectizeInput(
      "dotplot_genes",
      "Dotplot genes",
      choices = NULL,
      multiple = TRUE,
      options = list(maxOptions = 2000)
    )
  })

  output$dotplot_group_ui <- renderUI({
    obj <- seu()
    req(obj)
    meta_cols <- colnames(obj@meta.data)
    selectizeInput(
      "dotplot_group",
      "Dotplot group (metadata)",
      choices = meta_cols,
      selected = meta_cols[1],
      multiple = TRUE,
      options = list(maxOptions = 2000)
    )
  })

  observeEvent(list(seu(), input$assay), {
    obj <- seu()
    if (is.null(obj) || is.null(input$assay)) {
      dotplot_gene_choices(character())
      return()
    }
    feats <- get_assay_features(obj, input$assay)
    dotplot_gene_choices(feats)
    updateSelectizeInput(session, "dotplot_genes", choices = feats, selected = NULL, server = TRUE)
  })

  obj_filtered <- reactive({
    obj <- seu()
    req(obj)
    if (is.null(input$subset_col) || input$subset_col == "none") return(obj)
    req(input$subset_val)
    meta <- obj@meta.data
    if (!(input$subset_col %in% colnames(meta))) return(obj)
    keep <- rownames(meta)[as.character(meta[[input$subset_col]]) == input$subset_val]
    if (length(keep) == 0) return(obj)
    subset(obj, cells = keep)
  })

  umap_plot_obj <- reactive({
    obj <- obj_filtered()
    req(obj)
    req(input$reduction)
    req(input$assay)

    obj_use <- obj
    DefaultAssay(obj_use) <- input$assay

    if (!(input$reduction %in% Reductions(obj))) {
      validate("Selected reduction not found in object.")
    }

    split_by <- if (!is.null(input$split_by) && input$split_by != "none") input$split_by else NULL
    split_nrow <- if (!is.null(input$split_nrow) && input$split_nrow > 0) input$split_nrow else NULL
    split_ncol <- if (!is.null(input$split_ncol) && input$split_ncol > 0) input$split_ncol else NULL
    if (!is.null(split_by) && is.null(split_ncol) && !is.null(split_nrow)) {
      n_panels <- length(unique(obj_use@meta.data[[split_by]]))
      split_ncol <- ceiling(n_panels / split_nrow)
    }

    if (input$color_by == "metadata") {
      req(input$meta)
      meta_vec <- obj_use@meta.data[[input$meta]]
      if (is.numeric(meta_vec)) {
        p <- FeaturePlot(
          obj_use,
          reduction = input$reduction,
          features = input$meta,
          pt.size = input$pt_size,
          split.by = split_by,
          ncol = split_ncol
        )
        if (is.null(split_by)) {
          p <- p + ggtitle(input$meta)
        }
      } else {
        p <- DimPlot(
          obj_use,
          reduction = input$reduction,
          group.by = input$meta,
          label = input$label,
          pt.size = input$pt_size,
          split.by = split_by,
          ncol = split_ncol
        )
      }
    } else {
      req(input$gene)
      p <- FeaturePlot(
        obj_use,
        reduction = input$reduction,
        features = input$gene,
        pt.size = input$pt_size,
        split.by = split_by,
        ncol = split_ncol
      )
      if (is.null(split_by)) {
        p <- p + ggtitle(input$gene)
      }
    }

    if (!is.null(split_by) && !is.null(split_nrow) && requireNamespace("patchwork", quietly = TRUE)) {
      p <- p + patchwork::plot_layout(nrow = split_nrow, ncol = split_ncol)
    }

    p
  })

  output$umap_plot <- renderPlot({
    umap_plot_obj()
  }, width = function() input$plot_width, height = function() input$plot_height)

  output$violin_plot <- renderPlot({
    obj <- obj_filtered()
    req(obj)
    req(input$assay)

    if (is.null(input$split_by) || input$split_by == "none") {
      return(NULL)
    }

    obj_use <- obj
    DefaultAssay(obj_use) <- input$assay

    feature_name <- NULL
    if (input$color_by == "metadata") {
      req(input$meta)
      meta_vec <- obj_use@meta.data[[input$meta]]
      if (!is.numeric(meta_vec)) {
        return(NULL)
      }
      feature_name <- input$meta
    } else {
      req(input$gene)
      feature_name <- input$gene
    }

    VlnPlot(
      obj_use,
      features = feature_name,
      group.by = input$split_by,
      pt.size = 0
    ) +
      ggtitle(feature_name) +
      theme_classic() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, color = "black"),
        axis.title = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", size = 10),
        axis.ticks.length = unit(4, "pt")
      )
  }, width = function() input$violin_width, height = function() input$violin_height)

  dotplot_plot_obj <- reactive({
    obj <- obj_filtered()
    req(obj)
    req(input$assay)
    req(input$dotplot_genes)
    req(input$dotplot_group)
    if (length(input$dotplot_genes) == 0) return(NULL)
    if (length(input$dotplot_group) == 0) return(NULL)

    obj_use <- obj
    DefaultAssay(obj_use) <- input$assay

    scDotPlot(obj_use, features = input$dotplot_genes, group.by = input$dotplot_group)
  })

  output$dotplot_plot <- renderPlot({
    dotplot_plot_obj()
  }, width = function() input$dotplot_width, height = function() input$dotplot_height)

  output$heatmap_plot <- renderPlot({
    obj <- obj_filtered()
    req(obj)
    req(input$assay)
    req(input$dotplot_genes)
    req(input$dotplot_group)
    if (length(input$dotplot_genes) == 0) return(NULL)
    if (length(input$dotplot_group) == 0) return(NULL)

    obj_use <- obj
    DefaultAssay(obj_use) <- input$assay

    scHeatmap(
      obj_use,
      features = input$dotplot_genes,
      group.by = input$dotplot_group
    )
  }, width = function() input$heatmap_width, height = function() input$heatmap_height)

  output$download_heatmap <- downloadHandler(
    filename = function() {
      paste0("heatmap_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(input$show_heatmap)
      obj <- obj_filtered()
      req(obj)
      req(input$assay)
      req(input$dotplot_genes)
      req(input$dotplot_group)
      if (length(input$dotplot_genes) == 0) return(NULL)
      if (length(input$dotplot_group) == 0) return(NULL)

      obj_use <- obj
      DefaultAssay(obj_use) <- input$assay

      scHeatmap(
        obj_use,
        features = input$dotplot_genes,
        group.by = input$dotplot_group,
        save_path = file,
        width = px_to_in(input$heatmap_width),
        height = px_to_in(input$heatmap_height),
        units = "in"
      )
    }
  )

  output$download_umap <- downloadHandler(
    filename = function() {
      paste0("umap_tsne_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      p <- umap_plot_obj()
      req(p)
      ggsave(
        filename = file,
        plot = p,
        device = "pdf",
        width = px_to_in(input$plot_width),
        height = px_to_in(input$plot_height),
        units = "in"
      )
    }
  )

  output$download_dotplot <- downloadHandler(
    filename = function() {
      paste0("dotplot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      p <- dotplot_plot_obj()
      req(p)
      ggsave(
        filename = file,
        plot = p,
        device = "pdf",
        width = px_to_in(input$dotplot_width),
        height = px_to_in(input$dotplot_height),
        units = "in"
      )
    }
  )

  output$download_violin <- downloadHandler(
    filename = function() {
      paste0("violin_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      p <- isolate({
        obj <- obj_filtered()
        req(obj)
        req(input$assay)

        if (is.null(input$split_by) || input$split_by == "none") {
          return(NULL)
        }

        obj_use <- obj
        DefaultAssay(obj_use) <- input$assay

        feature_name <- NULL
        if (input$color_by == "metadata") {
          req(input$meta)
          meta_vec <- obj_use@meta.data[[input$meta]]
          if (!is.numeric(meta_vec)) {
            return(NULL)
          }
          feature_name <- input$meta
        } else {
          req(input$gene)
          feature_name <- input$gene
        }

        VlnPlot(
          obj_use,
          features = feature_name,
          group.by = input$split_by,
          pt.size = 0
        ) +
          ggtitle(feature_name) +
          theme_classic() +
          theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, color = "black"),
            axis.title = element_text(color = "black", size = 10),
            axis.text.y = element_text(color = "black", size = 10),
            axis.text.x = element_text(color = "black", size = 10),
            axis.ticks.length = unit(4, "pt")
          )
      })
      req(p)
      ggsave(
        filename = file,
        plot = p,
        device = "pdf",
        width = px_to_in(input$violin_width),
        height = px_to_in(input$violin_height),
        units = "in"
      )
    }
  )
}

shinyApp(ui, server)
