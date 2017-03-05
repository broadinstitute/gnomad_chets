library(shiny)
library(DT)
library(lazyeval)
library(shinythemes)
library(ggflags) # devtools::install_github('baptiste/ggflags')

if (!('shiny_data' %in% ls(globalenv()))) {
  source('../sampleqc.R')
  shiny_data = final_gnomad_meta()
  # shiny_data = european_shiny()
}

library(stringi)
format_label = function(x) {
  if (x %in% platform_names) {
    return(platform_names[[x]])
  } else {
    return(stri_trans_totitle(gsub('_', ' ', x)))
  }
}

ui = function(request) {
  fluidPage(theme = shinytheme('flatly'),
   titlePanel("gnomAD sample QC"),
   
   fluidRow(
     column(3,
      wellPanel(
        selectInput('data',
                    "PCA Data",
                    choices = list('Genotype' = 'genotype',
                                   'Europe' = 'europe',
                                   'Missingness' = 'missingness')),
        selectInput("plot_type",
                    "Plot Type",
                    choices = list('Scatter' = 'point',
                                   'Scatter (with flags!)' = 'flag',
                                   # '2D Density' = 'density_2d',
                                   'Density' = 'density',
                                   'Histogram' = 'histogram'
                                   )),
        radioButtons('source', 'Source data', c('All', 'ExAC', 'gnomAD'), inline = T),
        radioButtons('xtype', 'X data type', c('PC', 'Metric'), inline = T),
        conditionalPanel(
          condition = "input.xtype == 'PC'",
          sliderInput("xpc",
                      "x PC",
                      min = 1,
                      max = 10,
                      value = 1)
        ),
        conditionalPanel(
          condition = "input.xtype == 'Metric'",
          selectInput('xmetric',
                      'x Metric',
                      choices = names(shiny_data)[36:55])
        ),
        conditionalPanel(
          condition = "input.plot_type == 'point'",
          div(
            radioButtons('ytype', 'Y data type', c('PC', 'Metric'), inline = T),
            conditionalPanel(
              condition = "input.ytype == 'PC'",
              sliderInput("ypc",
                          "y PC",
                          min = 1,
                          max = 10,
                          value = 2)
            ),
            conditionalPanel(
              condition = "input.ytype == 'Metric'",
              selectInput('ymetric',
                          'y Metric',
                          choices = c(names(shiny_data)[36:55], 'freemix', 'pct_chimera'))
            )
          )
        ),
        conditionalPanel(
          condition = "input.plot_type == 'histogram'",
          sliderInput("bins",
                      "Bins",
                      min = 5,
                      max = 100,
                      value = 20)
        ),
        sliderInput("probability",
                    "Max population probability",
                    min = 0,
                    max = 100,
                    value = 100),
        conditionalPanel(
          condition = "input.plot_type == 'point'",
          sliderInput("alpha",
                    "Alpha",
                    min = 0,
                    max = 1,
                    value = 0.8)
        ),
        selectInput("color",
                    "Color by:",
                    choices = list('Overall platform' = 'overall_platform',
                                   'Predicted Population' = 'predicted_pop',
                                   'Known Population' = 'known_pop',
                                   'Internal/External' = 'ex_in',
                                   'Platform' = 'platform',
                                   'ExAC/gnomAD' = 'source',
                                   'PI' = 'pi',
                                   'Sex' = 'sex',
                                   'Detailed platform' = 'platform')),
        sliderInput('downsample_level', 'Downsample', min = 0.01, max = 1, value = 0.1),
        actionButton("refresh", "Refresh plot"),
        radioButtons("qc", "QC",
                     choices = c("All", "Release"), inline = T),
        radioButtons("permission", "Permission",
                  choices = c("All", unique(shiny_data$permission)), inline = T),
        radioButtons("external", "Internal or External",
                    choices = c("All", unique(shiny_data$ex_in)), inline = T)
      ),
      wellPanel(
        downloadButton('downloadPlot', 'Download Plot'),
        fluidRow(
          column(6,
                 numericInput('widthExport', 'Width (in)', 7)
          ), column(6,
                    numericInput('heightExport', 'Height (in)', 5)
          )
        ),
        bookmarkButton()
        # themeSelector()
      )
     ),
     column(9,
            fluidRow(
              column(4, offset=8, textOutput('samplesPlotted'))),
            # mainPanel(
             plotOutput("pcaPlot", brush = "plot_brush"),
             verbatimTextOutput('selectedText'),
            tabsetPanel(
              tabPanel("Projects", dataTableOutput('projectTable')),
              tabPanel("Samples", dataTableOutput('sampleTable'))
            )
            # )
     )
   )
)}

server <- shinyServer(function(input, output) {
  
  get_plot_data = reactive({
    print('Getting plot data...')
    out_data = shiny_data
    if (input$source != 'All') {
      out_data %<>% filter(source == input$source)
    }
    out_data$predicted_pop[out_data$probability < (input$probability/100 - 0.01)] = 'oth'
    if (input$data == 'missingness') {
      out_data = out_data %>% select(-c(pc1:pc10))
      names(out_data) = gsub('missingness_', '', names(out_data))
    }
    if (input$data == 'europe') {
      out_data = out_data %>% select(-c(pc1:pc20))
      names(out_data) = gsub('.eur', '', names(out_data))
    }
    if (input$color == 'platform') {
      out_data = subset(out_data, !grepl('|', platform, fixed=T))
    }
    if (input$external != "All") {
      out_data = subset(out_data, ex_in == input$external)
    }
    if (input$qc == 'Release') {
      out_data = subset(out_data, drop_status == 'keep')
    }
    if (input$permission != 'All') {
      out_data = subset(out_data, permission == input$permission)
    }
    out_data
  })

  output$samplesPlotted = renderText({
    data_points = round(nrow(get_plot_data())*input$downsample_level)
    paste0('Plotted: ', data_points, ' (', round(data_points/nrow(shiny_data)*100,1), '%)')
  })
  output$selectedText <- renderText({
    plot_data = get_plot_data()
    text = ''
    if (!is.null(input$plot_brush)) {
      e = input$plot_brush
      brush_text = ''
      if (input$plot_type == 'density') {
        brush_text = paste("Selected on plot: xmin =", round(e$xmin, 2), "xmax =", round(e$xmax, 2))
      } else {
        brush_text = paste("Selected on plot: xmin =", round(e$xmin, 4), "xmax =", round(e$xmax, 4),
             "ymin =", round(e$ymin, 4), "ymax =", round(e$ymax, 4))
      }
      text = paste(text, brush_text, '\n')
    }
    if (input$refresh != 0) {
      selected_data = check_selected_rows()
      if (!is.null(selected_data)) {
        text = paste(text, "Selected from table:", paste(selected_data[[2]], collapse=', '), "\n")
      }
    }
    text
  })
  
  output$projectTable = renderDataTable({
    print('Rendering table data...')
    get_project_table_data()
  })
  
  xmetric = reactive({
    ifelse (input$xtype == "PC", paste0('pc', input$xpc), input$xmetric)
  })
  
  ymetric = reactive({
    ifelse (input$ytype == "PC", paste0('pc', input$ypc), input$ymetric)
  })
  
  output$sampleTable = renderDataTable({
    print('Rendering table data...')
    columns = c("sample", "sample_name_in_vcf", "project_or_cohort",
                "permission", "ex_in", "pi", "description", 
                "gross_platform", "predicted_pop",
                "known_pop",
                xmetric(),
                ymetric()
    )
    get_sample_table_data() %>%
      select_(.dots=columns) 
  })
  
  check_selected_rows = eventReactive(input$refresh, {
    if (!is.null(input$sampleTable_rows_selected)) {
      table_data = get_sample_table_data()
      list('sample', table_data$sample[input$sampleTable_rows_selected])
    } else if (!is.null(input$projectTable_rows_selected)) {
      table_data = get_project_table_data()
      list('project', table_data$project[input$projectTable_rows_selected])
    } else {
      NULL
    }
  })
  
  get_project_table_data = reactive({
    columns = c(~n(), paste0('round(mean(', xmetric(), ', na.rm=T), 3)'), paste0('round(mean(', ymetric(), ', na.rm=T), 3)'))
    column_names = c('n', paste('mean', xmetric()), paste('yean', ymetric()))
    selected_data = get_sample_table_data()
    selected_data %>% group_by(project, ex_in, pi, description, gross_platform) %>% 
      summarize_(.dots = setNames(columns, column_names)) %>%
      arrange(desc(n))
  })
  
  get_sample_table_data = reactive({
    print('Getting table data...')
    selected_data = get_plot_data()
    if (!is.null(input$plot_brush)) {
      selected_data = brushedPoints(selected_data, input$plot_brush, xvar = xmetric(), yvar = ymetric())
      selected_data = subset(selected_data, !is.na(sample))
    }
    selected_data
  })

  
   output$pcaPlot <- renderPlot({
     withProgress(message = "Loading plot...", {generate_plot()})
   })
   
   generate_plot = reactive({
     print('Rendering plot data...')
     plot_data = get_plot_data()
     xpc = xmetric()
     ypc = ymetric()
     set.seed(42)
     plot_data$known_pop[is.na(plot_data$known_pop)] = 'unk'
     if (input$color == 'known_pop' & input$plot_type == 'flag') {
       plot_data %<>% filter(known_pop != 'unk')
     }
     p = ggplot(sample_frac(plot_data, input$downsample_level)) + theme_classic() + theme(text = element_text(size=16))
     if (input$color == 'overall_platform') {
       p = p + scale_color_manual(values=platform_colors, labels=platform_names, guide = guide_legend(title='Overall platform', reverse=TRUE)) + scale_fill_manual(values=platform_colors, labels=platform_names, guide = guide_legend(reverse=TRUE)) + aes(alpha = overall_platform) + scale_alpha_manual(values=platform_alphas)
     #} else if (input$color == 'known_pop' || input$color == 'predicted_pop') { # Commented for now to do European work
     } else if (input$color == 'predicted_pop') {
       p = p + scale_color_manual(values=pop_colors, labels=pop_names, guide = guide_legend(title='Population')) + scale_fill_manual(values=pop_colors, labels=pop_names, guide = guide_legend(title='Population'))
     }
     if (input$color == 'known_pop') {
       alphas = rep(0.9, length(unique(plot_data$known_pop)))
       names(alphas) = unique(plot_data$known_pop)
       alphas['unk'] = 0.1
       p = p + aes(alpha = known_pop) + scale_alpha_manual(values=alphas, guide = F)
       p = p + scale_color_manual(values=pop_colors, labels=pop_names, guide = guide_legend(title='Population'))
       p = p + scale_fill_manual(values=pop_colors, labels=pop_names, guide = guide_legend(title='Population'))
     }
     manual_alpha = input$color == 'known_pop'
     if (input$plot_type == 'flag') {
       p + aes_string(x = xpc, y = ypc) + aes(country=known_pop) + geom_flag() + scale_country(labels=sub_pop_names, guide=guide_legend(title='Population'))
     } else if (input$plot_type == 'point') {
       p = p + aes_string(x = xpc, y = ypc, col = input$color)
       if (input$refresh == 0) {
         # p + geom_point(alpha=input$alpha)
         if (manual_alpha) {
           p + geom_point()
         } else {
           p + geom_point(alpha=input$alpha)
         }
       } else {
         selected_data = check_selected_rows()
         if (is.null(selected_data)) {
           if (manual_alpha) {
             p + geom_point()
           } else {
             p + geom_point(alpha=input$alpha)
           }
         } else {
           type = selected_data[[1]]
           selected_data = selected_data[[2]]
           p + geom_point(alpha=0.2) +
             geom_point(aes_string(shape = type), filter_(plot_data, .dots=interp(~type %in% selected_data, type = as.name(type))), alpha=1, col='darkblue')
         }
       }
     } else if (input$plot_type == 'density') {
       if (input$refresh == 0) {
         p + geom_density(alpha=0.5) + aes_string(x = xpc, col = input$color, fill = input$color)
       } else {
         selected_data = check_selected_rows()
         if (is.null(selected_data)) {
           p + geom_density(alpha=0.5) + aes_string(x = xpc, col = input$color, fill = input$color)
         } else {
           type = selected_data[[1]]
           selected_data = selected_data[[2]]
           p + 
             geom_density(data=filter_(plot_data, .dots=interp(~!type %in% selected_data, type = as.name(type))), alpha=0.5) + 
             geom_density(data=filter_(plot_data, .dots=interp(~type %in% selected_data, type = as.name(type))), alpha=0.8, color='darkblue', fill='darkblue') +
             aes_string(x = xpc, col = input$color, fill = input$color)
         }
       }
     } else if (input$plot_type == 'histogram') {
       p + geom_histogram(position = "dodge", bins = input$bins) + aes_string(x = xpc, col = input$color, fill = input$color)
     } else if (input$plot_type == '2d_density') { # Not working yet
       clear = rgb(1, 1, 1, 0)
       smoothScatter(plot_data[,xpc], plot_data[,ypc], colramp=colorRampPalette(c(clear, clear)), col='white')
       l_ply(rev(names(platform_colors)), function(x) {
        color_ramp = colorRampPalette(c(clear, paste0(platform_colors[[x]], 'CC')), bias=1000, alpha = T)
        smoothScatter(subset(plot_data, overall_platform == x)[,xpc], subset(plot_data, overall_platform == x)[,ypc], colramp = color_ramp, col=platform_colors[[x]], add=T)
       })
     }
   })
   output$downloadPlot <- downloadHandler(
     filename = function() { paste0(input$data, '.png') },
     content = function(file) {
       ggsave(file, plot = generate_plot(), device = "png", height = as.numeric(input$heightExport), width = as.numeric(input$widthExport), units='in')
     }
   )
})

shinyApp(ui = ui, server = server, enableBookmarking = "server")