library(shiny)

# See above for the definitions of ui and server
ui <- fluidPage(
  
  # App title ----
  titlePanel("Hello Shiny!"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "alpha",
                  label = "alpha",
                  min = 0,
                  max = 1,
                  value = 0.1, step = 0.01),
      
      sliderInput(inputId = "beta",
                  label = "beta",
                  min = 0,
                  max = 10,
                  value = 6), 
      
      sliderInput(inputId = "gamma",
                  label = "gamma",
                  min = -1,
                  max = 1,
                  value = 0, step = 0.01), 
      
      sliderInput(inputId = "delta",
                  label = "delta",
                  min = 0,
                  max = 10,
                  value = 5, step = 0.01), 
      
      sliderInput(inputId = "cc",
                  label = "cc",
                  min = 0,
                  max = 100,
                  value = 5, step = 1)
    )
    
    ,
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "distPlot")
      
    )
  )
)
server <- server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  
  # parameterize model --------------------------------------------------------------------------------------------------- 
  source('~/Dropbox/projects/Competitive_HOI/code/sim_functions.R')
  
  output$distPlot <- renderPlot({

    alpha <- 0.1            # factors determining trade-off curve 
    beta <- 6
    gamma <- 0
    delta <- 5
    cc <- 5
    TO_pars <- list( alpha = alpha, beta = beta , gamma  = gamma, delta = delta, cc = cc)
    
    curve(TO_fun(x, TO_pars), 0, 10)
        
    alpha <- input$alpha            # factors determining trade-off curve 
    beta <- input$beta
    gamma <- input$gamma
    delta <- input$delta
    cc <- input$cc
    TO_pars <- list( alpha = alpha, beta = beta , gamma  = gamma, delta = delta, cc = cc, col = 'red')
    
    curve(TO_fun(x, TO_pars), add = T, 0, 10)
    
  })
  
}

shinyApp(ui = ui, server = server)