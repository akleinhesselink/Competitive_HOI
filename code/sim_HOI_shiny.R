library(shiny)

ui <- fluidPage(
  sidebarLayout( 
    sidebarPanel(
      sliderInput('soil', 'soil moisture', 0, 300, 200),
      br(),
      sliderInput('epsilon', 'soil drying rate', 0, 1, 0.01),
      br(),
      actionButton("run", "run simulation")
    ), 
    mainPanel(
      plotOutput("plot")
    )
  )
)

server <- function(input, output) {
  library(dplyr)
  library(tidyr)
  library(deSolve)
  library(ggplot2)
  library(stringr)
  
  source('sim_functions.R')
  source('figure_pars.R')
  # parameterize model --------------------------------------------------------------------------------------------------- 
  tiny <- .Machine$double.eps
  times <- 125             # length of simulation in days 
  #soil_m <- input$soil          # initial soil moisture (mm water in upper 500 mm of soil)
  pulse <- 0               # amount of water supplied per day in mm  
  rainy <- 10             # duration of rainy period 
  r <- c(4.2, 2.9, 2.3) # max uptake rates mm of water per g of plant per day
  K <- c(110, 30, 0.5)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
  m <- 0.09                # tissue respiration and loss rate g per g per day 
  q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
  epsilon <- 0.01         # rate of water evaporation and runoff mm per mm per day
  seedlings <- c(1, 0, 0)      # number of seedlings 
  seedling_mass <- c(0.005) # seed/seedling mass in g 
  conversion <- 0.1        # proportion live biomass converted to seed mass 
  R <- seq(0, 500, length.out = 1000)
  
  seeds <- c(1,1,1)
  #State <- c(soil_m, seeds*seedling_mass)
  parms <- list( r = r, K = K, m = m, p = c(rep(pulse, rainy), rep(0, times - rainy)), epsilon = epsilon, q = q)
  
  epsilon <- reactive({ 
      input$epsilon
  })
   
  State <- reactive({ 
    c(input$soil, seeds*seedling_mass)
  })

  dat <- eventReactive(input$run, {
    State <- State()
    parms$epsilon <- epsilon()
    ode(y=State, times = seq( 1, times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root )
  })

  output$plot <- renderPlot({
    plot_timeseries(dat(), parms, col = my_colors)
  })
}

shinyApp(ui = ui, server = server)


