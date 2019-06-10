rm(list = ls())
library(plotly)

load('output/pred_surfaces.rda')

make_3d_plot <- function(surface) { 
  
  p <- plot_ly(z = ~surface ) %>% 
    add_surface(
    contours = list(
      z = list(
        show=TRUE,
        usecolormap=TRUE,
        highlightcolor="#ff0000",
        project=list(z=TRUE)
      )
    )
  ) %>%
    layout(
      scene = list(
        camera=list(
          eye = list(x=1.87, y=0.88, z=-0.64)
        )
      )
    )
  return( p )
}

p1 <- make_3d_plot(pred_surface1)
p2 <- make_3d_plot(pred_surface2)
p3 <- make_3d_plot(pred_surface3)

p1 
p2
p3
