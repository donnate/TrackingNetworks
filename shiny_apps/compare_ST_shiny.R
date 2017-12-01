library(shiny)
path2dir='~/Dropbox/TrackingNetworkChanges/'    ##### Needs to be changed accordingly
setwd(path2dir)
source('./tests_synthetic_data/main_tests.R")
compare_spanning_trees_shiny<- function(seq_graphs) {
  
  shinyApp(
    
    ui = fluidPage(responsive = FALSE,
                   # Application title
                   titlePanel("Compare the different  Minimum Spanning Trres"),
                   
                   # Sidebar with a slider input for number of bins 
                   sidebarLayout(
                     
                     # Show a plot of the generated distribution
                     mainPanel(
                       plotOutput("distPlot1"),
                     )
                   )
    ),
    
    server = function(input, output, session) {
      
     
      output$distPlot1 <- renderPlot({
        dist=compare_spanning_trees(seq_graphs, go_plot=FALSE)
        plot(dist$mst1, vertex.size=5,main="Spanning Tree from Poly")
        plot(dist$mst2, vertex.size=5,main="Spanning Tree from Hamming")
        plot(dist$mst3, vertex.size=5,main="Spanning Tree from IM")
        plot(dist$mst4, vertex.size=5,main="Spanning Tree from HIM")
        plot(dist$mst5, vertex.size=5,main="Spanning Tree from ST Sim")
        plot(dist$mst6, vertex.size=5,main="Spanning Tree from Lasso Sim")
      })
      
      
      
    },
    
    options = list(height = 500)
  )
  
}
