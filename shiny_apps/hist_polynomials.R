library(shiny)
path2dir='~/Dropbox/TrackingNetworkChanges/'    ##### Needs to be changed accordingly
setwd(path2dir)
source("./tests_synthetic_data/main_tests.R")
hist_polynomials <- function() {
  
  shinyApp(
    
    ui = fluidPage(responsive = FALSE,
                   # Application title
                   titlePanel("Histograms for the Polynomial Distance"),
                   
                   # Sidebar with a slider input for number of bins 
                   sidebarLayout(
                     sidebarPanel(
                       selectInput("opts", label = "Graph type:",
                                   choices = c(0, 1, 2, 3,4), selected = 1),
                       selectInput("order_max", label = "Order of Polynomial:",
                                   choices=c(1,2,3,4,5,6,7),selected=3),
                       sliderInput("alpha", label = "Alpha:",
                                   min = 0.1, max = 1.5, value = 1, step = 0.1),
                       sliderInput("y_max", label = "Histogram Height adjustment:",
                                   min = 100, max = 700, value = 100, step = 50),
                       selectInput("B", label = "Bootstrap sample:",
                                   choices = c(50, 100, 500, 1000), selected = 500)
                       
                     ),
                     
                     # Show a plot of the generated distribution
                     mainPanel(
                       plotOutput("distPlot"),
                       plotOutput("ExPlot"),
                       tableOutput("summary")
                     )
                   )
    ),
    
    server = function(input, output, session) {
        
        test_bootstrap1<-reactive({
          args<-list(p=0.4,power=0.9,islands.n=3,islands.size=10,islands.pin=0.3,n.inter=3,K=3,block.sizes=c(10,10,10),pm=cbind( c(.4,0.1, .001), c(.1,0.2, .01),c(.001,0.01, .5)))
          N=30
          B=input$B
          test_bootstrap_shiny(N,input$alpha,input$order_max,args,input$opts,B)
        })
        
        
        
        output$distPlot <- renderPlot({
        hist(test_bootstrap1(),xlim=c(0,1),ylim=c(0,input$y_max),breaks=15,xlab="Distance",main="Histogram of the distance distribution", col = 'darkgray', border = 'white')
        })
        output$ExPlot <- renderPlot({
          args<-list(p=0.4,power=0.9,islands.n=3,islands.size=10,islands.pin=0.3,n.inter=3,K=3,block.sizes=c(10,10,10),pm=cbind( c(.4,0.1, .001), c(.1,0.2, .01),c(.001,0.01, .5)))
          N=30
          B=input$B
          A<-generate_realistic_adjacency(N,args,input$opts, verbose=FALSE)
          plot(A,main="Example of a graph")
        })
        output$summary <- renderTable({
          summary_statistics<-data.frame(c(mean(test_bootstrap1()),sd(test_bootstrap1())),row.names = c("mean","sd"),fix.empty.names = F)
        })
        
      },
    
    options = list(height = 500)
  )
}
