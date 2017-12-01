library(shiny)
path2dir='~/Dropbox/TrackingNetworkChanges/'    ##### Needs to be changed accordingly
setwd(path2dir)
source("./tests_synthetic_data/main_tests.R")
change_point_shiny<- function() {
  
  shinyApp(
    
    ui = fluidPage(responsive = FALSE,
                   # Application title
                   titlePanel("Change Point Detection at T=10"),
                   
                   # Sidebar with a slider input for number of bins 
                   sidebarLayout(
                     sidebarPanel(
                       selectInput("opts", label = "Graph type:",
                                   choices = c(0, 1, 2, 3,4), selected = 1),
                       sliderInput("m", label = "Nb of Edges changed at each time step:",
                                   min = 6, max = 10, value = 5, step = 2),
                       sliderInput("m2", label = "Nb of Edges changed at each time step after T=10:",
                                   min = 6, max = 20, value = 10, step = 2),
                       sliderInput("p_disp", label = "Probability of disappearance:",
                                   min = 0.01, max = 0.5, value = 0.1, step = 0.05),
                       sliderInput("p_disp2", label = "Probability of disappearance after T=10:",
                                   min = 0.01, max = 0.5, value = 0.15, step = 0.05),
                       sliderInput("p_creation", label = "Probability of creation:",
                                   min = 0.01, max = 0.2, value = 0.01, step = 0.01),
                       selectInput("order_max", label = "Order of Polynomial:",
                                   choices=c(1,2,3,4,5,6,7),selected=3),
                       sliderInput("alpha", label = "Alpha:",
                                   min = 0.1, max = 1.5, value = 1, step = 0.1),
                       sliderInput("y_max", label = "Plot Height adjustment:",
                                   min = 0.2, max = 100, value = 1, step = 1)
                       
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
      
      data<-reactive({
        N=30
        T=21
        arg_list<-list(power=0.9,islands.n=3,islands.size=10,islands.pin=0.3,n.inter=3,K=3,block.sizes=c(10,10,10),pm=cbind( c(.4,0.1, .001), c(.1,0.2, .01),c(.001,0.01, .5)))
        
        if(input$opts!=0){
          test=test_change_point_realistic(N,input$m,input$m2,input$p_disp,input$p_disp2,input$p_creation,input$opts, T=21,verbose=FALSE,loc=FALSE,compare_LASSO=TRUE,go_plot=FALSE,initial_message=FALSE)
        }
        else{
          test_change_point(N,p=0.4,prop=0.05,prop2=0.2,p2=0.4, T=21,verbose=FALSE,go_plot=FALSE,initial_message=FALSE)
        }
      })
      
      
      
      output$distPlot <- renderPlot({
        dist=data()$dist
        T=21
        plot(1:(T-1), 1-dist[2,], type='l',col=2, xlab="time", ylab="distance", ylim=c(0,input$y_max),main="distances with a change of regime at t=10")
        for (j in 3:7){
          points(1:(T-1), dist[j,], type='l',col=j)
        }
        legend(1,input$y_max,c("jaccard","Spanning Trees","Hamming","IM","HIM","Poly"),col=2:7,lty=1,lwd = 2,cex=0.8)
        
      })
      output$distPlot <- renderPlot({
        heatmap.2(data()$Sim_Lasso,Rowv=FALSE,Colv=FALSE,dendrogram = "none",main="Lasso Similarity")
      })
      
      output$ExPlot <- renderPlot({
        args<-list(p=0.4,power=0.9,islands.n=3,islands.size=10,islands.pin=0.3,n.inter=3,K=3,block.sizes=c(10,10,10),pm=cbind( c(.4,0.1, .001), c(.1,0.2, .01),c(.001,0.01, .5)))
        N=30
        A<-generate_realistic_adjacency(N,args,input$opts, verbose=FALSE)
        plot(A,main="Example of a graph")
      })
      output$summary <- renderTable({
        dist=data()$dist
      })

      
    },
    
    options = list(height = 500)
  )
  
}
