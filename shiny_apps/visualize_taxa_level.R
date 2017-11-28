library(shiny)
library(igraph)
library(plotly)
source("~/Dropbox/Distances/tests_synthetic_data/main_tests.R")
subjects=c(10003,10004,10005,10006,10008,10014,10017,10018,10020 ,10022,10023,10028,10031,10032,10034,10036,10039,10040,10043,10045,10046,10047,10101,19005,19007)
#load(file="~/Dropbox/Food_network/code/notes_extended/PNAS_graph_data/data_post_PNASvis_pb.RData")
load(file="/Users/cdonnat/Distances/Misc/histogram.RData")

visualize_results_taxa<- function() {
  
  shinyApp(
    
    ui = fluidPage(
      # Application title
      titlePanel("Visualization of the results for the PNAS dataset-taxa level "),
      
      # Sidebar with a slider input for number of bins 
      sidebarLayout(
        sidebarPanel(
          #helpText("Select a subject from the PNAS Data"),
          sliderInput("thres", label = "Threshold:",
                      min = 0.01, max = 0.7, value = 0.4, step = 0.05),
          selectInput("distance_type", label = "Distance_type",
                      choices=c("Poly","ST","HIM","Lasso"))
          
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
          plotOutput("treePlot"),
          plotlyOutput("PlotlyPlot"),
          plotOutput("HeatPlot")
        )
      )
    ),
    
    server = function(input, output, session) {
      
      Mydata<-reactive({
        distance_type=input$distance_type
        kept_taxa<-c()
        for (i in 1:length(histogram)){
          if (length(histogram[[i]]>0)){
            kept_taxa<-c(kept_taxa,filter_taxa(histogram[[i]],df,data_otu,input$thres))
          }
          
        }
        kept_taxa<-unique( names(kept_taxa))
        graph_seq=create_graph_per_time_point(histogram, df,data_otu,kept_taxa,input$thres)
        
        is_null=sapply(1:length(graph_seq),FUN=function(x){ return(is.null(graph_seq[[x]]) || sum(graph_seq[[x]])==0)})
        graph_seq_bis<-c()
        for(i in 1:length(is_null)){
          if (!is_null[i]){
            graph_seq_bis<-c(graph_seq_bis, graph_seq[i])
          }
        }
        labels_d=matrix("yellow", length(retained),1)
        labels_d[which(retained>0)]="red"
        distances=compare_distances_per_taxa(graph_seq_bis,distance_type=input$distance_type)
        
        taxa_graph=graph_from_adjacency_matrix(distances, mode = c("undirected"), weighted = TRUE)
        tree=mst(taxa_graph,weights = edge_attr(taxa_graph, "weight"))
        #plot(tree,vertex.color=colors_s_post)
        
        return(list(d=distances,tree=tree,labels_d=labels_d))
        # list(tree=tree,df=df,labels_s_post=labels_s_post,colors_s_post=colors_s_post,labels_s_term=labels_s_term,d=d,dist=dist,subject=s, addendum=addendum)
      })
      
      
      
      output$treePlot <- renderPlot({
        tree=Mydata()$tree
        #df=Mydata()$df
        colors_s_post=Mydata()$labels_d
        plot(tree,vertex.color=colors_s_post,main=paste("Minimum Spanning Tree", sep=" "))
        #plot(tree,vertex.color=colors_s_post)
        legend('topleft',legend=c("Post","Pre birth"),col=c("red","yellow"),pch=19)
      })
      
      output$HeatPlot <- renderPlot({
        d=Mydata()$d
        heatmap.2(as.matrix(d), Rowv=F, Colv=F,dendrogram ="none",main="Heatmap of the distances between Samples")
      })
      
      output$PlotlyPlot <- renderPlotly({
        plot_tsne_representation(Mydata()$d,input$distance_type,Mydata()$labels_d)
        
      })
      
      
    },
    
    options = list(height = 500)
  )
  
}
