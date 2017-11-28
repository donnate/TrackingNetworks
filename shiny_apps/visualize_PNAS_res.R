library(shiny)
library(igraph)
library(plotly)
source("~/Dropbox/Distances/tests_synthetic_data/main_tests.R")
subjects=c(10003,10004,10005,10006,10008,10014,10017,10018,10020 ,10022,10023,10028,10031,10032,10034,10036,10039,10040,10043,10045,10046,10047,10101,19005,19007)


visualize_results_PNAS<- function() {
  
  shinyApp(
    
    ui = fluidPage(
                   # Application title
                   titlePanel("Visualization of the results for the PNAS dataset "),
                   
                   # Sidebar with a slider input for number of bins 
                   sidebarLayout(
                     sidebarPanel(
                       #helpText("Select a subject from the PNAS Data"),
                       selectInput("subj", label = "Subject id:",
                                   choices=1:24)
                       
                     ),
                     
                     # Show a plot of the generated distribution
                     mainPanel(
                       plotOutput("treePlot"),
                       plotlyOutput("PlotlyPlot"),
                       plotOutput("HeatPlot"),
                       plotOutput("DistPlot"),
                       verbatimTextOutput("summary")
                     )
                   )
    ),
    
    server = function(input, output, session) {
      
      Mydata<-reactive({
        load(paste("~/Dropbox/Distances/Misc/PNAS_graph_data/graph_seq_complete_no_df",input$subj,".RData",sep=""))
        #load(file="~/Dropbox/Food_network/code/notes_extended/PNAS_graph_data/data_post_PNASvis_pb.RData")
        s=subjects[input$subj]
        #index_s<-df$SampleID[which(df$SubjectID ==s )]
        #       labels_s_post=df$PostPreg[which(df$SubjectID ==s )]
        #labels_s_post<-sapply(labels_s_post, FUN=function(x){
        # ifelse(x,1,2)
        #})
        #y=which(df$SubjectID ==s )
        #xx=y[1]
        #addendum=ifelse(df$Term[xx], "(Term Pregnancy)", ifelse(df$Preterm[xx], "(Preterm Pregnancy)", "(Marginal Pregnancy)"))
        #labels_s_term=ifelse(df$Term[xx], 3, ifelse(df$Preterm[xx], 1, 2))
        #subject_graph=graph_from_adjacency_matrix(d, mode = c("undirected"), weighted = TRUE)
        #tree=mst(subject_graph,weights = edge_attr(subject_graph, "weight"))
        #dist=triu(d,1)[which(triu(d,1)!=0)]
        #colors_s_post=sapply(which(df$SubjectID==s),FUN=function(x){
        # ifelse(df$Post[x], "red","yellow")
        #})
        tsne1=tsne(res$d, initial_config = NULL, k = 3, initial_dims = 30, perplexity = 100, max_iter = 1000, min_cost = 0, epoch_callback = NULL, whiten = TRUE,epoch=100)
        labels=sapply(res$colors_s_post, FUN=function(x){ ifelse(x=="yellow", "Pre-birth sample", "Post Birth Sample")})
        return(list(addendum=res$addendum,tree=res$tree,labels_s_post=res$labels_s_post,d=res$d,dist=res$dist,subject=res$subject,colors_s_post=res$colors_s_post,labels=labels,tsneD=tsne1))
        # list(tree=tree,df=df,labels_s_post=labels_s_post,colors_s_post=colors_s_post,labels_s_term=labels_s_term,d=d,dist=dist,subject=s, addendum=addendum)
      })
      
      
      
      output$treePlot <- renderPlot({
        tree=Mydata()$tree
        #df=Mydata()$df
        s=Mydata()$subject
        colors_s_post=Mydata()$colors_s_post
        plot(tree,vertex.color=colors_s_post,main=paste("Minimum Spanning Tree for subject", Mydata()$subject, Mydata()$addendum, sep=" "))
        #plot(tree,vertex.color=colors_s_post)
        legend('topleft',legend=c("Post","Pre birth"),col=c("red","yellow"),pch=19)
      })
      
      output$HeatPlot <- renderPlot({
        d=Mydata()$d
        heatmap.2(as.matrix(d), Rowv=F, Colv=F,dendrogram ="none",main="Heatmap of the distances between Samples")
      })
      
      output$PlotlyPlot <- renderPlotly({
        d <- Mydata()$d
        labels=Mydata()$labels
        tsne1=Mydata()$tsneD
        text(tsne1, labels=1:ncol(d),col="black")
        plot_ly(data.frame(tsne1), x = ~tsne1[,1], y = ~tsne1[,2], z = ~tsne1[,3],mode = 'text',text = 1:nrow(tsne1),color= labels,colors= c('#BF382A', '#0C4B8E'),
                     textfont = list(color = '#000000', size = 16)) %>%
          add_markers() %>%
          layout(scene = list(xaxis = list(title = 'tsne 1st axis'),
                              yaxis = list(title = 'tsne 2nd axis'),
                              zaxis = list(title = 'tsne 3rd axis')))
      
      })
      
      output$DistPlot <- renderPlot({
        
          hist(Mydata()$dist, main="Histogram of the distances")
      })
      output$summary <- renderPrint({
        dist <- Mydata()$dist
        summary(dist)
      })
      
      
    },
    
    options = list(height = 500)
  )
  
}
