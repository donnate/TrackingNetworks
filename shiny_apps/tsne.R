library(shiny)
library(igraph)
source("~/Dropbox/Distances/tests_synthetic_data/main_tests.R")
subjects=c(10003,10004,10005,10006,10008,10014,10017,10018,10020 ,10022,10023,10028,10031,10032,10034,10036,10039,10040,10043,10045,10046,10047,10101,19005,19007)
#load(file="~/Dropbox/Food_network/code/notes_extended/PNAS_graph_data/data_post_PNASvis_pb.RData")

library(tsne)
library(plotly)
#load(paste("~/Dropbox/Food_network/code/notes_extended/PNAS_graph_data/graph_seq_complete_no_df",input$subj,".RData",sep=""))

test_tsne=FALSE
if (test_tsne==TRUE){
  tsne1=tsne(res$d, initial_config = NULL, k = 3, initial_dims = 30, perplexity = 100, max_iter = 1000, min_cost = 0, epoch_callback = NULL, whiten = TRUE,epoch=100)
  plot(tsne1,pch=19, t='n')
  text(tsne1, labels=1:ncol(res$d),col="black")
  labels=sapply(res$colors_s_post, FUN=function(x){ ifelse(x=="yellow", "Pre-birth sample", "Post Birth Sample")})
  p <- plot_ly(data.frame(tsne1), x = ~tsne1[,1], y = ~tsne1[,2], z = ~tsne1[,3],mode = 'text',text = 1:ncol(res$d),color= labels,colors= c('#BF382A', '#0C4B8E'),
               textfont = list(color = '#000000', size = 16)) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'tsne 1st axis'),
                        yaxis = list(title = 'tsne 2nd axis'),
                        zaxis = list(title = 'tsne 3rd axis')))
  show(p)
}



plot_tsne_representation<-function(distances,distance_type,labels_d){
  tsne1=tsne(distances, initial_config = NULL, k = 3, initial_dims = 30, perplexity = 25, max_iter = 1000, min_cost = 0, epoch_callback = NULL, whiten = TRUE,epoch=100)
  plot(tsne1,pch=19, t='n')
  text(tsne1, labels=1:ncol(distances),col=labels_d)
  labels=sapply(labels_d, FUN=function(x){ ifelse(x=="yellow", "Pre-birth sample", "Post Birth Sample")})
  p <- plot_ly(data.frame(tsne1), x = ~tsne1[,1], y = ~tsne1[,2], z = ~tsne1[,3],mode = 'text',text = 1:ncol(distances),color= labels,colors= c('#BF382A', '#0C4B8E'),
               textfont = list(color = '#000000', size = 16)) %>%
    add_markers() %>%
    layout(title = distance_type,scene = list(xaxis = list(title = 'tsne 1st axis'),
                        yaxis = list(title = 'tsne 2nd axis'),
                        zaxis = list(title = 'tsne 3rd axis')))
  show(p)
}
