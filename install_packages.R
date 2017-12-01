list.of.packages <- c("glmnet", "sna", "ergm", "fpc", "ggplot2", "Rcpp","glmnet","gplots","igraph","LCA","Matrix","nettools",
                      "PBSmodelling","tsne","tictoc","plotly")


installed=installed.packages()
list.of.packages_bioclite=c("minerva","minet","dtw","WGCNA","impute","sincell")
for (pck in list.of.packages_bioclite){
  if (!(pck %in% installed )){
    source("https://bioconductor.org/biocLite.R")
    biocLite(pck)
    #install.packages(pck, repos='https://bioconductor.org/biocLite.R',quiet = TRUE)
  }
  require(pck, character.only=TRUE,quietly = TRUE)
}

for (pck in list.of.packages){
  if (!(pck %in% installed )){
    install.packages(pck, repos='http://cran.us.r-project.org',quiet = TRUE)
  }
  require(pck, character.only=TRUE,quietly = TRUE)
}




