list.of.packages <- c("glmnet", "sna", "ergm", "fpc", "ggplot2", "Rcpp","glmnet","gplots","igraph","LCA","Matrix","nettools",
    "PBSmodelling","tsne","tictoc","plotly","infotheo","RSQLite",
    "Kendall","MASS","XML","PBSmodelling")


installed=installed.packages()
for (pck in list.of.packages){
    if (!(pck %in% installed )){
        install.packages(pck, repos='http://cran.us.r-project.org',quiet = TRUE)
    }
    require(pck, character.only=TRUE,quietly = TRUE)
}



list.of.packages_bioclite=c("minerva","minet","dtw","WGCNA","impute","sincell")
options(useHTTPS=FALSE)
install.packages("BiocInstaller", repos="http://bioconductor.org/packages/3.3/bioc")
for (pck in list.of.packages_bioclite){
  if (!(pck %in% installed )){
    source("https://bioconductor.org/biocLite.R")
    biocLite(pck)
    #install.packages(pck, repos='https://bioconductor.org/biocLite.R',quiet = TRUE)
  }
  require(pck, character.only=TRUE,quietly = TRUE)
}





