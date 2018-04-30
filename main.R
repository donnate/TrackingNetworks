rm(list=ls()) ### clear environment for fresh start

dir=paste('~/Dropbox/TrackingNetworkChanges',sep='')  ### to be modified appropriately
setwd(dir)
library(RColorBrewer)
library(gplots)

source("./tests_synthetic_data/change_point_type_graph.R")
source("./tests_synthetic_data/tests_bootstrap.R")
source("./tests_synthetic_data/evolving_graphs.R")



source("./tests_synthetic_data/dynamics.R")
source("./distances/spanning_trees.R")
source("./distances/distances.R")

