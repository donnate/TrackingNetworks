# Distances


This folder gathers the R implementation of functions computing distances between two instances of graphs on the same (identified) nodes. The heatwave-based distances is the only one not contained in this folder, and was implemented with python -- for convenience purposes.

The different files are as follows
+ __compare_distances.R__ contains a function computing pairwise distances (for a few different distances) on a given sequence of graphs. This allows to compare the performance of the different distances (stability, sensitivity, ability to recognize a certain class of graphs, change point detection, etc.).
+ __distances.T__ contains the R implementation of most distances.
+ __spanning_trees.R__ contains the R implementation of the spanning tree distance.
