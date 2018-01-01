# TrackingNetworks


This provides the code for implementing different distances between a set of aligned graphs (that is, different instances on the same set of nodes, where each node is identified). This folder also provides code for testing and comparing the performance of these different distances in a variety of situations:
+ __Synthetic experiments__: we test the stability of the different distances and their ability to (a) recognize graphs belonging to the same class  and (b) detect change points in the evolution process or generation process.
+ __Real-life cases__: different distances were used to analyze real-life data: a microbiome study (the 2011 Relman antibiotics study), as well as a recipe network


The list of all required packages can be (found and ) downloaded by running __install_packages.R__.
Once the required packages have been installed, in order to load all the functions and firectly start working, simply run __main.R__.

Update: January 1st, 2018: code has been cleaned up in order to be more user friendly, and we are currently working on providing an Rmarkdown folder detailing how to use the different functions and displaying all the experiments shown in the paper.
