# TrackingNetworks: Tests on Synthetic Data


This folder contains all the code necessary to reprodice te synthetic experiments presented in the paper.
The folder is organized as follows:
+ __main_test.R__ contains
+ __test_classes_graphs.R__ contains
+ __test_edge_deletion.R__ contains
+ __dynamics.R__ contains all the functions for generating  random graphs and perturbations (to mimick an evolution process). The evolution process used here always follows the same guidelines:  a perturbed version of a given graph is created by  selecting a proportion (or number) of edges as "potentially modifiable edges". These edges are then either (a) rewired with probability p or (b) deleted with probability p_disp, so that p+p_disp<=1 We also allow for the possibility of "edge creation": each "non-existing" edge can be added to the new perturbed graph with probability p_creation.
+ __test_bootstrap.R__ contains functions for generating distributions of the distances between two graphs of the same type. This is useful to undertand the stability (spread, mean, etc) of a given distance and thus its potential sensitivity to change points.
+ __generated_graphs__ is a folder where the graph sequences produced in every experiment can be stored
+ __saved_data__ is another folder where temporary data can be stored (R session, etc.)
