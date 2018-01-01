# TrackingNetworks: Tests on Synthetic Data


This folder contains all the code necessary to reprodice te synthetic experiments presented in the paper.
The folder is organized as follows:
+ __dynamics.R__ contains all the functions for generating  random graphs and perturbations (to mimick an evolution process). The evolution process used here always follows the same guidelines:  a perturbed version of a given graph is created by  selecting a proportion (or number) of edges as "potentially modifiable edges". These edges are then either (a) rewired with probability p or (b) deleted with probability p_disp, so that p+p_disp<=1 We also allow for the possibility of "edge creation": each "non-existing" edge can be added to the new perturbed graph with probability p_creation.
+ __evolving_graphs.R__ contains functions fro tracking the evolution of a graph through time.
+ __change_point_type_graph.R__ contains functions for comparing independent graphs of the same type (same generation process). The goal is to assess if a change point in the generation mechanism is picked up by the different distances.
+ __test_edge_deletion.R__ contains
+ __test_bootstrap.R__ contains functions for generating distributions of the distances between two graphs of the same type. This is useful to undertand the stability (spread, mean, etc) of a given distance and thus its potential sensitivity to change points.
+ __generated_graphs__ is a folder where the graph sequences produced in every experiment can be stored
+ __saved_data__ is another folder where temporary data can be stored (R session, etc.)
