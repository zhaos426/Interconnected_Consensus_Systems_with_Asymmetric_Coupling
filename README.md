# Interconnected_Systems_with_Aysmmetric_Coupling
We provide here code to accompany paper: Asymmetric Coupling of Directed Networks Accelerates Diffusion and Consensus

The following python files are included:
```
network_properties.py -- contains code to initialize the networks in the paper and generate the basic properties of the networks, such as the number of nodes and layers, intralayer networkx graph, interlayer edge list, etc.
  Laplacians.py -- contains code to get Laplacians of the networks.
	eigenvalues.py -- contains code to generate the exact or predicted eigenvalues of the Laplacian matrix.
	strong_limit.py -- contains code to construct the strong limit theory used by predicting the eigenvalues of the Laplacian matrix.
	plot.py -- contains code to plot figures in the paper and supplement.
```


The following Jupyter notebooks are included:
```
	demo.ipny -- simple demo code to reproduce the findings in the paper.
	figure_1_and_2.ipynb -- contains code to reproduce figs 1 and 2 in the paper.
	figure_3_and_4.ipynb -- contains code to reproduce figs 3 and 4 in the paper.	
	supplement_fig_1.ipynb -- contains code to illustrate how the non-robust optima arises due to a bifurcation in which two eigenvalues collide to create a complex pair of eigenvalues.
	supplement_fig_2.ipynb -- contains code to visualize the social network in the paper.
```
