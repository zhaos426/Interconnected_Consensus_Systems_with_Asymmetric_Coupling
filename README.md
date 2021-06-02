# Interconnected_Systems_with_Aysmmetric_Coupling
We provide here code to accompany paper: Asymmetric Coupling of Directed Networks Accelerates Diffusion and Consensus


The following Jupyter notebooks are included:
```
demo.ipny -- simple demo code to reproduce the findings in the paper.
figure_1_and_2.ipynb -- contains code to reproduce figs 1 and 2 in the paper.
figure_3_and_4.ipynb -- contains code to reproduce figs 3 and 4 in the paper.	
supplement_fig_1.ipynb -- contains code to illustrate how the non-robust optima arises due to a bifurcation in which two eigenvalues collide to create a complex pair of eigenvalues.
supplement_fig_2.ipynb -- contains code to visualize the social network in the paper.
```

The following python files are included:
```
network_properties.py -- contains code to initialize the networks in the paper and generate the basic properties of the networks, such as the number of nodes and layers, intralayer networkx graph, interlayer edge list, etc.
Laplacians.py -- contains code to get Laplacians of the networks.
eigenvalues.py -- contains code to generate the exact or predicted eigenvalues of the Laplacian matrix.
strong_limit.py -- contains code to construct the strong limit theory used by predicting the eigenvalues of the Laplacian matrix.
plot.py -- contains code to plot figures in the paper and supplement.
```

The following real world network data file is included:
```
Krackhardt-High-Tech_Multiplex_Social:
	The multiplex social network consists of 3 kinds of relationships (Advice, Friendship and "Reports to") between managers of a high-tech company. In our paper, we only consider Advice layer.
	There are 21 nodes in each layer, labelled with integer ID between 1 and 21, with 312 connections.
	The multiplex is directed and unweighted, stored as edges list in the edges file.
```

### Reference
  [1] Taylor, D., Myers, S. A., Clauset, A., Porter, M. A., & Mucha, P. J. (2017). Eigenvector-based centrality measures for temporal networks. Multiscale Modeling & Simulation, 15(1), 537-574.
  
  [2] Gomez, S., Diaz-Guilera, A., Gomez-Gardenes, J., Perez-Vicente, C. J., Moreno, Y., & Arenas, A. (2013). Diffusion dynamics on multiplex networks. Physical review letters, 110(2), 028701.
  
  [3] Tejedor, A., Longjas, A., Foufoula-Georgiou, E., Georgiou, T. T., & Moreno, Y. (2018). Diffusion dynamics and optimal coupling in multiplex networks with directed layers. Physical Review X, 8(3), 031071.
  
  [4] D. Krackhardt, Social Networks9, 109 (1987).
