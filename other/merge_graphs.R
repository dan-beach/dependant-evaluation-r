
kidney <- read.graph('~/phd/data/slant_cancer/results/graph_plots/CAKI2_KIDNEY_graph_changing_weights.gml', format = 'gml')
breast <- read.graph('~/phd/data/slant_cancer/results/graph_plots/AU565_BREAST_graph_changing_weights.gml', format = 'gml')

V(kidney)$dependent
V(kidney)$dependent == V(breast)$dependent
sum(V(kidney)$dependent != V(breast)$dependent)

V(kidney)$b_dependent <- V(breast)$dependent

V(kidney)$diff <- V(kidney)$dependent == V(breast)$dependent

write.graph(kidney, '~/phd/data/slant_cancer/results/graph_plots/representative_graph.gml', format = 'gml')

