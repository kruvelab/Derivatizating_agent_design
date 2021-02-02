#install.packages('reticulate')
library(reticulate)
use_python("/usr/local/bin/python")

# reticulate::py_install("rdkit")
# reticulate::py_install("scipy")

source_python("matrix_to_smiles.py")
np <- import("numpy", convert = FALSE)
rdkit <- import("rdkit", conver = FALSE)


graph_to_SMILES = function(connected_graph) {
  nodes = connected_graph$atoms
  edges = connected_graph$adjacency_matrix
  SMILES = MolFromGraphs(c(t(nodes)), edges)
}