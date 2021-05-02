library(tidyverse)
library(rcdk)
library(igraph)

connect_substructures = function(SMILES1, SMILES2) {
  graph1_list = SMILES_to_graph(SMILES1)
  graph2_list = SMILES_to_graph(SMILES2)
  graph1 = set.vertex.attribute(graph1_list$graph, "name", value=as.character(1:length(graph1_list$atoms)))
  #while numbering the second graph we need to start from where the first graph finishes
  graph2 = set.vertex.attribute(graph2_list$graph, "name", 
                                value = as.character((length(graph1_list$atoms)+1):(length(graph1_list$atoms)+length(graph2_list$atoms))))
  union_graph = igraph::union(graph1, graph2)
  union_graph_matrix = as.matrix(as_adjacency_matrix(union_graph))
  union_graph_atoms = bind_rows(graph1_list$names_of_atoms, graph2_list$names_of_atoms)
  union_graph = graph_from_adjacency_matrix(union_graph_matrix)
  
  available_connections_first_SMILES = atoms_available_for_connection(graph1_list$names_of_atoms, 
                                                                      graph1_list$matrix_with_bond_count)
  available_connections_second_SMILES = atoms_available_for_connection(graph2_list$names_of_atoms, 
                                                                       graph2_list$matrix_with_bond_count)
  
  connections = c(sample(available_connections_first_SMILES, size = 1),
                  sample(available_connections_second_SMILES, size = 1)+length(graph1_list$atoms))
  union_graph = union_graph %>%
    add.edges(connections)
  union_graph_matrix = as.matrix(as_adjacency_matrix(union_graph))
  combined_molecule = list("graph" = union_graph, 
                           "adjacency_matrix" = union_graph_matrix, 
                           "atoms" = union_graph_atoms)
  return(combined_molecule)
}

#for adding amino acids we need a modification as only one position is available for connections in each amino acid
connect_substructures_AH = function(SMILES1, SMILES2) {
  graph1_list = SMILES_to_graph(SMILES1)
  graph2_list = SMILES_to_graph(SMILES2)
  graph1 = set.vertex.attribute(graph1_list$graph, "name", value=as.character(1:length(graph1_list$atoms)))
  #while numbering the second graph we need to start from where the first graph finishes
  graph2 = set.vertex.attribute(graph2_list$graph, "name", 
                                value = as.character((length(graph1_list$atoms)+1):(length(graph1_list$atoms)+length(graph2_list$atoms))))
  union_graph = igraph::union(graph1, graph2)
  union_graph_matrix = as.matrix(as_adjacency_matrix(union_graph))
  union_graph_atoms = bind_rows(graph1_list$names_of_atoms, graph2_list$names_of_atoms)
  union_graph = graph_from_adjacency_matrix(union_graph_matrix)
  
  available_connections_first_SMILES = atoms_available_for_connection(graph1_list$names_of_atoms, 
                                                                      graph1_list$matrix_with_bond_count)
  connections = c(sample(available_connections_first_SMILES, size = 1),
                  length(graph1_list$atoms)+1)
  union_graph = union_graph %>%
    add.edges(connections)
  union_graph_matrix = as.matrix(as_adjacency_matrix(union_graph))
  combined_molecule = list("graph" = union_graph, 
                           "adjacency_matrix" = union_graph_matrix, 
                           "atoms" = union_graph_atoms)
  return(combined_molecule)
}


SMILES_to_graph = function(SMILES) {
  m <- parse.smiles(SMILES)[[1]]
  matrix_with_bond_count <- get.connection.matrix(m)
  #matrix <- get.adjacency.matrix(m)
  atoms <- get.atoms(m)
  names_of_atoms <- tibble(atoms = character())
  for (i in 1:length(atoms)) {
    names_of_atoms <- names_of_atoms %>%
      add_row(atoms = get.symbol(atoms[[i]]))
  }
  graph = igraph::graph_from_adjacency_matrix(matrix_with_bond_count)
  graph = list("graph" = graph, 
               "atoms" = atoms, 
               "matrix_with_bond_count" = matrix_with_bond_count,
               "names_of_atoms" = names_of_atoms)
  return(graph)
}

atoms_available_for_connection = function(names_of_atoms, matrix_with_bond_count) {
  names_of_atoms = names_of_atoms %>%
    add_column(bond_count = colSums(matrix_with_bond_count))
  names_of_atoms = names_of_atoms %>%
    mutate(is_available_for_connection = case_when(
      atoms == "O" & bond_count < 2 ~ 1,
      atoms == "C" & bond_count < 4 ~ 1,
      atoms == "N" & bond_count < 3 ~ 1,
      atoms == "F" & bond_count < 1 ~ 1,
      atoms == "Cl" & bond_count < 1 ~ 1,
      atoms == "Br" & bond_count < 1 ~ 1,
      atoms == "I" & bond_count < 1 ~ 1,
      atoms == "S" & bond_count < 6 ~ 1,
      TRUE ~ 0
    ))
  return(which(names_of_atoms$is_available_for_connection == 1))
}


#base structure----
# SMILES1 <- "c1cc(O)ccc1(O)"
# SMILES2 <- "C(=O)O"
# test = connect_substructures(SMILES1, SMILES2)
# test$graph
# plot(test$graph)
