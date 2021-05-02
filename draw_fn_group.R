source("connecting_graphs.R")
source("graph_to_SMILES.R")
source("add_AH.R")


list_of_fn_groups = read_delim("SMILES_of_fingerprints_corrected.csv",
                               delim = ",",
                               col_names = TRUE)

list_of_fn_groups = list_of_fn_groups %>%
  filter(str_count(SMILES, pattern = "#") == 0) %>%
  mutate(probability = 1/dim(list_of_fn_groups)[1])

fn_groups = list_of_fn_groups$SMILES

new_structures = tibble(
  starting_SMILES = character(),
  add_SMILES = character(),
  SMILES = character(),
  this_depth = character()
)


add_fn_group = function(start_SMILES, depth, breadth, df) {
  print(start_SMILES)
  if (depth > 0) {
    df = make_breadth(start_SMILES, depth, breadth, df)
    return(df)
  }
  else {
    return(df)
  }
  return(df)
}

make_breadth = function(start_SMILES, depth, breadth, df) {
  for (i in 1:breadth) {
    n = sample(1:length(fn_groups), size = 1)
    fn_group = fn_groups[n]
    print(fn_group)
    new_structure = connect_substructures(start_SMILES, fn_group)
    new_structure = graph_to_SMILES(new_structure)
    df = df %>%
      bind_rows(
        c(starting_SMILES = start_SMILES,
          add_SMILES = fn_group,
          SMILES = new_structure,
          this_depth = depth))
    df = add_fn_group(new_structure, depth - 1, breadth, df)
  }
  return(df)
}

organize_df = function(df) {
  j = as.numeric(max(df$this_depth))
  i = j
  new_df = df %>%
    filter(this_depth == i) %>%
    select(-this_depth)
  
  new_starting = paste("starting_SMILES", i)
  new_yield = paste("starting_SMILES", i-1)
  
  new_df[[new_starting]] = new_df$starting_SMILES
  new_df[[new_yield]] = new_df$SMILES
  
  new_df = new_df %>%
    select(-starting_SMILES, -SMILES)
  for (i in (j-1):1) {
    new_starting = paste("starting_SMILES", i)
    new_add = paste("add_SMILES")
    new_yield = paste("starting_SMILES", i-1)
    df_sub = df %>%
      filter(this_depth == i) %>%
      select(-this_depth) 
    df_sub[[new_starting]] = df_sub$starting_SMILES
    df_sub[[new_add]] = df_sub$add_SMILES
    df_sub[[new_yield]] = df_sub$SMILES
    df_sub = df_sub %>%
      select(-starting_SMILES, -add_SMILES, -SMILES)
    
    new_df = new_df %>%
      left_join(df_sub)
  }
  return(new_df)
}



