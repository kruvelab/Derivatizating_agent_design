library(tidyverse)
source("draw_fn_group.R")
source("evaluation_function.R")


base_structure = "C"
MW_max = 300
search_width = 2

MonteCarlo = function(search_width, MW_max, base_structure) {
  new_structures = tibble(
    starting_SMILES = character(),
    add_SMILES = character(),
    SMILES = character(),
    this_depth = character()
  )
  collected_IE_values = tibble()
  added_fn_groups = tibble()
  # best_IE_value_previous = 0
  # best_IE_value = (IE_calculator(
  #   tibble(SMILES = Deriv_AH(base_structure)$Deriv_Gly_SMILES)
  # ))
  while (molecularmass(base_structure) < MW_max  ) { #&& (best_IE_value - best_IE_value_previous) > 0.2
    test = add_fn_group(base_structure, search_width, search_width, new_structures)
    test = organize_df(test)
    
    test = test %>%
      group_by(`starting_SMILES 0`) %>%
      mutate(SMILES_with_Phe = Deriv_AH(`starting_SMILES 0`)$Deriv_Phe_SMILES,
             SMILES_with_Gly = Deriv_AH(`starting_SMILES 0`)$Deriv_Gly_SMILES,
             SMILES_with_Arg = Deriv_AH(`starting_SMILES 0`)$Deriv_Arg_SMILES,
             SMILES_with_Glu = Deriv_AH(`starting_SMILES 0`)$Deriv_Glu_SMILES) %>%
      ungroup() 
    
    test_IE = calculate_IE_for_AH(test)
    collected_IE_values = collected_IE_values %>%
      bind_rows(test_IE)
    best_add_fn_group = evaluation_function(test_IE)
    base_structure = connect_substructures(base_structure, best_add_fn_group)
    base_structure = graph_to_SMILES(base_structure)
    # best_IE_value_previous = best_IE_value
    # best_IE_value = (IE_calculator(
    #   tibble(SMILES = Deriv_AH(base_structure)$Deriv_Gly_SMILES)
    # ))
    # print(best_IE_value)
    added_fn_groups = added_fn_groups %>%
      bind_rows(tibble(SMILES_added = best_add_fn_group,
                       SMILES_got = base_structure,
                       #logIE = best_IE_value
                       ))
  }
  write_delim(collected_IE_values,
              paste(paste(MW_max, search_width, sep = "_"),
                    make_clean_names(Sys.time()), 
                    ".csv", 
                    sep = ""),
              delim = ",")
  write_delim(added_fn_groups,
              paste(paste(MW_max, search_width, sep = "_"),
                    make_clean_names(Sys.time()),
                    "added_fn_grous", 
                    ".csv", 
                    sep = ""),
              delim = ",")
  
}


