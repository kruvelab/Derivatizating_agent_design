library(tidyverse)
source("draw_fn_group.R")
source("evaluation_function.R")

MonteCarlo = function(search_depth, search_width, MW_max, base_structure) {
  new_structures = tibble(
    starting_SMILES = character(),
    add_SMILES = character(),
    SMILES = character(),
    this_depth = character()
  )
  collected_IE_values = tibble()
  added_fn_groups = tibble()
  while (molecularmass(base_structure) < MW_max  ) { 
    test = add_fn_group(base_structure, search_depth, search_width, new_structures)
    test = organize_df(test)
    derivatized_test = tibble()
    for (smiles in test$`starting_SMILES 0`) {
      derivatized_test = derivatized_test %>%
        bind_rows(Deriv_AH(smiles))
    }
    test = test %>% 
      left_join(derivatized_test)
    test_IE = calculate_IE_for_AH(test)
    collected_IE_values = collected_IE_values %>%
      bind_rows(test_IE)
    test_IE = test_IE %>%
      na.omit()
    best_add_fn_group = evaluation_function(test_IE)
    print(best_add_fn_group)
    base_structure = connect_substructures(base_structure, best_add_fn_group)
    base_structure = graph_to_SMILES(base_structure)
    added_fn_groups = added_fn_groups %>%
      bind_rows(tibble(SMILES_added = best_add_fn_group,
                       SMILES_got = base_structure,
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


