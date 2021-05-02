
evaluation_function = function(IE_table) {
  evaluation_table =   IE_table %>%
    group_by(add_SMILES) %>%
    summarise(median_IE_Phe = median(IE_Phe),
              median_IE_Gly = median(IE_Gly),
              median_IE_Arg = median(IE_Arg),
              median_IE_Glu = median(IE_Glu)) %>%
    ungroup() %>%
    mutate(IE_Phe_rank = rank(median_IE_Phe),
           IE_Gly_rank = rank(median_IE_Gly),
           IE_Arg_rank = rank(median_IE_Arg),
           IE_Glu_rank = rank(median_IE_Glu)) %>% 
    mutate(score = IE_Phe_rank+3*IE_Gly_rank + 2*IE_Arg_rank + IE_Glu_rank) 
  best_fn_group = evaluation_table%>%
    filter(score == max(score)) %>%
    select(add_SMILES)
  return(best_fn_group$add_SMILES)
}

