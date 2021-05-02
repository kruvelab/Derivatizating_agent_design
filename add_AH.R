source("connecting_graphs.R")
source("graph_to_SMILES.R")
source("IE_calculator.R")

Phe = "COC(=O)NC(C(=O)O)c1ccccc1"
Gly = "COC(=O)NCC(=O)O"
Arg = "COC(=O)NC(CCCNC(=N)N)C(=O)O"
Glu = "COC(=O)NC(CCC(=O)O)C(=O)O"

SMILES_to_graph(Phe)
SMILES_to_graph(Gly)

Deriv_AH = function(DerivRegent_SMILES) {
  Deriv_Phe_SMILES = connect_substructures_AH(DerivRegent_SMILES, Phe)
  Deriv_Gly_SMILES = connect_substructures_AH(DerivRegent_SMILES, Gly)
  Deriv_Arg_SMILES = connect_substructures_AH(DerivRegent_SMILES, Arg)
  Deriv_Glu_SMILES = connect_substructures_AH(DerivRegent_SMILES, Glu)
  Deriv_Phe_SMILES = graph_to_SMILES(Deriv_Phe_SMILES)
  Deriv_Gly_SMILES = graph_to_SMILES(Deriv_Gly_SMILES)
  Deriv_Arg_SMILES = graph_to_SMILES(Deriv_Arg_SMILES)
  Deriv_Glu_SMILES = graph_to_SMILES(Deriv_Glu_SMILES)
  Deriv_AH_sturctures = list(
    "Deriv_Phe_SMILES" = Deriv_Phe_SMILES,
    "Deriv_Gly_SMILES" = Deriv_Gly_SMILES,
    "Deriv_Arg_SMILES" = Deriv_Arg_SMILES,
    "Deriv_Glu_SMILES" = Deriv_Glu_SMILES
  )
  return(Deriv_AH_sturctures)
}

calculate_IE_for_AH = function(df_with_structures) {
  
  SMILES_list_Phe = df_with_structures %>% 
    select(SMILES_with_Phe) %>% 
    na.omit() %>%
    rename(SMILES = SMILES_with_Phe) %>%
    unique()

  IE_Phe = IE_calculator_data_frame(SMILES_list_Phe)
  
  SMILES_list_Gly = df_with_structures %>% 
    select(SMILES_with_Gly) %>% 
    na.omit() %>%
    rename(SMILES = SMILES_with_Gly) %>%
    unique()
  
  IE_Gly = IE_calculator_data_frame(SMILES_list_Gly)
  
  SMILES_list_Arg = df_with_structures %>% 
    select(SMILES_with_Arg) %>% 
    na.omit() %>%
    rename(SMILES = SMILES_with_Arg) %>%
    unique()
  
  IE_Arg = IE_calculator_data_frame(SMILES_list_Arg)
  
  SMILES_list_Glu = df_with_structures %>% 
    select(SMILES_with_Glu) %>% 
    na.omit() %>%
    rename(SMILES = SMILES_with_Glu) %>%
    unique()
  
  IE_Glu = IE_calculator_data_frame(SMILES_list_Glu)
  
  df_with_structures = df_with_structures %>%
    left_join(IE_Phe %>%
                rename(SMILES_with_Phe = SMILES,
                       IE_Phe = logIE_pred)) %>%
    left_join(IE_Gly %>%
                rename(SMILES_with_Gly = SMILES,
                       IE_Gly = logIE_pred)) %>%
    left_join(IE_Arg %>%
                rename(SMILES_with_Arg = SMILES,
                       IE_Arg = logIE_pred)) %>%
    left_join(IE_Glu %>%
                rename(SMILES_with_Glu = SMILES,
                       IE_Glu = logIE_pred))
}
