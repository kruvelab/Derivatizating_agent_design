library(tidyverse)
source("PaDEL_descs_calculator.R")
source("connecting_graphs.R")
source("graph_to_SMILES.R")
source("compound_eluent.R")


list_of_fn_groups = read_delim("SMILES_of_fingerprints_corrected.csv",
                               delim = ",",
                               col_names = TRUE)
list_of_fn_groups = list_of_fn_groups %>%
  filter(str_count(SMILES, pattern = "#") == 0)
list_of_fn_groups = list_of_fn_groups$SMILES
regressor_pos = readRDS("regressor_pos.rds")

structures = tibble()

for (j in 1:50) {
  base_stucture = sample(list_of_fn_groups, size = 1)
  i = 1
  for (i in 1:5) {
    add_fn = sample(list_of_fn_groups, size = 1)
    base_stucture = connect_substructures(base_stucture, add_fn)
    base_stucture = graph_to_SMILES(base_stucture)
    structures = structures %>%
      bind_rows(tibble("SMILES" = base_stucture))
    print(base_stucture)
    i = i + 1
  }
}



descs = PaDEL(structures)
descs = descs %>%
  mutate(
    organic_modifier = "MeCN",
    organic = 80,
    pH.aq. = 2.7,
    NH4 = 0,
    viscosity =  viscosity(organic,organic_modifier),
    surface_tension = surfacetension(organic,organic_modifier),
    polarity_index = polarityindex(organic,organic_modifier)) 

IE_pred = descs %>% 
  mutate(logIE_pred_pos = 0)
prediction_pos =  predict(regressor_pos, newdata = IE_pred, predict.all = TRUE)
prediction_pos = prediction_pos$aggregate
IE_pred <- IE_pred %>%
  mutate(logIE_pred_pos = prediction_pos) %>%
  select(SMILES,logIE_pred_pos, everything())
