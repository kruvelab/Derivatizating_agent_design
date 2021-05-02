source("PaDEL_descs_calculator.R")
regressor_pos = readRDS("regressor_pos.rds")
source("compound_eluent.R")

IE_calculator = function(structure) {
  structure = tibble("SMILES" = structure)
  descs = PaDEL(structure)
  descs = descs %>%
    #we assume that compound elutes at 80% of MeCN and 0.1% formic acid solution is used
    mutate(
      organic_modifier = "MeCN",
      organic = 80,
      pH.aq. = 2.7,
      NH4 = 0,
      viscosity =  viscosity(organic,organic_modifier),
      surface_tension = surfacetension(organic,organic_modifier),
      polarity_index = polarityindex(organic,organic_modifier)) 
  
  IE_pred = descs %>% 
    mutate(logIE_pred = 0)
  prediction_pos =  predict(regressor_pos, newdata = IE_pred, predict.all = TRUE)
  prediction_pos = prediction_pos$aggregate
  IE_pred <- IE_pred %>%
    mutate(logIE_pred = prediction_pos) %>%
    select(SMILES,logIE_pred, everything())
  return(IE_pred$logIE_pred)
  
}


IE_calculator_data_frame = function(SMILES_list) {

  # exsisting_descs = read_delim("descs.csv",
  #                             delim = ",",
  #                             col_names = TRUE)
  descs = PaDEL(SMILES_list)
  # exsisting_descs = exsisting_descs %>%
  #   bind_rows(descs)
  # write_delim(exsisting_descs,
  #             "descs.csv",
  #             delim = ",")
  descs = descs %>%
    #we assume that compound elutes at 80% of MeCN and 0.1% formic acid solution is used
    mutate(
      organic_modifier = "MeCN",
      organic = 80,
      pH.aq. = 2.7,
      NH4 = 0,
      viscosity =  viscosity(organic,organic_modifier),
      surface_tension = surfacetension(organic,organic_modifier),
      polarity_index = polarityindex(organic,organic_modifier)) 
  
  IE_pred = descs %>% 
    mutate(logIE_pred = 0)
  prediction_pos =  predict(regressor_pos, newdata = IE_pred, predict.all = TRUE)
  prediction_pos = prediction_pos$aggregate
  IE_pred = IE_pred %>%
    mutate(logIE_pred = prediction_pos) %>%
    select(SMILES,logIE_pred)
  return(IE_pred)
}
