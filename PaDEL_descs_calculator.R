library(tidyverse)
library(stringr)
library(rjson)

PaDEL = function(SMILES_list) {
  descs = tibble()
  command = "java -jar descriptor-cli-0.1a-SNAPSHOT-all.jar"
  for (i in 1:dim(SMILES_list)[1]) {
    smiles = SMILES_list$SMILES[i]
    command_final = paste(command, smiles, sep =" ")
    descs_this_smiles = calculation_by_SMILES(smiles, command_final)
    descs = descs %>%
      bind_rows(descs_this_smiles)
  }
  # write_delim(descs,
  #             "descs.csv",
  #             delim = ",")
  return(descs)
}

calculation_by_SMILES = function(smiles, command_final) {
  tryCatch(
    {
      javaOutput = system(command_final, intern = TRUE)
      output_string = str_split(javaOutput, pattern = " ")
      descs_json = paste(output_string[1][[1]], output_string[2][[1]], output_string[3][[1]], output_string[4][[1]], sep = "")
      descs_this_smiles = fromJSON(descs_json)
      descs_this_smiles = data.frame(descs_this_smiles)
      descs_this_smiles = descs_this_smiles %>%
        mutate(SMILES = smiles)
      return(descs_this_smiles)
    },
    error = function(e) {
      return(tibble())
    }
  )
}
