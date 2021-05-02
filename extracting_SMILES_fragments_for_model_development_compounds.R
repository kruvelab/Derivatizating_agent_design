library(tidyverse)
library(stringr)
library(webchem)
library(rcdk)
setwd("C:/Users/annel/Nextcloud/minu asjad/Sissejuhatus andmeteadusesse/projekt/code")
source('fingerprints_hex_to_dataframe.R')

fn_CAS_to_CID = function(CAS) {
  CID_list = webchem::cts_convert(CAS, "cas", "pubchem cid", match = "first")
  CID = CID_list[[CAS]]
  return(CID)
}

#----Combining all neccessary PubChemCIDs----
setwd("C:/Users/annel/Nextcloud/mudeli script ja failid/andmed/retraining_190919/positive/PubChem fingerprints/with CID")

compounds_Gunda <- read_delim("compounds_Gunda_data.csv",
                               delim = ";",
                               col_names = TRUE)
compounds_Gunda = compounds_Gunda %>%
  rename(name = s) %>%
  mutate(X2 = case_when(
    X2 != "#N/A" ~ X2,
    TRUE ~ "0")) %>%
  mutate(X3 = case_when(
    X3 != NA ~ X3,
    TRUE ~ 0))

compounds_Gunda = compounds_Gunda %>%
  mutate(CID = case_when(
    X2 != 0 ~ as.double(X2),
    X3 != 0 ~ as.double(X3),
    TRUE ~ 0
  )) %>%
  filter(CID != 0)

compounds_Tingting <- read_delim("compounds_Tingting_data.csv",
                              delim = ";",
                              col_names = TRUE)


compounds_Tingting = compounds_Tingting %>%
  mutate(X2 = case_when(
    X2 != "#N/A" ~ X2,
    TRUE ~ "0")) %>%
  mutate(X3 = case_when(
    X3 != "#N/A" ~ X3,
    TRUE ~ "0"))

compounds_Tingting = compounds_Tingting %>%
  mutate(CID = case_when(
    X2 != 0 ~ as.double(X2),
    X3 != 0 ~ as.double(X3),
    TRUE ~ 0
  )) %>%
  filter(CID != 0)

compounds_Jaanus <- read_delim("compounds_Jaanus_data.csv",
                                 delim = ";",
                                 col_names = TRUE)

compounds_Jaanus = compounds_Jaanus %>%
  mutate(X2 = case_when(
    X2 != "#N/A" ~ X2,
    TRUE ~ "0")) %>%
  mutate(X3 = case_when(
    X3 != "#N/A" ~ X3,
    TRUE ~ "0"))

compounds_Jaanus = compounds_Jaanus %>%
  mutate(CID = case_when(
    X2 != 0 ~ as.double(X2),
    X3 != 0 ~ as.double(X3),
    TRUE ~ 0
  )) %>%
  filter(CID != 0)

compounds_ENTACT <- read_delim("compounds_ENTACT_data.csv",
                               delim = ",",
                               col_names = TRUE)
DTX_CAS <- read_delim("Dsstox_CAS_number_name.csv",
                               delim = ",",
                               col_names = TRUE) %>%
  rename(DTXSID = dsstox_substance_id)

compounds_ENTACT = compounds_ENTACT %>%
  left_join(DTX_CAS)

compounds_ENTACT = compounds_ENTACT %>%
  group_by(casrn) %>%
  mutate(PubChemCID = fn_CAS_to_CID(casrn)) %>%
  ungroup()

compounds_ENTACT = compounds_ENTACT %>%
  rename(CID = PubChemCID) %>%
  rename(name = preferred_name)

compounds_ENTACT = compounds_ENTACT %>%
  mutate(CID = as.double(CID))

compounds_Karin = read_delim("Karin_compounds.csv",
                             delim = ",",
                             col_names = TRUE)

compounds_Karin = compounds_Karin %>%
  rename(name = Name) %>%
  select(name, InChIKey, CAS_Nr) %>%
  unique()
  
compounds_Karin_CID = compounds_Karin %>%
  filter(CAS_Nr != "74223-64-6") %>%
  select(CAS_Nr) %>%
  group_by(CAS_Nr) %>%
  mutate(CID = fn_CAS_to_CID(CAS_Nr)) %>%
  ungroup()

compounds_Karin = compounds_Karin %>%
  left_join(compounds_Karin_CID)

compounds_Karin = compounds_Karin %>%
  mutate(CID = as.double(CID))

#collect all CID values that need to be looked up and 
compounds <- compounds_ENTACT %>% select(name, CID) %>%
  bind_rows(compounds_Gunda %>% select(name, CID)) %>%
  bind_rows(compounds_Jaanus %>% select(name, CID)) %>%
  bind_rows(compounds_Karin %>% select(name, CID)) %>%
  bind_rows(compounds_Tingting %>% select(name, CID)) %>%
  na.omit() %>%
  unique()

compounds = compounds %>%
  rename(PubChemCID = CID)

write_delim(compounds,
            "all_compounds_with_CID.csv",
            delim = ",")

compounds = compounds %>%
  select(PubChemCID)

setwd("C:/Users/annel/OneDrive - Kruvelab/Kruvelab/computational/PubChem")

data_path <- "C:/Users/annel/OneDrive - Kruvelab/Kruvelab/computational/PubChem/fingerprints"

files <- dir(data_path, pattern = "*.csv") # get file names
remove <- c("fingerprint_Compound_", ".sdf.csv") #idetify the part of filename that will be removed
Fingerprints_summary <- tibble() #here we will collect all matchies 

setwd(data_path)
for (filename in files) {
  if (length(Fingerprints_summary$PubChemCID) < length(compounds$PubChemCID)) {
    datafile <- read_delim(filename,
                           delim = "\t",
                           col_names = TRUE) %>%
      rename(PubChemCID = V2,
             SMILES = V3,
             Fingerprint = V4)
    print(filename)
    PubCID_range <- str_remove_all(filename, paste(remove, collapse = "|"))
    PubCID_range <- str_split(PubCID_range, "_")
    first_CID <- as.double(PubCID_range[[1]][1])
    last_CID <- as.double(PubCID_range[[1]][2])
    compounds_small <- compounds %>%
      filter(PubChemCID > first_CID & PubChemCID < last_CID) %>%
      left_join(datafile)
    Fingerprints_summary <- Fingerprints_summary %>%
      bind_rows(compounds_small)
  } else {
    break
  }
}
setwd("C:/Users/annel/Nextcloud/mudeli script ja failid/andmed/retraining_190919/positive/PubChem fingerprints/with CID")


write_delim(Fingerprints_summary,
            "Fingerprints_all.csv",
            delim = ",")
Fingerprints_summary <- Fingerprints_summary %>%
  na.omit()

Decoded_figenrprints <- tibble()
for (fingerprint in Fingerprints_summary$Fingerprint) {
  current_fingerprint <- fingerprint_hex_to_dataframe(fingerprint)
  print(current_fingerprint)
  Decoded_figenrprints <- Decoded_figenrprints %>%
    bind_rows(current_fingerprint)
}

fingerprints <- as_tibble(cbind(Fingerprints_summary, Decoded_figenrprints))

write_delim(fingerprints,
            "Fingerprints_all.csv",
            delim = ",")

fingerprints = read_delim("Fingerprints_all.csv",
                             delim = ",",
                             col_names = TRUE)

fingerprints = fingerprints %>%
  select(-SMILES, -X1, -SDFlineStart, -SDFlineEnd, -SDFID, -Fingerprint) %>%
  select(everything(), PubChemCID)

fingerprints = fingerprints %>%
  select_if(~mean(.) >0.01)

# 
# fingerprints = fingerprints %>%
#   select(-caret::nearZeroVar(fingerprints)) #removing fingerprints that do not change significantely between samples

write_delim(fingerprints,
            "fingerprints.csv",
            delim = ",")

list_of_SMILES = read_delim("PubChem_descs_names.csv",
                            delim = ",",
                            col_names = TRUE)

SMILES = tibble(Position = fingerprints %>% colnames())

SMILES = SMILES %>%
  filter(Position != "PubChemCID")

SMILES = SMILES %>%
  left_join(list_of_SMILES)

write_delim(SMILES,
            "SMILES_of_fingerprints.csv",
            delim = ",")  

SMILES = read_delim("SMILES_of_fingerprints_corrected.csv",
                    delim = ",",
                    col_names = TRUE)

SMILES$SMILES[30]

SMILES = SMILES %>%
  mutate(SMILES = str_replace_all(SMILES, 
                                  pattern = "~", 
                                  replacement = ""))

parse.smiles(SMILES$SMILES[30])[[1]]


write_delim(SMILES,
            "SMILES_of_fingerprints_corrected.csv",
            delim = ",")  

