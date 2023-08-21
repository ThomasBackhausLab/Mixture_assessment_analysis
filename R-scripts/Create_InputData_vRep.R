##############
### Standardtasks at the beginning
##############
rm(list=ls()) #remove ALL objects 
cat("\014") # clear console window prior to new run
gc() # clear memory

Sys.setenv(LANG = "en") #Let's keep stuff in English
Sys.setlocale("LC_ALL")
options(scipen = 999)
#############
# Global Variables
#############
RunAll <- TRUE
vers <- 1

#############
# Author and version
#############
# Author: Pedro A. Inostroza
# version: Repositoty

#############
# libs
#############
packages<- c("tidyverse", "readxl", "writexl")
lapply(packages, library, character.only= TRUE)

inwd  <- "/Users/xinope/M_GU/R/RProjects/FRAM_CECs_Chile2018/raw_data"
outwd <- "/Users/xinope/M_GU/R/RProjects/FRAM_CECs_Chile2018/outputs"
setwd("/Users/xinope/M_GU/R/RProjects/FRAM_CECs_Chile2018")

####################################
# load ecotox input data
####################################
# ECOTOX data set
ECOTOX<- read_delim(paste0(inwd,"/ECOTOX Summary_v12_Rev12_09_15_2022.csv"),locale=locale(decimal_mark=","),show_col_types= F)

# QSAR-VEGA predictions (chronic model only)
vega<- read_excel(paste0(inwd,"/QSAR_vega_predictions.xlsx"), col_names= TRUE, trim_ws= TRUE)

# QSAR-ECOSAR predictions 
ecosar<- read_delim(paste0(inwd,"/QSAR_ecosar_predictions.csv"),locale=locale(decimal_mark="."), show_col_types= F)

####################################
# load chemical input data
####################################
# List with molecular weight
Chemical_MW<- read_delim(paste0(inwd,"/CAS_Lookup_MW_for_ECOTOX_Output_v9_Rev_15 (2).csv"),locale=locale(decimal_mark="."),show_col_types= F, delim = "\t")

# List of quantified chemicals in the River Aconcagua
InputChems <- readRDS(paste0(inwd, "/MZQuant_Water_River_Chile_2018_vRepository.rds"))

# UFZ target list with MDL values
InputMdl<- read_excel(paste0(inwd, "/Target_ECs_WANA_2018_v4.xlsx"), col_names= TRUE, col_types= NULL) 

############
# Defined factor and levels
############
sites.order<- c("RS1", "RS2", "RS3", "T1", "T2", "T3", "R1", "R2", "R3")

############
# Removed symbols
############
symbols_remove<- c(">", "<", "=", ">=", "<=", "~")

####################################
# Preliminary quality check
####################################

### Measured Environmental Concentration Data

#-   Imported columns as character and thus measured environmental concentrations (MECs) need to be transformed to numeric.
#-   Concentrations converted to uM since effect data (US ECOTOX Database) is in uMolar
#-   MDL (Method Detection Limit) values were imported from UFZ list

InputChems[sites.order] <- sapply(InputChems[sites.order],as.numeric)

Chemical_MW<- Chemical_MW %>%
  mutate_at("mol_wt", as.numeric)

############
# re-arrange InputChems 
############

InputChems2<- InputChems%>% select(chemical_name, cas_number, RS1, RS2, RS3, T1, T2, T3, R1, R2, R3) %>%
  rename(CASRN = cas_number) %>%
  merge(InputMdl[c("CASRN", "MDL_ngL")], by = "CASRN", all.x = TRUE) %>%
  pivot_longer(cols = all_of(sites.order), names_to = "sites", values_to = "MEC_ngL") %>%
  rename(cas_number = CASRN) %>%
  mutate(conc_ngL = ifelse(MEC_ngL == 0, MDL_ngL, MEC_ngL)) %>%
  merge(Chemical_MW[c("cas_number", "mol_wt")], by = "cas_number", all.x = TRUE) %>%
  mutate(MDL_uM = (MDL_ngL/1000)/mol_wt)%>%
  mutate(MEC_uM = (MEC_ngL/1000)/mol_wt)%>%
  mutate(conc_uM = (conc_ngL/1000)/mol_wt)%>%
  select(cas_number, chemical_name, sites, MDL_uM, MEC_uM, conc_uM)

# select cas numbers for further filtering (Effect Data section)
cas_InputChems<- InputChems2$cas_number

####################################
# Experimental and In-Silico Effect Data
####################################

#US EPA ECOTOX Database was filtered as followed:
#  
#-   only freshwater exposure media
#-   only chronic values
#-   only chemicals present in InputChems file

# The definition of a chronic effect data was based on the duration time of the experiment according to 
# the Australian and New Zealand water directive (Table 1) (https://www.waterquality.gov.au/sites/default/files/documents/warne-wqg-derivation2018.pdf).

#-   Toxic Test = Chronic
#-   Fish and amphibians = adults/juveniles = \>= 21d
#-   Macroinvertebrates = adults/juveniles/larvae = \>= 14d
#-   Microalgae = \> 24h

#Relevant Endpoints

#-   Fish and amphibians = adults/juveniles = ALL
#-   Macroinvertebrates = adults/juveniles/larvae = ALL
#-   Microalgae = ALL


EC10_ecotox<- ECOTOX %>% select(cas_number,
                                latin_name,
                                species_group,
                                trophic_level,
                                media_type,
                                obs_duration_mean,
                                endpoint,
                                AC_Class,
                                conc1_mean_op_original,
                                EC10_Equiv_uM)%>%
  mutate_at("EC10_Equiv_uM", as.numeric) %>%
  mutate_at("obs_duration_mean", as.numeric) %>%
  filter(media_type == "FW" & AC_Class == "Chronic" & conc1_mean_op_original == "None") %>%
  filter(cas_number %in% cas_InputChems) %>%
  mutate(group = case_when(
    species_group == "Algae" & obs_duration_mean > 1  ~ "positive",
    species_group == "Fungi" & obs_duration_mean >1  ~ "positive",
    species_group == "Crustaceans" & obs_duration_mean >= 14  ~ "positive",
    species_group == "Invertebrates" & obs_duration_mean >=14  ~ "positive",
    species_group == "Molluscs" & obs_duration_mean >=14  ~ "positive",
    species_group == "Insects/Spiders" & obs_duration_mean >=14  ~ "positive",
    species_group == "Worms" & obs_duration_mean >=14  ~ "positive",
    species_group == "Fish" & obs_duration_mean >= 21  ~ "positive",
    species_group == "Amphibians" & obs_duration_mean >= 21  ~ "positive"))%>%
  filter(group == "positive") %>%
  filter(EC10_Equiv_uM > 0) %>%
  mutate(species_group = case_when(
    endsWith(trophic_level, "Producer")  ~ "alg_ecotox_uM",
    endsWith(trophic_level, "Primary Consumer")  ~ "crust_ecotox_uM",
    endsWith(trophic_level, "Secondary Consumer")  ~ "fish_ecotox_uM")) %>%
  drop_na(species_group) %>%
  group_by(cas_number, species_group) %>%
  summarise(ecotox_uM = exp(mean(log(EC10_Equiv_uM))))

EC10_vega <- vega %>%
  rename(cas_number = "Id") %>%
  rename(algae_vega = "Algae Chronic (NOEC) Toxicity model (IRFMN) - assessment") %>%
  rename(crust_vega = "Daphnia Magna Chronic (NOEC) Toxicity model (IRFMN) - assessment") %>%
  rename(fish_vega = "Fish Chronic (NOEC) Toxicity model (IRFMN) - assessment") %>%
  pivot_longer(cols= c("algae_vega",
                       "crust_vega",
                       "fish_vega"), names_to= "source", values_to= "conc") %>%
  separate(conc, c("conc", "score"), "mg/L") %>%
  mutate(score = str_replace_all(score, "[(]", "")) %>%
  mutate(score = str_replace_all(score, "[)]", "")) %>%
  mutate(vega_score = case_when(
    startsWith(score, " low")  ~ "Low",
    startsWith(score, " moderate")  ~ "Moderate",
    startsWith(score, " good")   ~ "Good",
    startsWith(score, " EXP") ~ "Experimental")) %>%
  drop_na(conc) %>%
  filter(conc >0) %>%
  mutate_at("conc", as.numeric) %>%
  filter(cas_number %in% cas_InputChems) %>%
  merge(Chemical_MW[c("cas_number", "mol_wt")], by = "cas_number", all.x = TRUE) %>%
  filter(mol_wt > 0)%>%
  # VEGA concentrations in mg/L
  # multiple per 1000 = ug/L
  mutate(vega_uM = (conc*1E3)/mol_wt) %>%
  mutate(species_group = case_when(
    startsWith(source, "algae")  ~ "alg_vega_uM",
    startsWith(source, "crust")  ~ "crust_vega_uM",
    startsWith(source, "fish")   ~ "fish_vega_uM"))%>%
  distinct(cas_number, species_group, .keep_all = TRUE)%>%
  select(cas_number, species_group, vega_uM)

# Warning message:
# Expected 2 pieces. Missing pieces filled with `NA` in 3 rows [1057, 1058, 1059]. # chemicals without CAS numbers

EC10_ecosar<- ecosar %>%
  select(original_CAS, `value (mg/L)`, model_organism)%>%
  rename(cas_number = "original_CAS")%>%
  rename(conc_mgL = "value (mg/L)")%>%
  rename(species_group = "model_organism")%>%
  filter(cas_number %in% cas_InputChems) %>%
  merge(Chemical_MW[c("cas_number", "mol_wt")], by = "cas_number", all.x = TRUE) %>%
  filter(mol_wt > 0)%>%
  # ECOSAR concentrations in mg/L
  # multiple per 1000 = ug/L
  mutate(ecosar_uM = (conc_mgL*1E3)/mol_wt) %>%
  mutate(species_group = case_when(
    startsWith(species_group, "algae")  ~ "alg_ecosar_uM",
    startsWith(species_group, "daphnia")  ~ "crust_ecosar_uM",
    startsWith(species_group, "fish")   ~ "fish_ecosar_uM"))%>%
  distinct(cas_number, species_group, .keep_all = TRUE)%>%
  select(cas_number, species_group, ecosar_uM)

####################################
# Consolidate and create InputData file for the mixture risk assessment script.
####################################

InputData <- InputChems2 %>% merge(EC10_ecotox[c("cas_number", "species_group", "ecotox_uM")], by = "cas_number", all.x = TRUE) %>%
  merge(EC10_vega[c("cas_number", "species_group", "vega_uM")], by = "cas_number", all.x = TRUE) %>%
  merge(EC10_ecosar[c("cas_number", "species_group", "ecosar_uM")], by = "cas_number", all.x = TRUE) %>%
  pivot_wider(names_from = species_group, values_from = ecosar_uM) %>% select(!"NA") %>%
  pivot_wider(names_from = species_group.x, values_from = ecotox_uM) %>% select(!"NA") %>%
  pivot_wider(names_from = species_group.y, values_from = vega_uM) %>% select(!"NA") %>%
  mutate_all(~replace(., is.na(.), -7777))%>%
  select(sites, cas_number, chemical_name, MEC_uM, MDL_uM, conc_uM,
         alg_ecotox_uM, crust_ecotox_uM, fish_ecotox_uM,
         alg_ecosar_uM, crust_ecosar_uM, fish_ecosar_uM,
         alg_vega_uM, crust_vega_uM, fish_vega_uM)

####################################
# Export InputData
####################################

# rds format
saveRDS(InputData, paste0(outwd, "/Summary_Table_Chile2018_vRepository", ".rds"))
# tab separated
write_tsv(InputData, paste0(outwd, "/Summary_Table_Chile2018_vRepository", ".tsv"))
# csv
write_delim(InputData, paste0(outwd, "/Summary_Table_Chile2018_vRepository", ".csv"), delim = "\t")



