library(haven)
library(data.table) 
library(tidyverse)
library(lubridate)
start_time <- Sys.time()
conflicted::conflict_prefer("filter", "dplyr")
`%+%` <- function(x, y)  mapply(sum, x, y, MoreArgs = list(na.rm = TRUE))

df <- read_sas("data/OSFP_F/table_osfp_visites.sas7bdat")  

# Vérification car pas supposé perdre de patients dans l'opérations 
nombre_ligne_init <- nrow(df)

# transformation des colonnes contenant date en type date  
df <- df %>% mutate_if( str_detect(colnames(.),"date"),as_date) %>%  
  # ajout d'une variable instrumental sévérité du sas  
  mutate(CHA_id_sexe = as.numeric(CHA_id_sexe == 1)) %>% # recodé le sex entre 0 et 1 
  mutate(INS_IAH = coalesce(PSG_IAH_heure, PSG_nbrDesatO2),
         PPC_nomVentilateur = as.numeric(PPC_nomVentilateur == "")) %>% # transformation des noms de ventilateur en booléen car je ne sais pas quoi faire d'un facteur a 5 k modalité
  select(-PSG_IAH_heure,-PSG_nbrDesatO2) %>% 
  group_by(id_patient) %>%  
  arrange(num_visite) %>%  
  mutate(INS_visite_Max = max(num_visite)) %>%  # ajout d'une variable instrumentale nombres max de visite  
  fill(INS_IAH,CHA_taille) %>%   # rmplacé les manquant de sévérité et tailles par la derniere valeurs  
  mutate(CHA_BMI = ifelse(is.na(CHA_BMI) , (CHA_poids / (CHA_taille/100)^2),CHA_BMI)) %>% # caclcule des BMI  
  mutate(INS_tmp_entre_rdv = CHA_age - dplyr::lag(CHA_age)) %>% # ajout d'une variable instrumentale temps entre les rdv  
  mutate(INS_tmp_entre_rdv = ifelse(is.na(INS_tmp_entre_rdv),0,INS_tmp_entre_rdv)) %>% 
  mutate(CHA_id_sexe = ifelse((is.na(CHA_id_sexe) & SYM_troubleErection == 1),1,CHA_id_sexe)) %>% # les patient avec des troubles de l'erection sont des hommes
  fill(CHA_id_sexe) %>% # remplace les manquant du sex par la derniere valeurs 
  mutate(SOM_estimDureeMoyDeSom_Cor = as.numeric(SOM_estimDureeMoyDeSom_Cor, units="secs")/3600) %>%  # transformation du temps de sommeil en durée en heure
  mutate(SOM_estimDureeMoyDeSom_Cor = ifelse(SOM_estimDureeMoyDeSom_Cor >= 14,NA,SOM_estimDureeMoyDeSom_Cor)) %>% 
  
  arrange(id_patient) %>% 
  ungroup() %>% 
  # ajout d'une correction sur ronflements sans s car sinon ne mathc pas échelles
  rename( "SYM_ronflement" = "SYM_ronflements")

# Nouveaux pré traitement groupement des effets indésirable de la PPC

vec_effet_indez_ppc <- c("PPC_aerophagie",
"PPC_boucheSeche",
"PPC_cephalee",
"PPC_genePsychologique",
"PPC_lesionsCutanees",
"PPC_nezBouche",
"PPC_sensationEtouffement",
"PPC_toleranceEntourage",
"PPC_yeuxIrrites")



df <- df %>% mutate(INS_PPC_effet_ind = reduce(select(.,all_of(vec_effet_indez_ppc)),
                                  `%+%`)) %>% select(-all_of(vec_effet_indez_ppc))

# remove column echelle et autre column corrélé

df <- df %>% select(everything(),-contains("GDS"),GDS_gazDuSang,-contains("EFR"),EFR_spirometrie,-ends_with("_echelle"))


# Choix de retirer des variables pas ou peu informative
df <- df %>% select(-c(SOM_id_autrePathoSommeil,
                        ATC_obesite,
                        CHA_EESP_Medecin,
                        CHA_EESP_Patient,
                        CHA_id_profession,
                        PPC_id_TypePPC,
                        TTT_activitePhysique,
                        TTT_conseilsHygienoDietetiques),
                    -SOM_estimDureeMoyDeSom_Cor, # retirer car pause problème (durée inf à l'obs et trop grande  + vlaur a  0)

                    # retrai%>% taille poids comprise dans le BMI
                    -CHA_poids,
                    -CHA_taille,

                    )


# Variables ATCD  remplace NA par le max du patient
# version DPLY posan%>% des problèmes de performances


df <- df %>% 
  group_by(id_patient) %>%
  mutate_at(vars(contains("ATC_")),~ifelse(
    is.numeric(.),
    ifelse(sum(!is.na(.))!=0,max(.,na.rm = TRUE),NA),
    .)) %>% 
  ungroup()  %>% 
  replace_na(list(ATC_id_typeDiabete = 0))






 unique_mesure <- setDT(df)[, lapply(.SD,
                                          function(x) length(unique(na.omit(x)))), 
                                 by = id_patient] %>%
   .[, lapply(.SD, max)] %>% 
   data.frame()

df <- df %>% tibble()



colonne_name_mesure_unique <- unique_mesure %>%
  select_if(. <= 1) %>% 
  colnames() %>% 
  gsub("_nb_valeur","",.)

#select_if(function(col) n_distinct(col, na.rm = TRUE) > 1)


colonne_renseigne_une_fois <-  df %>% 
  select(id_patient,num_visite,all_of(colonne_name_mesure_unique)) %>% 
  # filter_at(vars(unique_mesure %>% select_if(. <= 1) %>% 
  #                  colnames() %>%
  #                  gsub("_nb_valeur","",.)),any_vars(!is.na(.))) %>% #  vérification que les autres  visite exclusivment composé de NA 
  fill() %>% 
  rename_all(function(x) paste0("Mesure_unique_", x)) %>%
  rename(id_patient = 1,num_visite = 2)


ncol_df_temp <- ncol(df)
df <- colonne_renseigne_une_fois %>% 
  inner_join(df,by = c("id_patient" = "id_patient","num_visite" = "num_visite"))   %>% 
  select(-all_of(colonne_name_mesure_unique))

if (ncol_df_temp != ncol(df) ) {
  stop('Nombres de colones différent après récupérations des variables uniques')}
rm(ncol_df_temp)

# recehrche des variables avec une variance null


df <- df %>% 
  ungroup() %>% 
  filter(num_visite == 1) %>% 
  select( -num_visite ) %>% 
  select(-contains("Mesure_unique_") ) %>% 
  rename_all(function(x) paste0("Vis_init_", x)) %>% 
  select_if(!(colSums(is.na(.)) == nrow(.)) ) %>% # remove column sans valeurs renseigné a la premireres visites 
  inner_join(df,by = c("Vis_init_id_patient" = "id_patient")) %>% 
  rename(id_patient = 1)

# Vérification car pas supposé perdre de patients dans l'opérations 
nombre_ligne_final <- nrow(df)


if (nombre_ligne_init != nombre_ligne_final ) {
  stop('Nombres de ligne différent après pré traitement')}

rm(nombre_ligne_final,nombre_ligne_init)
 
colonnes_ss_var_une_val <-  df %>%  # liste des colonnes n'ayant q'une valeur différente de na 
  select_if(function(col) n_distinct(col,na.rm = TRUE) <= 1) %>% colnames() 
 
df <-  df %>% 
  select_if(function(col) n_distinct(col,na.rm = TRUE) > 1)  # retire les variables visite init sans valeurs toute = a 0


# Missing in diabete to 0 
df <- df %>% 
  replace_na(list(Vis_init_ATC_id_typeDiabete = 0,
                  ATC_id_typeDiabete = 0)) 


saveRDS(df, file = paste0("data/genere/data_prett",".rds"))

end_time <- Sys.time()
end_time - start_time