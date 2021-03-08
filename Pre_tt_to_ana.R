

# packages 
library(data.table)
library(tidyverse)
library(lubridate)
library(caret)
# fonction 
source("Source_papier_prett.R",
       encoding = "UTF-8")
# var 

Id <- "id_patient"
name_expo <- paste0(quote(INS_obs_categ))
name_outcome <- paste0(quote(SYM_echelleEpworth))
visite_seuil <- 2
# data



Data_post_pre_tt <- read_rds("data/genere/data_prett.rds") %>% 
  select(-Mesure_unique_INS_visite_Max) %>% # remove instrumental var
  filter(num_visite != 1)   # remove first visit because include in visit-init var



# Select number of visit 
df_EPW <- Data_post_pre_tt %>%
  filter(num_visite <= visite_seuil) %>% 
  # Select patients with SAS
  filter(!is.na(Vis_init_INS_IAH),
         Vis_init_INS_IAH >= 5) %>% #INS_IAH 
  # REmove patients without adherence
  group_by(id_patient)  %>% 
  mutate(INS_observance = sum(is.na(PPC_observanceMoy_finale))) %>%
  filter(INS_observance == 0 ) %>% 
  select(-INS_observance) %>% # retrait de la variables instrumentale
  # REmove patients without Epworth 
  group_by(id_patient) %>%
  mutate(naEPW = sum(is.na(SYM_echelleEpworth)|is.na(Vis_init_SYM_echelleEpworth))) %>% # utile car permet de scale sur plus de rendez-vous
  filter(naEPW == 0) %>% 
  select(-naEPW) %>% 
  ungroup()  


# Adherence to categorical
# Quantile obs
quantile_grp <- df_EPW %>% 
  group_by(id_patient)   %>% 
  ungroup()  %>% 
  do(data.frame(t(round(quantile(df_EPW$PPC_observanceMoy_finale, probs = c(0.25,0.5,0.75)),0))
  )) %>% 
  t %>% 
  as.vector() %>% c(0,.)

# Check because of round 
if(any(diff(round(quantile(df_EPW$PPC_observanceMoy_finale, probs = c(0.25,0.5,0.75)),0)) == 0)) {stop("groupe d'observance superposé à cause de l'arrondie")}


df_obs_categ <- df_EPW %>% 
  mutate(INS_obs_categ = cut(PPC_observanceMoy_finale, 
                             breaks=c(quantile_grp, Inf),
                             include.lowest = TRUE,
                             labels=c("1","2","3","4"))) %>% 
  mutate(INS_obs_categ = as.numeric(INS_obs_categ)) 

# Table temporaire pour sortie
cut_cluster <- df_obs_categ %>% 
  group_by(INS_obs_categ) %>% 
  summarise(min = min(PPC_observanceMoy_finale), max = max(PPC_observanceMoy_finale),.groups = "drop" ) %>% 
  arrange(max)

vec_cut <- c(0,cut_cluster$max)
# Save cut for result printing
saveRDS(vec_cut, file = paste0("data/genere/obs_cut",visite_seuil,".rds"))
df_obs_categ <- df_obs_categ %>% select(-PPC_observanceMoy_finale)
# Gestion des dates et ID

# Date column
date_ddn <- df_obs_categ %>% 
  ungroup() %>% 
  select(which(sapply(.,is.Date)),contains("date"),contains("ddn")#,contains("INS_tmp_entre_rdv")
  ) %>% 
  colnames()
# Id column
id_idv <- df_obs_categ %>% 
  ungroup() %>% 
  select(contains("id_"),contains("idv")) %>% 
  select(-id_patient,
         -Mesure_unique_ATC_id_typeDiabete,
         -PPC_id_typeMasque,
         -Mesure_unique_CHA_id_sexe) %>% 
  colnames()


# remove date and idand transformation to numeric
df_var_select_man <- df_obs_categ  %>%   
  ungroup() %>% 
  select(-one_of(date_ddn),
         -one_of(id_idv)) %>% 
  mutate_all(as.numeric) %>% 
  select_if(!(colSums(is.na(.)) == nrow(.))) # remove char col 



# Remove num visit in case there is only one
if(visite_seuil == 2) df_var_select_man <-  df_var_select_man %>% select(-contains("num_visite"))

# Remove duplicate col
check_duplic <- distinc_col(df_var_select_man)

col_ident <- check_duplic$colonne_suprime

# Gestion des cas ou variables visit_init et var temps dep sont identique
temp_tt_colonne_dupli <- lapply(seq_along(col_ident), function(x) {
  court <- names(which.min(sapply(col_ident[[x]], nchar)))
  restant <- str_remove(col_ident[x] %>% unlist(), court) 
  if (all(restant %in% c("Vis_init_", ""))) {
    court
  }
}) 

duplicate_fact_init <- temp_tt_colonne_dupli %>% compact %>% unlist


# Create list of var who are identical and need to be manually check
duplicate_verif_manu <-   col_ident[lapply(temp_tt_colonne_dupli, function(x) is.null(x) ) %>% unlist()] 

if (!is_empty(duplicate_verif_manu))  {stop("vérif manuel a faire")
  } else df_var_select_no_dup <- check_duplic$df


# Remove var with more than 60 % missing Val

nb_val_manq_par_var <- compte_na_par_var_par_grp(df_var_select_no_dup,Id,colnames(df_var_select_no_dup)) 
seuil_na <- 0.6
frac_val_manq_par_var <- nb_val_manq_par_var %>% 
  mutate_at(vars(-presence_na),list(~(./sum(.)))) %>% 
  filter(presence_na == (visite_seuil - 1)) %>% select(-presence_na)
trop_manquant <- frac_val_manq_par_var %>% select_if(. > 1-seuil_na ) %>% colnames()

corelation_var_df <- df_var_select_no_dup %>% select(-all_of(trop_manquant))

# colinear varaible

seuil_cor <- 0.7
mat_cor <- cor(corelation_var_df %>% select(-all_of(Id))  ,use =  "pairwise.complete.obs")

var_cor <- which(abs(mat_cor) > seuil_cor & (row(mat_cor) != col(mat_cor)) ,arr.ind = TRUE) %>% 
  as_tibble() %>% 
  pivot_wider(names_from = row, values_from = col,values_fn = list(col = list)) %>%
  t.df()


var_correlle <- lapply(unique(
  lapply(
    lapply(split(var_cor, row.names(var_cor)), unlist),
    function(x) sort(
      as.numeric(
        unique(x))
    ))), 
  function(x) colnames(mat_cor)[x])


# Gestion des cas ou variables visit_init et var temps dep corrélé
Var_corr_init_temp_dep <- lapply(var_correlle, function(x) {
  if (length(x) == 2){
    court <- names(which.min(sapply(x, nchar)))
    long <- names(which.max(sapply(x, nchar)))
    different_part <- str_remove(long, court)
    if (different_part == "Vis_init_") {
      result <- court
    }
  }}
) %>% compact %>% unlist

Var_corr_no_init_temp_dep <- lapply(var_correlle, function(x) {
  if (length(x) == 2){
    court <- names(which.min(sapply(x, nchar)))
    long <- names(which.max(sapply(x, nchar)))
    different_part <- str_remove(long, court)
    if (different_part != "Vis_init_") {
      result <- TRUE#c(court,long)
    } else {result <- FALSE}
  } else if (length(x) != 2) {
    result <- TRUE
  }
}
) %>% 
  # compact %>%
  unlist


# list of correlated variables
Var_cor_select_manu <- var_correlle[Var_corr_no_init_temp_dep]  

# variable choose from list above 
Select_manu_cor <- c("Vis_init_CHA_BMI")



var_cor_retire <- c(Var_corr_init_temp_dep, 
                    select_manuel_verif(Var_cor_select_manu,Select_manu_cor))



df_few_pat <- corelation_var_df %>% select(-all_of(var_cor_retire))


# Remove var with to few moda
df_pre_impute <- df_few_pat %>% select(-all_of(nearZeroVar(.,freqCut = 80/20,names = TRUE)))



df_post_imput <- impute_si_changement(df_pre_impute,"data/genere/data_impute.rds",reimputation = FALSE)



# test 
## pas de NA 
if (df_post_imput %>% summarise_all(~sum(is.na(.))) %>% rowSums(.) != 0 ) {
  stop('Row avec Na')}
# test pas de tb erectile chez les femmes

if("Vis_init_SYM_troubleErection" %in% colnames(df_post_imput)){
  if ((df_post_imput  %>% filter(Mesure_unique_CHA_id_sexe == 0) %>% summarise(sum(Vis_init_SYM_troubleErection))) != 0  ) {
    stop('trouble erectile chez les femmes')}
}
saveRDS(df_pre_impute, file = paste0("data/genere/data_pre_impute_nb_vis_",visite_seuil,".rds"))
saveRDS(df_post_imput, file = paste0("data/genere/data_prett_nb_vis_",visite_seuil,".rds"))


