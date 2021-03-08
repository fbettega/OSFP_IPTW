# packages 
library(tidyverse)
library(nnet)
library(caret)
library(smd)
library(kableExtra)
# fonction 
source("Source_papier_prett.R",
       encoding = "UTF-8")
# var 

Id <- "id_patient"
name_expo <- paste0(quote(INS_obs_categ))
name_outcome <- paste0(quote(SYM_echelleEpworth))
visite_seuil <- 2
# data
df_modele_pds <- read_rds("data/genere/data_prett_nb_vis_2.rds") 



# perform model variables selection
select_var  <-  model_selection(data_frame = df_modele_pds,name_expo , name_outcome,rerun = FALSE)


Var_model_outcome <- select_var$var_select$model_outcome





# function that compute outcome estimation with linear regression and IPWRA 
# and compute weight  IPTW  for IPW
function_res_papier <- function(data, index,params_fun ){
  df_fun <- data[index,]

  
  Id_patient_analyse <- df_fun %>% select(all_of(c(params_fun$Id_data,params_fun$expo)))

  res_function_poids <- do.call(calcule_pds_stage,c(list(donne = df_fun),
                                                    params_fun))
  
  

  poids_strap <- res_function_poids$res_intermediaire$poids$poids_tronc$poids_trunc %>% 
    unlist(use.names = FALSE)

  
  donnee_pre_mod_fin <- df_fun %>% mutate_at(params_fun$expo,list(~relevel(., ref = "4"
  ))) %>% select(-all_of(params_fun$Id_data))
  
  model_pred_outcomes <- glm( as.formula(paste0(params$out_come, " ~ .")),
                              data = donnee_pre_mod_fin,
                              family = gaussian(), 
                              weights = poids_strap
  )
  
  
  model_pred_outcomes_non_pond <- glm(as.formula(paste0(params$out_come, " ~ .")),
                                       data = donnee_pre_mod_fin,
                                       family = gaussian()
  )
  
  df_pred_ss_pond <-  predict(model_pred_outcomes_non_pond,newdata = donnee_pre_mod_fin) 
  
  unique_expo <- df_fun %>% select(all_of(params_fun$expo)) %>% unlist(use.names = FALSE)

  expo_fact_fun <- params_fun$expo

  df_pred <- donnee_pre_mod_fin %>% select(-all_of(params_fun$expo))  %>% #data.table::
    merge(data.frame(temp_name = unique(unique_expo))) %>% rename(!!expo_fact_fun:=temp_name)
  #print(df_pred)
  treatment_effect_df <- predict(model_pred_outcomes,newdata = df_pred) %>% matrix(nrow = nrow(df_fun) ,ncol = length(unique(unique_expo)))
  
  colnames(treatment_effect_df)  <- paste0("treatment_grp_", unique(unique_expo))
  
  return(list(var_selec = Var_model_outcome,
              poids = poids_strap,
              reg_pond= coef(model_pred_outcomes),
              Treatment_effect = treatment_effect_df,
              df_id_grp_tt = Id_patient_analyse,
              result_func = res_function_poids$df,
              reg_non_pond = coef(model_pred_outcomes_non_pond),
              estim_reg_lin = df_pred_ss_pond,
              ind = index
  ))
}


params <- list(expo = name_expo,
               covar = Var_model_outcome, #var_mod_pds ,
               Id_data = Id ,
               weighting_function = multinomial_IPTW,
               out_come = name_outcome,
               percentile_tronc = 0.01
               )
donnee_boot <- select_var$data_frame$df_analyse %>% 
  mutate_at( typages_function(select_var$data_frame$df_analyse,10)$colonne_type$facteur,
             factor) # %>% 
  #as.data.frame() 


if("OK" =="OK"){
  boot_strap_full_process <- boot_strap_fun(donnee_boot,function_res_papier,1000,params)
  saveRDS(boot_strap_full_process, file = paste0("data/genere/bott_strap_full_process.rds"))
}

  
  
boot_strap_full_process <- readRDS("data/genere/bott_strap_full_process.rds")


###############################################################################
# boot strap analysis 
# truncated weight 
boot_strap_pds_trunc <- lapply(boot_strap_full_process[6,],function(x){
  Trucature_pds_function(x$Weight,c(0,1,5,10,25,50)/100) %>%   summarise_all(list(~mean(.),
                                                                                  ~min(.),
                                                                                  ~max(.))) %>%
    t.df() %>% 
    separate(key, 
             into = c("truncations",
                      "fun"),
             sep = "_") %>% 
    pivot_wider(names_from = fun ,
                values_from = col_1)  %>% 
    rename_all(function(x) c("Truncations","Mean","Minimum","Maximum"))
}) %>% bind_rows() 
summarry_boot_strap_pds_trunc <- boot_strap_pds_trunc %>% 
  mutate(Truncations=factor(Truncations,
                            levels = c("(0; 1)",
                                       "(0.01; 0.99)",
                                       "(0.05; 0.95)",
                                       "(0.1; 0.9)",
                                       "(0.25; 0.75)",
                                       "(0.5; 0.5)")) ) %>% 
  group_by(Truncations) %>% 
  summarise_at(vars(-group_cols()),~paste0(
    sprintf("%.01f",round(mean(.),1)) ,
    "(",sprintf("%.01f",round(quantile(.,0.05),1)),"; ",
    sprintf("%.01f",round(quantile(.,0.95),1)), ")")
  )  

# Weight summary among bootstrap
descriptif_pds_boot <- lapply(boot_strap_full_process[6,],function(x){
  x %>% group_by(exposition)  %>% 
    summarise("mean_W" = mean(Weight) ,
              "min_W" = min(Weight),
              "max_W" = max(Weight),
              "mean_SW" = mean(SWeight),
              "min_SW" = min(SWeight),
              "max_SW" = max(SWeight),
              .groups = "drop"
    )
}) %>% bind_rows() 


summarydescriptif_pds_boot <- descriptif_pds_boot %>% 
  group_by(exposition) %>% 
  #mutate_all(list(~as.character())) %>% 
  summarise_at(vars(-group_cols()),~paste0(sprintf("%.01f",round(mean(.),1)) , "(",sprintf("%.01f",round(quantile(.,0.05),1)),"; ",sprintf("%.01f",round(quantile(.,0.95),1)), ")")) 


# Bootstrap comparaison of treatment effect estimation

boot_strap_coef_reg_non_pond <- apply(boot_strap_full_process,2,function(x){
  weight_fun <- x[2]$poids
  id_boot_strap <- x[5]$df_id_grp_tt %>% select(-INS_obs_categ)
  df_diff_simple <- id_boot_strap %>% left_join(donnee_boot , by = c("id_patient"))# %>% select(-id_patient)
  
  
  rbind(x[7]$reg_non_pond %>% as.data.frame() %>% rename("value" = 1) %>% rownames_to_column() %>% filter(str_detect(rowname,"INS_obs_categ")) %>% t.df(pivot = "rowname") %>% mutate(key = "reg_pond"),
        
        x[3]$reg_pond %>% as.data.frame() %>% rename("value" = 1) %>% rownames_to_column() %>% filter(str_detect(rowname,"INS_obs_categ")) %>% t.df(pivot = "rowname") %>% mutate(key = "reg_non_pond"),
        
        df_diff_simple %>% mutate(weight = weight_fun ) %>% 
          group_by(INS_obs_categ) %>% summarise(normal_mean = mean(SYM_echelleEpworth), IPTW = weighted.mean(x = SYM_echelleEpworth, w = weight)) %>% t.df("INS_obs_categ") %>% mutate_at(vars(-key),~(. - `4`)) %>% select(-`4`) %>% rename_at(vars(-key),~paste0("INS_obs_categ",.))
  )
  
}) %>% bind_rows() 





df_4meth <- boot_strap_coef_reg_non_pond %>% pivot_longer(-key) %>% group_by(key,name) %>% summarise(moy = mean(value),ci_sup = quantile(value,0.975),ci_inf = quantile(value,0.025) ,.groups = "drop") %>%
  mutate(key = str_replace(key,"normal_mean", "Simple mean"),
         key = str_replace(key,"reg_non_pond", "Simple regression"),
         key = str_replace(key,"reg_pond", "IPWRA")) %>%
  mutate(name = str_remove(name,"INS_obs_categ")) %>% 
  mutate(key = factor(key,levels = c("Simple mean","IPTW","Simple regression","IPWRA"))) 

##############################################################################



list_id_boot <- lapply(boot_strap_full_process[5,],function(x){
  
  x$id_patient[x$id_patient %notin% donnee_boot$id_patient]
  
  
})



##############################################################################
# Bootstrap Variables balance before and after weighting



nb_grp <- donnee_boot %>% select(all_of(name_expo)) %>% 
  unique() %>% 
  unlist(use.names = FALSE) %>% 
  sort


df_smd <- select_var$data_frame$df_analyse %>% select(all_of(c(Id,name_expo,Var_model_outcome)))#one_hot_fb(select_var$data_frame$df_analyse ,list_factor =  c(  "PPC_id_typeMasque"  ))
descriptif_pds_boot <- apply(boot_strap_full_process,2,function(i){

  weight_fun <- i[2]$poids
  id_boot_strap <- i[5]$df_id_grp_tt %>% select(-INS_obs_categ)
  
  
  
  df_func <- id_boot_strap %>% left_join(df_smd , by = c("id_patient"))# %>% select(-id_patient)
  
  mat_var <- df_func %>% select(-id_patient,-all_of(c(name_expo))) %>% as.matrix()
  
  
  
  df_max_func <- data.frame(sm_max_non_pondere = apply(smd(mat_var,df_func$INS_obs_categ),2, function(x)abs(mean(x))),
                       sm_max_pondere = apply(smd(mat_var,df_func$INS_obs_categ,weight_fun),2,function(x)abs(mean(x)))) #
})


df_max <- apply(simplify2array(lapply(descriptif_pds_boot,as.matrix)), 1:2, mean) %>% as.data.frame() %>% 
  rownames_to_column()  


##################################################################################
# Print paper result
vec_cut <- read_rds("data/genere/obs_cut2.rds")

# patient per group and our range 
for (i in seq_along(vec_cut[-1])) {
  cat("\\item", table(df_modele_pds$INS_obs_categ)[i]," ( ", round(table(df_modele_pds$INS_obs_categ)[i]*100/sum(table(df_modele_pds$INS_obs_categ)),1) , "\\% )"," patients with an adherence between ",roud_hour(vec_cut)[i], " and ", roud_hour(vec_cut)[i + 1] ,'\n')
}

groupe_obs_table_latex <-  sapply(seq_along(vec_cut[-1]),function(i) (paste0(vec_cut[i],"-", vec_cut[i + 1]," h")))




alpha_tableau_resume <- 0.05
var_instrumental_name  <- c("INS_obs_categ" = "adherence groups", "Vis_init_INS_IAH" = "apnea hypopnea index", "INS_PPC_effet_ind" = "number of ADR types under CPAP", "INS_tmp_entre_rdv" = "duration since diagnosis (year)")




table_rename <- one_hot_fb(select_var$data_frame$df_analyse ,list_factor =  c(  "PPC_id_typeMasque"  )) %>% 
  select(-Id) %>%  rename_variables(var_instrum_name = var_instrumental_name) 



# nb var model fin 

paste0("variables at diagnosis : " , paste0(str_remove_all(firstlow(table_rename$complete_table$Label[table_rename$complete_table$var_name %in% Var_model_outcome[str_detect(Var_model_outcome,
                                                                                                                                                                             "^Vis_init_|^Mesure_unique_")]]),"Diagnosis |\\(.+\\)"), collapse = ", "))

paste0("and the variables under CPAP treatment : " , paste0(str_replace(str_remove_all(firstlow(table_rename$complete_table$Label[table_rename$complete_table$var_name %in% Var_model_outcome[!str_detect(Var_model_outcome,"^Vis_init_|^Mesure_unique_")]]),"\\(.+\\)"),"ADR","adverse drug reaction"), collapse = ", "))

select_var$data_frame$df_analyse  %>% summarise(epworth = mean(SYM_echelleEpworth), 
                             sd = sd(SYM_echelleEpworth), 
                             base_epworth =  mean(Vis_init_SYM_echelleEpworth),
                             Sd_base = sd(Vis_init_SYM_echelleEpworth),
                             diff_epworth = median(Vis_init_SYM_echelleEpworth - SYM_echelleEpworth ), 
                             SD_diff = sd(Vis_init_SYM_echelleEpworth - SYM_echelleEpworth )
                            ) %>%
  arrondie_df(2)


select_var$data_frame$df_analyse %>%  group_by(INS_obs_categ) %>% 
  summarise(epworth = mean(SYM_echelleEpworth),
            SD = sd(SYM_echelleEpworth),
            base_epworth =  mean(Vis_init_SYM_echelleEpworth),
            Sd_base = sd(Vis_init_SYM_echelleEpworth),
            diff_epworth = mean(Vis_init_SYM_echelleEpworth - SYM_echelleEpworth ), 
            SD_diff = sd(Vis_init_SYM_echelleEpworth - SYM_echelleEpworth ),
  ) %>%
  arrondie_df(2)


quantile(select_var$data_frame$df_analyse$Vis_init_SYM_echelleEpworth,c(0.25,0.75))


# number of patient with severe osa 
sum(select_var$data_frame$df_analyse $Vis_init_INS_IAH >= 30)
round(mean(select_var$data_frame$df_analyse $Vis_init_INS_IAH >= 30) * 100,1)


# Final result 
df_4meth %>% mutate(ci = paste0("(",round(ci_inf,2),";",round(ci_sup,2),")"),moy = round(moy,2)) %>% select(-ci_sup ,- ci_inf )%>% pivot_wider(names_from = key,values_from = c(moy,ci))



##################################################################################
# Figure and table generation

# Netoyage du fichier tex
path_latex <- "figure_papier/table_iptw.tex"

sub_file_latex_preambule <- {"\\documentclass[../main.tex]{subfiles}
\\begin{document}"}



write(sub_file_latex_preambule, path_latex, append=FALSE)
# génération du code pour le papier inutile dans le rapport 
summarydescriptif_pds_boot  %>% select(exposition,ends_with("_W")) %>% rename(`Adherence group` = exposition,
                                                                              Mean = mean_W,
                                                                              Minimum = min_W,
                                                                              Maximum = max_W) %>% 
  mutate(`Adherence group` =  groupe_obs_table_latex) %>% 
  kable(format = "latex",
        align = "c",
        booktabs = TRUE,
        escape = FALSE,
        label = "Unstabilized_table",
        linesep = "",
        caption = "Distribution of unstabilized weights"
  ) %>%
  kable_styling(latex_options =  c( "striped","HOLD_position", "scale_down","repeat_header"))%>%
  add_footnote(c("Number are express in mean (5th percentile; 95th percentile) of bootstrap iterations"),
               notation ="symbol") %>% ecriture_fichier(path_latex)
write("\\clearpage", path_latex, append=TRUE)  


# troncature des poids 
summarry_boot_strap_pds_trunc %>%
  arrondie_df(1) %>% 
  kable(format = "latex",
        align = "c",
        booktabs = TRUE,
        escape = FALSE,
        label = "weight_truncations",
        linesep = "",
        caption = "Unstabilized weight truncations"
  ) %>%
  kable_styling(latex_options =  c( "striped","HOLD_position", "scale_down","repeat_header")) %>%   add_footnote(c("Number are express in mean ( 5th percentile; 95th percentile) of bootstrap iterations"),
                                                                                                                 notation ="symbol") %>% 
  ecriture_fichier(path_latex)
write("\\clearpage", path_latex, append=TRUE)     



# in paper we reduct number of variables in tables 
Variables_core_paper <- c("INS_obs_categ","SYM_echelleEpworth","Mesure_unique_CHA_id_sexe","Vis_init_CHA_age","Vis_init_CHA_BMI","Vis_init_INS_IAH","PPC_IRAH","Vis_init_FDR_fumeur",
                          "Vis_init_SYM_echelleDepression","SYM_echelleDepression","Vis_init_SYM_echellePichot","SYM_echellePichot",
                          "Vis_init_SYM_cephaleesMatinales","SYM_cephaleesMatinales","Vis_init_SYM_fatigueMatinale","SYM_fatigueMatinale",
                          "Mesure_unique_ATC_diabete","INS_PPC_effet_ind","INS_tmp_entre_rdv")
foot_note <- list(CPAP = list(var = c("INS_PPC_effet_ind","PPC_IRAH"),foot_notes =  "CPAP : Continuous Positive Airway Pressure"), 
                  ADR = list(var = c("INS_PPC_effet_ind"),foot_notes =  "ADR : adverse drug reaction"))


table__resume_rapport_non_rename <- select_var$data_frame$df_analyse %>%  
  select(any_of(Variables_core_paper)) 


foot_note_var_presence <- lapply(foot_note, function(x)
  x$foot_notes [any(x$var %in% colnames( table__resume_rapport_non_rename))]
) %>% unlist() 
# Tables avec les variables renommer pour jointure

table__resume_rapport <- table__resume_rapport_non_rename %>%  rename_variables(var_instrum_name = var_instrumental_name)  %>% .$table_rename
table_resume_rapport_latex_ss_escape <- table_resume_latex(df_one_hot_encode_fonc =  table__resume_rapport,
                                                           name_expo_fonc = "adherence groups",
                                                           nom_grp  = "adherence grp",
                                                           p_val = TRUE,
                                                           alpha = alpha_tableau_resume)

write("\\clearpage", path_latex, append=TRUE)  
header_landscape <- "\\newgeometry{a4paper,left=1in,right=1in,top=1in,bottom=1in,nohead}" # Il faut changer les marges des pages avant de landscape pour éviter un saut de pages indésiré
write(header_landscape, path_latex, append=TRUE)
table_resume_rapport_latex_ss_escape %>%
  rename_all(~c("all adherence grp",paste0(groupe_obs_table_latex, " (",seq_along(groupe_obs_table_latex),")" ))) %>% 
  # Marquage modèle de poids
  rownames_to_column() %>% 
  column_to_rownames() %>% 
  kable(format = "latex",
        align = "c",
        booktabs = TRUE,
        escape = FALSE,
        label = "Table_resume_pop",
        linesep = "",
        caption = "Table of variables for each adherence group"
  ) %>%
  kable_styling(latex_options =  c( "striped","HOLD_position",
                                    "repeat_header"),
  ) %>%
  add_footnote(c("t-test was performed for the quantitative variables and a Pearson's Chi squared test for the categorical variables after application of a Bonferroni correction for multiple testing.",
                 "1,2,3,4 numbers in subscript refers to columns statistically different at the 5% threshold. e.g. 1 means that there is a statistically significant difference between the 0-4h adherence group and the 0-4h adherence group (1) for the variable in question.",
                 foot_note_var_presence
                 
                 #"* after index corresponding to variables with too few people to perform test"
  ),threeparttable= TRUE,
  notation ="symbol") %>% 
  landscape() %>% 
  ecriture_fichier(path_latex)
end_landscape <- "\\restoregeometry % Restore the global document page margins"
write(end_landscape, path_latex, append=TRUE)
write("\\clearpage", path_latex, append=TRUE)  




df_box_plot_fig2 <- select_var$data_frame$df_analyse %>%  
  rename_variables(var_instrum_name = var_instrumental_name)  %>% 
  .$table_rename %>% 
  select(`adherence groups`,`diagnosis epworth sleepiness scale`,`epworth sleepiness scale`) %>% pivot_longer(-`adherence groups`) %>% mutate(`adherence groups` = as.factor(as.character(`adherence groups`)), name = factor(name,levels = c("diagnosis epworth sleepiness scale","epworth sleepiness scale")))  


ggplot(data = df_box_plot_fig2, aes(x = name ,y = value, fill = `adherence groups` )) + geom_boxplot() +
  
  ggsci::scale_fill_lancet(labels =groupe_obs_table_latex)  + labs(y = "Epworth score",
                                    x = paste0(""
                                    ),
                                    fill = "CPAP adherence group") 
ggsave("figure_papier/figure_2.jpeg")


# SD bootstrap 

df_smd_plot <- df_max %>%  filter(rowname %in% Var_model_outcome)
 # filter(rowname %in% Variables_core_paper)
df_plot_max <- table_rename$complete_table %>% 
  select(var_name,Label)  %>% 
  right_join(df_smd_plot ,by = c("var_name" = "rowname")) %>% 
  select(-var_name) %>% mutate(Label = factor(Label,levels = .$Label[order(abs(df_smd_plot$sm_max_non_pondere))])) %>% 
  pivot_longer(-Label) %>% 
  mutate(name = str_replace_all(name,"sm_max_non_pondere","SMD without Weighting"),`Weighting` = str_replace_all(name,"sm_max_pondere","SMD with Weighting")) %>% select(-name)

ggplot(data = df_plot_max,
       mapping = aes(x = Label, y = value, group = Weighting, color = Weighting)) +
  # geom_line() +
  geom_point() +
  
  geom_hline(yintercept = 0.1, color = "black", size = 0.1) +
  scale_y_continuous(breaks = sort(c(seq(0, 1, length.out=5), 0.1))) +
  coord_flip() +
  theme_bw() + theme(legend.key = element_blank()) + 
  ggplot2::labs(title = "Standardized mean differences \n before and after weighting",y = "SMD",
                x = "") 

ggsave("figure_papier/figure_3.jpeg")



pd <- position_dodge(0.1) # move them .05 to the left and right
ggplot(df_4meth %>% mutate(name = groupe_obs_table_latex[as.numeric(name)]) ,
       aes(x=name,y = moy )) + 
  geom_errorbar(aes(ymin=ci_inf , ymax=ci_sup ), colour="black", width=.1, position=pd ) +
  geom_point(position=pd , size=3) + facet_wrap(~ key) +
  labs(x = "Adherence groups",y = "Variation of Epworth score")+
  theme(strip.text  = element_text(size = 12 , face = "bold" ))
ggsave("figure_papier/figure_4.jpeg")




write("\\end{document}", path_latex, append=TRUE)


