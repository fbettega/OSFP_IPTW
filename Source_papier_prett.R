#########################################################################################
# retire les colonnes iodentiques
distinc_col <- function(df_func,return_list_col_supr = TRUE){
  column_unique_df <- df_func %>% t.df() %>% distinct_at(vars(-key),.keep_all = TRUE)
  df_col_unique <- df_func %>% select(all_of(column_unique_df$key))
  
  if (return_list_col_supr) {
    col_supr <- colnames(df_func)[colnames(df_func) %notin% column_unique_df$key]
    df_col_supr <- df_func %>% select(all_of(col_supr))
    
    
    liste_col_supr <- column_comparator(df_col_supr,df_col_unique)
    
    res <- list(df = df_col_unique, colonne_suprime =  liste_col_supr)
  } else {res <- df_col_unique}
  return(res)
}

#########################################################################################
# liste les  colonnes iodentiques
column_comparator <- function(col_a_compar,df_compar_func){ # recherche les collonnes identique et les regroupe sous la forme de nom de collones supprimé : nom de colonnes restant matchant 
  # Attention une colonnes restante étant égale a 2 coçlonnes supprimé génère donc deux liste
  liste_colum_identique <-   lapply(seq_along(df_compar_func),function(x) { 
    col_a_compar_temp <- df_compar_func[,x] %>% unlist( use.names = FALSE)
    column_compared <- (col_a_compar_temp == col_a_compar ) | 
      (is.na(col_a_compar_temp) & is.na(col_a_compar ))
    matching_col <- c(names(df_compar_func[,x]),names(which(apply(column_compared,2,all))))
  })
  liste_colum_identique  <- lapply(liste_colum_identique, function(x) x[length(x) > 1]) %>% compact()
  return(liste_colum_identique)
}
##########################################################################################

# count number of NA per variables avec les individus groupé par une variables
compte_na_par_var_par_grp <- function(df,group_col,colonnes){
  df_NA_var_fonc <- df %>% select(all_of(group_col),all_of(colonnes)) %>% setDT 
  nb_val_manq_par_var <- df_NA_var_fonc %>% 
    .[, lapply(.SD,  function(x) sum(is.na(x))), group_col] %>% 
    select(-all_of(group_col)) %>% 
    gather(name,presence_na) %>%  # reshape datset
    count(name, presence_na) %>%  # count combinations
    pivot_wider(names_from = name,
                values_from = n,
                values_fill = list(n = 0))
  return(nb_val_manq_par_var)
}

##########################################################################################
# Fonction usuful when i need to hard code choose beetween to colinear variables check that i choose all nedd varaibles and no more 
select_manuel_verif <- function(fuc_list_manu,choose_elem){
  # cette fonction ne sert que de vérifications 
  if(length(fuc_list_manu) != length(choose_elem)) {
    print("number of element choose differ with list length")
    if(length(fuc_list_manu) > length(choose_elem)) {
      stop("not enougth manual choosen elements")
    } else { stop("to many manual choosen elements") }
  } else if(length(fuc_list_manu) == length(choose_elem)){ 
    # vérif que tout les éléments corresponde bien a ceux dans lequels manuellement choisir
    bool_in_list <- sapply(seq_along(choose_elem), function(x) choose_elem[x] %in% fuc_list_manu[[x]])
    if (all(bool_in_list)) {
      # selection de tout les éléements qui ne sont pas celui manuel
      res <- sapply(seq_along(choose_elem), function(x) fuc_list_manu[[x]][fuc_list_manu[[x]]  %notin% choose_elem[x]]) %>% 
        compact %>% unlist
    } else {
      stop(paste0("choose elements ",choose_elem[!bool_in_list], " not in list"))
    } 
  }
  return(res)
}
##########################################################################################

# transposé dataframe
t.df <- function(df,pivot=NULL){
  if (is.null(pivot)){
    pivot <- "row_id"
    df <- df %>% mutate(row_id=paste0("col_",1:nrow(df) ))
  }
  res <- df %>% pivot_longer(cols = -!!pivot,"key","value") %>%
    pivot_wider(names_from = !!pivot,values_from = value)
  return(res) 
}
##########################################################################################
`%notin%` <- Negate(`%in%`)

##########################################################################################

##Imputation si changment de données##
####################################
#####################################################
impute_si_changement <- function(data_frame_a_verif,path,reimputation = FALSE){
  df_impute_ex <- tryCatch(read_rds(path),  
                           error = function(e) data.frame())  
  colanmes_manq_func <- data_frame_a_verif %>%
    select_if(data_frame_a_verif %>% 
                summarise_all(list(~sum(is.na(.)))) != 0) %>% 
    colnames()

  if ((all(sort(colnames(data_frame_a_verif)) == sort(colnames(df_impute_ex)) ) & 
       (nrow(data_frame_a_verif) ==  nrow(df_impute_ex))) & !reimputation) {
    res <- df_impute_ex
  } else if (any(sort(colnames(data_frame_a_verif)) != sort(colnames(df_impute_ex))) | reimputation | nrow(df_impute_ex) == 0 | nrow(df_impute_ex) != nrow(data_frame_a_verif))  {
    cat("nécessité de réimputer les données 
      ça va être long + ou - une heure")
      imputed_Data_func <- mice::mice(data_frame_a_verif,
                                      m = 10,
                                      maxit = 50,
                                      method = 'pmm',
                                      printFlag = FALSE)
    saveRDS(imputed_Data_func, file = paste0("data/genere/imputation_object.rds"))
    if (!is.null(imputed_Data_func$loggedEvents)) { 
      print("imputation avec warning")
      print(imputed_Data_func$loggedEvents)}
    res <- mice::complete(imputed_Data_func)
    saveRDS(res, file = paste0("data/genere/data_impute.rds"))
  } else {stop("Condition non remplie problème fonction")}
  return(res)
}

#########################################################################################

# selection des variables en univarié


# a améliorer pour rendre généraliste
model_selection <- function(data_frame ,name_expo_fonc , name_outcome_fonc, seuil_pval = c(0.2,0.2),  #vecteur de seuil  de pval considéré comme significatifs first outcome second weight
                            rerun = TRUE) {
  
  
  if (rerun){
    
    df_model_select_var <- data_frame %>% group_by_at(name_expo_fonc) %>% slice_sample(prop= 0.2) %>% ungroup %>% as.data.frame()
    # other patients
    df_analyse_final <- data_frame %>% anti_join(df_model_select_var , by = "id_patient") 
    
    df_select_var_base <- df_model_select_var %>% select(-id_patient)
    var_to_select_from <- df_select_var_base %>% select(-{{name_expo_fonc}},{{name_outcome_fonc}}) %>% colnames()
    
    model_outcome <- df_select_var_base %>% 
      select(-all_of(name_expo_fonc)) %>% 
      gather(measure, value, -all_of(name_outcome_fonc)) %>%
      mutate(value = as.numeric(value)) %>%
      group_by(measure) %>% 
      nest() %>%
      ungroup() %>% 
      mutate(fit = map(data, ~ glm(paste0(name_outcome_fonc,"~ value"), data = .x), family = gaussian(link = "identity"),na.action = na.omit),
             tidied = map(fit, broom::tidy)) %>%
      unnest(tidied) %>% 
      filter(term != "(Intercept)") %>% # on elnlève l'intercept
      select(-data,-fit,-term)  # retire les donnée
    
    
    
    
    
    P_val_out_come_colnames <-  model_outcome %>%
      select(measure,contains("p.value")) %>% 
      filter(p.value < seuil_pval[1]) %>% 
      select(measure) %>% 
      unlist(use.names = FALSE) %>% 
      sort
    
    
    Var_model_outcome <- var_to_select_from[var_to_select_from %in% P_val_out_come_colnames]
    res <- list(data_frame = list(df_analyse = df_analyse_final, 
                                  df_var_selec = df_model_select_var ),
                var_select = list(model_outcome = Var_model_outcome))
    
    saveRDS(res, file = paste0("data/genere/select_var.rds"))
    
    
  } else {
    res <- readRDS("data/genere/select_var.rds")
  }
  
  
  return(res)
}
##########################################################################################
# fonction maison pour rechercher le type des données 

## objectifs suivant différence entre int et float
typages_function <- function(df,nb_moda_max_fact = NULL ){
  if (is.null(nb_moda_max_fact)) {
    df %>% 
      summarise_all(list(~n_distinct(na.omit(.)))) %>%
      t.df() %>% filter(col_1 > 2) %>% 
      arrange(col_1) %>% 
      print
    stop('Si la liste des facteurs n\'est pas fournis le nombre de modalité à partir
           duquel un facteur doit etre considéré comme un numéric avec `nb_moda_max_fact = `,
           \n pour vous aidez dans le choix du nombre de modalité la liste des variables
           avec plus de deux modalité différente est présenté au  dessus')}
  else {
    temp_moda_par_var <- df %>% summarise_all(list(~n_distinct(na.omit(.)))) %>% 
      t.df() %>% 
      mutate(binaire = col_1 == 2, numeric = col_1 >= nb_moda_max_fact,
             multinomial = (col_1 < nb_moda_max_fact & col_1 > 2)) %>% 
      arrange(col_1) 
    list_factor <- temp_moda_par_var %>% 
      filter(multinomial) %>%
      select(key) %>% 
      unlist(use.names = FALSE)
    liste_booleen <- temp_moda_par_var %>% 
      filter(binaire) %>%
      select(key) %>% 
      unlist(use.names = FALSE)
    liste_numeric <- temp_moda_par_var %>% 
      filter(numeric) %>%
      select(key) %>% 
      unlist(use.names = FALSE)
    res <- list(colonne_type = list(facteur = list_factor,
                                    booleen = liste_booleen,
                                    numerique = liste_numeric),
                data_frame_tot = temp_moda_par_var)
  }
  return(res)
}
##########################################################################################
# fonction maison pour le one hot encoding

one_hot_fb <- function(df, nb_moda_max_fact = NULL, list_factor = NULL){
  if (is.null(list_factor)) {
    list_factor <- typages_function(df, nb_moda_max_fact)$colonne_type$facteur
    cat("Le nombre maximum de modalité par facteur est de ",
        nb_moda_max_fact,
        "\n pour supprimer ce warning utilisez `list_factor = ", "c( ",
        paste0("\"",list_factor,"\"",collapse = " , "),
        " )` \n au lieu de `nb_moda_max_fact = ", nb_moda_max_fact,"`")
  }
  
  
  if (all(is.na(list_factor))) { # cas ou il n'y aucun facteur
    res <- df
  } else {
    df <- df %>% mutate_at( list_factor,as.factor)
    dmy <- dummyVars(paste0(" ~ ", paste0(list_factor,collapse = " + ")), data = df)
    trsf <- data.frame(predict(dmy, newdata = df))
    res <- df %>% select(-all_of(list_factor)) %>% cbind(trsf) # ajout all_of retirer si bug 
  }
  return(res)
  # reste a ajouter une partie qui renomme les varibles mieux 
}


# fonction servant à la troncature des poids
#########################################################################################
fun_trunc <- function(x,.probs) {
  pmin(pmax(x, quantile(x, probs = .probs)), 
       quantile(x, probs = 1-.probs))}



##function Truncature des poids
## prend un vecteur de poids et le trunc selon un vecteur de percentile
#########################################################################################
Trucature_pds_function <- function(vect_poids,percentile_tronc) {
  
  res <-  map(percentile_tronc,function(x) fun_trunc(vect_poids,x)) %>% 
    as.data.frame() %>% 
    structure(names = paste0("(",
                             as.character(percentile_tronc),
                             "; ",
                             as.character(1 - percentile_tronc),
                             ")" ) )
  return(res)
}

##function générant la table de base du modèle de poids pré insertion des résultat##
#########################################################################################

Table_base_result_weighting <- function(donne,expo,Id_data,PS_table){
  exposure <- donne[,as.character(expo)]
  res <- data.frame(ID = Id_data,exposition = exposure , PS = NA) %>% 
    setNames(c("ID","exposition", "PS")) 
  res <- res %>% 
    group_by(exposition) %>%
    mutate(n = n(),numerator = n / nrow(.)) %>% select(-n) # numerator for stabilisation
  
  # Probability of being expose to observed outcomes
  res$PS <- sapply(1:nrow(PS_table), 
                   function(x) PS_table[x,as.character(as.numeric(res$exposition[x]))]) 
  
  res <- res %>% mutate(SWeight = numerator/PS, Weight = 1/PS) 
  return(res)
}

##Ponderation via regression multinomial
#########################################################################################
# Premiere fonction du genre le but est de prendre une matrice de donne 
#et retourne une probabilité d'appartenance à chacun des groupes
# ainsi que le modele en question 
# destiné à etre utilsier deans calcules pds_stage
multinomial_IPTW <- function(donne,expo,covar,Id_data){
  # mutlinomial regression to compute probability to belong in each group knowing confunders
  mod1 <-   multinom(formula =   eval(parse(text =
                                              paste("as.numeric(",expo,")", # Exposure
                                                    paste0("~ ",paste0(covar, collapse = " + ")), # all confunders additive models but can easily be change with another kind of model
                                                    sep = "")
  )),data = donne, na.action = na.fail , # vérif d'erreur
  trace = FALSE)
  
  # Predict using previously train models 
  proba_tps_inv <- cbind(ID = Id_data,predict(mod1, type = "probs") %>% 
                           as.data.frame() %>% 
                           rename("1" = names(.)[1]))
  return(list(matrice_proba = proba_tps_inv , model_ponderation = mod1))
}


##function principale du stage##
#########################################################################################
# Note perso le temps d'execution et très majoritairent liée à la regression  > 90 %
##function principale du stage##
#########################################################################################
# Note perso le temps d'execution et très majoritairent liée à la regression  > 90 %
calcule_pds_stage <- function(donne, 
                              weighting_function,
                              supplementary_weighting_function_params = NULL, # add here a name list of named parameters
                              Id_data = NULL,
                              expo,
                              covar,
                              out_come,
                              percentile_tronc = c(0,1,5,10,25,50)/100 ,
                              talbe_and_plot = TRUE){
  # fetch function parameters
  tempcall <- match.call()
  
  
  # TODO Séparer le calcules des poids par IPTW dans une autre fonction
  # TODO globalement cette fonction devrait devenir uniquement des appel à fonction de type 
  # TODO prevoir des tests pour les différentes fonctions 
  # Calcules de poids avec méthode X puis on déroules
  # La méthode X doit prendre en entrée une table avec 
  # Probablement à réflechir est a delete
  
  exposure <- donne[,as.character(tempcall$expo)]
  ID <- donne[,as.character(tempcall$Id_data)]
  
  weighting_function_params <- list(donne = donne,
                                    expo = tempcall$expo,
                                    covar = covar ,
                                    Id_data = ID)
  if (!is.null(supplementary_weighting_function_params)){
    weighting_function_params <- append(weighting_function_params,supplementary_weighting_function_params)
  }

  
  ponderation <- do.call(weighting_function,weighting_function_params)

  
  proba_tps_inv <- ponderation$matrice_proba %>% as.data.frame()
  
  
  
  
  
  res <- Table_base_result_weighting(donne,tempcall$expo,ID,proba_tps_inv)
  # troncature des poids 
  poid_trunc_df <- Trucature_pds_function(res$Weight,percentile_tronc)
  
  poid_trunc_stab_df <- Trucature_pds_function(res$SWeight,percentile_tronc) 
  
  ### graphique et autres table explicative
  

  return(list(df = res,
              res_intermediaire = list(
                poids = list(
                  poids_tronc = list(poids_trunc = poid_trunc_df,
                                     poids_trunc_stab = poid_trunc_stab_df)),
                regression_modele_pds = list(
                  regression_temps_ind = ponderation$model_ponderation,
                  data_frame_coef = proba_tps_inv)
                )
  )
  )
}

#########################################################################################

#########################################################################################
# Boot strap other a function
boot_strap_fun <- function(df,fonction, nb_iter,params_fun){
  replicate(nb_iter,fonction(df,sample(1:nrow(df), size = nrow(df), replace = TRUE),params_fun))
}

###############################################################################
# RESULT FUNCTION
##########################################################################################
# convertie les heures  avec virgule en heure et minute
roud_hour <- function(decimal_hour){
  heure <- floor(decimal_hour)
  minutes <- round(60 * (decimal_hour - floor(decimal_hour)), 0)
  res <- sprintf("%02d h %02d min", heure, minutes)
  return(res)
}

#
################################################################################
#to lower first char
firstlow <- function(x) {
  substr(x, 1, 1) <- tolower(substr(x, 1, 1))
  x
}

##########################################################################################
# fait des arrondie + notation scientifique pour les nombres a virugules dans les DF en 
# conservant les integer telquel
arrondie_df <- function(df_func,digit_to_round = 3){
  res <-   df_func %>%  
    rownames_to_column() %>% 
    mutate_if(is.numeric,
              # ancien avec probablement un problème dans l'ordre des ifelse
              # ~ifelse(.%%1==0,as.character(round(.,0)),ifelse((. > 10^3| 1/abs(.) > 10^3 ),
              #         formatC(., format = "e", digits = 1), # tentative de gestion des nombres 
              #         as.character(round(.,3))
              #         )
              # )
              ~ifelse(.%%1==0 & . < 10^3,as.character(round(.,0)),ifelse((. > 10^3| 1/abs(.) > 10^3 ),
                                                                         formatC(., format = "e", digits = 1), # tentative de gestion des nombres 
                                                                         as.character(round(.,digit_to_round))
              )
              )
    ) %>% 
    column_to_rownames() 
  return(res)
}



##########################################################################################
# Fonction qui sert a rename les colonnes depuis un fichier
rename_variables <- function(table_func, var_instrum_name = NULL, path_to_var_lab = "data/List_of_variables.xls")  {
  # récup fichier avec les noms 
  label_var <- readxl::read_excel(path_to_var_lab) %>% select(Variable,Label) %>% mutate(Label =   unlist(lapply(strsplit(Label, " "),function(x) 
    paste0(ifelse(str_count(x,"[A-Z]") == 1,firstlow(x),x),collapse = " ")
  )))

  lapply(strsplit(label_var$Label, " "),function(x) 
    paste0(ifelse(str_count(x,"[A-Z]") == 1,firstlow(x),x),collapse = " ")
  )
  #On ajoute le fait que c'est le delta epworth
  # script variable name
  actual_name <- data.frame(var_name = colnames(table_func)) %>% 
    mutate(suffix_factor = str_extract(var_name,"\\.+\\d$"), # récupération du cas ou il y a eu one hot encode
           preffix_all = str_extract(var_name,"^Vis_init_|^Mesure_unique_"), # récupération de mes 3 preffix crée
           var_ins = str_extract(var_name,"INS_"), # variables instrumentale que tu devras nommer toi meme
           real_name = str_remove_all(var_name,"INS_|^Vis_init_|^Mesure_unique_|\\.+\\d$") # on enlève tout ce qu'on a détecté
    )
  join_table <- actual_name %>% left_join(label_var, by = c("real_name" = "Variable")) # joionture avec les noms
  # gestions des variables instrumentales
  name_need_supply <- join_table %>% 
    filter(!is.na(var_ins)) %>% 
    select(var_name) %>% 
    distinct() %>% 
    unlist(use.names = FALSE)
  if(is.null(var_instrum_name)){
    cat("You must supply name for variable you create \n variables you must name are \n")
    print(name_need_supply) 
    cat("\n for that ` var_instrum_name = c(\"variable_name\" = \"new name\")`\n")
    stop()
  } else if (length(name_need_supply) != length(var_instrum_name)) { 
    # heler pour les variables manquante ou en trop
    out_temp <- if_else(condition = length(name_need_supply) > length(var_instrum_name),
                        true = paste0("Not enought name provide, missing name for \n ",
                                      paste(name_need_supply[name_need_supply%notin%names(var_instrum_name)], collapse = " ,")),
                        false = paste0("to many name provide, don't need name for \n ",
                                       paste(var_instrum_name[names(var_instrum_name)%notin%name_need_supply], collapse = " ,"))
    )
    cat(out_temp)
    stop()
  } else {
    instrum_name <- data.frame(var_name = names(var_instrum_name) , Label = var_instrum_name) 
    complete_name_table <- join_table %>% 
      left_join(instrum_name,by = "var_name") %>% 
      mutate(Label = coalesce(Label.x, Label.y)) %>% 
      select(-contains("Label.")) %>% 
      mutate(Label = ifelse(!is.na(suffix_factor),
                            paste0(Label, suffix_factor),
                            Label),# on remet les indices de facteur
             Label = ifelse((!is.na(preffix_all)&preffix_all == "Vis_init_"),
                            paste0("diagnosis ",Label),
                            Label)
      )
    short_name_table <- complete_name_table$var_name
    names(short_name_table) <- complete_name_table$Label 
    res <- table_func %>% rename(all_of(short_name_table))
  }
  return(list(table_rename = res,
              complete_table = complete_name_table ))
}


ci_fb <- function(con_int,x)   {
  #qnorm(con_int)*(sd(x)/sqrt(length(x)))
  qt(con_int,df=length(x)-1)*(sd(x)/sqrt(length(x)))
}



#########################################################################################
# Fonction rapide pour écrire dans un fichier
ecriture_fichier <- function(table_latex,path_to_file) {
  write(table_latex, path_to_file,append = TRUE)
}




##########################################################################################
# table descriptive des variables

# df_one_hot_encode_fonc généré avec one_hot_fb
# Version prévue pour échapper les caracter latex
table_resume_latex <- function(df_one_hot_encode_fonc,name_expo_fonc,nom_grp = "clusters", p_val = FALSE, arrondie = TRUE, alpha = 0.05,ponderation = NULL) {
  
  # Typage des variables a réusmer en booleen ou numeric
  # Car normalement df pré_one hot encode
  
  table_bool_var <-   df_one_hot_encode_fonc %>% 
    select(-all_of(name_expo_fonc)) %>% 
    summarise_all(list(~n_distinct(., na.rm = TRUE))) %>% 
    t.df %>% 
    rename(booleen = col_1) %>% 
    {temp_verif_bool <<- .} %>% 
    mutate(booleen = booleen <= 2)
  
  if (any(temp_verif_bool$booleen < 2)) {
    print(temp_verif_bool$key[temp_verif_bool$booleen < 2])
    stop("moins de deux valeurs distinct pour une variables")}
  
  
  if (!is.null(ponderation)){
    df_one_hot_encode_fonc  <- df_one_hot_encode_fonc %>% mutate_at(vars(-name_expo_fonc),funs(.*ponderation)) #mutate_at(vars(-name_expo_fonc),funs(if(is_whole(.)){round(.*ponderation,0)} else{.*ponderation})) 
  }
  
  
  name_all_grp <- paste0("all ",nom_grp)
  
  all_cluster_descript_var <- df_one_hot_encode_fonc %>% 
    recup_var_table_res(name_expo_fonc) %>% 
    #{ifelse(arrondie ,arrondie_df(.), . )} %>% print %>% 
    
    mutate(nb_NA = ifelse(nb_NA == 0,"",paste0("NA:", nb_NA )))  
  
  nb_grp <- df_one_hot_encode_fonc %>% select(all_of(name_expo_fonc)) %>% 
    unique() %>% 
    unlist(use.names = FALSE) %>% 
    sort
  
  group_cluster_descript_var <- lapply(nb_grp, function(x) {
    col_name <- paste0(nom_grp,"_",x)
    res <- df_one_hot_encode_fonc  %>%
      filter(!!sym(name_expo_fonc) == x) %>% 
      recup_var_table_res(name_expo_fonc)  %>% 
      mutate(nb_NA = ifelse(nb_NA == 0,"",paste0("NA:", nb_NA ))) %>% 
      inner_join(table_bool_var,by = "key") %>% 
      mutate({{col_name}} := ifelse( booleen
                                     , paste0(round(n,0),"(",round(pourcent,1), "%)", nb_NA),
                                     paste0(round(moy,1) ,"(",
                                            round(std,1), ")", nb_NA))) %>%
      select(all_of(col_name))
    return(res)
  }
  )
  
  table_res <- all_cluster_descript_var %>%  
    inner_join(table_bool_var,by = "key") %>% 
    mutate( {{name_all_grp}} := ifelse(booleen
                                       , paste0(round(n,0),"(",round(pourcent,1), "%)",  nb_NA),
                                       paste0(round(moy,1) ,"(",
                                              round(std,1), ")", nb_NA) ) )  %>% 
    select(key,all_of(name_all_grp)) %>% 
    cbind(bind_cols(group_cluster_descript_var)) %>% 
    data.frame(., row.names = 1) 
  # rename_at(vars(contains(".grp")),funs(str_replace(.,"\\."," ")))
  
  
  table_res <- table_res %>% rename_all(list(~str_replace_all(.,"\\."," ")))
  # si choix de calculer les p-val 
  if (p_val) {
    # prevoir un groupe de plus pour le toutes les catégorie
    nb_group_pval <-  as.character(c(as.numeric(nb_grp),max(as.numeric(nb_grp)) + 1 ))
    # création de toutes les combinaisons de groupe a tester
    combin_grp <- nb_group_pval %>%  combn(2)
    # Création du groupe supplémentaire tout les groupes en dupliquant la dataframe avec un  groupes de plus
    df_pval <- df_one_hot_encode_fonc %>% 
      mutate(!!sym(name_expo_fonc) := max(as.numeric(nb_group_pval))) %>% 
      rbind(df_one_hot_encode_fonc)
    
    non_boolean_var <- table_bool_var %>% filter(!booleen) %>% select(key) %>% unlist(use.names = FALSE)
    boolean_var <- table_bool_var %>% filter(booleen) %>% select(key) %>% unlist(use.names = FALSE)
    # création  de la table avec p-value pour chaque combin et rename chaque colonnes a_b 
    combin_ttest_pval <- apply(combin_grp, 2, function(x)
      df_pval %>%
        select(sym(name_expo_fonc),all_of(non_boolean_var)) %>% 
        #summarise_at(vars(-(sym(name_expo_fonc))),list(~t.test(.[!!sym(name_expo_fonc) == x[1]], .[!!sym(name_expo_fonc)  ==  x[2]])$p.value)) %>% 
        summarise_at(vars(-(sym(name_expo_fonc))),list(~t.test(.[!!sym(name_expo_fonc) == x[1]], .[!!sym(name_expo_fonc)  ==  x[2]])$p.value)) %>% 
        t.df %>% 
        rename_at("col_1",list( ~paste0(x[1],"_",x[2])))
    )
    if (length(boolean_var) != 0){
      combin_chisq_pval <- apply(combin_grp, 2, function(x)
        df_pval %>%
          select(sym(name_expo_fonc),all_of(boolean_var)) %>% 
          dplyr::summarise_at(vars(-sym(name_expo_fonc)),
                              list(~ifelse(sum(.[!!sym(name_expo_fonc) == x[1]],na.rm = TRUE) < 8| sum(.[!!sym(name_expo_fonc)  ==  x[2]],na.rm = TRUE) < 8,
                                           NA,
                                           prop.test(
                                             x = c(sum(.[!!sym(name_expo_fonc) == x[1]],na.rm = TRUE), sum(.[!!sym(name_expo_fonc)  ==  x[2]],na.rm = TRUE)), # compute number of success
                                             n = c(sum(!is.na(.[!!sym(name_expo_fonc) == x[1]])), sum(!is.na(.[!!sym(name_expo_fonc)  ==  x[2]])))
                                           )$p.value))
          ) %>% 
          t.df %>% 
          rename_at("col_1",list( ~paste0(x[1],"_",x[2])))
      )
      combin_total_pval <- mapply(rbind,combin_chisq_pval,combin_ttest_pval,SIMPLIFY=FALSE)
      
    } else combin_total_pval <-combin_ttest_pval
    
    
    
    
    # transformation de la p-value en booléen en avec comme seuil le alpha définis en appliquant une correction de bonneferonni en considérant le nombre total de test nb ligne X nb collones
    result_pval <- bind_cols(combin_total_pval) %>%
      rename(key = key...1) %>% 
      select(key,contains("_")) %>%
      mutate_at(vars(-key),list(~(. < (alpha / ncol(combin_grp))
      ))
      ) %>%  # hypothèse et correction de bonneferonnie 
      mutate_at(vars(-key), function(x) {
        x_var <- rlang::enquo(x)
        ifelse(x , rlang::quo_name(x_var), "non") # remplacement des p-val non signif par une chaine spécifique
      }) %>% 
      mutate_at(vars(-key), function(x) {
        x_var <- rlang::enquo(x)
        ifelse(is.na(x) , paste0(rlang::quo_name(x_var),"*"),x) # remplacement des p-val non signif par une chaine spécifique
      })
    
    # REcherche avec une simili boucle des p-val signif pour chaque colonnnes 
    # on en lève la chaine spécifique de non corrélation 
    df_pval_final <- lapply(nb_group_pval, function(x) {
      result_pval %>% select(key,contains(x)) %>% 
        mutate_at(vars(contains(x)),list(~str_remove_all(.,paste(c("_",x), collapse = "|")))) %>% 
        unite(!!sym(paste0(nom_grp,"_",x)) ,contains(x),sep = ",")
    }
    ) %>% bind_cols() %>% 
      rename(key = key...1) %>% 
      select(key,contains(nom_grp)) %>% 
      rename(!!sym(name_all_grp) := paste0(nom_grp,"_",max(as.numeric(nb_group_pval)))) %>%   
      mutate_all(list(~str_remove_all(.,"non,|,non"))) %>%  
      mutate_all(list(~str_remove_all(.,"non"))) %>% 
      mutate_all(list(~str_remove_all(.,paste0(",",as.character(length(nb_grp) +1),"|",as.character(length(nb_grp) +1)) ))) # nouveau symbole pour all adherence grou^p
    
    
    
    if(df_pval_final %>%  transmute_at(vars(-key),list(~str_detect(.,"non"))) %>% as.matrix() %>% any) {
      stop("il reste des p-val non traité")}
    
    # Gestion des tables latex pour que les différences statisquement significative soit en subscript
    # en échappant les underscore
    table_res_pval <-  table_res %>% 
      rownames_to_column() %>%
      pivot_longer(-rowname,values_to = "valeur") %>%
      inner_join((df_pval_final %>% pivot_longer(-key,values_to = "pvalue") ), by = c("rowname" = "key", "name" = "name")) %>% 
      mutate(combin = paste0(valeur,"\\textsubscript{",pvalue, "}")) %>% 
      select(rowname,name,combin) %>% 
      pivot_wider(names_from = name, values_from = combin) %>% 
      column_to_rownames()   
    
    table_res_pval <- table_res_pval %>% 
      rownames_to_column() %>% 
      mutate_at(vars(-rowname),list(~str_replace_all(.,"%","\\\\%"))) %>% 
      column_to_rownames() %>% 
      #mutate_all(funs(str_replace_all(.,"%","\\\\%"))) %>% 
      select(all_of(name_all_grp),sort(tidyselect::peek_vars()))
    
    
    rownames(table_res_pval) <- str_replace_all(rownames(table_res_pval),"_","\\\\_")
    colnames(table_res_pval) <- str_replace_all(colnames(table_res_pval),"_","\\\\_")
    
    
  } else {table_res_pval <- table_res %>% 
    select(all_of(name_all_grp),sort(tidyselect::peek_vars()))  }
  tbl_nb_ind_grp <-  table(df_one_hot_encode_fonc[,name_expo_fonc])
  nb_pat_par_grp <- c(nrow(df_one_hot_encode_fonc),
                      tbl_nb_ind_grp[order(as.numeric(names(tbl_nb_ind_grp)))] )
  
  nb_pat_par_grp <-  paste0(nb_pat_par_grp , " (" , round(nb_pat_par_grp*100/nrow(df_one_hot_encode_fonc),1) ," \\%)")
  
  
  res <- rbind(`Number of patient` = nb_pat_par_grp,table_res_pval)
  return(res)
}

##########################################################################################
# function récupérant les variables pour table descriptive des variables

# df_one_hot_encode_fonc généré avec one_hot_fb

recup_var_table_res <- function(df_one_hot_encode_fonc,name_expo_fonc){
  res <- df_one_hot_encode_fonc %>% 
    select(-all_of(name_expo_fonc)) %>% 
    mutate_if(is.factor, ~as.numeric(as.character(.))) %>% # points litigieux j'utilise cette méthode pour convertir mes facteur booleen en numeric a surveilllé a l'avenir
    summarise_all(list(fonc_moy = ~mean(.,na.rm = TRUE),
                       fonc_std =~sd(.,na.rm = TRUE),
                       fonc_med = ~median(.,na.rm = TRUE),
                       fonc_quart1 = ~quantile(.,0.25,na.rm = TRUE),
                       fonc_quart2 = ~quantile(.,0.75,na.rm = TRUE),
                       fonc_n = ~sum(.,na.rm = TRUE), 
                       fonc_pourcent = ~mean(.,na.rm = TRUE)*100,
                       fonc_nb_NA = ~sum(is.na(.))
    )
    ) %>% 
    pivot_longer(cols =  everything(),
                 names_to = c(".value", "level"),
                 names_pattern = "(.*)_fonc_(.*)") %>% 
    t.df(.,"level") 
  return(res)
  }
