#Code for: Variation in anthelmintic responses are driven by 
#genetic differences among diverse C. elegans wild strains 
#Updated 2022 Nov 09

#load packages
require(tidyverse)
devtools::install_github("AndersenLab/easyXpress")
require(easyXpress)
install.packages("drc")
install.packages("ddpcr")
# install.packages("RColorBrewer")
install.packages("cowplot")
require(cowplot)
install.packages("ggbeeswarm")
install.packages("ggrepel")
install.packages('ggnewscale')
install.packages("magrittr")
install.packages('lmtest')
install.packages("sandwich")
#install.packages(sommer)
install.packages("kableExtra")
install.packages("beepr")
install.packages("futile.logger")

# remove.packages("RColorBrewer")

library(tidyverse)
library(easyXpress)
library(drc)
library(ddpcr)
library(RColorBrewer)
library(cowplot)
library(ggbeeswarm)
library(ggrepel)
library(ggnewscale)
library(magrittr)
library(lmtest)
library(sandwich)
library(sommer)
library(kableExtra)
library(beepr)
library(boot)
library(lme4)
library(tidyverse)
library(futile.logger)
library(sommer)
library(data.table)

# If having issues installing package sommer may need to specify package version
remove.packages("sommer") 
install.packages("remotes") #install to load a particular version of package
library(remotes)
install_version("sommer", type = "source", "4.1.3") #version used in tx dr analysis

wget https://cran.r-project.org/src/contrib/Archive/sommer/sommer_4.1.3.tar.gz
packageurl <- "https://cran.r-project.org/src/contrib/Archive/sommer/sommer_4.1.3.tar.gz"
install.packages(packageurl, repos = NULL, type="source")
  
install.packages('sommer', type = "binary")
library(sommer)


#set working directory
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

#load full dataset 50 pixels 
#load("all_processed_assays_drugspecificMDHD_50.RData") #50 pixel cutoff 

#load full dataset 30 pixels 
load("all_processed_assays_drugspecificMDHD_30.RData") #30 pixel cutoff 


#load genotype matrix
geno.matrix <- data.table::fread("data/Genotype_Matrix.tsv") %>% 
  as.data.frame()

genos <- data.table::fread("data/Genotype_Matrix.tsv") %>% 
  as.data.frame()

#load strain list
strain.list <- data.table::fread("data/strain_list.tsv", header = FALSE)

#load anthelmintic drugs and classes
anth.classes <- data.table::fread("data/anth.classes.csv")

#rename drug to anthelmintic
anth.classes <- anth.classes %>% 
  rename(Anthelmintic = drug)

#reorder anth classes to plot by major drug groups first 
anth.classes <- anth.classes %>% 
  dplyr::mutate(big_class = factor(big_class, levels = c("Benzimidazoles", "Macrocyclic_lactones","Nicotinic_acetylcholine_receptor_agonists", 
                                                         "Nicotinic_acetylcholine_receptor_antagonists", "Crystal_protein", "Cyclooctadepsipeptides",
                                                         "Schistosomicides", "Other")))

#load C. elegans strain data
strain.data <- data.table::fread(file = "data/CelegansStrainData.tsv")


#############
# FUNCTIONS #
#############

rood <- brewer.pal(n = 5, name = 'Reds')
blauw <- brewer.pal(n = 5, name = 'Blues')
groen <-brewer.pal(n = 5, name = 'Greens')


#all strain colors from lab assignments (except CX11314, does not have assignment yet)
strain_colors       <- c("#0000ff", "#FFA500", "#5F9EA0","#a021f0", "#AA4A44", "#67697C")
names(strain_colors) <- c("CB4856","N2", "DL238",  "JU775"  ,"CX11314",  "MY16")


concatenate.assays <- function(x){
  assay <- get(x)
  assay.df <- assay$summarized_processed %>%
    dplyr::mutate(Metadata_Plate = gsub(x = Metadata_Plate,
                                        pattern = "p",
                                        replacement = ""),
                  well_censor_reason = as.character(well_censor_reason),
                  notes = as.character(notes))
  
  nested.doses <- assay$summarized_processed %>%
    dplyr::select(drug, concentration_um) %>%
    dplyr::distinct() %>%
    dplyr::group_by(drug) %>%
    tidyr::nest()
  dose.rank.list <- list()
  for(i in 1:length(nested.doses$drug)){
    dose.rank.list[[i]] <- nested.doses$data[[i]] %>%
      dplyr::mutate(dose.rank = 1:nrow(nested.doses$data[[i]]),
                    drug = nested.doses$drug[[i]])
  }
  dose.ranks.df <- dose.rank.list %>%
    Reduce(rbind,.)
  assay.df %>%
    dplyr::full_join(.,dose.ranks.df)
  
  
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = na.rm,
                         ...)
  H <- 1.5 * stats::IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

collapseDoses <- function(nested.raw.data){
  if(length(unique(nested.raw.data$concentration_um)) == 1){
    return(nested.raw.data)
  } else {
    nested.raw.data$concentration_um <- unique(nested.raw.data$concentration_um)[1]
    return(nested.raw.data)
  }
}

cleanData <- function(raw.data){
  
  nested.raw.data <- raw.data %>%
    dplyr::group_by(dose.rank) %>%
    tidyr::nest()
  
  raw.data <- purrr::map(nested.raw.data$data, collapseDoses) %>%
    Reduce(rbind,.)
  
  if(length(levels(as.factor(raw.data$Metadata_Experiment))) > 1){
    
    # Remove Low n and Censored Wells (changing for small filter test, used to be 5 to 40)
    n.drug.df <- raw.data %>%
      dplyr::filter(n < 40,
                    n > 3,
                    is.na(well_censor))
    
    # Statistical Outlier Removal
    outlier.removed.drug.df <- n.drug.df %>%
      dplyr::group_by(concentration_um, strain, drug) %>%
      dplyr::mutate(median_wormlength_um = remove_outliers(median_wormlength_um)) %>%
      dplyr::filter(!is.na(median_wormlength_um))
    
    # Complete dose representation
    complete.doses <- outlier.removed.drug.df %>%
      dplyr::ungroup() %>%
      dplyr::select(Metadata_Experiment, concentration_um) %>%
      dplyr::distinct() %>%
      dplyr::group_by(concentration_um) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::filter(n > length(unique(outlier.removed.drug.df$Metadata_Experiment))*0.8) #remove measurements not represented in at least 80% (0.8) of independent assays, changing for small filter test 
    
    outlier.removed.drug.df  %>%
      dplyr::filter(concentration_um %in% complete.doses$concentration_um)
    
    # Regression | Normalize by regressing variation from assay and technical replicate effects
    median.reg <- lm(data = outlier.removed.drug.df, formula = median_wormlength_um ~ Metadata_Experiment + bleach)
    mean.reg <- lm(data = outlier.removed.drug.df, formula = mean_wormlength_um ~ Metadata_Experiment + bleach)
    outlier.removed.drug.df <- cbind(outlier.removed.drug.df, median.reg$residuals, mean.reg$residuals)
    colnames(outlier.removed.drug.df) <- c(colnames(outlier.removed.drug.df)[-c(length(colnames(outlier.removed.drug.df))-1,length(colnames(outlier.removed.drug.df)))],
                                           "median_wormlength_um_reg","mean_wormlength_um_reg")
    
    
    # Control Delta
    if(!0 %in% outlier.removed.drug.df$concentration_um){
      full_df_proc <-  outlier.removed.drug.df
    } else {
      control_values <- outlier.removed.drug.df %>%
        dplyr::filter(concentration_um == 0) %>% # filter to control wells
        dplyr::group_by(strain, bleach, Metadata_Experiment, drug) %>%
        dplyr::mutate(control_pheno = mean(median_wormlength_um_reg)) %>% # get mean control value for each trait and strain
        dplyr::distinct(control_pheno, bleach, strain, Metadata_Experiment) # make it neat
      
      full_df_proc <-  outlier.removed.drug.df %>%
        dplyr::ungroup() %>%
        dplyr::left_join(., control_values) %>% # join control values for each trait
        dplyr::mutate(median_wormlength_um_reg = median_wormlength_um_reg - control_pheno)
    }
    
    n.raw <- raw.data %>%
      dplyr::group_by(Metadata_Experiment, bleach, drug, strain) %>%
      dplyr::summarise(raw = n())
    # red in data retention  
    n.cleaned <- n.drug.df %>%
      dplyr::group_by(Metadata_Experiment, bleach, drug, strain) %>%
      dplyr::summarise(cleaned = n())
    
    # blue in data retention
    n.outlier.removed <- outlier.removed.drug.df %>%
      dplyr::group_by(Metadata_Experiment, bleach, drug, strain) %>%
      dplyr::summarise(outlier.removed = n())
    
    filtering.summary <- n.raw %>%
      dplyr::full_join(., n.cleaned) %>%
      dplyr::full_join(., n.outlier.removed)
    filtering.summary[is.na(filtering.summary)] <- 0
    
    drug.filtering.summary <- filtering.summary %>%
      dplyr::mutate(pct.retained = (outlier.removed/raw)*100,
                    pct.outlier = ((cleaned-outlier.removed)/raw)*100,
                    pct.cleaned = ((raw-cleaned)/raw)*100,
                    check = pct.retained + pct.outlier + pct.cleaned) %>%
      tidyr::pivot_longer(cols = -c(Metadata_Experiment, bleach, drug, strain), names_to = "stat")
    
  } else {
    
    # Remove Low n and Censored Wells
    n.drug.df <- raw.data %>%
      dplyr::filter(n < 40,
                    n > 3,
                    is.na(well_censor))
    
    
    # Outlier Removal
    outlier.removed.drug.df <- n.drug.df %>%
      dplyr::group_by(concentration_um, strain, drug) %>%
      # dplyr::group_by(concentration_um) %>%
      dplyr::mutate(median_wormlength_um = remove_outliers(median_wormlength_um)) %>%
      dplyr::filter(!is.na(median_wormlength_um))
    
    complete.doses <- outlier.removed.drug.df %>%
      dplyr::ungroup() %>%
      dplyr::select(Metadata_Experiment, concentration_um) %>%
      dplyr::distinct() %>%
      dplyr::group_by(concentration_um) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::filter(n > length(unique(outlier.removed.drug.df$Metadata_Experiment))*0.8)
    
    outlier.removed.drug.df  %<>%
      dplyr::filter(concentration_um %in% complete.doses$concentration_um)
    
    # Regression
    median.reg <- lm(data = outlier.removed.drug.df, formula = median_wormlength_um ~ bleach)
    mean.reg <- lm(data = outlier.removed.drug.df, formula = mean_wormlength_um ~ bleach)
    outlier.removed.drug.df <- cbind(outlier.removed.drug.df, median.reg$residuals, mean.reg$residuals)
    colnames(outlier.removed.drug.df) <- c(colnames(outlier.removed.drug.df)[-c(length(colnames(outlier.removed.drug.df))-1,length(colnames(outlier.removed.drug.df)))],
                                           "median_wormlength_um_reg","mean_wormlength_um_reg")
    
    # Control Delta
    if(!0 %in% outlier.removed.drug.df$concentration_um){
      full_df_proc <-  outlier.removed.drug.df
    } else {
      control_values <- outlier.removed.drug.df %>%
        dplyr::filter(concentration_um == 0) %>% # filter to control wells
        dplyr::group_by(strain, bleach, Metadata_Experiment, drug) %>%
        dplyr::mutate(control_pheno = mean(median_wormlength_um_reg)) %>% # get mean control value for each trait and strain
        dplyr::distinct(control_pheno, bleach, strain, Metadata_Experiment) # make it neat
      
      full_df_proc <-  outlier.removed.drug.df %>%
        dplyr::ungroup() %>%
        dplyr::left_join(., control_values) %>% # join control values for each trait
        dplyr::mutate(median_wormlength_um_reg = median_wormlength_um_reg - control_pheno)
    }
    
    n.raw <- raw.data %>%
      dplyr::group_by(Metadata_Experiment, bleach, drug, strain) %>%
      dplyr::summarise(raw = n())
    
    n.cleaned <- n.drug.df %>%
      dplyr::group_by(Metadata_Experiment, bleach, drug, strain) %>%
      dplyr::summarise(cleaned = n())
    
    n.outlier.removed <- outlier.removed.drug.df %>%
      dplyr::group_by(Metadata_Experiment, bleach, drug, strain) %>%
      dplyr::summarise(outlier.removed = n())
    
    filtering.summary <- n.raw %>%
      dplyr::full_join(., n.cleaned) %>%
      dplyr::full_join(., n.outlier.removed)
    filtering.summary[is.na(filtering.summary)] <- 0
    
    drug.filtering.summary <- filtering.summary %>%
      dplyr::mutate(pct.retained = (outlier.removed/raw)*100,
                    pct.outlier = ((cleaned-outlier.removed)/raw)*100,
                    pct.cleaned = ((raw-cleaned)/raw)*100,
                    check = pct.retained + pct.outlier + pct.cleaned) %>%
      tidyr::pivot_longer(cols = -c(Metadata_Experiment, bleach, drug, strain), names_to = "stat")
  }
  
  full_df_proc %<>%
    dplyr::ungroup() %>%
    dplyr::group_by(concentration_um, strain, drug) %>%
    dplyr::mutate(median_wormlength_um_reg = remove_outliers(median_wormlength_um_reg))
  
  return(list(full_df_proc,drug.filtering.summary))
}


#################################
# dose-response model inference #
#################################

# general four-parameter log-logistic model function 
# fixed = is used for fixing parameters at given values 

# b = slope of DR curve 
# c = upper asymptote of DR curve 
# d = lower asymptote of DR curve 
# e = specified effective dose 

safe_LL4_strain <- purrr::safely(.f = ~ drc::drm(data = ., 
                                                 formula = value ~ concentration_um, strain,
                                                 pmodels=list(~strain-1,  ~1, ~1, ~strain-1),
                                                 fct = LL.4(fixed=c(NA, -600, NA, NA))),
                                 otherwise = "Unable to optimize model")


#No strain in model
safe_LL4_no_strain <- purrr::safely(.f = ~ drc::drm(data = ., #model fitting function 
                                                    formula = value ~ concentration_um, #no strain
                                                    pmodels=list(~1,  ~1, ~1, ~1),
                                                    fct = LL.4(fixed=c(NA, -600, NA, NA))),
                                    otherwise = "Unable to optimize model")




safe_ED <- function(fit){
  if(is.character(fit$result)){print("Unable to optimize model")}
  else {
    drc::ED(object = fit$result, c(10, 50, 90), interval = "delta", display = FALSE) #extract strain specific EC10, 50, 90 values
  }
}


safe_slope_strain <- function(fit){
  if(is.character(fit$result)){print("Unable to extract slope")}
  else {
    #slopes<- data.frame(fit$result$coefficients[1:8]) %>% 8 strains for Sam? 
      slopes<- data.frame(fit$result$coefficients[1:6]) %>%
       # dplyr::rename(slope = `fit.result.coefficients.1.8.`) %>%8 strains for Sam? 
      dplyr::rename(slope = `fit.result.coefficients.1.6.`) %>%
      dplyr::mutate(strain = rownames(.),
                    strain = gsub(strain, pattern = "b:strain", replacement = ""))
    rownames(slopes) <- NULL
    return(slopes)
  }
}

#data used in lower asymp
# View(data [10:17])

safe_lower_asym <- function(fit){
  if(is.character(fit$result)){print("Unable to extract slope")}
  else {
    
    asyms <- data.frame(fit$result$coefficients[10:17]) %>% #data[10:17] = q90_worm_length, max_worm_length, cv_worm_length, n, metadata, metaexperiment, magnification, assay type
      dplyr::rename(lower.asymp = `fit.result.coefficients.10.17.`) %>%
      dplyr::mutate(strain = rownames(.),
                    strain = gsub(strain, pattern = "d:strain", replacement = ""))
    rownames(asyms) <- NULL
    return(asyms)
  }
}


safe_predict <- function(data, fit){
  if(is.character(fit$result)){print("Unable to optimize model")}
  else{
    expand.grid(conc = exp(seq(log(max(data$concentration_um)),
                               log(min(data$concentration_um[data$concentration_um!=min(data$concentration_um)])),
                               length = 100))) %>%
      dplyr::bind_cols(., tibble::as.tibble(predict(fit$result, newdata=., interval="confidence")))
  }
  
}


fit_DRC_models <- function(data, traits = "median_wormlength_um_reg", ...){
  # extract concentration variable name

  conc_var_name <- str_subset(names(data), pattern = "concentration_um")
  
  
  ##### old code #####
  # complete.doses <- data %>%
  #   dplyr::ungroup() %>%
  #   dplyr::select(strain, drug, concentration_um) %>%
  #   dplyr::distinct() %>%
  #   dplyr::group_by(concentration_um) %>%
  #   dplyr::summarise(n()) %>%
  #   dplyr::filter(`n()` == 8)
  #
  
  # # gather trait data, nest by grouping variables and trait, fit drc, find ED10, ED50, ED90, and make predictions.
  
  # WORKING
  # LL4_fit_data <- data %>%
  #   dplyr::select(Metadata_Plate, Metadata_Well, Metadata_Experiment, bleach, drug, strain,
  #                 contains("wormlength"), concentration_um := !!glue::glue("{conc_var_name}")) %>%
  #   tidyr::gather(trait, value, -drug, -strain, -concentration_um, -Metadata_Plate, -Metadata_Well, -Metadata_Experiment, -bleach) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::select(drug, strain, trait, value, concentration_um) %>% # reorder for drc::drm
  #   dplyr::filter(trait %in% traits) %>% # filter to traits provided in arguments
  #   # dplyr::group_by(drug, strain, trait, Metadata_Experiment) %>%
  #   dplyr::group_by(drug, strain, trait) %>% # should eventually be "..." passed from function + trait
  #   tidyr::nest() %>%
  #   mutate(fit = purrr::map(data, safe_LL4),
  #          slope = purrr::map(fit, safe_slope),
  #          lower_asym = purrr::map(fit, safe_lower_asym),
  #          ED = purrr::map(fit, safe_ED),
  #          predictions = purrr::map2(data, fit, safe_predict),
  #          model = "LL.4")
  #####
  
# View(data)  
  
  LL4_fit_data <- data %>%
    dplyr::select(Metadata_Plate, Metadata_Well, Metadata_Experiment, bleach, drug, strain, 
                  contains("wormlength"), concentration_um := !!glue::glue("{conc_var_name}")) %>%
    tidyr::gather(trait, value, -drug, -strain, -concentration_um, -Metadata_Plate, -Metadata_Well, -Metadata_Experiment, -bleach) %>%
    dplyr::ungroup() %>%
    dplyr::select(drug, strain, trait, value, concentration_um) %>% # reorder for drc::drm
    dplyr::filter(trait %in% traits) %>% # filter to traits provided in arguments
    # dplyr::group_by(drug, strain, trait, Metadata_Experiment) %>%
    dplyr::group_by(drug, trait) %>% # should eventually be "..." passed from function + trait
    tidyr::nest() %>%
    mutate(fit = purrr::map(data, safe_LL4_strain),
           slope = purrr::map(fit, safe_slope_strain),
           lower_asym = purrr::map(fit, safe_lower_asym),
           ED = purrr::map(fit, safe_ED),
           model = "LL.4")
  return(LL4_fit_data)
}


logticks <- function(datavar,type) {
  minimum <- 1/10^abs(floor(log10(min(datavar, na.rm=TRUE))))
  maximum <- 1*10^abs(floor(log10(max(datavar, na.rm=TRUE)))+1)
  multiple <- floor(log10(maximum/minimum))
  
  yourtickvector <- c()
  
  if (type=="breaks") {
    
    yourtickvector <- c(minimum)
    
    for (x in seq(0,multiple)) {
      
      andadd <- seq(minimum*10^x,minimum*10^(x+1),minimum*10^x)[-1]
      
      yourtickvector <- c(yourtickvector,andadd)
      
    }
  } else if (type=="labels") {
    
    for (x in seq(0,multiple)) {
      
      andadd <- c(minimum*10^x,rep("",8))
      
      yourtickvector <- c(yourtickvector,andadd)
      
    }
    
    yourtickvector <- c(yourtickvector,minimum*10^multiple)
    
  }
  
  return(yourtickvector)
  
}


dose.check.plot <- function(dat){
  strain_colors       <- c("#0000ff", "#FFA500", "#5F9EA0","#a021f0", "#AA4A44", "#67697C")
  names(strain_colors) <- c("CB4856","N2", "DL238",  "JU775"  ,"CX11314",  "MY16")
  dose.check <- ggplot(dat) +
    theme_bw(base_size = 8) +
    # geom_jitter(mapping = aes(x = strain, y = median_wormlength_um_reg, colour = Metadata_Experiment),size = 1) +
    geom_quasirandom(mapping = aes(x = strain, y = median_wormlength_um_reg, colour = Metadata_Experiment), size = 1) +
    geom_boxplot(mapping = aes(x = strain, y = median_wormlength_um_reg),
                 outlier.shape = NA, alpha = 0.5) +
    facet_wrap(.~round(concentration_um,2), scales = "free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "top") +
    labs(fill = "Strain",
         y = "Relative Worm Length (um)",
         x = "Genotype") + #Previously "Strain"
    theme(plot.title = element_text(hjust = 0.5), legend.box="vertical") +
    ggplot2::theme(legend.position = "none") + #see if this adds legends aos 
    # ylim(c(-810,200)) +
    #scale_colour_brewer(palette = "Dark2") +
    ggtitle(unique(dat$drug))
  
  ggsave(dose.check, filename = paste0("output/", paste(paste(strsplit(unique(dat$drug),split = " ")[[1]], collapse = "_"), "dose.check.plot", sep = "_"),".png"),
         width = 10, height = 6)
}


data.retention <- function(data){
  
  A <- data %>%
    dplyr::filter(stat %in% c("pct.retained","pct.outlier","pct.cleaned")) %>%
    dplyr::mutate(bleach = paste("Bleach", bleach, sep = " ")) %>%
    ggplot(., mapping = aes(x = strain, y = value, fill = stat)) +
    theme_bw() +
    geom_col() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "top",
          panel.grid = element_blank()) +
    scale_fill_brewer(palette = "Set1", name = "Cleaning Category") +
    facet_wrap(.~Metadata_Experiment+bleach) +
    ggtitle("Data Retention Cleaning: Before") +
    labs(x = "Genotype", y = "Experimental Retention (% of Expected Wells)")
  
  data.retention.df <- data %>%
    tidyr::pivot_wider(names_from = stat, values_from = value) %>%
    dplyr::filter(pct.cleaned < 35) %>%
    dplyr::group_by(Metadata_Experiment, bleach) %>%
    dplyr::summarise(n.strains = n()) %>%
    dplyr::filter(n.strains > 4)
  
  
  B <-  data.retention.df %>%
    dplyr::left_join(.,data) %>%
    dplyr::filter(stat %in% c("pct.retained","pct.outlier","pct.cleaned")) %>%
    # tidyr::pivot_longer(cols = -c(Metadata_Experiment, bleach), names_to = "stat") %>%
    ggplot(., mapping = aes(x = strain, y = value, fill = stat)) +
    theme_bw() +
    geom_col() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "top",
          panel.grid = element_blank()) +
    scale_fill_brewer(palette = "Set1", name = "Cleaning Category") +
    facet_wrap(.~Metadata_Experiment+bleach) +
    ggtitle("Data Retention Cleaning: After") +
    labs(x = "Genotype", y = "")
  
  AB <- cowplot::plot_grid(A + theme(legend.position = "none"),
                           B + theme(legend.position = "none"))
  
  title <- ggdraw() + draw_label(unique(data$drug), fontface = "bold")
  data.retention.plot <- cowplot::plot_grid(title, AB, ncol = 1,
                                            rel_heights = c(0.1, 1))
  ggsave(data.retention.plot + theme(plot.background = element_rect(fill = "white")), filename = paste0("output/", paste(paste(strsplit(unique(data$drug),split = " ")[[1]], collapse = "_"), "data.retention.plot", sep = "_"),".png"),
         width = 10, height = 6)
  
  return(data.retention.df[,1:2])
  
}


# se.DRC <- function(x, dodge.width = 0.03){ #originally dodge.width 0.15
#   strain_colors       <- c("#0000ff", "#FFA500", "#5F9EA0","#a021f0", "#AA4A44", "#67697C")
#   names(strain_colors) <- c("CB4856","N2", "DL238",  "JU775"  ,"CX11314",  "MY16")
#   
#   norm <- x %>%
#     dplyr::group_by(strain) %>%
#     dplyr::filter(concentration_um == 0) %>%
#     dplyr::summarize(mean_ctrl = mean(mean_wormlength_um_reg)) %>%
#     dplyr::ungroup() %>%
#     dplyr::select(strain, mean_ctrl) %>% ##added trait
#     dplyr::full_join(x, by = c("strain")) %>% ## "trait" instead of "condition"
#     dplyr::distinct() %>%
#     dplyr::mutate(norm_pheno = mean_wormlength_um_reg - mean_ctrl)
#   
#   dose_avg <- norm %>%
#     dplyr::group_by(strain, concentration_um) %>%
#     dplyr::summarise(avg = mean(norm_pheno), st_dev = sd(norm_pheno))
#   
#   strains <- dose_avg %>%
#     dplyr::mutate(log.concentration = round(log(concentration_um), 2))
#   
#   # conc.range <- range(strains[which(strains$log.concentration != -Inf),]$concentration_um)[2] - range(strains[which(strains$log.concentration != -Inf),]$concentration_um)[1]
  # dodge.width = conc.range/650
 
  
  #ORIGINALLY WHAT I WAS USING FOR SE.DRC 
  # SE.DRC <- strains[which(strains$log.concentration != -Inf),] %>%
  #   ggplot(.,mapping = aes(x=concentration_um, y=avg)) +
  #   theme_bw(base_size = 18) +
  #   geom_line(mapping = aes(group = strain, colour = strain),
  #             position = position_dodge2(width = dodge.width)) +
  #   geom_pointrange(mapping = aes(ymin=avg-st_dev,
  #                                 ymax = avg+st_dev,
  #                                 colour = strain),
  #                   position = position_dodge2(width = dodge.width),
  #                   size = 0.4) +
  #   
  #   scale_colour_manual(values = strain_colors, name = "Strain") +
  #   scale_y_continuous("Normalized Worm Length (µm)") +
  #   scale_x_log10() +
  #   theme(panel.grid = element_blank(),
  #         legend.position = "none",
  #         axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  #   labs(x="Concentration (μM)") +
  #   ggtitle(unique(x$drug)) + 
  #   ggplot2::theme(legend.position = "none") #see if this adds strain legends aos 
  # ggsave(SE.DRC, filename = paste0("output/", paste(paste(strsplit(unique(x$drug),split = " ")[[1]], collapse = "_"), "SE.DRC.plot", sep = "_"),".png"),
  #        width = 6, height = 6)
  # print(SE.DRC)
  
#}

se.DRC <- function(x, dodge.width = 0.03, point.size = 0.4){
  
  strain_colors       <- c("#0000ff", "#FFA500", "#5F9EA0","#a021f0", "#AA4A44", "#67697C")
  names(strain_colors) <- c("CB4856","N2", "DL238",  "JU775"  ,"CX11314",  "MY16")
  
  norm <- x %>%
    dplyr::group_by(strain) %>%
    dplyr::filter(concentration_um == 0) %>%
    dplyr::summarize(mean_ctrl = mean(mean_wormlength_um_reg)) %>%
    dplyr::ungroup() %>%
    dplyr::select(strain, mean_ctrl) %>% ##added trait
    dplyr::full_join(x, by = c("strain")) %>% ## "trait" instead of "condition"
    dplyr::distinct() %>%
    dplyr::mutate(norm_pheno = mean_wormlength_um_reg - mean_ctrl) 
  
  dose_avg <- norm %>%
    dplyr::group_by(strain, concentration_um) %>%
    dplyr::summarise(avg = mean(norm_pheno), st_dev = sd(norm_pheno))
  
  strains <- dose_avg %>%
    dplyr::mutate(log.concentration = round(log(concentration_um), 2))
  
  # conc.range <- range(strains[which(strains$log.concentration != -Inf),]$concentration_um)[2] - range(strains[which(strains$log.concentration != -Inf),]$concentration_um)[1] 
  # dodge.width = conc.range/650
  
  #BZs make bottom asymptote -500
  #MLs make bottom asymptote -600 
  #NAChAs make bottom asymptote -500 (need to remove names for this so doesn't interfere with DRC)
  SE.DRC <- strains[which(strains$log.concentration != -Inf),] %>%
    ggplot(.,mapping = aes(x=concentration_um, y=avg)) +
    theme_bw(base_size = 18) + 
    geom_line(mapping = aes(group = strain, colour = strain), 
              position = position_dodge2(width = dodge.width), size = 0.15) +
    geom_pointrange(mapping = aes(ymin=avg-st_dev, 
                                  ymax = avg+st_dev,
                                  colour = strain), 
                    position = position_dodge2(width = dodge.width),
                    size = point.size) + 
    
    scale_colour_manual(values = strain_colors, name = "Strain") +
    scale_y_continuous("Normalized Worm Length (µm)", limits = c(-600,100), breaks = seq(-600, 0, 100)) + 
    scale_x_log10() +
    # theme(panel.grid = element_blank(), #use old theme to get rid of legend
    #         legend.position = "none",
    #        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
    theme(panel.grid = element_blank(),
      legend.position = c(0.15,0.25),
      legend.background = element_blank(),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      axis.text = element_text(colour = "black")) +
    labs(x="Concentration (μM)") + 
    ggtitle(unique(x$drug))
  ggsave(SE.DRC, filename = paste0("output/", paste(paste(strsplit(unique(x$drug),split = " ")[[1]], collapse = "_"), "SE.DRC.plot", sep = "_"),".png"), 
         width = 6, height = 6)
  print(SE.DRC)
  
}


#using cry5b [[11]] as an example to figure out where we need to make adjustments in dose.response.summary 

#data = tx.nested.dose.responses$data[[11]]
#Anthelmintic = tx.nested.dose.responses$drug[[11]]

dose.response.summary <- function(data, Anthelmintic){
  
  # Data cleaning
  raw.data <- data %>%
    dplyr::mutate(drug = Anthelmintic)
  cleaned.data <- suppressMessages(cleanData(raw.data))
  
  
  retained.bleaches <- data.retention(cleaned.data[[2]])
  dose.check.plot(dat = cleaned.data[[1]])
  
  
  # are any doses represented by fewer than 10 total replicates per strain?
  cleaned.filtered.data <- retained.bleaches %>%
    dplyr::left_join(.,cleaned.data[[1]]) %>%
    dplyr::group_by(strain, concentration_um) %>%
    dplyr::count() %>%
    dplyr::arrange(n) %>%
    dplyr::filter(n > 10) %>% 
    dplyr::select(-n) %>%
    dplyr::left_join(.,cleaned.data[[1]])
  
  se.DRC(cleaned.filtered.data)
  
  #suppressMessages evaluates its expression in a context that ignores all ‘simple’ diagnostic messages.
  drug.drm <- suppressMessages(fit_DRC_models(data = cleaned.filtered.data,
                                              traits = "median_wormlength_um_reg"))
  
  
  # Population-wide LOAELs
  population.wide.LOAEL <- aov(data = cleaned.filtered.data, formula = median_wormlength_um_reg ~ as.factor(concentration_um))
  LOAEL.tukey <- TukeyHSD(population.wide.LOAEL)
  
  LOAEL.tests <- rownames(LOAEL.tukey$`as.factor(concentration_um)`) %>%
    data.frame() %>%
    tidyr::separate(data = ., col = `.`, c("dose1","dose2"), sep = "-") %>%
    dplyr::filter(dose2 == 0) %>%
    rownames() %>% as.numeric()
  
  adj.Ps <- LOAEL.tukey$`as.factor(concentration_um)`[LOAEL.tests,] %>%
    data.frame()
  
  doses.adj.Ps <- rownames(LOAEL.tukey$`as.factor(concentration_um)`) %>%
    data.frame() %>%
    tidyr::separate(data = ., col = `.`, c("dose1","dose2"), sep = "-") %>%
    dplyr::filter(dose2 == 0) %>%
    cbind(.,adj.Ps)
  rownames(doses.adj.Ps) <- NULL
  
  
  LOAEL <-  doses.adj.Ps %>%
    dplyr::filter(p.adj < 0.05 & diff < 0) %>%
    # dplyr::filter(p.adj < 0.05) %>% do you want to filter on significantly lower responses from control if hormetic?
    dplyr::mutate(tx = Anthelmintic) %>%
    dplyr::rename(LOAEL = dose1) %>%
    dplyr::select(-dose2)
  rownames(LOAEL) <- NULL
  overall.tx.LOAEL <- LOAEL[which(LOAEL$LOAEL == min(LOAEL$LOAEL)),]
  
  
  # Strain-specific LOAELs
  LOAEL.strain.nested <- cleaned.filtered.data %>%
    dplyr::group_by(strain) %>%
    tidyr::nest()
  strain.LOAELs <- list()
  for(i in 1:length(LOAEL.strain.nested$strain)){
    strain.LOAEL <- aov(data = LOAEL.strain.nested$data[[i]], formula = median_wormlength_um_reg ~ as.factor(concentration_um))
    strain.LOAEL.tukey <- TukeyHSD(strain.LOAEL)
    
    strain.LOAEL.tests <- rownames(strain.LOAEL.tukey$`as.factor(concentration_um)`) %>%
      data.frame() %>%
      tidyr::separate(data = ., col = `.`, c("dose1","dose2"), sep = "-") %>%
      dplyr::filter(dose2 == 0) %>%
      rownames() %>% as.numeric()
    
    strain.adj.Ps <- strain.LOAEL.tukey$`as.factor(concentration_um)`[strain.LOAEL.tests,] %>%
      data.frame() %>%
      dplyr::select(p.adj)
    
    strain.doses.adj.Ps <- rownames(strain.LOAEL.tukey$`as.factor(concentration_um)`) %>%
      data.frame() %>%
      tidyr::separate(data = ., col = `.`, c("dose1","dose2"), sep = "-") %>%
      dplyr::filter(dose2 == 0) %>%
      cbind(.,strain.adj.Ps)
    rownames(strain.doses.adj.Ps) <- NULL
    
    
    strain.LOAEL.df <-  strain.doses.adj.Ps %>%
      dplyr::filter(p.adj < 0.05) %>%
      dplyr::mutate(tx = Anthelmintic) %>%
      dplyr::rename(LOAEL = dose1) %>%
      dplyr::select(-dose2) %>%
      dplyr::mutate(strain = LOAEL.strain.nested$strain[[i]])
    rownames(strain.LOAEL.df) <- NULL
    strain.LOAELs[[i]] <- strain.LOAEL.df
  }
  
  
  strain.tx.LOAEL <- strain.LOAELs %>%
    Reduce(rbind,.)
  
    
    ED.table.df <- data.frame(drug.drm$ED)
    ED.table.df$params <- rownames(ED.table.df)
    rownames(ED.table.df) <- NULL
    ED.df <- ED.table.df %>%
      tidyr::separate(params, c("x","strain","metric"),sep = ":") %>%
      dplyr::select(-x) %>%
      dplyr::mutate(metric = paste0("EC",metric)) %>%
      dplyr::mutate(drug = Anthelmintic)
    
    
    
    # Dose-response statistics
    strain.model.fit <- safe_LL4_strain(drug.drm$data[[1]])
    no.strain.model.fit <- safe_LL4_no_strain(drug.drm$data[[1]])
    summarized.model.fit <- summary(strain.model.fit$result)
    strain.coefs <- data.frame(summarized.model.fit$coefficients)
    strain.coefs$params <- rownames(strain.coefs)
    rownames(strain.coefs) <- NULL
    slope.df <- strain.coefs %>%
      tidyr::separate(params, c("parameter","strain"), sep = ":") %>%
      dplyr::mutate(strain = gsub(strain, pattern = "strain", replacement = "")) %>%
      dplyr::filter(parameter == "b") %>%
      dplyr::rename(metric = parameter) %>%
      dplyr::mutate(drug = Anthelmintic) %>%
      dplyr::select(Estimate, Std..Error, metric, strain, drug, p.value)
    
    # not super important
    intercept <- strain.coefs %>%
      tidyr::separate(params, c("parameter","strain"), sep = ":") %>%
      # dplyr::mutate(strain = gsub(strain, pattern = "strain", replacement = "")) %>%
      dplyr::filter(parameter == "c") %>%
      dplyr::select(Estimate) %>%
      as.numeric()
    
    # not super important
    asymp.df <- strain.coefs %>%
      tidyr::separate(params, c("parameter","strain"), sep = ":") %>%
      dplyr::mutate(strain = gsub(strain, pattern = "strain", replacement = "")) %>%
      dplyr::filter(parameter == "d") %>%
      dplyr::rename(metric = parameter) %>%
      dplyr::mutate(drug = Anthelmintic,
                    int = intercept,
                    lower.asymp = Estimate + int) %>%
      dplyr::select(Estimate, Std..Error, metric, lower.asymp, strain, drug, p.value)
    
    
    dose.response.parameters.df <- ED.df %>%
      dplyr::full_join(.,slope.df) %>%
      dplyr::full_join(.,asymp.df)
    
    
    drc.anova <- anova(no.strain.model.fit$result, strain.model.fit$result)
    
    
    EDcomps.df <- suppressMessages(data.frame(drc::EDcomp(strain.model.fit$result, c(10,10))))
    EDcomps.df$comps <- rownames(EDcomps.df)
    rownames(EDcomps.df) <- NULL
    ed.comps <- EDcomps.df %>%
      dplyr::mutate(drug = Anthelmintic)
    
    
    Slope.comps.df <- suppressMessages(data.frame(drc::compParm(strain.model.fit$result, "b")))
    Slope.comps.df$comps <- rownames(Slope.comps.df)
    rownames(Slope.comps.df) <- NULL
    slope.comps <- Slope.comps.df %>%
      dplyr::mutate(drug = Anthelmintic)
    
    
    
    return(list(overall.tx.LOAEL, strain.tx.LOAEL,
                dose.response.parameters.df,
                summarized.model.fit, drc.anova, ed.comps, slope.comps,
                cleaned.filtered.data))
  }
  

dose.response.plots.only <- function(data, Anthelmintic, SE.drc.dodge.width = 0.01, SE.drc.point.size = 0.4){
  raw.data <- data %>%
    dplyr::mutate(drug = Anthelmintic)
  cleaned.data <- suppressMessages(cleanData(raw.data))
  
  
  retained.bleaches <- data.retention(cleaned.data[[2]])
  dose.check.plot(dat = cleaned.data[[1]])
  
  cleaned.filtered.data <- retained.bleaches %>%
    dplyr::left_join(.,cleaned.data[[1]]) %>%
    dplyr::group_by(strain, concentration_um) %>%
    dplyr::count() %>%
    dplyr::arrange(n) %>%
    dplyr::filter(n > 5) %>%
    dplyr::select(-n) %>%
    dplyr::left_join(.,cleaned.data[[1]])
  
  se.DRC(cleaned.filtered.data, dodge.width = SE.drc.dodge.width, point.size = SE.drc.point.size)
}

make.reps <- function(x){
  x$replicate <- seq(1,nrow(x))
  x %>%
    dplyr::rename(value = median_wormlength_um_reg)
}

###########################
####### HERITABILITY ######
###########################

#Heritability code from latest github relase as of 2022/08/31
#lab herit calc: https://github.com/AndersenLab/calc_heritability/blob/main/bin/20210716_H2_script.R

H2 <- function(d){
  strain.fit <- lme4::lmer(data = d, formula = value ~ 1 + (1|strain))
  variances <- lme4::VarCorr(x = strain.fit)
  A <- as.data.frame(variances)
  Vg <- A$vcov[1]
  Ve <- A$vcov[2]
  H2 <- Vg/(Vg+Ve)
  return(H2)
}

h2 <- function(d, geno_matrix){
  pheno_strains <- unique(d$strain)
  A <- sommer::A.mat(t(geno_matrix[,colnames(geno_matrix) %in% pheno_strains]))
  
  df_y <- d %>%
    dplyr::arrange(strain) %>%
    dplyr::select(strain, value) %>%
    dplyr::mutate(strain = as.character(strain)) %>% 
    dplyr::mutate(strain = factor(strain, levels = unique(.$strain))) 
  h2_res <- sommer::mmer(value~1,
                         random=~vs(strain,Gu=A), 
                         data=df_y)
  h2 <- vpredict(h2_res, h2 ~ (V1) / ( V1+V2))[[1]][1]
  h2_SE <- vpredict(h2_res, h2 ~ (V1) / ( V1+V2))[[2]][1]
  
  h2_df <- data.frame(h2) %>%
    dplyr::mutate(h2.upper = h2 + h2_SE,
                  h2.lower = h2 - h2_SE)
  return(h2_df)
}

H2.bootstrapping.calc <- function(d, nreps = 100, boot = T){
  
  if(boot == T){
    # Broad-sense Heritability
    H2.point <- H2(d = d)
    h2.point <- h2(d = d, geno_matrix = genos)
    H2.boots <- list()
    for(i in 1:nreps) {
      if(i %% 10 == 0){
        print(paste0((i/nreps)*100,"%"))
      }
      #################################
      # Bootstrap within strain ##
      #################################
      nested <- d %>%
        dplyr::group_by(strain) %>%
        tidyr::nest()
      boot.strain <- list()
      for(j in 1:length(nested$strain)){
        boot.strain[[j]] <- nested$data[[j]][sample(seq(1:nrow(nested$data[[j]])),replace = T),] %>%
          dplyr::mutate(strain = nested$strain[[j]])
      }
      boot <- boot.strain %>%
        Reduce(rbind,.)
      
      ##################################
      ## Bootstrap agnostic of strain ##
      ##################################
      # boot <- d[sample(seq(1:nrow(d)), replace = T),]
      
      check <- boot %>%
        dplyr::group_by(strain) %>%
        dplyr::summarise(n())
      if(1 %in% check$`n()`){
        print("Only 1 Strain Sampled in Bootstrap - Skipping")
        next
      }
      # Broad-Sense Heritability
      H2.est <- H2(d = boot)
      H2.boots[i] <- H2.est
    }
    
    H2.boots.vec <- unlist(H2.boots)
    H2.quantiles <- c(quantile(H2.boots.vec, probs = seq(0,1,0.05)))
    H2.CI <- data.frame(H2.point, 
                        as.numeric(H2.quantiles[2]), 
                        as.numeric(H2.quantiles[21])) %>%
      `colnames<-`(c("H2.Point.Estimate","H2.5.perc","H2.95.perc"))
    
    
    
    return(list(H2.CI,H2.boots.vec,h2.point))
    
  } else {
    
    H2.point <- H2(d = d)
    # h2.point <- h2(d = d)
    H2.CI <- data.frame(H2.point, 
                        NA, 
                        NA) %>%
      `colnames<-`(c("H2.Point.Estimate","H2.5.perc","H2.95.perc"))
    return(H2.CI)
  }
  
}
drugHeritability <- function(rep.data, concentration_um){
  
  # Broad-sense Heritability
  no.reps <- rep.data %>%
    dplyr::group_by(strain) %>%
    dplyr::summarise(n = n())
  
  
  if(1 %in% no.reps$n){
    H2.est <- data.frame("NA", "NA", "NA") %>%
      dplyr::mutate(drug = unique(rep.data$drug),
                    concentration = concentration_um) %>%
      `colnames<-`(c("H2", "lower", "upper", "drug", "concentration_um"))
    
    h2.est <- data.frame("NA", "NA", "NA") %>%
      dplyr::mutate(drug = unique(rep.data$drug),
                    concentration = concentration_um) %>%
      `colnames<-`(c("h2", "upper", "lower", "drug", "concentration_um"))
    
    list(H2.est, h2.est)
    
  } else {
    
    herits <- suppressMessages(H2.bootstrapping.calc(d = rep.data, boot = T))
    
    H2.est <- as.data.frame(herits[[1]]) %>%
      dplyr::mutate(drug = unique(rep.data$drug),
                    concentration = concentration_um) %>%
      `colnames<-`(c("H2", "lower", "upper", "drug", "concentration_um"))
    
    H2.replicates <- as.data.frame(herits[[2]]) %>%
      dplyr::mutate(drug = unique(rep.data$drug),
                    concentration = concentration_um) %>%
      `colnames<-`(c("H2.rep", "drug", "concentration_um"))
    
    h2.est <- as.data.frame(herits[[3]]) %>%
      dplyr::mutate(drug = unique(rep.data$drug),
                    concentration = concentration_um) %>%
      `colnames<-`(c("h2", "upper", "lower", "drug", "concentration_um"))
    
    list(H2.est,H2.replicates,h2.est)
  }
}

heritability.calculation <- function(data, Anthelmintic){
  # Data cleaning
  raw.data <- data %>%
    dplyr::mutate(drug = Anthelmintic)
  cleaned.data <- suppressMessages(cleanData(raw.data))
  retained.bleaches <- data.retention(cleaned.data[[2]])
  cleaned.filtered.data <- retained.bleaches %>%
    dplyr::left_join(.,cleaned.data[[1]])
  
  dose.nested <- cleaned.filtered.data %>%
    dplyr::group_by(concentration_um) %>%
    tidyr::nest() 
  
  dose.nested$data <- purrr::map(dose.nested$data, make.reps)
  dose.data.nested.reps <- dose.nested %>%
    tidyr::unnest(cols = c(data)) %>%
    dplyr::group_by(concentration_um) %>%
    tidyr::nest()
  
  
  herits <- purrr::map2(dose.data.nested.reps$data, 
                        dose.data.nested.reps$concentration_um,
                        drugHeritability)
  
}

plot_HH <- function(herit.data){
  ggplot(herit.data, mapping = aes(x = H2, y = h2, colour = log(concentration_um))) + 
    theme_bw(base_size = 11) +
    geom_abline(slope = 1, alpha = 0.5, linetype = 3, size = 0.25) +
    # geom_smooth(method = "lm", se = F, colour = "black", linetype = 1, size = 0.5) +
    # stat_smooth(method="lm_right", fullrange=TRUE, colour = "black", se = F, size = 0.5, linetype = 1,alpha = 0.5) + 
    geom_errorbar(aes(ymin = h2.lower, ymax = h2.upper), size = 0.2) +
    geom_errorbarh(aes(xmin = H2.lower, xmax = H2.upper), size = 0.2) + 
    facet_grid(. ~ drug) + 
    scale_colour_gradient2(low = "darkblue", mid = "violet", high = "darkred", midpoint = 2, 
                           name = expression(log[10](Concentration (µM)))) + 
    theme(panel.grid = element_blank(),
          axis.title.y = element_text(size = 8.5, angle = 0, vjust = 0.5),
          axis.title.x = element_text(size = 8.5),
          axis.text = element_text(size = 7), 
          strip.text = element_text(size = 7),
          legend.position="bottom",
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7)) + 
    scale_x_continuous(breaks = seq(0,1,0.5), limits = c(0,1)) + 
    scale_y_continuous(breaks = seq(0,1,0.5), limits = c(0,1)) + 
    labs(x = bquote(italic(H^2)),
         y = bquote(italic(h^2)))
}

gather.herit.ranges <- function(all.herits){
  print(all.herits[[1]][[3]]$drug)
  dose.herits <- list()
  for(i in 1:length(all.herits)){
    if(length(all.herits[[i]]) == 3){
      H2.estimate <- all.herits[[i]][[1]] %>%
        dplyr::rename(H2.upper = upper, H2.lower = lower)
      
      h2.estimate <- all.herits[[i]][[3]] %>%
        dplyr::rename(h2.upper = upper, h2.lower = lower)
      
      dose.herits[[i]] <- H2.estimate %>%
        dplyr::full_join(h2.estimate,.) %>%
        suppressMessages() %>%
        dplyr::select(drug, concentration_um, everything())
    } else {
      H2.estimate <- all.herits[[i]][[1]] %>%
        dplyr::rename(H2.upper = upper, H2.lower = lower)
      
      h2.estimate <- all.herits[[i]][[2]] %>%
        dplyr::rename(h2.upper = upper, h2.lower = lower)
      
      dose.herits[[i]] <- H2.estimate %>%
        dplyr::full_join(h2.estimate,.) %>%
        suppressMessages() %>%
        dplyr::select(drug, concentration_um, everything())
    }
    
  }
  dose.herits %>% Reduce(rbind,.)
}

h2.comps.to.individual.strains <- function(EC10, this.strain){
  strain.EC10.comp <- EC10 %>%
    dplyr::filter(drug %in% h2.comp.table$Anthelmintic) %>%
    dplyr::select(drug, Estimate) %>%
    dplyr::rename(Anthelmintic = drug) %>%
    dplyr::left_join(., h2.comp.table %>%
                       dplyr::select(Anthelmintic, `Top Heritable Dose`))
  strain.h2.comp.r2 <- summary(lm(data = strain.EC10.comp, formula = log10(`Top Heritable Dose`) ~ log10(Estimate)))
  lm.params <- list(r2 = format(strain.h2.comp.r2$r.squared, digits = 3),
                    pval = format(as.numeric(pf(strain.h2.comp.r2$fstatistic[1],
                                                strain.h2.comp.r2$fstatistic[2],
                                                strain.h2.comp.r2$fstatistic[3],
                                                lower.tail=FALSE)), digits=3)) %>% 
    data.frame() %>%
    dplyr::mutate(strain = this.strain)
  
}


############
############
# ANALYSIS #
############
############

#############
#############
# 30 pixels #
#############
#############

#combine all assays into list
#assays have MDHD removed from problematic drugs, outlier assays removed, Worm_Length > 30 pixels removed (smallest anticipated L1 size)

anthDR.drugMDHD30.data.list <- as.list(ls(pattern = "30"))

View(anthDR.drugMDHD30.data.list) #check that list contains 10 assays 

anthDR.drugMDHD30 <- purrr::map(anthDR.drugMDHD30.data.list, concatenate.assays) %>%
  Reduce(rbind,.) %>%
  dplyr::mutate(drug = str_to_sentence(drug)) %>%
  dplyr::mutate(drug = if_else(drug == "Abamectin", true = "Abamectin", false = drug),
                drug = if_else(drug == "Albendazole", true = "Albendazole", false = drug),
                drug = if_else(drug == "Benomyl", true = "Benomyl", false = drug),
                drug = if_else(drug == "Closantel", true = "Closantel", false = drug),
                drug = if_else(drug == "Cry5b", true = "Cry5B", false = drug),
                drug = if_else(drug == "Dec", true = "Diethylcarbamazine", false = drug),
                drug = if_else(drug == "Derquantel", true = "Derquantel", false = drug),
                drug = if_else(drug == "Doramectin", true = "Doramectin", false = drug),
                drug = if_else(drug == "Emodepside", true = "Emodepside", false = drug),
                drug = if_else(drug == "Eprinomectin", true = "Eprinomectin", false = drug),
                drug = if_else(drug == "Fenbendazole", true = "Fenbendazole", false = drug),
                drug = if_else(drug == "Ivermectin", true = "Ivermectin", false = drug),
                drug = if_else(drug == "Levamisole", true = "Levamisole", false = drug),
                drug = if_else(drug == "Mebendazole", true = "Mebendazole", false = drug),
                drug = if_else(drug == "Milbemycin", true = "Milbemycin", false = drug),
                drug = if_else(drug == "Monepantel_ly33414916", true = "Monepantel LY33414916", false = drug),
                drug = if_else(drug == "Monepantel_ly3348298", true = "Monepantel LY3348298", false = drug),
                drug = if_else(drug == "Morantel", true = "Morantel", false = drug),
                drug = if_else(drug == "Moxidectin", true = "Moxidectin", false = drug),
                drug = if_else(drug == "Niridazole", true = "Niridazole", false = drug),
                drug = if_else(drug == "Oxamniquine", true = "Oxamniquine", false = drug),
                drug = if_else(drug == "Piperazine", true = "Piperazine", false = drug),
                drug = if_else(drug == "Pyrantel_cytrate", true = "Pyrantel", false = drug),
                drug = if_else(drug == "Pzq", true = "Praziquantel", false = drug),
                #drug = if_else(drug == "ricobendazole", true = "Ricobendazole", false = drug),
                drug = if_else(drug == "Selamectin", true = "Selamectin", false = drug),
                drug = if_else(drug == "Thiabendazole", true = "Thiabendazole", false = drug))
anthDR.drugMDHD30[which(anthDR.drugMDHD30$strain == "PD1074"),]$strain <- "N2"


anthDR.drugMDHD30 %>%
  dplyr::distinct(food_type, drug) %>% #identify distinct food sources
  dplyr::group_by(food_type, drug) %>% #group by distinct food sources
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = food_type, values_from = n) %>%
  dplyr::filter(!is.na(`15hHB101_20200813`)) #one food source used in assays


anthDR.drugMDHD30 <- anthDR.drugMDHD30 %>% 
  filter(!is.na(drug)) #remove rows with NA value in drug
 

control.wells.variance <- anthDR.drugMDHD30%>%
  dplyr::filter(concentration_um == 0) %>%
  dplyr::group_by(Metadata_Experiment, bleach) %>%
  dplyr::summarise(mean.n = mean(n),
                   sd.n = sd(n),
                   cv.n = (sd.n/mean.n),
                   mean.median_wormlength_um = mean(median_wormlength_um),
                   sd.median_wormlength_um = sd(median_wormlength_um),
                   cv.worm.length = sd(median_wormlength_um)/mean(median_wormlength_um))


cv.n.cutoff = 1 #not cutting off cv 

annotated.control.wells.variance <- control.wells.variance %>%
  dplyr::mutate(censor = if_else(condition = cv.n > cv.n.cutoff, true = "CENSOR", false = "NO CENSOR")) %>%
  dplyr::select(Metadata_Experiment, bleach, censor, cv.n, cv.worm.length) %>%
  dplyr::rename(censor.cv.n = cv.n, censor.cvWL = cv.worm.length)

#full data
save(anthDR.drugMDHD30, file = "output/master.anthDR.drugMDHD30.Rdata")

tx.nested.dose.responses.MDHD30 <- anthDR.drugMDHD30%>%
  dplyr::full_join(., annotated.control.wells.variance) %>%
  dplyr::filter(censor == "NO CENSOR") %>%
  dplyr::group_by(drug) %>%
  tidyr::nest() #creates a list-column of df

dose.response.summaries.drugspecificMDHD30 <- purrr::map2(tx.nested.dose.responses.MDHD30$data,
                                                          tx.nested.dose.responses.MDHD30$drug,
                                                          dose.response.summary)

tx.nested.dose.responses.MDHD30$drug

# #############
# #############
# # 50 pixels #
# #############
# #############
# 
# 
# #combine all assays into list
# #assays have MDHD removed from problematic drugs, outlier assays removed, Worm_Length > 30 pixels removed (smallest anticipated L1 size)
# 
# anthDR.drugMDHD50.data.list <- as.list(ls(pattern = "50"))
# 
# View(anthDR.drugMDHD50.data.list) #check that list contains 10 assays 
# 
# # remove(anthDR.drugMDHD30.data.list)
# library(dplyr)
# anthDR.drugMDHD50 <- purrr::map(anthDR.drugMDHD50.data.list, concatenate.assays) %>%
#   Reduce(rbind,.) %>%
#   dplyr::mutate(drug = str_to_sentence(drug)) %>%
#   dplyr::mutate(drug = if_else(drug == "Abamectin", true = "Abamectin", false = drug),
#                 drug = if_else(drug == "Albendazole", true = "Albendazole", false = drug),
#                 drug = if_else(drug == "Benomyl", true = "Benomyl", false = drug),
#                 drug = if_else(drug == "Closantel", true = "Closantel", false = drug),
#                 drug = if_else(drug == "Cry5b", true = "Cry5B", false = drug),
#                 drug = if_else(drug == "Dec", true = "Diethylcarbamazine", false = drug),
#                 drug = if_else(drug == "Derquantel", true = "Derquantel", false = drug),
#                 drug = if_else(drug == "Doramectin", true = "Doramectin", false = drug),
#                 drug = if_else(drug == "Emodepside", true = "Emodepside", false = drug),
#                 drug = if_else(drug == "Eprinomectin", true = "Eprinomectin", false = drug),
#                 drug = if_else(drug == "Fenbendazole", true = "Fenbendazole", false = drug),
#                 drug = if_else(drug == "Ivermectin", true = "Ivermectin", false = drug),
#                 drug = if_else(drug == "Levamisole", true = "Levamisole", false = drug),
#                 drug = if_else(drug == "Mebendazole", true = "Mebendazole", false = drug),
#                 drug = if_else(drug == "Milbemycin", true = "Milbemycin", false = drug),
#                 drug = if_else(drug == "Monepantel_LY33414916", true = "Monepantel LY33414916", false = drug),
#                 drug = if_else(drug == "Monepantel_LY3348298", true = "Monepantel LY3348298", false = drug),
#                 drug = if_else(drug == "Morantel", true = "Morantel", false = drug),
#                 drug = if_else(drug == "Moxidectin", true = "Moxidectin", false = drug),
#                 drug = if_else(drug == "Niridazole", true = "Niridazole", false = drug),
#                 drug = if_else(drug == "Oxamniquine", true = "Oxamniquine", false = drug),
#                 drug = if_else(drug == "Piperazine", true = "Piperazine", false = drug),
#                 drug = if_else(drug == "Pyrantel_cytrate", true = "Pyrantel citrate", false = drug),
#                 drug = if_else(drug == "Pzq", true = "Praziquantel", false = drug),
#                 #drug = if_else(drug == "ricobendazole", true = "Ricobendazole", false = drug),
#                 drug = if_else(drug == "Selamectin", true = "Selamectin", false = drug),
#                 drug = if_else(drug == "Thiabendazole", true = "Thiabendazole", false = drug))
# anthDR.drugMDHD50[which(anthDR.drugMDHD50$strain == "PD1074"),]$strain <- "N2"
# 
# 
# anthDR.drugMDHD30 %>%
#   dplyr::distinct(food_type, drug) %>% #identify distinct food sources
#   dplyr::group_by(food_type, drug) %>% #group by distinct food sources
#   dplyr::summarise(n = n()) %>%
#   tidyr::pivot_wider(names_from = food_type, values_from = n) %>%
#   dplyr::filter(!is.na(`15hHB101_20200813`)) #one food source used in assays
# 
# 
# anthDR.drugMDHD30 <- anthDR.drugMDHD30 %>% 
#   filter(!is.na(drug)) #remove rows with NA value in drug
# 
# 
# 
# control.wells.variance <- anthDR.drugMDHD30%>%
#   dplyr::filter(concentration_um == 0) %>%
#   dplyr::group_by(Metadata_Experiment, bleach) %>%
#   dplyr::summarise(mean.n = mean(n),
#                    sd.n = sd(n),
#                    cv.n = (sd.n/mean.n),
#                    mean.median_wormlength_um = mean(median_wormlength_um),
#                    sd.median_wormlength_um = sd(median_wormlength_um),
#                    cv.worm.length = sd(median_wormlength_um)/mean(median_wormlength_um))
# 
# 
# cv.n.cutoff = 1 #not cutting off cv 
# 
# annotated.control.wells.variance <- control.wells.variance %>%
#   dplyr::mutate(censor = if_else(condition = cv.n > cv.n.cutoff, true = "CENSOR", false = "NO CENSOR")) %>%
#   dplyr::select(Metadata_Experiment, bleach, censor, cv.n, cv.worm.length) %>%
#   dplyr::rename(censor.cv.n = cv.n, censor.cvWL = cv.worm.length)
# 
# 
# cv.well.n.block.figure <- annotated.control.wells.variance %>%
#   ggplot(., mapping = aes(x = Metadata_Experiment, y = censor.cv.n, colour = censor)) +
#   theme_bw(base_size = 9) +
#   geom_hline(yintercept = cv.n.cutoff, linetype = 3) +
#   geom_jitter() +
#   scale_color_manual(values = c("red","black"), name = "Data Censored") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         legend.position = "top",
#         panel.grid = element_blank()) +
#   labs(y = "Well Titer Coefficient of Variation",
#        x = "Assay")
# 
# cv.n.cv.WL.figure <- annotated.control.wells.variance %>%
#   ggplot(., mapping = aes(y = censor.cv.n, x = censor.cvWL, color = as.factor(censor))) +
#   theme_bw(base_size = 9) +
#   geom_hline(yintercept = cv.n.cutoff, linetype = 3) +
#   geom_point() +
#   theme(panel.grid = element_blank(),
#         axis.title.y = element_blank()) +
#   scale_color_manual(values = c("red","black"), name = "Data Censored") +
#   labs(y = "Well Titer Coefficient of Variation",
#        x = "Median Worm Length Coefficient of Variation")
# 
# cv.cleaning.supp.fig <- cowplot::plot_grid(cv.well.n.block.figure + theme(legend.position = "none"),
#                                            cv.n.cv.WL.figure, align = "h", labels = "AUTO")
# 
# ggsave(cv.cleaning.supp.fig + theme(plot.background = element_rect(fill = "white",colour = NA)),
#        filename = "figures/testcheck.png", width = 7.5, height = 3.5)
# 
# #full data
# save(anthDR.drugMDHD50, file = "data/anthDR.drugMDHD50.Rdata")
# 
# tx.nested.dose.responses.MDHD30 <- anthDR.drugMDHD30%>%
#   dplyr::full_join(., annotated.control.wells.variance) %>%
#   dplyr::filter(censor == "NO CENSOR") %>%
#   dplyr::group_by(drug) %>%
#   tidyr::nest() #creates a list-column of df

# #aos | need a cleaned dataset of all assays that is not nested to test heritability functions from lab 
# tx.nested.dose.responses_nonest.noMDHD <- dose.responses.df.noMDHD %>%
#   dplyr::full_join(., annotated.control.wells.variance) %>%
#   dplyr::filter(censor == "NO CENSOR") %>%
#   dplyr::group_by(drug, concentration_um, strain)

#tx.nested.dose.responses_nonest.noMDHD$drug

#cleaned data not nested 
#View(tx.nested.dose.responses_nonest)

# tx.nested.dose.responses is the cleaned data set. No models applied yet. 
#tx.nested.dose.responses$drug #list of 24 drugs 

#View(tx.nested.dose.responses$data[1][1]) #24 lists for each drug 

#view drug order 
#tx.nested.dose.responses$drug

# dose.response.summaries.drugspecificMDHD30 <- purrr::map2(tx.nested.dose.responses.MDHD30$data,
#                                                           tx.nested.dose.responses.MDHD30$drug,
#                                                           dose.response.summary)
# 
# 
# dose.response.summaries.drugspecificMDHD50 <- purrr::map2(tx.nested.dose.responses.MDHD50$data,
#                                                           tx.nested.dose.responses.MDHD50$drug,
#                                                           dose.response.summary)
# 
###############
###############
# RESUME HERE #
###############
###############

#individual dataframes to validate model selection for 30 pixels 

#always double check drug order before saving dataframes
# saving lists of cleaned.filtered.data for each drug and plotting to see how cleaning looks and test models
moxidectin.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30 [[1]][[8]]
write.csv(moxidectin.cleaned.filtered.df, "data/moxidectin.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

milbemycin.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30 [[2]][[8]]
write.csv(milbemycin.cleaned.filtered.df, "data/milbemycin.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

emodepside.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30 [[3]][[8]]
write.csv(emodepside.cleaned.filtered.df, "data/emodepside.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

closantel.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30 [[4]][[8]]
write.csv(closantel.cleaned.filtered.df, "data/closantel.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

benomyl.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30 [[5]][[8]]
write.csv(benomyl.cleaned.filtered.df, "data/benomyl.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

morantel.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30[[6]][[8]]
write.csv(morantel.cleaned.filtered.df, "data/morantel.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

derquantel.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30 [[7]][[8]]
write.csv(derquantel.cleaned.filtered.df, "data/derquantel.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

niridazole.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30 [[8]][[8]]
write.csv(niridazole.cleaned.filtered.df, "data/niridazole.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

dec.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30 [[9]][[8]]
write.csv(dec.cleaned.filtered.df, "data/dec.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

cry5b.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30 [[10]][[8]]
write.csv(cry5b.cleaned.filtered.df, "data/cry5b.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

pzq.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30 [[11]][[8]]
write.csv(pzq.cleaned.filtered.df, "data/pzq.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

piperazine.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30 [[12]][[8]]
write.csv(piperazine.cleaned.filtered.df, "data/piperazine.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

oxamniquine.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30 [[13]][[8]]
write.csv(oxamniquine.cleaned.filtered.df, "data/oxamniquine.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

monepantel_LY3348298.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30 [[14]][[8]]
write.csv(monepantel_LY3348298.cleaned.filtered.df, "data/monepantel_LY3348298.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

monepantel_LY33414916.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30[[15]][[8]]
write.csv(monepantel_LY33414916.cleaned.filtered.df, "data/monepantel_LY33414916.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

pyrantel_citrate.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30[[16]][[8]]
write.csv(pyrantel_citrate.cleaned.filtered.df, "data/pyrantel_citrate.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

eprinomectin.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30[[17]][[8]]
write.csv(eprinomectin.cleaned.filtered.df, "data/eprinomectin.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

selamectin.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30[[18]][[8]]
write.csv(selamectin.cleaned.filtered.df, "data/selamectin.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

fenbendazole.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30[[19]][[8]]
write.csv(fenbendazole.cleaned.filtered.df, "data/fenbendazole.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

mebendazole.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30[[20]][[8]]
write.csv(mebendazole.cleaned.filtered.df, "data/mebendazole.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

thiabendazole.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30[[21]][[8]]
write.csv(thiabendazole.cleaned.filtered.df, "data/thiabendazole.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

albendazole.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30[[22]][[8]]
write.csv(albendazole.cleaned.filtered.df, "data/albendazole.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

levamisole.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30[[23]][[8]]
write.csv(levamisole.cleaned.filtered.df, "data/levamisole.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

abamectin.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30 [[24]][[8]]
write.csv(abamectin.cleaned.filtered.df, "data/abamectin.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

ivermectin.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30[[25]][[8]]
write.csv(ivermectin.cleaned.filtered.df, "data/ivermectin.cleaned.filtered.df.30pixels.csv", row.names = FALSE)

doramectin.cleaned.filtered.df <- dose.response.summaries.drugspecificMDHD30[[26]][[8]]
write.csv(doramectin.cleaned.filtered.df, "data/doramectin.cleaned.filtered.df.30pixels.csv", row.names = FALSE)



tr.dose.response.summaries <- purrr::transpose(dose.response.summaries.drugspecificMDHD30) 


dose.response.parameter.summaries <- purrr::map(.x = tr.dose.response.summaries,
                                                .f = function(x){Reduce(rbind,x)})

View(dose.response.parameter.summaries)

###########
## LOAEL ##
###########
#LOAEL Not reported in manuscript. 

# overall.LOAEL.summary.df <- dose.response.parameter.summaries[[1]] %>%
#   # Reduce(rbind,.) %>%
#   dplyr::mutate(LOAEL = as.numeric(LOAEL))
# 
# strain.LOAEL.summary.df <- dose.response.parameter.summaries[[2]] %>%
#   # Reduce(rbind,.) %>%
#   dplyr::select(tx, LOAEL, strain) %>%
#   dplyr::group_by(tx, strain) %>%
#   dplyr::summarise(min(LOAEL)) %>%
#   dplyr::mutate(`min(LOAEL)` = as.numeric(`min(LOAEL)`)) %>%
#   dplyr::ungroup() %>%
#   tidyr::pivot_wider(values_from = `min(LOAEL)`, names_from = strain) %>%
#   dplyr::rename(Anthelmintic = tx) #%>% 
#  # dplyr::rename("N2" = "PD1074")
# 
# LOAEL.summary <- overall.LOAEL.summary.df %>%
#   dplyr::select(tx, LOAEL) %>%
#   dplyr::rename(Anthelmintic = tx,
#                 `Population-wide LOAEL` = LOAEL) %>%
#   dplyr::full_join(.,strain.LOAEL.summary.df)

# write.csv(LOAEL.summary, file = "output/LOAEL.summary.csv", quote = F, row.names = F) #aos checking loeal steps

#rename "drug" in anth.classes to "Anthelmintic"
# anth.classes <- anth.classes %>%
#   dplyr::rename(Anthelmintic = drug)

#create table with population wide and strain specific LOAEL
# LOAEL.table <- LOAEL.summary %>%
#   #dplyr::filter(!Anthelmintic %in% c("Triclabendazole","Ricobendazole")) %>%
#   dplyr::full_join(.,anth.classes) %>% #join by Anthelmintic
#   dplyr::arrange(big_class) %>%
#   dplyr::select(big_class, everything(), -diluent, -class) %>%
#   dplyr::rename(`Anthelmintic Class` = big_class)

#Table S2. LOAEL summary table
# write.csv(LOAEL.table, file = "output/table1_LOAEL_table_30pixelfilter.csv", row.names = F)


################################
### DOSE RESPONSE PARAMETERS ###
################################

#dose.responses.df = all easyXpress data with conc assays 

max.doses <- anthDR.drugMDHD30 %>% #input = dose.responses.df
  dplyr::select(drug, concentration_um) %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(max.dose = max(concentration_um))

all.DR.params <- dose.response.parameter.summaries[[3]]

View(all.DR.params)

##########
## EC10 ##
##########

EC10 <- dose.response.parameter.summaries[[3]] %>%
  dplyr::filter(metric == "EC10") %>%
  dplyr::mutate(null.model = is.na(Std..Error)) %>%
  dplyr::left_join(.,anth.classes %>%
                    dplyr::rename(drug = Anthelmintic))%>%
  dplyr::left_join(.,max.doses) %>%
  dplyr::mutate(above.max.conc = Estimate >= max.dose) %>%
  dplyr::mutate(big_class = gsub(big_class, pattern = "_", replacement = " "),
                flag = if_else(null.model == TRUE | Lower < 0,
                              true = TRUE, false = FALSE))

View(EC10)

above.max.concs <- EC10  %>%
  dplyr::group_by(drug, above.max.conc) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(above.max.conc == TRUE)

View(above.max.concs) 

flagged <- EC10  %>%
  dplyr::group_by(drug, null.model) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(null.model == TRUE)

EC10.unflagged <- EC10 %>%
  dplyr::filter(null.model == FALSE,
                above.max.conc == FALSE)

EC10.filtered <- EC10.unflagged %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n == 6) %>% #6 strains
  dplyr::left_join(.,EC10.unflagged) %>%
  dplyr::select(-n)

View(EC10.filtered)

##########
## EC50 ##
##########

EC50 <- dose.response.parameter.summaries[[3]] %>%
  #Reduce(rbind,.) %>%
  dplyr::filter(metric == "EC50") %>%
  dplyr::mutate(null.model = is.na(Std..Error)) %>%
  dplyr::left_join(.,anth.classes %>%
                     dplyr::rename(drug = Anthelmintic)) %>%
  dplyr::left_join(.,max.doses) %>%
  dplyr::mutate(above.max.conc = Estimate >= max.dose) %>%
  dplyr::mutate(big_class = gsub(big_class, pattern = "_", replacement = " "),
                flag = if_else(null.model == TRUE | above.max.conc == TRUE | Lower < 0, true = TRUE, false = FALSE))


flagged <- EC50 %>%
  dplyr::group_by(drug, flag) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = flag, values_from = n) %>%
  dplyr::filter(!is.na(`TRUE`)) %>%
  dplyr::left_join(.,EC50) %>%
  dplyr::ungroup() %>%
  dplyr::select(drug, strain, null.model, above.max.conc, flag) %>%
  data.frame()

above.max.concs <- EC50 %>%
  dplyr::group_by(drug, above.max.conc) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = above.max.conc, values_from = n) %>%
  dplyr::filter(!is.na(`TRUE`))

EC50.unflagged <- EC50 %>%
  dplyr::filter(flag == FALSE)

EC50.filtered <- EC50.unflagged %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n == 6) %>%
  dplyr::left_join(.,EC50.unflagged) %>%
  dplyr::select(-n)

##########
## EC90 ##
##########

EC90 <- dose.response.parameter.summaries[[3]] %>%
  # Reduce(rbind,.) %>%
  dplyr::filter(metric == "EC90") %>%
  dplyr::mutate(null.model = is.na(Std..Error)) %>%
  dplyr::left_join(.,anth.classes %>%
                     dplyr::rename(drug = Anthelmintic)) %>%
  dplyr::left_join(.,max.doses) %>%
  dplyr::mutate(above.max.conc = Estimate >= max.dose) %>%
  dplyr::mutate(big_class = gsub(big_class, pattern = "_", replacement = " "),
                flag = if_else(null.model == TRUE | above.max.conc == TRUE | Lower < 0, true = TRUE, false = FALSE))

flagged <- EC90 %>%
  dplyr::group_by(drug, flag) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = flag, values_from = n) %>%
  dplyr::filter(!is.na(`TRUE`)) %>%
  dplyr::left_join(.,EC90) %>%
  dplyr::ungroup() %>%
  dplyr::select(drug, strain, null.model, above.max.conc, flag) %>%
  data.frame()

above.max.concs <- EC90 %>%
  dplyr::group_by(drug, above.max.conc) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = above.max.conc, values_from = n) %>%
  dplyr::filter(!is.na(`TRUE`))

EC90.unflagged <- EC90 %>%
  dplyr::filter(flag == FALSE)

EC90.filtered <- EC90.unflagged %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n == 6) %>%
  dplyr::left_join(.,EC90.unflagged) %>%
  dplyr::select(-n)


reduced.EC10.filtered <- EC10.filtered %>%
  dplyr::select(drug, strain, metric, Estimate, Std..Error)
reduced.EC50.filtered <- EC50.filtered %>%
  dplyr::select(drug, strain, metric, Estimate, Std..Error)
reduced.EC90.filtered <- EC90.filtered %>%
  dplyr::select(drug, strain, metric, Estimate, Std..Error)

reduced.endpoint <- reduced.EC10.filtered %>%
  dplyr::full_join(., reduced.EC50.filtered) %>%
  dplyr::full_join(., reduced.EC90.filtered) %>%
  dplyr::arrange(metric, drug) %>%
  dplyr::rename(estimate = Estimate,
                SE = Std..Error)

write.csv(reduced.endpoint, "output/endpoints_30pixelfilter.csv", row.names = F)

 
View(EC10.filtered)

######################
### SLOPE (b) PLOT ###
######################

# Complete slope plot figure. Not in manuscript but available to look at all slopes across drugs.
#AOS added this section to compare 50 vs 30 pixels 
# filtered.slope.data <- all.DR.params %>% 
#   dplyr::filter(metric == "b") #make a dataframe that contains only slopes 
# 
# complete.slope.plot <- filtered.slope.data %>%
#   dplyr::mutate(drug = if_else(drug == "Dec", true = "Diethylcarbamazine", false = drug),
#                 drug = if_else(drug == "Monepantel_ly33414916", true = "Monepantel LY33414916", false = drug),
#                 drug = if_else(drug == "Monepantel_ly3348298", true = "Monepantel LY3348298", false = drug),
#                 drug = if_else(drug == "Pyrantel_cytrate", true = "Pyrantel citrate", false = drug),
#                 drug = if_else(drug == "Pzq", true = "Praziquantel", false = drug)) %>%
#   ggplot(., mapping = aes(x = strain, y = Estimate, ymin = Lower, ymax = Upper,
#                           color = strain)) +
#   theme_bw(base_size = 11) +
#   geom_pointrange(position = position_quasirandom(groupOnX = FALSE),
#                   size = 0.5) +
#   scale_color_manual(values = strain_colors, name = "Strain") +
# facet_wrap(drug ~., scales = "free", nrow = 7) +
#   theme(strip.text = element_text(angle = 0, size = 10.5),
#         axis.text = element_text(color = "black"),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank(),
#         legend.position = "top") +
#   labs(y = "Slope Estimate",
#        x = "Strain")
# 
# complete.slope.plot
# 
# ggsave(complete.slope.plot + theme(plot.background = element_rect(fill = "white",colour = NA)),
#        filename = "output/complete.slope.plot_30pixels.png", width = 11, height = 8.5)



#################
### EC10 PLOT ###
#################

EC10.filtered[which(EC10.filtered$strain == "PD1074"),]$strain <- "N2"


# ###########################################
# ### COMPLETE EC10 PLOT - BENZIMIDAZOLES ###
# ###########################################
# complete.EC10.plot.BZ <- EC10.filtered %>%
#   dplyr::filter(!drug %in% c("triclabendazole","zts_00507244", "zts_00507244",
#                                "zts_0007437", "zts_7027178","zts_7027187", "Abamectin",
#                              "Closantel","Cry5B", "Diethylcarbamazine", "Derquantel","Doramectin", "Emodepside",
#                              "Eprinomectin", "Ivermectin", "Levamisole", "Milbemycin", "Monepantel LY33414916", 
#                              "Monepantel LY3348298", "Morantel", "Moxidectin", "Niridazole", "Oxamniquine", "Piperazine",
#                              "Pyrantel", "Praziquantel", "Selamectin" )) %>% # to make plot with just BZs represented
#   dplyr::mutate(drug = if_else(drug == "albendazole", true = "Albendazole", false = drug),
#                 drug = if_else(drug == "benomyl", true = "Benomyl", false = drug),
#                 drug = if_else(drug == "fenbendazole", true = "Fenbendazole", false = drug),
#                 drug = if_else(drug == "mebendazole", true = "Mebendazole", false = drug),
#                 drug = if_else(drug == "thiabendazole", true = "Thiabendazole", false = drug)) %>% 
#   ggplot(., mapping = aes(x = strain, y = Estimate, ymin = Lower, ymax = Upper,
#                           color = strain)) + 
#   theme_bw(base_size = 11) + #text size
#   geom_pointrange(position = position_quasirandom(groupOnX = FALSE),
#                   size = 0.5) +
#   scale_color_manual(values = strain_colors, name = "Strain") + 
#   facet_wrap(drug ~., scales = "free", nrow = 1) +
#   theme(strip.text = element_text(angle = 0, size = 9),
#         axis.text = element_text(color = "black"),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank(),
#         legend.position = "top") + 
#   labs(y = "EC10 Estimate (µM)")
# 
# complete.EC10.plot.BZ
# 
# ggsave(complete.EC10.plot.BZ + theme(plot.background = element_rect(fill = "white",colour = NA)),
#        filename = "output/complete.EC10.plot_30pixels_Bzs.png", width = 6, height = 2.23)
# 
# View(EC10.filtered)
# 
# 

# #################################################
# ### COMPLETE EC10 PLOT - MACROCYCLIC LACTONES ###
# #################################################
# complete.EC10.plot.ML <- EC10.filtered %>%
#   dplyr::filter(!drug %in% c("triclabendazole","zts_00507244", "zts_00507244",
#                              "zts_0007437", "zts_7027178","zts_7027187",
#                              "Closantel","Cry5B", "Diethylcarbamazine", "Derquantel", "Emodepside",
#                               "Levamisole", "Monepantel LY33414916","Albendazole", "Benomyl", "Fenbendazole",
#                              "Monepantel LY3348298", "Morantel",  "Niridazole", "Oxamniquine", "Piperazine",
#                              "Pyrantel", "Praziquantel", "Mebendazole","Thiabendazole" )) %>% # to make plot with just BZs represented
#   dplyr::mutate(drug = if_else(drug == "abamectin", true = "Abamectin", false = drug),
#                 drug = if_else(drug == "doramectin", true = "Doramectin", false = drug),
#                 drug = if_else(drug == "eprinomectin", true = "Eprinomectin", false = drug),
#                 drug = if_else(drug == "ivermectin", true = "Ivermectin", false = drug),
#                 drug = if_else(drug == "milbemycin", true = "Milbemycin", false = drug),
#                 drug = if_else(drug == "moxidectin", true = "Moxidectin", false = drug),
#                 drug = if_else(drug == "selamectin", true = "Selamectin", false = drug)) %>% 
#   ggplot(., mapping = aes(x = strain, y = Estimate, ymin = Lower, ymax = Upper,
#                           color = strain)) + 
#   theme_bw(base_size = 11) + 
#   geom_pointrange(position = position_quasirandom(groupOnX = FALSE),
#                   size = 0.5) +
#   scale_color_manual(values = strain_colors, name = "Strain") + 
#   facet_wrap(drug ~., scales = "free", nrow = 1) + 
#   theme(strip.text = element_text(angle = 0, size = 9),
#         axis.text = element_text(color = "black"),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank(),
#         legend.position = "top") + 
#   labs(y = "EC10 Estimate (µM)",
#        x = "Strain")
# 
# complete.EC10.plot.ML
# 
# ggsave(complete.EC10.plot.ML+ theme(plot.background = element_rect(fill = "white",colour = NA)),
#        filename = "output/complete.EC10.plot_30pixels_MLs.png", width = 9, height = 2.24)
# 
# View(EC10.filtered)



# ###########################################
# ### COMPLETE EC10 PLOT - NACHAS ###
# ###########################################
# complete.EC10.plot.NACHA <- EC10.filtered %>%
#   dplyr::filter(!drug %in% c("triclabendazole","zts_00507244", "zts_00507244",
#                              "zts_0007437", "zts_7027178","zts_7027187",
#                              "Closantel","Cry5B", "Diethylcarbamazine", "Emodepside",
#                             "Monepantel LY33414916","Albendazole", "Benomyl", "Fenbendazole",
#                              "Monepantel LY3348298", "Niridazole", "Oxamniquine", "Piperazine",
#                              "Praziquantel", "Mebendazole","Thiabendazole", "Abamectin",
#                             "Doramectin","Eprinomectin","Ivermectin","Milbemycin","Moxidectin","Selamectin")) %>% # to make plot with just BZs represented
#   dplyr::mutate(drug = if_else(drug == "derquantel", true = "Derquantel", false = drug),
#                 drug = if_else(drug == "levamisole", true = "Levamisole", false = drug),
#                 drug = if_else(drug == "morantel", true = "Morantel", false = drug),
#                 drug = if_else(drug == "pyrantel_cytrate", true = "Pyrantel", false = drug)) %>% 
#   ggplot(., mapping = aes(x = strain, y = Estimate, ymin = Lower, ymax = Upper,
#                           color = strain)) + 
#   theme_bw(base_size = 11) + 
#   geom_pointrange(position = position_quasirandom(groupOnX = FALSE),
#                   size = 0.5) +
#   scale_color_manual(values = strain_colors, name = "Strain") + 
#   facet_wrap(drug ~., scales = "free", nrow = 1) + 
#   theme(strip.text = element_text(angle = 0, size = 9),
#         axis.text = element_text(color = "black"),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank(),
#         legend.position = "top") + 
#   labs(y = "EC10 Estimate (µM)",
#        x = "Strain")
# 
# complete.EC10.plot.NACHA
# 
# ggsave(complete.EC10.plot.NACHA + theme(plot.background = element_rect(fill = "white",colour = NA)),
#        filename = "output/complete.EC10.plot_30pixels_NACHAs.png", width = 6, height = 2.23)
# 
# View(EC10.filtered)


#EC10 estimates table 

EC10.estimates.table <- EC10.filtered %>%
  dplyr::rename(Anthelmintic = drug,
                EC10 = Estimate,
                SE = Std..Error) %>%
  dplyr::select(Anthelmintic, EC10, SE, strain) %>%
  dplyr::mutate(EC10 = round(EC10,4),
                SE = round(SE,4)) %>% 
  tidyr::unite("EC10 (µM)",EC10:SE, 
               sep = " ± ") %>%
  tidyr::pivot_wider(names_from = strain, values_from = `EC10 (µM)`)

View(EC10.estimates.table)

#Supplemental Table 1
write.csv(EC10.estimates.table, "output/EC10.estimates.table.30pixels.csv", row.names = F, quote = F)

View(dose.response.parameter.summaries)

n.EC.comp.tests <- dose.response.parameter.summaries[[6]] %>%
  dplyr::filter(drug %in% EC10.filtered$drug) %>%
  tidyr::separate(comps, c("strains","fracs"), sep = ":") %>%
  tidyr::separate(strains, c("strain1","strain2"), sep = "/") %>%
  dplyr::filter(strain2 == "N2" | strain1 == "N2") %>%
  nrow()

View(n.EC.comp.tests)

BF <- 0.05/n.EC.comp.tests

options(scipen = 999999)

View(EC10.filtered)

View(dose.response.parameter.summaries[[6]])

EC10.relative.potency <- dose.response.parameter.summaries[[6]] %>%
  tidyr::separate(comps, c("strains","fracs"), sep = ":") %>%
  tidyr::separate(strains, c("strain1","strain2"), sep = "/") %>%
  dplyr::filter(drug %in% unique(EC10.filtered$drug),
                strain2 == "N2" | strain1 == "N2") %>%
  dplyr::mutate(N2.normed.estimate = if_else(strain1 == "N2", true = 1/Estimate, false = Estimate)) %>%
  dplyr::mutate(N2.normed.diff = N2.normed.estimate-1,
                focal.strain = if_else(strain1 == "N2", true = strain2, false = strain1)) %>%
  dplyr::left_join(.,anth.classes %>%
                     dplyr::rename(drug = Anthelmintic)) %>%
  dplyr::mutate(big_class = gsub(big_class, pattern = "_", replacement = " ")) %>%
  dplyr::mutate(start = 1,
                sig = if_else(condition = p.value < BF, true = "SIG", false = "NONSIG"))

View(EC10.relative.potency)

#Supplementary Table 4 
write.csv(EC10.relative.potency, "output/EC10.relative.potency.table_2.csv")

  
EC10.relative.potency.figure <- EC10.relative.potency %>%
    ggplot(., mapping = aes(y = drug, x = N2.normed.diff, 
                          xmin = N2.normed.diff-Std..Error, 
                          xmax = N2.normed.diff+Std..Error,
                          color = focal.strain,
                          alpha = sig)) + 
  theme_bw(base_size = 11) +
  geom_vline(xintercept = 0, linetype = 3, colour = "orange") + 
  geom_pointrange(position = position_dodge(width = 0.2),
                  size = 0.5) + 
  scale_color_manual(values = strain_colors[c(1,3:6)], name = "Strain") +
  scale_alpha_manual(values = c(0.15,1), guide = "none") +
  facet_grid(big_class~., scales = "free", space = "free") + 
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_text(color = "black"),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "none") + 
  labs(x = "Resistance Compared to Reference (N2)")


fig.3.legend <- cowplot::get_legend(complete.slope.plot)


EC10.relative.potency.figure

complete.EC10.figure <- cowplot::plot_grid(fig.3.legend, EC10.relative.potency.figure,  rel_heights = c(2, 20), ncol = 1)
ggsave(complete.EC10.figure, filename = "output/Supp.Fig.3.EC10.relative.potency.ALL.png", width = 12, height = 8)


complete.EC10.plot

EC10.fig.2 <- cowplot::plot_grid(complete.EC10.plot +  
                                   theme(strip.background = element_blank(),strip.text.y = element_blank(),
                                         plot.margin = unit(c(0.5,0.5,0.25,0.5), "cm"),
                                         legend.position = "none"),
                                 EC10.relative.potency.figure + 
                                   theme(plot.margin = unit(c(0.5,0.5,0.25,0.5), "cm"),
                                         legend.position = "none"), 
                                 rel_widths = c(1.2,1), labels = "AUTO")

complete.EC10.fig.2.plot <- cowplot::plot_grid(EC10.fig.2, 
                                               fig.2.legend, 
                                               rel_heights = c(20,2), ncol = 1)

complete.EC10.fig.2.plot 

ggsave(complete.EC10.fig.2.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), filename = "output/complete.EC10.fig.2.plot.30pixels.below2.png", width = 16, height = 8)


## EC10 Relative Potency Tests ##


#slope code from most recent github update - 2022/08/31 
slopes.filtered <- EC10.filtered %>%
  dplyr::select(strain, drug) %>%
  dplyr::left_join(.,all.DR.params) %>%
  dplyr::left_join(.,anth.classes %>%
                     dplyr::rename(drug = Anthelmintic)) %>%
  dplyr::filter(metric == "b") %>%
  dplyr::mutate(big_class = gsub(big_class, pattern = "_", replacement = " "))


slope.estimates.table <- slopes.filtered %>%
  dplyr::rename(Anthelmintic = drug,
                Slope = Estimate,
                SE = Std..Error) %>%
  dplyr::select(Anthelmintic, Slope, SE, strain) %>%
  dplyr::mutate(Slope = round(Slope,2),
                 SE = round(SE,2)) %>% 
  #               Anthelmintic = if_else(Anthelmintic == "2,4-D", true = "2_4-D", false = Anthelmintic)) %>% 
  tidyr::unite("Slope",Slope:SE, 
               sep = " ± ") %>%
  tidyr::pivot_wider(names_from = strain, values_from = Slope)
write.csv(slope.estimates.table, "output/slope.estimates.table.csv", row.names = F, quote = F)

n.slope.comp.tests <- dose.response.parameter.summaries[[7]] %>%
  tidyr::separate(comps, c("strain1","strain2"), sep = "/") %>%
  dplyr::mutate(strain1 = gsub(strain1, pattern = "strain", replacement = "")) %>%
  dplyr::mutate(strain2 = gsub(strain2, pattern = "strain", replacement = "")) %>% 
  dplyr::filter(drug %in% unique(EC10.filtered$drug),
                strain2 == "N2" | strain1 == "N2") %>%
  dplyr::mutate(N2.normed.estimate = if_else(strain1 == "N2", true = 1/Estimate, false = Estimate)) %>% 
  nrow()

BF <- 0.05/n.slope.comp.tests

options(scipen = 999999)

relative.slope <- dose.response.parameter.summaries[[7]] %>%
  tidyr::separate(comps, c("strain1","strain2"), sep = "/") %>%
  dplyr::mutate(strain1 = gsub(strain1, pattern = "strain", replacement = "")) %>%
  dplyr::mutate(strain2 = gsub(strain2, pattern = "strain", replacement = "")) %>% 
  dplyr::filter(drug %in% unique(EC10.filtered$drug),
                strain2 == "N2" | strain1 == "N2") %>%
  dplyr::mutate(N2.normed.estimate = if_else(strain1 == "N2", true = 1/Estimate, false = Estimate)) %>%
  dplyr::mutate(N2.normed.diff = N2.normed.estimate-1,
                focal.strain = if_else(strain1 == "N2", true = strain2, false = strain1)) %>%
  dplyr::left_join(.,anth.classes %>%
                     dplyr::rename(drug = Anthelmintic)) %>%
  dplyr::mutate(big_class = gsub(big_class, pattern = "_", replacement = " ")) %>%
  dplyr::mutate(start = 1,
                sig = if_else(condition = p.value < BF, true = "SIG", false = "NONSIG"))

## slope Relative Potency Tests ##
slope.relative.potency.supp.table <- dose.response.parameter.summaries[[7]] %>%
  tidyr::separate(comps, c("strain1","strain2"), sep = "/") %>%
  dplyr::mutate(strain1 = gsub(strain1,pattern = "strain",replacement = ""),
                strain2 = gsub(strain2,pattern = "strain",replacement = "")) %>%
  dplyr::mutate(`p < 0.05` = if_else(condition = p.value < 0.05, true = "T", false = "F"),
                p.value = round(p.value,6),
                p.value = if_else(p.value < 0.000001, true =  "<< 0.00001", false = as.character(p.value))) %>%
  dplyr::rename(`Relative Potency Estimate` = Estimate,
                SE = Std..Error,
                Anthelmintic = drug) %>%
  dplyr::select(Anthelmintic, strain1, strain2, everything()) %>%
  data.frame()

#supplementary table 4 
write.csv(slope.relative.potency.supp.table, "output/slope.relative.potency.supp.table.csv", row.names = F, quote = F)


###################################################################
##### Figure 3. A + B SLOPE AND EC10 PLOTTING - BENZIMIDAZOLES ####
###################################################################
complete.EC10.plot.BZ <- EC10.filtered %>%
  dplyr::filter(!drug %in% c("Abamectin",
                             "Closantel","Cry5B", "Diethylcarbamazine", "Derquantel","Doramectin", "Emodepside",
                             "Eprinomectin", "Ivermectin", "Levamisole", "Milbemycin", "Monepantel LY33414916", 
                             "Monepantel LY3348298", "Morantel", "Moxidectin", "Niridazole", "Oxamniquine", "Piperazine",
                             "Pyrantel", "Praziquantel", "Selamectin" )) %>% # to make plot with just BZs represented
  dplyr::mutate(drug = if_else(drug == "albendazole", true = "Albendazole", false = drug),
                drug = if_else(drug == "benomyl", true = "Benomyl", false = drug),
                drug = if_else(drug == "fenbendazole", true = "Fenbendazole", false = drug),
                drug = if_else(drug == "mebendazole", true = "Mebendazole", false = drug),
                drug = if_else(drug == "thiabendazole", true = "Thiabendazole", false = drug)) %>% 
  ggplot(., mapping = aes(x = strain, y = Estimate, ymin = Lower, ymax = Upper,
                          color = strain)) + 
  theme_bw(base_size = 11) + #text size
  geom_pointrange(position = position_quasirandom(groupOnX = FALSE),
                  size = 0.5) +
  scale_color_manual(values = strain_colors, name = "Strain") + 
  facet_wrap(drug ~., scales = "free", nrow = 1) + 
  theme(strip.text = element_text(angle = 0, size = 9),
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") + 
  labs(y = "EC10 Estimate (µM)",
       x = " ")

#BZ EC10 plot with no legend
complete.EC10.plot.BZ

#Legend for BZ EC10 and slope plot 
fig.3.legend_BZ <- cowplot::get_legend(complete.EC10.plot.BZ)

# ggsave(complete.EC10.plot.BZ + theme(plot.background = element_rect(fill = "white",colour = NA)),
#        filename = "output/complete.EC10.plot_30pixels_Bzs.png", width = 6, height = 2.23)

complete.slope.plot.BZ <- slopes.filtered %>%
  dplyr::filter(!drug %in% c("Abamectin",
                             "Closantel","Cry5B", "Diethylcarbamazine", "Derquantel","Doramectin", "Emodepside",
                             "Eprinomectin", "Ivermectin", "Levamisole", "Milbemycin", "Monepantel LY33414916", 
                             "Monepantel LY3348298", "Morantel", "Moxidectin", "Niridazole", "Oxamniquine", "Piperazine",
                             "Pyrantel", "Praziquantel", "Selamectin" )) %>% # to make plot with just BZs represented
  dplyr::mutate(drug = if_else(drug == "albendazole", true = "Albendazole", false = drug),
                drug = if_else(drug == "benomyl", true = "Benomyl", false = drug),
                drug = if_else(drug == "fenbendazole", true = "Fenbendazole", false = drug),
                drug = if_else(drug == "mebendazole", true = "Mebendazole", false = drug),
                drug = if_else(drug == "thiabendazole", true = "Thiabendazole", false = drug)) %>% 
  ggplot(., mapping = aes(x = strain, y = Estimate, ymin = Estimate-Std..Error, ymax = Estimate+Std..Error,
                          color = strain)) + 
  theme_bw(base_size = 11) + 
  geom_pointrange(position = position_quasirandom(groupOnX = FALSE),
                  size = 0.5) +
  scale_color_manual(values = strain_colors, name = "Strain") + 
  facet_wrap(drug ~., scales = "free", nrow = 1) + 
  theme(strip.text = element_text(angle = 0, size = 9),
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(), legend.position="none") +
  labs(y = "Estimated Slope",
       x = "Strain")

complete.slope.plot.BZ

summary.plot.EC10.slope.BZ <- cowplot::plot_grid(complete.EC10.plot.BZ, complete.slope.plot.BZ, ncol = 1, align="v", rel_heights = c(0.5,0.5))

summary.plot.EC10.slope.BZ

BZ.LEGEND.EC10.SLOPE.PLOT <- cowplot::plot_grid(fig.3.legend_BZ, summary.plot.EC10.slope.BZ,ncol = 1, rel_heights = c(0.1,                                                                                                   0.4,                                                                                             0.5))

BZ.LEGEND.EC10.SLOPE.PLOT

#Figure 3 A + B 
ggsave(BZ.LEGEND.EC10.SLOPE.PLOT + theme(plot.background = element_rect(fill = "white",colour = NA)), 
        filename = "output/BZ.LEGEND.EC10.SLOPE.PLOT.png", width = 8.75, height = 4) #8.75 max height for PLOS

#These steps keep A and B relatively the same height - width is difficult because of digit differences on y axis 
# BZ.EC10 <- cowplot::plot_grid(complete.EC10.plot.BZ + theme(legend.position = "none"), rel_widths = c(0.5))
# 
# BZ.SLOPE <- cowplot::plot_grid(complete.slope.plot.BZ  + 
#                           theme(legend.position = "none"), rel_widths = c(0.5))
# 
# BZ.LEGEND.EC10.SLOPE.PLOT <- cowplot::plot_grid(fig.3.legend_BZ, BZ.EC10, BZ.SLOPE,ncol = 1, rel_heights = c(0.1, 
#                                                                                                  0.4, 
#                                                                                                  0.4))
#                                                 
# BZ.LEGEND.EC10.SLOPE.PLOT
# 
# ggsave(BZ.LEGEND.EC10.SLOPE.PLOT + theme(plot.background = element_rect(fill = "white",colour = NA)), 
#        filename = "manuscript_figures/fig.3.BZ.EC10.SLOPE.png", width = 7, height = 4)


######################################################################
### Figure 5. A + B SLOPE AND EC10 PLOTTING - MACROCYCLIC LACTONES ###
######################################################################
complete.EC10.plot.ML <- EC10.filtered %>%
  dplyr::filter(!drug %in% c("triclabendazole","zts_00507244", "zts_00507244",
                             "zts_0007437", "zts_7027178","zts_7027187",
                             "Closantel","Cry5B", "Diethylcarbamazine", "Derquantel", "Emodepside",
                             "Levamisole", "Monepantel LY33414916","Albendazole", "Benomyl", "Fenbendazole",
                             "Monepantel LY3348298", "Morantel",  "Niridazole", "Oxamniquine", "Piperazine",
                             "Pyrantel citrate", "Praziquantel", "Mebendazole","Thiabendazole" )) %>% # to make plot with just BZs represented
  dplyr::mutate(drug = if_else(drug == "abamectin", true = "Abamectin", false = drug),
                drug = if_else(drug == "doramectin", true = "Doramectin", false = drug),
                drug = if_else(drug == "eprinomectin", true = "Eprinomectin", false = drug),
                drug = if_else(drug == "ivermectin", true = "Ivermectin", false = drug),
                drug = if_else(drug == "milbemycin", true = "Milbemycin", false = drug),
                drug = if_else(drug == "moxidectin", true = "Moxidectin", false = drug),
                drug = if_else(drug == "selamectin", true = "Selamectin", false = drug)) %>% 
  ggplot(., mapping = aes(x = strain, y = Estimate, ymin = Lower, ymax = Upper,
                          color = strain, group=drug)) + 
  theme_bw(base_size = 11) + 
  geom_pointrange(position = position_quasirandom(groupOnX = FALSE),
                  size = 0.5) +
  scale_color_manual(values = strain_colors, name = "Strain") + 
  facet_wrap(drug ~., scales = "free", nrow = 1) + 
  theme(strip.text = element_text(angle = 0, size = 9),
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") + 
  scale_y_continuous(breaks = scales::breaks_extended(n=5)) +
  labs(y = "EC10 Estimate (µM)",
       x = " ")

complete.EC10.plot.ML

# Save ML EC10 plot
# ggsave(complete.EC10.plot.ML+ theme(plot.background = element_rect(fill = "white",colour = NA)), 
#        filename = "output/complete.EC10.plot_30pixels_MLs.png", width = 9, height = 2.24)

complete.slope.plot.ML <- slopes.filtered %>%
  dplyr::filter(!drug %in% c("triclabendazole","zts_00507244", "zts_00507244",
                             "zts_0007437", "zts_7027178","zts_7027187",
                             "Closantel","Cry5B", "Diethylcarbamazine", "Derquantel", "Emodepside",
                             "Levamisole", "Monepantel LY33414916","Albendazole", "Benomyl", "Fenbendazole",
                             "Monepantel LY3348298", "Morantel",  "Niridazole", "Oxamniquine", "Piperazine",
                             "Pyrantel citrate", "Praziquantel", "Mebendazole","Thiabendazole" )) %>% # to make plot with just BZs represented
  dplyr::mutate(drug = if_else(drug == "abamectin", true = "Abamectin", false = drug),
                drug = if_else(drug == "doramectin", true = "Doramectin", false = drug),
                drug = if_else(drug == "eprinomectin", true = "Eprinomectin", false = drug),
                drug = if_else(drug == "ivermectin", true = "Ivermectin", false = drug),
                drug = if_else(drug == "milbemycin", true = "Milbemycin", false = drug),
                drug = if_else(drug == "moxidectin", true = "Moxidectin", false = drug),
                drug = if_else(drug == "selamectin", true = "Selamectin", false = drug)) %>% 
  ggplot(., mapping = aes(x = strain, y = Estimate, ymin = Estimate-Std..Error, ymax = Estimate+Std..Error,
                          color = strain)) + 
  theme_bw(base_size = 11) + 
  geom_pointrange(position = position_quasirandom(groupOnX = FALSE),
                  size = 0.5) +
  scale_color_manual(values = strain_colors, name = "Strain") + 
  facet_wrap(drug ~., scales = "free", nrow = 1) + 
  theme(strip.text = element_text(angle = 0, size = 9),
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(), legend.position="none") +
  labs(y = "Estimated Slope",
       x = "Strain")

 summary.plot.EC10.slope.ML <- cowplot::plot_grid(complete.EC10.plot.ML, complete.slope.plot.ML,  ncol = 1, align="v",rel_heights = c(0.5,0.5))
 summary.plot.EC10.slope.ML 

ML.LEGEND.EC10.SLOPE.PLOT <- cowplot::plot_grid(fig.3.legend_BZ,  summary.plot.EC10.slope.ML ,ncol = 1, rel_heights = c(0.1,                                                                                                   0.4, 
                                                                                                                      0.5))

ggsave(ML.LEGEND.EC10.SLOPE.PLOT + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "output/summary.plot.EC10.slope.ML.png", width = 8.75, height = 4) #8.75 max heigh PLOS

#These steps keep A and B relatively the same height - width is difficult because of digit differences on y axis 
# ML.EC10 <- cowplot::plot_grid(complete.EC10.plot.ML + theme(legend.position = "none"), rel_widths = c(0.5))
# 
# ML.SLOPE <- cowplot::plot_grid(complete.slope.plot.ML  + 
#                                  theme(legend.position = "none"), rel_widths = c(0.5))
# 
# ML.LEGEND.EC10.SLOPE.PLOT <- cowplot::plot_grid(fig.3.legend_BZ, ML.EC10, ML.SLOPE,ncol = 1, rel_heights = c(0.1, 
#                                                                                                              0.4, 
#                                                                                                              0.4))
# 
# ML.LEGEND.EC10.SLOPE.PLOT
# 
# ggsave(ML.LEGEND.EC10.SLOPE.PLOT + theme(plot.background = element_rect(fill = "white",colour = NA)), 
#        filename = "manuscript_figures/fig.5.ML.EC10.SLOPE.png", width = 7, height = 4)



########################################################
### Figure 7. A + B SLOPE AND EC10 PLOTTING - NAChAs ###
########################################################
complete.EC10.plot.NAChAs <- EC10.filtered %>%
  dplyr::filter(!drug %in% c("Closantel","Cry5B", "Diethylcarbamazine", "Emodepside",
                             "Monepantel LY33414916","Albendazole", "Benomyl", "Fenbendazole",
                             "Monepantel LY3348298", "Niridazole", "Oxamniquine", "Piperazine",
                             "Praziquantel", "Mebendazole","Thiabendazole", "Abamectin",
                             "Doramectin","Eprinomectin","Ivermectin","Milbemycin","Moxidectin","Selamectin", "Derquantel")) %>% # to make plot with just BZs represented
  dplyr::mutate(drug = if_else(drug == "levamisole", true = "Levamisole", false = drug),
                drug = if_else(drug == "morantel", true = "Morantel", false = drug),
                drug = if_else(drug == "pyrantel", true = "Pyrantel", false = drug)) %>% 
  ggplot(., mapping = aes(x = strain, y = Estimate, ymin = Lower, ymax = Upper,
                          color = strain, group=drug)) + 
  theme_bw(base_size = 11) + 
  geom_pointrange(position = position_quasirandom(groupOnX = FALSE),
                  size = 0.5) +
  scale_color_manual(values = strain_colors, name = "Strain") + 
  facet_wrap(drug ~., scales = "free", nrow = 1) + 
  theme(strip.text = element_text(angle = 0, size = 9),
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top") + 
  scale_y_continuous(breaks = scales::breaks_extended(n=5)) +
  labs(y = "EC10 Estimate (µM)",
       x = " ")


complete.slope.plot.NAChAs <- slopes.filtered %>%
  dplyr::filter(!drug %in% c("Closantel","Cry5B", "Diethylcarbamazine", "Emodepside",
                             "Monepantel LY33414916","Albendazole", "Benomyl", "Fenbendazole",
                             "Monepantel LY3348298", "Niridazole", "Oxamniquine", "Piperazine",
                             "Praziquantel", "Mebendazole","Thiabendazole", "Abamectin",
                             "Doramectin","Eprinomectin","Ivermectin","Milbemycin","Moxidectin","Selamectin", "Derquantel")) %>% # to make plot with just BZs represented
  dplyr::mutate(drug = if_else(drug == "levamisole", true = "Levamisole", false = drug),
                drug = if_else(drug == "morantel", true = "Morantel", false = drug),
                drug = if_else(drug == "pyrantel_cytrate", true = "Pyrantel citrate", false = drug)) %>% 
  ggplot(., mapping = aes(x = strain, y = Estimate, ymin = Estimate-Std..Error, ymax = Estimate+Std..Error,
                          color = strain)) + 
  theme_bw(base_size = 11) + 
  geom_pointrange(position = position_quasirandom(groupOnX = FALSE),
                  size = 0.5) +
  scale_color_manual(values = strain_colors, name = "Strain") + 
  facet_wrap(drug ~., scales = "free", nrow = 1) + 
  theme(strip.text = element_text(angle = 0, size = 9),
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(), legend.position="none") +
  labs(y = "Estimated Slope",
       x = "Strain")

summary.plot.EC10.slope.NAChAs <- cowplot::plot_grid(complete.EC10.plot.NAChAs, complete.slope.plot.NAChAs , ncol = 1, align="v",rel_heights = c(0.6,0.5))
summary.plot.EC10.slope.NAChAs

NACHAs.LEGEND.EC10.SLOPE.PLOT <- cowplot::plot_grid(fig.3.legend_BZ,  summary.plot.EC10.slope.NAChAs ,ncol = 1, rel_heights = c(0.1,                                                                                                   0.4, 
                                                                                                                        0.6))
ggsave(NACHAs.LEGEND.EC10.SLOPE.PLOT + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "output/summary.plot.EC10.slope.NAChAs.png", width = 8.75, height = 4) #8.75 max height PLOS

ggsave(summary.plot.EC10.slope.NAChAs,
       filename = "output/summary.plot.EC10.slope.NAChAs.png", width = 8.75, height = 6) #8.75 max height PLOS



###########################################
### COMPLETE EC10 PLOT - NACHAS ###
###########################################
complete.EC10.plot.NACHA <- EC10.filtered %>%
  dplyr::filter(!drug %in% c("Closantel","Cry5B", "Diethylcarbamazine", "Emodepside",
                             "Monepantel LY33414916","Albendazole", "Benomyl", "Fenbendazole",
                             "Monepantel LY3348298", "Niridazole", "Oxamniquine", "Piperazine",
                             "Praziquantel", "Mebendazole","Thiabendazole", "Abamectin",
                             "Doramectin","Eprinomectin","Ivermectin","Milbemycin","Moxidectin","Selamectin")) %>% # to make plot with just BZs represented
  dplyr::mutate(drug = if_else(drug == "derquantel", true = "Derquantel", false = drug),
                drug = if_else(drug == "levamisole", true = "Levamisole", false = drug),
                drug = if_else(drug == "morantel", true = "Morantel", false = drug),
                drug = if_else(drug == "pyrantel_cytrate", true = "Pyrantel citrate", false = drug)) %>% 
  ggplot(., mapping = aes(x = strain, y = Estimate, ymin = Lower, ymax = Upper,
                          color = strain)) + 
  theme_bw(base_size = 11) + 
  geom_pointrange(position = position_quasirandom(groupOnX = FALSE),
                  size = 0.5) +
  scale_color_manual(values = strain_colors, name = "Strain") + 
  facet_wrap(drug ~., scales = "free", nrow = 1) + 
  theme(strip.text = element_text(angle = 0, size = 9),
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top") + 
  labs(y = "EC10 Estimate (µM)",
       x = "Strain")

complete.EC10.plot.NACHA

ggsave(complete.EC10.plot.NACHA + theme(plot.background = element_rect(fill = "white",colour = NA)),
       filename = "output/complete.EC10.plot_30pixels_NACHAs.png", width = 6, height = 2.23)

View(EC10.filtered)



########################################################
### Figure 8. SLOPE AND EC10 PLOT COMPARED TO N2 #######
########################################################
#Relative EC10 figure compared to N2
EC10.relative.potency$big_class <- factor(EC10.relative.potency$big_class, levels = c("Benzimidazoles", "Macrocyclic lactones","Nicotinic acetylcholine receptor agonists", 
                                                                        "Nicotinic acetylcholine receptor antagonists", "Crystal protein", "Cyclooctadepsipeptides",
                                                                        "Schistosomicides", "Other"))



EC10.relative.potency.figure <- EC10.relative.potency %>%
  dplyr::filter(!(drug=="Albendazole" & focal.strain=="MY16")) %>% #in supp fig. - on different scale
  dplyr::filter(!(drug=="Fenbendazole" & focal.strain=="MY16")) %>% #in supp fig. - on different scale
  dplyr::filter(!(drug=="Mebendazole" & focal.strain=="MY16")) %>% #in supp fig. - on different scale
  dplyr::filter(!(drug=="Benomyl" & focal.strain=="MY16")) %>% #in supp fig. - on different scale
  dplyr::filter(!(drug=="Levamisole" & focal.strain=="DL238")) %>% #in supp fig. - on different scale
  dplyr::filter(!(drug=="Levamisole" & focal.strain=="CB4856")) %>% #in supp fig. - on different scale
  dplyr::mutate(drug = gsub(drug, pattern = "Pyrantel citrate", replacement = "Pyrantel")) %>% 
  ggplot(., mapping = aes(y = drug, x = N2.normed.diff, 
                          xmin = N2.normed.diff-Std..Error, 
                          xmax = N2.normed.diff+Std..Error,
                          color = focal.strain,
                          alpha = sig)) + 
  scale_y_discrete(limits = rev) + #places drug names in alphabetical order by class
   theme_bw(base_size = 11) +
  geom_vline(xintercept = 0, linetype = 5, colour = "orange") + 
  geom_pointrange(position = position_dodge(width = 0.2),
                  size = 0.5) + 
  scale_color_manual(values = strain_colors[c(1,3:6)], name = "Strain") +
  scale_alpha_manual(values = c(0.15,1), guide = "none") +
  facet_grid(big_class~., scales = "free", space = "free") + 
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_text(color = "black"),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "none") + 
  labs(x = "Resistance Compared to Reference (N2)") 
  
EC10.relative.potency.figure

#Organize drug classes in plotting order 
relative.slope$big_class <- factor(relative.slope$big_class, levels = c("Benzimidazoles", "Macrocyclic lactones","Nicotinic acetylcholine receptor agonists", 
                                                                        "Nicotinic acetylcholine receptor antagonists", "Crystal protein", "Cyclooctadepsipeptides",
                                                                        "Schistosomicides", "Other"))

View(relative.slope)
#Relative slope figure compared to N2
relative.slope.figure <- relative.slope %>%
  dplyr::filter(!(drug=="Milbemycin" & focal.strain=="DL238")) %>% #in supp fig. - on different scale
  dplyr::filter(!(drug=="Albendazole" & focal.strain=="MY16")) %>% #in supp fig. - on different scale
  dplyr::mutate(drug = gsub(drug, pattern = "Pyrantel citrate", replacement = "Pyrantel")) %>% 
   ggplot(., mapping = aes(y = drug, x = N2.normed.diff, 
                          xmin = N2.normed.diff-Std..Error, 
                          xmax = N2.normed.diff+Std..Error,
                          color = focal.strain,
                          alpha = sig)) + 
  scale_y_discrete(limits = rev) +
  theme_bw(base_size = 11) +
  geom_vline(xintercept = 0, linetype = 5, colour = "orange") + 
  geom_pointrange(position = position_dodge(width = 0.2),
                  size = 0.5) + 
  scale_color_manual(values = strain_colors[c(1,3:6)], name = "Strain") +
  scale_alpha_manual(values = c(0.15,1), guide = "none") +
  facet_grid(big_class~., scales = "free", space = "free") + #places by big class in alphabetical order
  theme(axis.text.x = element_text(color = "black"),
        #axis.text.y = element_text(color = "black", size = 11),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), #added
        axis.ticks.y = element_blank(), #added
        panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "none") + 
  labs(x = "Relative Slope Compared to Reference (N2)")
relative.slope.figure

# ggsave(relative.slope.figure, filename = "output/relative.slope.figure.alldrugs.png", width = 16, height = 8)
# ggsave(EC10.relative.potency.figure, filename = "output/EC10.relative.potency.ALL.png", width = 12, height = 8)

figure8 <- cowplot::plot_grid(EC10.relative.potency.figure +  
                                   theme(strip.background = element_blank(),strip.text.y = element_blank(),
                                         plot.margin = unit(c(0.5,0.5,0.25,0.5), "cm"),
                                         legend.position = "none"),
                              relative.slope.figure  + 
                                   theme(plot.margin = unit(c(0.5,0.5,0.25,0.5), "cm"),
                                         legend.position = "none"), 
                                 rel_widths = c(1.2,1.4), labels = "AUTO")

complete.figure8 <- cowplot::plot_grid(fig.3.legend_BZ, figure8, rel_heights = c(2,20), ncol = 1)
complete.figure8
ggsave(complete.figure8 + theme(plot.background = element_rect(fill = "white",colour = NA)), filename = "output/figure8_ec10_and_slope_compared_n2.png", width = 12, height = 6)
#too large for PLOS - rescaled in biorender 


########################################################
### Supplemental Figure 3.EC10 PLOT COMPARED TO N2 #####
########################################################

#Relative EC10 figure compared to N2
supp.fig.3.EC10.relative.potency.figure <- EC10.relative.potency %>%
  # dplyr::filter(!(drug=="Albendazole" & focal.strain=="MY16")) %>% #in supp fig. - on different scale
  # dplyr::filter(!(drug=="Fenbendazole" & focal.strain=="MY16")) %>% #in supp fig. - on different scale
  # dplyr::filter(!(drug=="Mebendazole" & focal.strain=="MY16")) %>% #in supp fig. - on different scale
  # dplyr::filter(!(drug=="Benomyl" & focal.strain=="MY16")) %>% #in supp fig. - on different scale
  # dplyr::filter(!(drug=="Levamisole" & focal.strain=="DL238")) %>% #in supp fig. - on different scale
  # dplyr::filter(!(drug=="Levamisole" & focal.strain=="CB4856")) %>% #in supp fig. - on different scale
  dplyr::mutate(drug = gsub(drug, pattern = "Pyrantel citrate", replacement = "Pyrantel")) %>% 
    ggplot(., mapping = aes(y = drug, x = N2.normed.diff, 
                          xmin = N2.normed.diff-Std..Error, 
                          xmax = N2.normed.diff+Std..Error,
                          color = focal.strain,
                          alpha = sig)) + 
  scale_y_discrete(limits = rev) +
  theme_bw(base_size = 11) +
  geom_vline(xintercept = 0, linetype = 5, colour = "orange") + 
  geom_pointrange(position = position_dodge(width = 0.2),
                  size = 0.5) + 
  scale_color_manual(values = strain_colors[c(1,3:6)], name = "Strain") +
  scale_alpha_manual(values = c(0.15,1), guide = "none") +
  facet_grid(big_class~., scales = "free", space = "free") + 
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_text(color = "black"),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "none") + 
  labs(x = "Resistance Compared to Reference (N2)")

supp.fig.3.EC10.relative.potency.figure

supp.fig3.4.legend <- cowplot::get_legend(complete.EC10.plot)

supp.fig3.EC10 <- cowplot::plot_grid(supp.fig.3.EC10.relative.potency.figure +  
                                theme(strip.background = element_blank(),strip.text.y = element_blank(),
                                      plot.margin = unit(c(0.5,0.5,0.25,0.5), "cm"),
                                      legend.position = "none"),
                              rel_widths = c(1.2))

complete.supp.fig.3 <- cowplot::plot_grid(fig.3.legend_BZ, supp.fig.3.EC10.relative.potency.figure, rel_heights = c(2,20), ncol = 1)

complete.supp.fig.3

ggsave(complete.supp.fig.3 + theme(plot.background = element_rect(fill = "white",colour = NA)), filename = "output/complete.supp.fig.3.png", width = 7.5, height = 5)


########################################################
### Supplemental Figure 4.SLOPE PLOT COMPARED TO N2 ####
########################################################
#Relative slope figure compared to N2
supp.fig.4.relative.slope.figure <- relative.slope %>%
  # dplyr::filter(!(drug=="Milbemycin" & focal.strain=="DL238")) %>% #in supp fig. - on different scale
  # dplyr::filter(!(drug=="Albendazole" & focal.strain=="MY16")) %>% #in supp fig. - on different scale
  dplyr::mutate(drug = gsub(drug, pattern = "Pyrantel citrate", replacement = "Pyrantel")) %>% 
  ggplot(., mapping = aes(y = drug, x = N2.normed.diff, 
                          xmin = N2.normed.diff-Std..Error, 
                          xmax = N2.normed.diff+Std..Error,
                          color = focal.strain,
                          alpha = sig)) + 
  scale_y_discrete(limits = rev) +
    theme_bw(base_size = 11) +
  geom_vline(xintercept = 0, linetype = 5, colour = "orange") + 
  geom_pointrange(position = position_dodge(width = 0.2),
                  size = 0.5) + 
  scale_color_manual(values = strain_colors[c(1,3:6)], name = "Strain") +
  scale_alpha_manual(values = c(0.15,1), guide = "none") +
  facet_grid(big_class~., scales = "free", space = "free") + 
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.y = element_blank(),
        # axis.text.y = element_blank(), #added
        # axis.ticks.y = element_blank(), #added
        panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "none") + 
  labs(x = "Relative Slope Compared to Reference (N2)")
supp.fig.4.relative.slope.figure


supp.fig3.4.legend <- cowplot::get_legend(complete.EC10.plot)

supp.fig4.slope <- cowplot::plot_grid(supp.fig.4.relative.slope.figure+  
                                       theme(strip.background = element_blank(),strip.text.y = element_blank(),
                                             plot.margin = unit(c(0.5,0.5,0.25,0.5), "cm"),
                                             legend.position = "none"),
                                     rel_widths = c(1.2))

complete.supp.fig.4 <- cowplot::plot_grid(fig.3.legend_BZ, supp.fig.4.relative.slope.figure, rel_heights = c(2,20), ncol = 1)

complete.supp.fig.4

ggsave(complete.supp.fig.4 + theme(plot.background = element_rect(fill = "white",colour = NA)), filename = "output/complete.supp.fig.4.png", width = 7.5, height = 5)


slope.aov.nested <- slopes.filtered %>%
  #dplyr::filter(!drug %in% c("Deltamethrin","Malathion")) %>%
  dplyr::select(drug, Estimate, strain, big_class) %>%
  dplyr::group_by(big_class) %>%
  tidyr::nest()

View(slope.aov.nested)

slope.aov <- function(estimate, class){
  if(length(unique(estimate$drug)) > 1){
    big.class.aov <- aov(estimate, formula = Estimate ~ strain + drug)
    aov.df <- data.frame(summary(big.class.aov)[[1]][[5]][1],
                         summary(big.class.aov)[[1]][[5]][2])
    colnames(aov.df) <- c("strain.p","tox.p")
    aov.df$class <- class
    return(list(aov.df, TukeyHSD(big.class.aov)))
  } else {
    return("Only one Anthelmintic!")
  }
  
}

slope.aov.list <- purrr::map2(slope.aov.nested$data,
                              slope.aov.nested$big_class,
                              slope.aov)

slope.aov.list.tr <- slope.aov.list %>%
  purrr::keep(., is.list) %>%
  purrr::transpose()

View(slope.aov.list.tr)
# Overall anova results for slope\
Reduce(rbind,slope.aov.list.tr[[1]])

# Macrocyclic lactones -  Tukey's HSD Results
slope.aov.list.tr[[2]][[1]][[2]] %>%
  data.frame() %>%
  dplyr::filter(p.adj < 0.05)

# Benzimidazoles - Tukey's HSD Results
slope.aov.list.tr[[2]][[2]][[2]] %>%
  data.frame() %>%
  dplyr::filter(p.adj < 0.05)

# Nicotinic acetylcholine receptors agonists -  Tukey's HSD Results
slope.aov.list.tr[[2]][[3]][[2]] %>%
  data.frame() %>%
  dplyr::filter(p.adj < 0.05)

# Nicotinic acetylcholine receptors antagonists -  Tukey's HSD Results
slope.aov.list.tr[[2]][[4]][[2]] %>%
  data.frame() %>%
  dplyr::filter(p.adj < 0.05)


## Slope Relative Potency Tests ##
dose.response.parameter.summaries[[7]] %>%
  tidyr::separate(comps, c("strain1","strain2"), sep = "/") %>%
  dplyr::mutate(strain1 = gsub(strain1, pattern = "strain", replacement = "")) %>%
  dplyr::mutate(strain2 = gsub(strain2, pattern = "strain", replacement = "")) %>% 
  dplyr::filter(drug %in% unique(EC10.filtered$drug),
                strain2 == "N2" | strain1 == "N2")

sig.EC.comps <- dose.response.parameter.summaries[[7]] %>%
  tidyr::separate(comps, c("strain1","strain2"), sep = "/") %>%
  dplyr::mutate(strain1 = gsub(strain1, pattern = "strain", replacement = "")) %>%
  dplyr::mutate(strain2 = gsub(strain2, pattern = "strain", replacement = "")) %>% 
  dplyr::filter(drug %in% unique(EC10.filtered$drug),
                strain2 == "N2" | strain1 == "N2") %>%
  dplyr::mutate(N2.normed.estimate = if_else(strain1 == "N2", true = 1/Estimate, false = Estimate)) %>%
  dplyr::mutate(N2.normed.diff = N2.normed.estimate-1,
                focal.strain = if_else(strain1 == "N2", true = strain2, false = strain1)) %>%
  dplyr::left_join(.,anth.classes %>%
                     dplyr::rename(drug = Anthelmintic)) %>%
  dplyr::filter(p.value < BF)

nrow(sig.EC.comps)


sig.slope.comps <- relative.slope %>%
  dplyr::filter(sig == "SIG")

nrow(sig.slope.comps)
sig.slope.comps %>%
  dplyr::group_by(big_class) %>%
  dplyr::count()

sig.slope.comps %>%
  dplyr::group_by(focal.strain) %>%
  dplyr::count()


slope.comp.sig <- sig.slope.comps %>%
  dplyr::mutate(sens = if_else(N2.normed.diff > 0, "steep", "shallow")) %>%
  dplyr::group_by(focal.strain, sens) %>%
  dplyr::count() %>%
  tidyr::pivot_wider(names_from = sens, values_from = n)


resistance <- slope.comp.sig %>%
  dplyr::select(focal.strain, shallow) %>%
  dplyr::mutate(not.shallow = 21-shallow)

strain.resistance.fisher.reps <- list()
for(i in 1:1000){
  rep <- fisher.test(resistance[,2:3], simulate.p.value = T)
  strain.resistance.fisher.reps[[i]] <- rep$p.value
}
c(mean(unlist(strain.resistance.fisher.reps)),sd(unlist(strain.resistance.fisher.reps)))


sensitive <- slope.comp.sig %>%
  dplyr::select(focal.strain, steep) %>%
  dplyr::mutate(not.steep = 21-steep)
strain.sensitive.fisher.reps <- list()
for(i in 1:1000){
  rep <- fisher.test(sensitive[,2:3], simulate.p.value = T)
  strain.sensitive.fisher.reps[[i]] <- rep$p.value
}
c(mean(unlist(strain.sensitive.fisher.reps)),sd(unlist(strain.sensitive.fisher.reps)))


# ##############################
# ###     HERITABILITY       ###
# ##############################
all.heritabilities <- purrr::map2(tx.nested.dose.responses.MDHD30$data,
                                  tx.nested.dose.responses.MDHD30$drug,
                                  heritability.calculation)


View(tx.nested.dose.responses.MDHD30$data)

View(all.heritabilities)
View(NARAntagonist.HH$data)

summarized.heritabilities <- purrr::map(all.heritabilities, gather.herit.ranges)

summarized.heritabilities.df <- Reduce(rbind,summarized.heritabilities) %>%
  dplyr::filter(h2 != "NA",
                drug %in% unique(EC10.filtered$drug),
                !drug %in% c("Deltamethrin","Malathion")) %>% 
  dplyr::mutate(h2 = as.numeric(h2),
                H2 = as.numeric(H2),
                h2.upper= as.numeric(h2.upper),
                h2.lower= as.numeric(h2.lower),
                H2.lower= as.numeric(H2.lower),
                H2.upper= as.numeric(H2.upper)) %>%
  dplyr::left_join(.,anth.classes %>%
                     dplyr::rename(drug = Anthelmintic)) %>%
  dplyr::mutate(big_class = gsub(big_class, pattern = "_", replacement = " "))


View(summarized.heritabilities.df)

# h2 for levamisole, pyrantel, morantel
NARagonists.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Nicotinic acetylcholine receptors agonists") %>%
  dplyr::filter(!drug %in% c("Derquantel")) %>% 
  plot_HH(.)

write.csv(NARagonists.HH$data, "NARagonists.HH.csv")

View(NARagonists.HH$data)
#h2 for just derquantel
NARagonists.HH_derquantel <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Nicotinic acetylcholine receptors agonists") %>%
  dplyr::filter(!drug %in% c("Morantel", "Levamisole", "Pyrantel citrate")) %>% 
  plot_HH(.)

write.csv(NARagonists.HH_derquantel$data, "NARagonists.HH_derquantel.csv")

BZ.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Benzimidazoles" ) %>%
  plot_HH(.)

write.csv(BZ.HH$data, "BZ.HH.csv")

BZ.HH.NoMY16 <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Benzimidazoles" ) %>%
  dplyr::filter(!strain %in% c("MY16")) %>% 
  plot_HH(.)


Cyclooctadepsipeptides.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Cyclooctadepsipeptides") %>%
  plot_HH(.)

write.csv(Cyclooctadepsipeptides.HH$data, "Cyclooctadepsipeptides.HH.csv")


Diethylcarbamazine.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Diethylcarbamazine") %>%
  plot_HH(.)

write.csv(Diethylcarbamazine.HH$data, "Diethylcarbamazine.HH.csv")


Crystal_protein.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Crystal protein") %>%
  plot_HH(.)

write.csv(Crystal_protein.HH$data, "Crystal_protein.HH.csv")


ML.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Macrocyclic lactones") %>%
  plot_HH(.)

write.csv(ML.HH$data, "ML.HH.csv")


NARAntagonists.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Nicotinic acetylcholine receptors antagonists") %>%
  plot_HH(.)

write.csv(NARAntagonists.HH $data, "NARAntagonists.HH .csv")


Pyrazinoisoquinolines.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Pyrazinoisoquinolines") %>%
  plot_HH(.)
write.csv(Pyrazinoisoquinolines.HH$data, "Pyrazinoisoquinolines.HH.csv")


Quinolines.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Quinolines") %>%
  plot_HH(.)
write.csv(Quinolines.HH$data, "Quinolines.HH.csv")


Salicylanilides.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Salicylanilides") %>%
  plot_HH(.)
write.csv(Salicylanilides.HH$data, "Salicylanilides.HH .csv")


Schistosomicides.HH <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Schistosomicides") %>%
  plot_HH(.)
write.csv(Schistosomicides.HH $data, "Schistosomicides.HH.csv")


fig.4.legend <- cowplot::get_legend(BZ.HH)
A <- cowplot::plot_grid(Crystal_protein.HH + theme(legend.position = "none"),
                        Cyclooctadepsipeptides.HH + theme(legend.position = "none"),
                        #Quinolines.HH + theme(legend.position = "none"),
                        labels = c("A", 'B'), rel_widths = c(0.35,0.35))


B <- cowplot::plot_grid(Diethylcarbamazine.HH + theme(legend.position = "none"),
                        Salicylanilides.HH + theme(legend.position = "none"),
                        labels = c("C", 'D'), rel_widths = c(0.45, 0.45))


C <- cowplot::plot_grid(Schistosomicides.HH + theme(legend.position = "none"),
                        # Pyrazinoisoquinolines.HH + theme(legend.position = "none"),
                        labels = c("E"), ncol = 2, rel_widths = c(0.35))

D <- cowplot::plot_grid(NARagonists.HH + theme(legend.position = "none"),
                        labels = c("C"), rel_widths = c(0.35, 0.35, 0.35))

F <- cowplot::plot_grid(NARagonists.HH_derquantel + theme(legend.position = "none"),
                        labels = c("C"), rel_widths = c(0.35))


E<- cowplot::plot_grid(NARAntagonists.HH + theme(legend.position = "none"),
                       labels = c("G"), rel_widths = c(0.45, 0.45))

F <- cowplot::plot_grid(BZ.HH + theme(legend.position = "none") +
                          facet_wrap(facets = drug ~., nrow = 2),
                        labels = c("C"), ncol = 1, rel_widths  = c(0.35))


G <- cowplot::plot_grid(ML.HH + theme(legend.position = "none") +
                          facet_wrap(facets = drug ~., nrow = 2),
                        labels = c("C"), ncol = 1)


# Benzimidazoles Heritability plot 
fig.BZ.H2.legend <- cowplot::get_legend(BZ.HH)
F <- cowplot::plot_grid(BZ.HH + theme(legend.position = "none") +
                          facet_wrap(facets = drug ~., nrow = 2),
                        labels = c("C"), ncol = 1, rel_widths  = c(0.35))

heritability.summary.plot <- cowplot::plot_grid(F, fig.BZ.H2.legend, ncol = 1, rel_heights = c(1, 0.35))
heritability.summary.plot
ggsave(heritability.summary.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "output/fig.H2.only_BZs_NOMY16.png", width = 7, height = 7)

#past sizing for BZ plot 
#width = 8.75, height = 5

# Macrocyclic lactones Heritability plot 
fig.ML.H2.legend <- cowplot::get_legend(ML.HH)
F <- cowplot::plot_grid(ML.HH + theme(legend.position = "none") +
                          facet_wrap(facets = drug ~., nrow = 2),
                        labels = c("C"), ncol = 1, rel_widths  = c(0.35))

heritability.summary.plot <- cowplot::plot_grid(F, fig.ML.H2.legend, ncol = 1, rel_heights = c(1, 0.35))
heritability.summary.plot
ggsave(heritability.summary.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "output/fig.H2.only_MLs_3.png", width = 7, height = 7)


# NAChAs Heritability plot 
fig.NAChAs.H2.legend <- cowplot::get_legend(NARagonists.HH)
F <- cowplot::plot_grid(NARagonists.HH+ theme(legend.position = "none") +
                          facet_wrap(facets = drug ~., nrow = 1),
                        labels = c("C"), ncol = 1, rel_widths  = c(0.35))

heritability.summary.plot <- cowplot::plot_grid(F, fig.NAChAs.H2.legend, ncol = 1, rel_heights = c(1, 0.35))
heritability.summary.plot
ggsave(heritability.summary.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "output/fig.H2.only_NAChAs_noderquantel.png", width = 7, height = 3.5)


# Derquantel Heritability plot 
fig.derquantel.H2.legend <- cowplot::get_legend(NARagonists.HH_derquantel)
F <- cowplot::plot_grid(NARagonists.HH_derquantel + theme(legend.position = "none"), rel_widths  = c(0.35))

heritability.summary.plot <- cowplot::plot_grid(F, fig.NAChAs.H2.legend, ncol = 1, rel_heights = c(1, 0.35))
heritability.summary.plot
ggsave(heritability.summary.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "output/fig.H2.only_derquantel.png", width = 3, height = 3.5)

# Closantel Heritability plot 
fig.4.legend <- cowplot::get_legend(Salicylanilides.HH)
G <- cowplot::plot_grid(Salicylanilides.HH  + theme(legend.position = "none"), rel_widths = c(0.35))

heritability.summary.plot <- cowplot::plot_grid(G, fig.4.legend, ncol = 1, rel_heights = c(1, 0.35))
heritability.summary.plot
ggsave(heritability.summary.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "output/fig.H2.only_closantel2.png", width = 3, height = 3.5)


# Emodepside Heritability plot 
fig.emoh2.legend <- cowplot::get_legend(Cyclooctadepsipeptides.HH)
G <- cowplot::plot_grid(Cyclooctadepsipeptides.HH  + theme(legend.position = "none"),
                        rel_widths = c(0.35))

heritability.summary.plot <- cowplot::plot_grid(G, fig.emoh2.legend, ncol = 1, rel_heights = c(1, 0.35))
ggsave(heritability.summary.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "output/fig.H2.only_emodepside.png", width = 3, height = 3.5)

# Diethylcarbamazine Heritability plot 
fig.diethylh2.legend <- cowplot::get_legend(Diethylcarbamazine.HH)
G <- cowplot::plot_grid(Diethylcarbamazine.HH  + theme(legend.position = "none"),
                        rel_widths = c(0.35))

heritability.summary.plot <- cowplot::plot_grid(G, fig.diethylh2.legend, ncol = 1, rel_heights = c(1, 0.35))
ggsave(heritability.summary.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "output/fig.H2.only_diethylcarbamazine_2.png", width = 3, height = 3.5)

# Cry5b Heritability plot 
fig.cry5bh2.legend <- cowplot::get_legend(Crystal_protein.HH)
G <- cowplot::plot_grid(Crystal_protein.HH + theme(legend.position = "none"),
                        rel_widths = c(0.35))

heritability.summary.plot <- cowplot::plot_grid(G, fig.cry5bh2.legend, ncol = 1, rel_heights = c(1, 0.35))
ggsave(heritability.summary.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "output/fig.H2.only_cry5b_2.png", width = 3, height = 3.5)

# # Oxamniquine Heritability plot - No values for this drug. No H2. 
# fig.oxamniquine.h2.legend <- cowplot::get_legend(Quinolines.HH)
# G <- cowplot::plot_grid(Quinolines.HH + theme(legend.position = "none"),
#                         labels = c("B"), rel_widths = c(0.35))
# 
# heritability.summary.plot <- cowplot::plot_grid(G, fig.oxamniquine.h2.legend, ncol = 1, rel_heights = c(1, 0.35))
# ggsave(heritability.summary.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), 
#        filename = "output/fig.H2.only_oxamniquine.png", width = 2.5, height = 2.5)


# Niridazole Heritability plot 
fig.nirdizaole.h2.legend <- cowplot::get_legend(Schistosomicides.HH)
G <- cowplot::plot_grid(Schistosomicides.HH + theme(legend.position = "none"),
                        rel_widths = c(0.35))

heritability.summary.plot <- cowplot::plot_grid(G, fig.nirdizaole.h2.legend, ncol = 1, rel_heights = c(1, 0.35))
ggsave(heritability.summary.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "output/fig.H2.only_niridzole_2.png", width = 3, height = 3.5)


# Monepantel drugs Heritability plot 
#seperate monepantel drugs 
NARAntagonists.HH.Mon916 <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Nicotinic acetylcholine receptors antagonists") %>%
  dplyr::filter(!drug == "Monepantel LY33414916") %>% 
  plot_HH(.)

NARAntagonists.HH.Mon298 <- summarized.heritabilities.df %>%
  dplyr::filter(big_class == "Nicotinic acetylcholine receptors antagonists") %>%
  dplyr::filter(!drug == "Monepantel LY3348298") %>% 
  plot_HH(.)

fig.monepantel.h2.legend <- cowplot::get_legend(NARAntagonists.HH)

G <- cowplot::plot_grid(NARAntagonists.HH.Mon916 + theme(legend.position = "none"),
                        rel_widths = c(0.35))

H <- cowplot::plot_grid(NARAntagonists.HH.Mon298 + theme(legend.position = "none"),
                        rel_widths = c(0.35))

heritability.summary.plot <- cowplot::plot_grid(G, fig.monepantel.h2.legend, ncol = 1, rel_heights = c(1, 0.35))

ggsave(heritability.summary.plot + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "output/fig.NARAntagonists.HH.Mon916.png", width = 3, height = 3.5)


###############################
### Supplementary Figure 1  ###
###############################

nested.sum.herits <- summarized.heritabilities.df %>%
  dplyr::group_by(drug) %>%
  tidyr::nest()
ranked.herits <- purrr::map2(nested.sum.herits$data,
                             nested.sum.herits$drug, 
                             function(x,y){
                               x %>%
                                 dplyr::mutate(rank = seq(1:nrow(x))) %>%
                                 dplyr::mutate(drug = y)}) %>%
  Reduce(rbind,.)
max.H2 <- ranked.herits %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(H2 = max(H2)) %>%
  dplyr::left_join(., ranked.herits) %>%
  dplyr::select(drug, concentration_um, H2, rank) %>%
  dplyr::arrange(rank) %>%
  data.frame()
max.H2
max.h2 <- ranked.herits %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(h2 = max(h2)) %>%
  dplyr::left_join(., ranked.herits) %>%
  dplyr::select(drug, concentration_um, h2, rank) %>%
  dplyr::arrange(rank) %>%
  data.frame()
max.h2

max.h2 %>%
  dplyr::group_by(rank) %>%
  dplyr::summarise(n = n())

avg.herits <- ranked.herits %>%
  dplyr::filter(concentration_um != 0) %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(mean(H2), 
                   sd(H2),
                   mean(h2),
                   sd(h2))
avg.herits %>%
  dplyr::filter(`mean(H2)` %in% range(avg.herits$`mean(H2)`))

avg.herits %>%
  dplyr::filter(`mean(h2)` %in% range(avg.herits$`mean(h2)`))



EC.10.diff.H2 <- EC10.filtered %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(mean.EC10 = mean(Estimate)) %>%
  dplyr::full_join(ranked.herits,.) %>%
  dplyr::mutate(herit.EC10.diff = abs(concentration_um - mean.EC10))

herit.at.EC10.closest.dose <- EC.10.diff.H2 %>%
  dplyr::select(drug, concentration_um, rank, mean.EC10, herit.EC10.diff) %>%
  dplyr::filter(!is.na(mean.EC10)) %>%
  dplyr::group_by(drug) %>% 
  dplyr::summarise(herit.EC10.diff = min(herit.EC10.diff)) %>%
  dplyr::left_join(.,EC.10.diff.H2) %>%
  dplyr::select(drug, concentration_um, H2, h2) %>%
  dplyr::left_join(., max.H2 %>%
                     dplyr::rename(max.H2 = H2,
                                   max.H2.conc = concentration_um) %>%
                     dplyr::select(-rank)) %>%
  dplyr::left_join(., max.h2 %>%
                     dplyr::rename(max.h2 = h2,
                                   max.h2.conc = concentration_um))

h2.comp.table <- herit.at.EC10.closest.dose %>%
  dplyr::filter(drug %in% unique(EC10.filtered$drug)) %>% 
  dplyr::select(-max.h2.conc) %>%
  dplyr::mutate(`EC10 Dose - Top Heritable Dose` =  concentration_um - max.H2.conc) %>%
  dplyr::select(drug, concentration_um,max.H2.conc, `EC10 Dose - Top Heritable Dose`, H2, max.H2, h2, max.h2, rank) %>%
  dplyr::rename(Anthelmintic = drug,
                `Closest Dosage to EC10` = concentration_um,
                `Top Heritable Dose` = max.H2.conc,
                `H2 EC10` = H2,
                `Top H2` = max.H2,
                `h2 EC10` = h2,
                `Top h2` = max.h2,
                `Top Heritable Dose Rank` = rank)
colnames(h2.comp.table)

strain.nested.EC10 <- EC10 %>%
  dplyr::filter(drug %in% unique(EC10.filtered$drug)) %>% 
  dplyr::select(drug, strain, Estimate) %>%
  dplyr::group_by(strain) %>%
  tidyr::nest()



# 
purrr::map2(strain.nested.EC10$data,
            strain.nested.EC10$strain,
            h2.comps.to.individual.strains)


View(h2.comp.table)
write.csv(h2.comp.table, "output/h2.comp.table.csv", row.names = T)

h2.comp.table.2 <- h2.comp.table %>% 
  dplyr::filter(!Anthelmintic %in% c("Cry5B", "Levamisole", "Morantel"))
#excluding drugs that have a value of 0 for "Closest Dosage to EC10"
#Drugs with value 0 = Cry5B, Levamisole, Morantel

View(h2.comp.table.2)

h2.comp.r2 <- summary(lm(data = h2.comp.table.2, formula = log10(`Top Heritable Dose`) ~ log10(`Closest Dosage to EC10`)))

l <- list(r2 = format(h2.comp.r2$r.squared, digits = 3),
          pval = format(as.numeric(pf(h2.comp.r2$fstatistic[1],
                                      h2.comp.r2$fstatistic[2],
                                      h2.comp.r2$fstatistic[3],
                                      lower.tail=FALSE)), digits=3, scientific = TRUE))


eq <- substitute(italic(r)^2 == r2*","~~italic(p) == pval,l)

eqstr <- as.character(as.expression(eq))


h2.comp.figure <- ggplot() +
  theme_bw(base_size = 11) + 
  geom_smooth(h2.comp.table, mapping = aes(x = log10(`Closest Dosage to EC10`), 
                                           y = log10(`Top Heritable Dose`)),
              method = "lm", colour = "darkred") + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_jitter(h2.comp.table, mapping = aes(x = log10(`Closest Dosage to EC10`), 
                                           y = log10(`Top Heritable Dose`)),
              width = 0.06, height = 0.06) + # some values are overplotted on log scale
  
  annotate(geom = "text", x = -1, y = 3, label = eqstr, parse = TRUE, colour = "darkred", size = 3) + 
  theme(panel.grid = element_blank()) + 
  labs(x = expression(log[10](`Closest Dosage to EC10 `(µM))),
       y = expression(log[10](`Most Heritable Dosage `(µM))))
h2.comp.figure

ggsave(h2.comp.figure + theme(plot.background = element_rect(fill = "white",colour = NA)), 
       filename = "output/fig.ec10_h2compfigure.png", width = 3.5, height = 3.5)


pivoted.top.herits <- max.H2 %>%
  tidyr::pivot_longer(cols = H2, names_to = "measurement") %>%
  dplyr::full_join(., max.h2 %>%
                     tidyr::pivot_longer(cols = h2, names_to = "measurement"))

pivoted.herits.at.EC10 <- herit.at.EC10.closest.dose %>%
  dplyr::rename(H2.EC10 = H2, 
                h2.EC10 = h2) %>%
  tidyr::pivot_longer(cols = c(H2.EC10,h2.EC10), names_to = "measurement")

top.v.EC10.herit <- pivoted.top.herits %>%
  dplyr::full_join(., pivoted.herits.at.EC10)

top.v.EC10.herit %>%
  dplyr::arrange(drug) %>%
  dplyr::filter(measurement %in% c("h2", "h2.EC10")) %>%
  dplyr::select(-rank) %>%
  dplyr::group_by(drug) 


#################################
###    All DR Curve Figures   ###
#################################
#supp.DR.dodge.widths <- c(rep(0.15, 12),0.04, 0.08,rep(0.04, 2),0.09, 0.04, 0.02,rep(0.10, 2),0.04,0.02,rep(0.04,2))
supp.DR.dodge.widths <- 0.01
DR.theme <- ggplot2::theme(axis.text = element_text(size = 7),
                           axis.title = element_text(size = 7), 
                           title = element_text(size = 7))
for(i in 1:length(tx.nested.dose.responses.MDHD30$data)){
  supp.DR <- dose.response.plots.only(tx.nested.dose.responses.MDHD30$data[[i]],
                                      tx.nested.dose.responses.MDHD30$drug[[i]],
                                      supp.DR.dodge.widths[[i]])
  #SE.drc.point.size = 0.05) 
  assign(paste0(gsub(tx.nested.dose.responses.MDHD30$drug[[i]],
                     pattern = " ",
                     replacement = "_"),"_SUPP_PANEL"), 
         value = supp.DR + 
           DR.theme + 
           ggplot2::theme(legend.position = "none"))
  ggsave(supp.DR, filename = paste0("output/fig.", i + 1,".png"), width = 6, height = 6)
}

fig.2 <- cowplot::plot_grid(plotlist = mget(list(ls(pattern = "SUPP_PANEL"))[[1]]), nrow = 5, cols = 5)
ggsave(fig.2,filename = "output/fig.all.drcurves.png", height = 10, width = 10)

