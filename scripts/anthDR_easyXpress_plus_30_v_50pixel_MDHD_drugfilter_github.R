## EasyXpress steps taken from: http://127.0.0.1:25554/library/easyXpress/doc/easyXpress.html

#install.packages("devtools")
#devtools::install_github("AndersenLab/easyXpress")
require(tidyverse)
require(easyXpress)
install.packages("drc")
install.packages("knitr")
install.packages("RCurl")
install.packages("RColorBrewer")
install.packages("cowplot")
install.packages("ggbeeswarm")
install.packages("ggrepel")
install.packages("ggnewscale")
install.packages("scales")

library(tidyverse)
library(dplyr)
library(easyXpress)
library(drc)
library(knitr)
library(RCurl)
library(RColorBrewer)
library(cowplot)
library(ggbeeswarm)
library(ggrepel)
library(ggnewscale)
library(scales)

#############################
#                           #
#          Assay A          #
#                           #
#############################

#Reading in the data
# Define experimental directory and file name
datafileA <- "20210205_assayA_Analysis-20220710.RData"

# Define processed image directory
# proc_img_dir_assayA (set wd)

rawA_50 <- easyXpress::readXpress(filedir = dirsA, rdafile = datafileA, design = TRUE) %>%
   dplyr::filter(Worm_Length > 50) #remove objects smaller than 50 pixels = 165 um 

rawA_30 <- easyXpress::readXpress(filedir = dirsA, rdafile = datafileA, design = TRUE) %>%
  dplyr::filter(Worm_Length > 30) #remove objects smaller than 30 pixels = 100 um 

View(rawA_50) #ensure Worm_Length > 50 filtered properly

#drugs in assay A: moxidectin, abamectin, milbemycin, emodepside

#perform model selection where MDHD is removed on a drug by drug basis | 50 pixels
model_selectedA_50_drugspecificMDHD <- easyXpress::modelSelection(rawA_50) %>% #perform model selection
  dplyr::filter(drug != 'abamectin') #remove abamectin from assay A
edge_flaggedA_50 <- easyXpress::edgeFlag(model_selectedA_50_drugspecificMDHD, radius=825, center_x=1024, center_y=1024)
raw_flaggedA_50 <- easyXpress::setFlags(edge_flaggedA_50, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayA_drugMDHD_50 <- easyXpress::process(raw_flaggedA_50, Metadata_Plate, Metadata_Well)

save(processed_assayA_drugMDHD_50, file = "/data/processed_assayA_drugMDHD_50.RData")


#perform model selection where MDHD is removed on a drug by drug basis | 30 pixels
model_selectedA_30_drugspecificMDHD <- easyXpress::modelSelection(rawA_30) %>% #perform model selection
  dplyr::filter(drug != 'abamectin') #remove abamectin from assay A 
edge_flaggedA_30 <- easyXpress::edgeFlag(model_selectedA_30_drugspecificMDHD, radius=825, center_x=1024, center_y=1024)
raw_flaggedA_30 <- easyXpress::setFlags(edge_flaggedA_30, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayA_drugMDHD_30 <- easyXpress::process(raw_flaggedA_30, Metadata_Plate, Metadata_Well)

save(processed_assayA_drugMDHD_30, file = "/data/processed_assayA_drugMDHD_30.RData")


#############################
#                           #
#          Assay B          #
#                           #
#############################

#Reading in the data
# Define experimental directory and file name
# dirsB - set wd 
datafileB <- "20210212_assayB_Analysis-20220710.RData"

#drugs in assay B and max doses:
# benomyl	100
# closantel	100
# cry5b	33.00
# dec	660
# derquantel	165.00
# morantel	330.00
# niridazole	330.00
# oxamniquine	145.00
# piperazine	330
# pzq	1100

rawB_30 <- easyXpress::readXpress(filedir = dirsB, rdafile = datafileB, design = TRUE) %>%
  dplyr::filter(Worm_Length > 30) #remove objects smaller than 30 pixels = 100 um 

rawB_50 <- easyXpress::readXpress(filedir = dirsB, rdafile = datafileB, design = TRUE) %>%
  dplyr::filter(Worm_Length > 50) #remove objects smaller than 50 pixels = 165 um 

#Filter by 30 
#perform model selection where MDHD is removed on an assay by assay selection 
model_selectedB_30_drugspecificMDHD_BENOMYL <- easyXpress::modelSelection(rawB_30)%>% 
  dplyr::filter(drug == 'benomyl' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_30_drugspecificMDHD_CRY5B <- easyXpress::modelSelection(rawB_30)%>% 
  dplyr::filter(drug == 'cry5b' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_30_drugspecificMDHD_DEC <- easyXpress::modelSelection(rawB_30)%>% 
  dplyr::filter(drug == 'dec' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_30_drugspecificMDHD_MORANTEL <- easyXpress::modelSelection(rawB_30)%>% 
  dplyr::filter(drug == 'morantel' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_30_drugspecificMDHD_NIRIDAZOLE <- easyXpress::modelSelection(rawB_30)%>% 
  dplyr::filter(drug == 'niridazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_30_drugspecificMDHD_OXAMNIQUINE <- easyXpress::modelSelection(rawB_30)%>% 
  dplyr::filter(drug == 'oxamniquine' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_30_drugspecificMDHD_PIPERAZINE <- easyXpress::modelSelection(rawB_30)%>% 
  dplyr::filter(drug == 'piperazine' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_30_drugspecificMDHD_PZQ <- easyXpress::modelSelection(rawB_30)%>% 
  dplyr::filter(drug == 'pzq' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_30_drugspecificMDHD_CLOSANTEL <- easyXpress::modelSelection(rawB_30)%>% 
  dplyr::filter(drug == 'closantel')

model_selectedB_30_drugspecificMDHD_DERQUANTEL <- easyXpress::modelSelection(rawB_30)%>% 
  dplyr::filter(drug == 'derquantel')

model_selectedB_30_drugspecificMDHD <- model_selectedB_30_drugspecificMDHD_BENOMYL %>% 
  full_join(model_selectedB_30_drugspecificMDHD_CRY5B, by = NULL) %>% 
  full_join(model_selectedB_30_drugspecificMDHD_DEC, by = NULL) %>% 
  full_join(model_selectedB_30_drugspecificMDHD_MORANTEL, by = NULL) %>% 
  full_join(model_selectedB_30_drugspecificMDHD_NIRIDAZOLE, by = NULL) %>% 
  full_join(model_selectedB_30_drugspecificMDHD_OXAMNIQUINE, by = NULL) %>% 
  full_join(model_selectedB_30_drugspecificMDHD_PIPERAZINE, by = NULL) %>% 
  full_join(model_selectedB_30_drugspecificMDHD_PZQ, by =NULL) %>% 
  full_join(model_selectedB_30_drugspecificMDHD_CLOSANTEL, by = NULL) %>% 
  full_join(model_selectedB_30_drugspecificMDHD_DERQUANTEL, by = NULL)

edge_flaggedB_drugspecific_30 <- easyXpress::edgeFlag(model_selectedB_30_drugspecificMDHD, radius=825, center_x=1024, center_y=1024)
raw_flaggedB_drugspecific_30 <- easyXpress::setFlags(edge_flaggedB_drugspecific_30, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayB_30_drugspecificMDHD <- easyXpress::process(raw_flaggedB_drugspecific_30, Metadata_Plate, Metadata_Well)

save(processed_assayB_30_drugspecificMDHD, file = "/data/processed_assayB_30_drugspecificMDHD.RData")

# 50 pixels
#perform model selection where MDHD is removed on an assay by assay selection 
model_selectedB_50_drugspecificMDHD_BENOMYL <- easyXpress::modelSelection(rawB_50)%>% 
  dplyr::filter(drug == 'benomyl' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_50_drugspecificMDHD_CRY5B <- easyXpress::modelSelection(rawB_50)%>% 
  dplyr::filter(drug == 'cry5b' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_50_drugspecificMDHD_DEC <- easyXpress::modelSelection(rawB_50)%>% 
  dplyr::filter(drug == 'dec' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_50_drugspecificMDHD_MORANTEL <- easyXpress::modelSelection(rawB_50)%>% 
  dplyr::filter(drug == 'morantel' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_50_drugspecificMDHD_NIRIDAZOLE <- easyXpress::modelSelection(rawB_50)%>% 
  dplyr::filter(drug == 'niridazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_50_drugspecificMDHD_OXAMNIQUINE <- easyXpress::modelSelection(rawB_50)%>% 
  dplyr::filter(drug == 'oxamniquine' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_50_drugspecificMDHD_PIPERAZINE <- easyXpress::modelSelection(rawB_50)%>% 
  dplyr::filter(drug == 'piperazine' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_50_drugspecificMDHD_PZQ <- easyXpress::modelSelection(rawB_50)%>% 
  dplyr::filter(drug == 'pzq' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedB_50_drugspecificMDHD_CLOSANTEL <- easyXpress::modelSelection(rawB_50)%>% 
  dplyr::filter(drug == 'closantel')

model_selectedB_50_drugspecificMDHD_DERQUANTEL <- easyXpress::modelSelection(rawB_50)%>% 
  dplyr::filter(drug == 'derquantel')

model_selectedB_50_drugspecificMDHD <- model_selectedB_50_drugspecificMDHD_BENOMYL %>% 
  full_join(model_selectedB_50_drugspecificMDHD_CRY5B, by = NULL) %>% 
  full_join(model_selectedB_50_drugspecificMDHD_DEC, by = NULL) %>% 
  full_join(model_selectedB_50_drugspecificMDHD_MORANTEL, by = NULL) %>% 
  full_join(model_selectedB_50_drugspecificMDHD_NIRIDAZOLE, by = NULL) %>% 
  full_join(model_selectedB_50_drugspecificMDHD_OXAMNIQUINE, by = NULL) %>% 
  full_join(model_selectedB_50_drugspecificMDHD_PIPERAZINE, by = NULL) %>% 
  full_join(model_selectedB_50_drugspecificMDHD_PZQ, by =NULL) %>% 
  full_join(model_selectedB_50_drugspecificMDHD_CLOSANTEL, by = NULL) %>% 
  full_join(model_selectedB_50_drugspecificMDHD_DERQUANTEL, by = NULL)

edge_flaggedB_drugspecific_50 <- easyXpress::edgeFlag(model_selectedB_50_drugspecificMDHD, radius=825, center_x=1024, center_y=1024)
raw_flaggedB_drugspecific_50 <- easyXpress::setFlags(edge_flaggedB_drugspecific_50, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayB_50_drugspecificMDHD <- easyXpress::process(raw_flaggedB_drugspecific_50, Metadata_Plate, Metadata_Well)

save(processed_assayB_50_drugspecificMDHD, file = "/data/processed_assayB_50_drugspecificMDHD.RData")



#Steps to seperate out derquantel for image overlays to see what model selection looks like at highest doses 

# model_selectedB_30_drugspecificMDHD_DERQUANTEL <- easyXpress::modelSelection(rawB_30)%>% 
#   dplyr::filter(drug == 'derquantel')
# edge_flaggedB_derquantel_30 <- easyXpress::edgeFlag(model_selectedB_30_drugspecificMDHD_DERQUANTEL, radius=825, center_x=1024, center_y=1024)
# raw_flaggedB_derquantel_30 <- easyXpress::setFlags(edge_flaggedB_derquantel_30 , cluster_flag = TRUE, well_edge_flag = TRUE)
# processed_assayB_30_derquantel_30 <- easyXpress::process(raw_flaggedB_derquantel_30, Metadata_Plate, Metadata_Well)
# save(processed_assayB_30_derquantel_30, file = "/data/processed_assayB_30_derquantel_30.RData")
# 
# model_selectedB_50_drugspecificMDHD_DERQUANTEL <- easyXpress::modelSelection(rawB_50)%>% 
#   dplyr::filter(drug == 'derquantel')
# edge_flaggedB_derquantel_50 <- easyXpress::edgeFlag(model_selectedB_50_drugspecificMDHD_DERQUANTEL, radius=825, center_x=1024, center_y=1024)
# raw_flaggedB_derquantel_50 <- easyXpress::setFlags(edge_flaggedB_derquantel_50 , cluster_flag = TRUE, well_edge_flag = TRUE)
# processed_assayB_50_derquantel_50 <- easyXpress::process(raw_flaggedB_derquantel_50, Metadata_Plate, Metadata_Well)
# save(processed_assayB_50_derquantel_50, file = "/data/validation_30pixels/processed_assayB_50_derquantel_50.RData")

# # Saving processed_data list element to new variable
# proc_data_derqunatel_30 <- processed_assayB_30_derquantel_30[[2]]
# proc_data_derqunatel_50  <- processed_assayB_50_derquantel_50[[2]]
# 


# 
# #highest dose concentration 
# derquantel_assayB_bleach1_p019_H03_30 <- easyXpress::viewWell(proc_data_derqunatel_30, proc_img_dir_derquantel, "p019", "H03", boxplot = TRUE) 
# ggsave(derquantel_assayB_bleach1_p019_H03_30 , file="derquantel_assayB_bleach1_p019_H03_30 .jpg", width=20, height=10, dpi=300)
# 
# derquantel_assayB_bleach2_p023_H01_30 <- easyXpress::viewWell(proc_data_derqunatel_30, proc_img_dir_derquantel, "p023", "H01", boxplot = TRUE) 
# ggsave(derquantel_assayB_bleach2_p023_H01_30 , file="derquantel_assayB_bleach2_p023_H01_30.jpg", width=20, height=10, dpi=300)
# 
# derquantel_assayB_bleach3_p024_H04_30 <- easyXpress::viewWell(proc_data_derqunatel_30, proc_img_dir_derquantel, "p024", "H04", boxplot = TRUE) 
# ggsave(derquantel_assayB_bleach3_p024_H04_30, file="derquantel_assayB_bleach3_p024_H04_30.jpg", width=20, height=10, dpi=300)
# 
# #highest dose concentration 
# derquantel_assayB_bleach1_p019_H03_50 <- easyXpress::viewWell(proc_data_derqunatel_50, proc_img_dir_derquantel, "p019", "H03", boxplot = TRUE) 
# ggsave(derquantel_assayB_bleach1_p019_H03_50 , file="derquantel_assayB_bleach1_p019_H03_50 .jpg", width=20, height=10, dpi=500)
# 
# derquantel_assayB_bleach2_p023_H01_50 <- easyXpress::viewWell(proc_data_derqunatel_50, proc_img_dir_derquantel, "p023", "H01", boxplot = TRUE) 
# ggsave(derquantel_assayB_bleach2_p023_H01_50 , file="derquantel_assayB_bleach2_p023_H01_50.jpg", width=20, height=10, dpi=500)
# 
# derquantel_assayB_bleach3_p024_H04_50 <- easyXpress::viewWell(proc_data_derqunatel_50, proc_img_dir_derquantel, "p024", "H04", boxplot = TRUE) 
# ggsave(derquantel_assayB_bleach3_p024_H04_50, file="derquantel_assayB_bleach3_p024_H04_50.jpg", width=20, height=10, dpi=500)



#############################
#                           #
#          Assay C          #
#                           #
#############################

#Reading in the data
# Define experimental directory and file name

datafileC <- "20210219_assayC_Analysis-20220710.RData"

rawC_30 <- easyXpress::readXpress(filedir = dirsC, rdafile = datafileC, design = TRUE) %>%
  dplyr::filter(Worm_Length > 30) #remove objects smaller than 30 pixels = 100 um 

rawC_50 <- easyXpress::readXpress(filedir = dirsC, rdafile = datafileC, design = TRUE) %>%
  dplyr::filter(Worm_Length > 50) #remove objects smaller than 50 pixels = 165 um 

# Assay C drugs 
# eprinomectin (keep MDHD)
# levamisole (keep MDHD)
# monepantel_LY33414916 (remove MDHD)
# monepantel_LY3348298 (keep MDHD)
# pyrantel_cytrate (remove MDHD)
# selamectin (keep MDHD)


#pixel filter 30 
model_selectedC_30_drugspecificMDHD_EPRINOMECTIN <- easyXpress::modelSelection(rawC_30)%>% 
  dplyr::filter(drug == 'eprinomectin')

model_selectedC_30_drugspecificMDHD_LEVAMISOLE <- easyXpress::modelSelection(rawC_30)%>% 
  dplyr::filter(drug == 'levamisole ')

model_selectedC_30_drugspecificMDHD_MONEPANTEL_LY33414916 <- easyXpress::modelSelection(rawC_30)%>% 
  dplyr::filter(drug == 'monepantel_LY33414916')

model_selectedC_30_drugspecificMDHD_MONEPANTEL_LY3348298 <- easyXpress::modelSelection(rawC_30)%>% 
  dplyr::filter(drug == 'monepantel_LY3348298')

model_selectedC_30_drugspecificMDHD_PYRANTELCITRATE <- easyXpress::modelSelection(rawC_30)%>% 
  dplyr::filter(drug == 'pyrantel_cytrate' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedC_30_drugspecificMDHD_SELAMECTIN <- easyXpress::modelSelection(rawC_30)%>% 
  dplyr::filter(drug == 'selamectin')

model_selectedC_30_drugspecificMDHD  <- model_selectedC_30_drugspecificMDHD_EPRINOMECTIN  %>% 
  full_join(model_selectedC_30_drugspecificMDHD_LEVAMISOLE, by = NULL) %>% 
  full_join(model_selectedC_30_drugspecificMDHD_MONEPANTEL_LY33414916, by = NULL) %>% 
  full_join(model_selectedC_30_drugspecificMDHD_MONEPANTEL_LY3348298, by = NULL) %>% 
  full_join(model_selectedC_30_drugspecificMDHD_PYRANTELCITRATE, by = NULL) %>% 
  full_join(model_selectedC_30_drugspecificMDHD_SELAMECTIN, by = NULL) 

edge_flaggedC_30_drugspecific <- easyXpress::edgeFlag(model_selectedC_30_drugspecificMDHD  , radius=825, center_x=1024, center_y=1024)
raw_flaggedC_30_drugspecific <- easyXpress::setFlags(edge_flaggedC_30_drugspecific , cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayC_30_drugspecificMDHD <- easyXpress::process(raw_flaggedC_30_drugspecific, Metadata_Plate, Metadata_Well)

save(processed_assayC_30_drugspecificMDHD, file = "/data/processed_assayC_30_drugspecificMDHD.RData")


#pixel filter 50 
model_selectedC_50_drugspecificMDHD_EPRINOMECTIN <- easyXpress::modelSelection(rawC_50)%>% 
  dplyr::filter(drug == 'eprinomectin')

model_selectedC_50_drugspecificMDHD_LEVAMISOLE <- easyXpress::modelSelection(rawC_50)%>% 
  dplyr::filter(drug == 'levamisole ')

model_selectedC_50_drugspecificMDHD_MONEPANTEL_LY33414916 <- easyXpress::modelSelection(rawC_50)%>% 
  dplyr::filter(drug == 'monepantel_LY33414916')

model_selectedC_50_drugspecificMDHD_MONEPANTEL_LY3348298 <- easyXpress::modelSelection(rawC_50)%>% 
  dplyr::filter(drug == 'monepantel_LY3348298')

model_selectedC_50_drugspecificMDHD_PYRANTELCITRATE <- easyXpress::modelSelection(rawC_50)%>% 
  dplyr::filter(drug == 'pyrantel_cytrate' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedC_50_drugspecificMDHD_SELAMECTIN <- easyXpress::modelSelection(rawC_50)%>% 
  dplyr::filter(drug == 'selamectin')

model_selectedC_50_drugspecificMDHD  <- model_selectedC_50_drugspecificMDHD_EPRINOMECTIN  %>% 
  full_join(model_selectedC_50_drugspecificMDHD_LEVAMISOLE, by = NULL) %>% 
  full_join(model_selectedC_50_drugspecificMDHD_MONEPANTEL_LY33414916, by = NULL) %>% 
  full_join(model_selectedC_50_drugspecificMDHD_MONEPANTEL_LY3348298, by = NULL) %>% 
  full_join(model_selectedC_50_drugspecificMDHD_PYRANTELCITRATE, by = NULL) %>% 
  full_join(model_selectedC_50_drugspecificMDHD_SELAMECTIN, by = NULL) 

edge_flaggedC_50_drugspecific <- easyXpress::edgeFlag(model_selectedC_50_drugspecificMDHD  , radius=825, center_x=1024, center_y=1024)
raw_flaggedC_50_drugspecific <- easyXpress::setFlags(edge_flaggedC_50_drugspecific , cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayC_50_drugspecificMDHD <- easyXpress::process(raw_flaggedC_50_drugspecific, Metadata_Plate, Metadata_Well)

save(processed_assayC_50_drugspecificMDHD, file = "/data/processed_assayC_50_drugspecificMDHD.RData")


 # Code for doing model overlays for selamectin to look at highest doses
model_selectedC_30_drugspecificMDHD_SELAMECTIN <- easyXpress::modelSelection(rawC_30)%>%
  dplyr::filter(drug == 'selamectin')
edge_flaggedC_selamectin_30 <- easyXpress::edgeFlag(model_selectedC_30_drugspecificMDHD_SELAMECTIN, radius=825, center_x=1024, center_y=1024)
raw_flaggedC_selamectin_30 <- easyXpress::setFlags(edge_flaggedC_selamectin_30, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayC_selamectin_30 <- easyXpress::process(raw_flaggedC_selamectin_30, Metadata_Plate, Metadata_Well)
save(processed_assayC_selamectin_30, file = "/data/validation_30pixels/processed_assayC_selamectin_30.RData")


model_selectedC_50_drugspecificMDHD_SELAMECTIN <- easyXpress::modelSelection(rawC_50)%>%
  dplyr::filter(drug == 'selamectin')
edge_flaggedC_selamectin_50 <- easyXpress::edgeFlag(model_selectedC_50_drugspecificMDHD_SELAMECTIN, radius=825, center_x=1024, center_y=1024)
raw_flaggedC_selamectin_50 <- easyXpress::setFlags(edge_flaggedC_selamectin_50, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayC_selamectin_50 <- easyXpress::process(raw_flaggedC_selamectin_50, Metadata_Plate, Metadata_Well)
save(processed_assayC_selamectin_50, file = "/data/validation_30pixels/processed_assayC_selamectin_50.RData")

# Saving processed_data list element to new variable
proc_data_selamectin_50 <- processed_assayC_selamectin_50[[2]]
proc_data_selamectin_30 <- processed_assayC_selamectin_30[[2]]
# Define processed image directory
proc_img_dir_selamectin <- '/projects/20210219_assayC/Analysis-20220710/processed_images'

#highest dose concentration
selamectin_assayC_bleach1_p037_H03_50 <- easyXpress::viewWell(proc_data_selamectin_50, proc_img_dir_selamectin, "p037", "H03", boxplot = TRUE)
ggsave(selamectin_assayC_bleach1_p037_H03_50, file="selamectin_assayC_bleach1_p037_H03_50.jpg", width=20, height=10, dpi=300)

selamectin_assayC_bleach1_p039_A04_50 <- easyXpress::viewWell(proc_data_selamectin_50, proc_img_dir_selamectin, "p039", "H04", boxplot = TRUE)
ggsave(selamectin_assayC_bleach1_p039_A04_50, file="selamectin_assayC_bleach1_p039_A04_50.jpg", width=20, height=10, dpi=300)

selamectin_assayC_bleach3_p041_A01_50 <- easyXpress::viewWell(proc_data_selamectin_50, proc_img_dir_selamectin, "p041", "H01", boxplot = TRUE)
ggsave(selamectin_assayC_bleach3_p041_A01_50, file="selamectin_assayC_bleach3_p041_A01_50.jpg", width=20, height=10, dpi=300)

#highest dose concentration
selamectin_assayC_bleach1_p037_H03_30 <- easyXpress::viewWell(proc_data_selamectin_30, proc_img_dir_selamectin, "p037", "H03", boxplot = TRUE)
ggsave(selamectin_assayC_bleach1_p037_H03_30, file="selamectin_assayC_bleach1_p037_H03_30.jpg", width=20, height=10, dpi=300)

selamectin_assayC_bleach1_p039_A04_30 <- easyXpress::viewWell(proc_data_selamectin_30, proc_img_dir_selamectin, "p039", "H04", boxplot = TRUE)
ggsave(selamectin_assayC_bleach1_p039_A04_30, file="selamectin_assayC_bleach1_p039_A04_30.jpg", width=20, height=10, dpi=300)

selamectin_assayC_bleach3_p041_A01_30 <- easyXpress::viewWell(proc_data_selamectin_30, proc_img_dir_selamectin, "p041", "H01", boxplot = TRUE)
ggsave(selamectin_assayC_bleach3_p041_A01_30, file="selamectin_assayC_bleach3_p041_A01_30.jpg", width=20, height=10, dpi=300)
# 


#############################
#                           #
#          Assay D          #
#                           #
#############################

#Reading in the data
# Define experimental directory and file name
datafileD <- "20210226_assayD_Analysis-20220710.RData"

rawD_30 <- easyXpress::readXpress(filedir = dirsD, rdafile = datafileD, design = TRUE) %>%
  dplyr::filter(Worm_Length > 30)  #remove objects smaller than 30 pixels = 100 um 

rawD_50 <- easyXpress::readXpress(filedir = dirsD, rdafile = datafileD, design = TRUE) %>%
  dplyr::filter(Worm_Length > 50) #remove objects smaller than 50 pixels = 165 um 

# Assay D drugs
# albendazole (remove MDHD)
# benomyl (remove MDHD)
# closantel (keep MDHD)
# dec (remove MDHD)
# fenbendazole (remove MDHD)
# mebendazole (remove MDHD)
# niridazole (remove MDHD)
# pzq (remove MDHD)
# thiabendazole (remove MDHD)

#perform model selection where MDHD is removed on an assay by assay selection 
# 30 pixels
model_selectedD_30_drugspecificMDHD_albendazole <- easyXpress::modelSelection(rawD_30)%>% 
  dplyr::filter(drug == 'albendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_30_drugspecificMDHD_benomyl <- easyXpress::modelSelection(rawD_30)%>% 
  dplyr::filter(drug == 'benomyl' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_30_drugspecificMDHD_closantel <- easyXpress::modelSelection(rawD_30)%>% 
  dplyr::filter(drug == 'closantel')

model_selectedD_30_drugspecificMDHD_dec <- easyXpress::modelSelection(rawD_30)%>% 
  dplyr::filter(drug == 'dec' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_30_drugspecificMDHD_fenbendazole <- easyXpress::modelSelection(rawD_30)%>% 
  dplyr::filter(drug == 'fenbendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_30_drugspecificMDHD_mebendazole <- easyXpress::modelSelection(rawD_30)%>% 
  dplyr::filter(drug == 'mebendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_30_drugspecificMDHD_niridazole <- easyXpress::modelSelection(rawD_30)%>% 
  dplyr::filter(drug == 'niridazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_30_drugspecificMDHD_pzq <- easyXpress::modelSelection(rawD_30)%>% 
  dplyr::filter(drug == 'pzq' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_30_drugspecificMDHD_thiabendazole <- easyXpress::modelSelection(rawD_30)%>% 
  dplyr::filter(drug == 'thiabendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_30_drugspecificMDHD <- model_selectedD_30_drugspecificMDHD_albendazole %>% 
  full_join(model_selectedD_30_drugspecificMDHD_benomyl, by = NULL) %>% 
  full_join(model_selectedD_30_drugspecificMDHD_closantel, by = NULL) %>% 
  full_join(model_selectedD_30_drugspecificMDHD_dec, by = NULL) %>% 
  full_join(model_selectedD_30_drugspecificMDHD_fenbendazole, by = NULL) %>% 
  full_join(model_selectedD_30_drugspecificMDHD_mebendazole, by = NULL) %>% 
  full_join(model_selectedD_30_drugspecificMDHD_niridazole, by = NULL) %>% 
  full_join(model_selectedD_30_drugspecificMDHD_pzq, by = NULL) %>% 
  full_join(model_selectedD_30_drugspecificMDHD_thiabendazole, by = NULL)

edge_flaggedD_30_drugspecific <- easyXpress::edgeFlag(model_selectedD_30_drugspecificMDHD, radius=825, center_x=1024, center_y=1024)
raw_flaggedD_30_drugspecific <- easyXpress::setFlags(edge_flaggedD_30_drugspecific, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayD_30_drugspecificMDHD <- easyXpress::process(raw_flaggedD_30_drugspecific, Metadata_Plate, Metadata_Well)

save(processed_assayD_30_drugspecificMDHD, file = "/data/processed_assayD_30_drugspecificMDHD.RData")


#50 pixels 
model_selectedD_50_drugspecificMDHD_albendazole <- easyXpress::modelSelection(rawD_50)%>% 
  dplyr::filter(drug == 'albendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_50_drugspecificMDHD_benomyl <- easyXpress::modelSelection(rawD_50)%>% 
  dplyr::filter(drug == 'benomyl' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_50_drugspecificMDHD_closantel <- easyXpress::modelSelection(rawD_50)%>% 
  dplyr::filter(drug == 'closantel')

model_selectedD_50_drugspecificMDHD_dec <- easyXpress::modelSelection(rawD_50)%>% 
  dplyr::filter(drug == 'dec' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_50_drugspecificMDHD_fenbendazole <- easyXpress::modelSelection(rawD_50)%>% 
  dplyr::filter(drug == 'fenbendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_50_drugspecificMDHD_mebendazole <- easyXpress::modelSelection(rawD_50)%>% 
  dplyr::filter(drug == 'mebendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_50_drugspecificMDHD_niridazole <- easyXpress::modelSelection(rawD_50)%>% 
  dplyr::filter(drug == 'niridazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_50_drugspecificMDHD_pzq <- easyXpress::modelSelection(rawD_50)%>% 
  dplyr::filter(drug == 'pzq' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_50_drugspecificMDHD_thiabendazole <- easyXpress::modelSelection(rawD_50)%>% 
  dplyr::filter(drug == 'thiabendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedD_50_drugspecificMDHD <- model_selectedD_50_drugspecificMDHD_albendazole %>% 
  full_join(model_selectedD_50_drugspecificMDHD_benomyl, by = NULL) %>% 
  full_join(model_selectedD_50_drugspecificMDHD_closantel, by = NULL) %>% 
  full_join(model_selectedD_50_drugspecificMDHD_dec, by = NULL) %>% 
  full_join(model_selectedD_50_drugspecificMDHD_fenbendazole, by = NULL) %>% 
  full_join(model_selectedD_50_drugspecificMDHD_mebendazole, by = NULL) %>% 
  full_join(model_selectedD_50_drugspecificMDHD_niridazole, by = NULL) %>% 
  full_join(model_selectedD_50_drugspecificMDHD_pzq, by = NULL) %>% 
  full_join(model_selectedD_50_drugspecificMDHD_thiabendazole, by = NULL)

edge_flaggedD_50_drugspecific <- easyXpress::edgeFlag(model_selectedD_50_drugspecificMDHD, radius=825, center_x=1024, center_y=1024)
raw_flaggedD_50_drugspecific <- easyXpress::setFlags(edge_flaggedD_50_drugspecific, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayD_50_drugspecificMDHD <- easyXpress::process(raw_flaggedD_50_drugspecific, Metadata_Plate, Metadata_Well)

save(processed_assayD_50_drugspecificMDHD, file = "/data/processed_assayD_50_drugspecificMDHD.RData")


#to test thiabendazole 50 vs 30 pixels looking at image overlays at highest doses 
# 
# #30 pixels
# edge_flaggedD_30_thiabendazole <- easyXpress::edgeFlag(model_selectedD_30_drugspecificMDHD_thiabendazole, radius=825, center_x=1024, center_y=1024)
# raw_flaggedD_30_thiabendazole <- easyXpress::setFlags(edge_flaggedD_30_thiabendazole, cluster_flag = TRUE, well_edge_flag = TRUE)
# processed_assayD_30_thiabendazole <- easyXpress::process(raw_flaggedD_30_thiabendazole , Metadata_Plate, Metadata_Well)
# 
# save(processed_assayD_30_thiabendazole, file = "//data/validation_30pixels/processed_assayD_30_thiabendazole.RData")
# 
# 
# #50 pixels
# edge_flaggedD_50_thiabendazole <- easyXpress::edgeFlag(model_selectedD_50_drugspecificMDHD_thiabendazole, radius=825, center_x=1024, center_y=1024)
# raw_flaggedD_50_thiabendazole <- easyXpress::setFlags(edge_flaggedD_50_thiabendazole, cluster_flag = TRUE, well_edge_flag = TRUE)
# processed_assayD_50_thiabendazole <- easyXpress::process(raw_flaggedD_50_thiabendazole , Metadata_Plate, Metadata_Well)
# 
# save(processed_assayD_50_thiabendazole, file = "/data/validation_30pixels/processed_assayD_50_thiabendazole.RData")
# 
# #view high dose wells to compare what the models are retaining 
# 
# # Saving processed_data list element to new variable
# proc_data_thiabendazole_50 <- processed_assayD_50_thiabendazole[[2]]
# 
# # Define processed image directory
# proc_img_dir_thiabendazole <- '/projects/20210226_assayD/Analysis-20220710/processed_images'
# 
# #highest dose concentration 
# thiabendazole_assayD_bleach1_p013_H09_50 <- easyXpress::viewWell(proc_data_thiabendazole_50 , proc_img_dir_thiabendazole, "p013", "H09", boxplot = TRUE) 
# ggsave(thiabendazole_assayD_bleach1_p013_H09_50, file="thiabendazole_assayD_bleach1_p013_H09_50.jpg", width=20, height=10, dpi=300)
# 
# thiabendazole_assayD_bleach2_p016_H10_50 <- easyXpress::viewWell(proc_data_thiabendazole_50 , proc_img_dir_thiabendazole, "p016", "H10", boxplot = TRUE) 
# ggsave(thiabendazole_assayD_bleach2_p016_H10_50, file="thiabendazole_assayD_bleach2_p016_H10_50.jpg", width=20, height=10, dpi=300)
# 
# thiabendazole_assayD_bleach3_p017_H07_50 <- easyXpress::viewWell(proc_data_thiabendazole_50 , proc_img_dir_thiabendazole, "p017", "H07", boxplot = TRUE) 
# ggsave(thiabendazole_assayD_bleach3_p017_H07_50 , file="thiabendazole_assayD_bleach3_p017_H07_50.jpg", width=20, height=10, dpi=300)
# 
# 
# proc_data_thiabendazole_30 <- processed_assayD_30_thiabendazole[[2]]
# 
# # Define processed image directory
# proc_img_dir_thiabendazole <- '/projects/20210226_assayD/Analysis-20220710/processed_images'
# 
# #highest dose concentration 
# thiabendazole_assayD_bleach1_p013_H09_30 <- easyXpress::viewWell(proc_data_thiabendazole_30 , proc_img_dir_thiabendazole, "p013", "H09", boxplot = TRUE) 
# ggsave(thiabendazole_assayD_bleach1_p013_H09_30, file="thiabendazole_assayD_bleach1_p013_H09_30.jpg", width=20, height=10, dpi=300)
# 
# thiabendazole_assayD_bleach2_p016_H10_30 <- easyXpress::viewWell(proc_data_thiabendazole_30 , proc_img_dir_thiabendazole, "p016", "H10", boxplot = TRUE) 
# ggsave(thiabendazole_assayD_bleach2_p016_H10_30, file="thiabendazole_assayD_bleach2_p016_H10_30.jpg", width=20, height=10, dpi=300)
# 
# thiabendazole_assayD_bleach3_p017_H07_30 <- easyXpress::viewWell(proc_data_thiabendazole_30 , proc_img_dir_thiabendazole, "p017", "H07", boxplot = TRUE) 
# ggsave(thiabendazole_assayD_bleach3_p017_H07_30 , file="thiabendazole_assayD_bleach3_p017_H07_30.jpg", width=20, height=10, dpi=300)
# 
# #CLOSEANTEL 30 vs 50 pixel test 
# 
# 
# #50 pixels
# edge_flaggedD_50_closantel <- easyXpress::edgeFlag(model_selectedD_50_drugspecificMDHD_closantel, radius=825, center_x=1024, center_y=1024)
# raw_flaggedD_50_closantel <- easyXpress::setFlags(edge_flaggedD_50_closantel, cluster_flag = TRUE, well_edge_flag = TRUE)
# processed_assayD_50_closantel <- easyXpress::process(raw_flaggedD_50_closantel, Metadata_Plate, Metadata_Well)
# 
# save(processed_assayD_50_closantel, file = "/data/validation_30pixels/processed_assayD_50_closantel.RData")
# 
# #30 pixels
# edge_flaggedD_30_closantel <- easyXpress::edgeFlag(model_selectedD_30_drugspecificMDHD_closantel, radius=825, center_x=1024, center_y=1024)
# raw_flaggedD_30_closantel <- easyXpress::setFlags(edge_flaggedD_30_closantel, cluster_flag = TRUE, well_edge_flag = TRUE)
# processed_assayD_30_closantel <- easyXpress::process(raw_flaggedD_30_closantel, Metadata_Plate, Metadata_Well)
# 
# save(processed_assayD_30_closantel, file = "/data/validation_30pixels/processed_assayD_30_closantel.RData")
# 
# # Saving processed_data list element to new variable
# proc_data_closantel_30 <- processed_assayD_30_closantel[[2]]
# 
# # Define processed image directory
# proc_img_dir_closantel <- '/projects/20210226_assayD/Analysis-20220710/processed_images'
# 
# #highest dose concentration 
# closantel_assayD_bleach1_p025_H09_30 <- easyXpress::viewWell(proc_data_closantel_30, proc_img_dir_closantel, "p025", "H09", boxplot = TRUE) 
# ggsave(closantel_assayD_bleach1_p025_H09_30 , file="closantel_assayD_bleach1_p025_H09_30.jpg", width=20, height=10, dpi=300)
# 
# closantel_assayD_bleach1_p027_H10_30 <- easyXpress::viewWell(proc_data_closantel_30, proc_img_dir_closantel, "p027", "H10", boxplot = TRUE) 
# ggsave(closantel_assayD_bleach1_p027_H10_30 , file="closantel_assayD_bleach1_p027_H10_30 .jpg", width=20, height=10, dpi=300)
# 
# closantel_assayD_bleach3_p029_H07_30 <- easyXpress::viewWell(proc_data_closantel_30, proc_img_dir_closantel, "p029", "H07", boxplot = TRUE) 
# ggsave(closantel_assayD_bleach3_p029_H07_30 , file="closantel_assayD_bleach3_p029_H07_30.jpg", width=20, height=10, dpi=300)
# 
# 
# #highest dose concentration 
# proc_data_closantel_50 <- processed_assayD_50_closantel[[2]]
# closantel_assayD_bleach1_p025_H09_50 <- easyXpress::viewWell(proc_data_closantel_50, proc_img_dir_closantel, "p025", "H09", boxplot = TRUE) 
# ggsave(closantel_assayD_bleach1_p025_H09_50 , file="closantel_assayD_bleach1_p025_H09_50.jpg", width=20, height=10, dpi=500)
# 
# closantel_assayD_bleach1_p027_H10_50 <- easyXpress::viewWell(proc_data_closantel_50, proc_img_dir_closantel, "p027", "H10", boxplot = TRUE) 
# ggsave(closantel_assayD_bleach1_p027_H10_50 , file="closantel_assayD_bleach1_p027_H10_50 .jpg", width=20, height=10, dpi=500)
# 
# closantel_assayD_bleach3_p029_H07_50 <- easyXpress::viewWell(proc_data_closantel_50, proc_img_dir_closantel, "p029", "H07", boxplot = TRUE) 
# ggsave(closantel_assayD_bleach3_p029_H07_50 , file="closantel_assayD_bleach3_p029_H07_50.jpg", width=20, height=10, dpi=500)



#############################
#                           #
#          Assay E          #
#                           #
#############################

# Define experimental directory and file name
datafileE <- "20210312_assayE_Analysis-20220710.RData"

rawE_50<- easyXpress::readXpress(filedir = dirsE, rdafile = datafileE, design = TRUE) %>%
  dplyr::filter(Worm_Length > 50) #remove objects smaller than 50 pixels = 165 um 

rawE_30 <- easyXpress::readXpress(filedir = dirsE, rdafile = datafileE, design = TRUE) %>%
  dplyr::filter(Worm_Length > 30) #remove objects smaller than 30 pixels = 100 um 

# Assay E drugs 
# abamectin (keep MDHD)
# benomyl (remove MDHD)
# closantel (keep MDHD)
# cry5b (remove MDHD)
# derquantel (keep MDHD)
# emodepside (keep MDHD)
# levamisole (keep MDHD)
# milbemycin (keep MDHD)
# moxidectin (remove MDHD)

# 30 pixels 
model_selectedE_30_drugspecificMDHD_abamectin <- easyXpress::modelSelection(rawE_30)%>% 
  dplyr::filter(drug == 'abamectin')

model_selectedE_30_drugspecificMDHD_benomyl <- easyXpress::modelSelection(rawE_30)%>% 
  dplyr::filter(drug == 'benomyl' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedE_30_drugspecificMDHD_closantel <- easyXpress::modelSelection(rawE_30)%>% 
  dplyr::filter(drug == 'closantel')

model_selectedE_30_drugspecificMDHD_cry5b <- easyXpress::modelSelection(rawE_30)%>% 
  dplyr::filter(drug == 'cry5b' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedE_30_drugspecificMDHD_derquantel <- easyXpress::modelSelection(rawE_30)%>% 
  dplyr::filter(drug == 'derquantel')

model_selectedE_30_drugspecificMDHD_emodepside <- easyXpress::modelSelection(rawE_30)%>% 
  dplyr::filter(drug == 'emodepside')

model_selectedE_30_drugspecificMDHD_levamisole <- easyXpress::modelSelection(rawE_30)%>% 
  dplyr::filter(drug == 'levamisole')

model_selectedE_30_drugspecificMDHD_milbemycin <- easyXpress::modelSelection(rawE_30)%>% 
  dplyr::filter(drug == 'milbemycin')

model_selectedE_30_drugspecificMDHD_moxidectin <- easyXpress::modelSelection(rawE_30)%>% 
  dplyr::filter(drug == 'moxidectin')

model_selectedE_30_drugspecificMDHD <- model_selectedE_30_drugspecificMDHD_abamectin  %>% 
  full_join(model_selectedE_30_drugspecificMDHD_benomyl, by = NULL) %>% 
  full_join(model_selectedE_30_drugspecificMDHD_closantel, by = NULL) %>% 
  full_join(model_selectedE_30_drugspecificMDHD_cry5b, by = NULL) %>% 
  full_join(model_selectedE_30_drugspecificMDHD_derquantel, by = NULL) %>% 
  full_join(model_selectedE_30_drugspecificMDHD_emodepside, by = NULL) %>% 
  full_join(model_selectedE_30_drugspecificMDHD_levamisole, by = NULL) %>% 
  full_join(model_selectedE_30_drugspecificMDHD_milbemycin, by = NULL) %>% 
  full_join(model_selectedE_30_drugspecificMDHD_moxidectin, by = NULL)

edge_flaggedE_30_drugspecific <- easyXpress::edgeFlag(model_selectedE_30_drugspecificMDHD, radius=825, center_x=1024, center_y=1024)
raw_flaggedE3_30_drugspecific <- easyXpress::setFlags(edge_flaggedE_30_drugspecific, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayE_30_drugspecificMDHD <- easyXpress::process(raw_flaggedE3_30_drugspecific, Metadata_Plate, Metadata_Well)

save(processed_assayE_30_drugspecificMDHD, file = "/data/processed_assayE_30_drugspecificMDHD.RData")


# 50 pixels 
model_selectedE_50_drugspecificMDHD_abamectin <- easyXpress::modelSelection(rawE_50)%>% 
  dplyr::filter(drug == 'abamectin')

model_selectedE_50_drugspecificMDHD_benomyl <- easyXpress::modelSelection(rawE_50)%>% 
  dplyr::filter(drug == 'benomyl' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedE_50_drugspecificMDHD_closantel <- easyXpress::modelSelection(rawE_50)%>% 
  dplyr::filter(drug == 'closantel')

model_selectedE_50_drugspecificMDHD_cry5b <- easyXpress::modelSelection(rawE_50)%>% 
  dplyr::filter(drug == 'cry5b' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedE_50_drugspecificMDHD_derquantel <- easyXpress::modelSelection(rawE_50)%>% 
  dplyr::filter(drug == 'derquantel')

model_selectedE_50_drugspecificMDHD_emodepside <- easyXpress::modelSelection(rawE_50)%>% 
  dplyr::filter(drug == 'emodepside')

model_selectedE_50_drugspecificMDHD_levamisole <- easyXpress::modelSelection(rawE_50)%>% 
  dplyr::filter(drug == 'levamisole')

model_selectedE_50_drugspecificMDHD_milbemycin <- easyXpress::modelSelection(rawE_50)%>% 
  dplyr::filter(drug == 'milbemycin')

model_selectedE_50_drugspecificMDHD_moxidectin <- easyXpress::modelSelection(rawE_50)%>% 
  dplyr::filter(drug == 'moxidectin')

model_selectedE_50_drugspecificMDHD <- model_selectedE_50_drugspecificMDHD_abamectin  %>% 
  full_join(model_selectedE_50_drugspecificMDHD_benomyl, by = NULL) %>% 
  full_join(model_selectedE_50_drugspecificMDHD_closantel, by = NULL) %>% 
  full_join(model_selectedE_50_drugspecificMDHD_cry5b, by = NULL) %>% 
  full_join(model_selectedE_50_drugspecificMDHD_derquantel, by = NULL) %>% 
  full_join(model_selectedE_50_drugspecificMDHD_emodepside, by = NULL) %>% 
  full_join(model_selectedE_50_drugspecificMDHD_levamisole, by = NULL) %>% 
  full_join(model_selectedE_50_drugspecificMDHD_milbemycin, by = NULL) %>% 
  full_join(model_selectedE_50_drugspecificMDHD_moxidectin, by = NULL)

edge_flaggedE_50_drugspecific <- easyXpress::edgeFlag(model_selectedE_50_drugspecificMDHD, radius=825, center_x=1050, center_y=1050)
raw_flaggedE3_50_drugspecific <- easyXpress::setFlags(edge_flaggedE_50_drugspecific, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayE_50_drugspecificMDHD <- easyXpress::process(raw_flaggedE3_50_drugspecific, Metadata_Plate, Metadata_Well)

save(processed_assayE_50_drugspecificMDHD, file = "/data/processed_assayE_50_drugspecificMDHD.RData")


# to test thiabendazole 50 vs 30 pixels looking at image overlays at highest doses

#30 pixels abamecint 
edge_flaggedD_30_abamectin <- easyXpress::edgeFlag(model_selectedE_30_drugspecificMDHD_abamectin, radius=825, center_x=1024, center_y=1024)
raw_flaggedD_30_abamectin <- easyXpress::setFlags(edge_flaggedD_30_abamectin, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayD_30_abamectin <- easyXpress::process(raw_flaggedD_30_abamectin , Metadata_Plate, Metadata_Well)

save(processed_assayD_30_abamectin, file = "/data/validation_30pixels/processed_assayD_30_abamectin.RData")


#50 pixels
edge_flaggedD_50_abamectin <- easyXpress::edgeFlag(model_selectedE_50_drugspecificMDHD_abamectin, radius=825, center_x=1024, center_y=1024)
raw_flaggedD_50_abamectin <- easyXpress::setFlags(edge_flaggedD_50_abamectin, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayD_50_abamectin <- easyXpress::process(raw_flaggedD_50_abamectin , Metadata_Plate, Metadata_Well)

save(processed_assayD_50_abamectin, file = "/data/validation_30pixels/processed_assayD_50_abamectin.RData")

#view high dose wells to compare what the models are retaining

# Saving processed_data list element to new variable
proc_data_abamectin_50 <- processed_assayD_50_abamectin[[2]]
proc_data_abamectin_30 <- processed_assayD_30_abamectin[[2]]

# Define processed image directory
proc_img_dir_abamectin <- '/projects/20210312_assayE/Analysis-20220710/processed_images'

#highest dose concentration
#bleach 1, 50 pixels 
abamectin_assayE_bleach1_p050_H09_50 <- easyXpress::viewWell(proc_data_abamectin_50 , proc_img_dir_abamectin, "p050", "H09", boxplot = TRUE)
ggsave(abamectin_assayE_bleach1_p050_H09_50, file="abamectin_assayE_bleach1_p050_H09_50.jpg", width=20, height=10, dpi=300)

#bleach 2, 50 pixels 
abamectin_assayE_bleach1_p051_H010_50 <- easyXpress::viewWell(proc_data_abamectin_50 , proc_img_dir_abamectin, "p051", "H10", boxplot = TRUE)
ggsave(abamectin_assayE_bleach1_p051_H010_50, file="abamectin_assayE_bleach1_p051_H010_50.jpg", width=20, height=10, dpi=300)

#bleach 3 , 50 pixels 
abamectin_assayE_bleach1_p053_H07_50 <- easyXpress::viewWell(proc_data_abamectin_50 , proc_img_dir_abamectin, "p053", "H07", boxplot = TRUE)
ggsave(abamectin_assayE_bleach1_p053_H07_50, file="abamectin_assayE_bleach1_p053_H07_50.jpg", width=20, height=10, dpi=300)


#bleach 1, 30 pixels 
abamectin_assayE_bleach1_p030_H09_30 <- easyXpress::viewWell(proc_data_abamectin_30 , proc_img_dir_abamectin, "p030", "H09", boxplot = TRUE)
ggsave(abamectin_assayE_bleach1_p030_H09_30, file="abamectin_assayE_bleach1_p030_H09_30.jpg", width=20, height=10, dpi=300)

#bleach 2, 30 pixels 
abamectin_assayE_bleach1_p051_H010_30 <- easyXpress::viewWell(proc_data_abamectin_30 , proc_img_dir_abamectin, "p051", "H10", boxplot = TRUE)
ggsave(abamectin_assayE_bleach1_p051_H010_30, file="abamectin_assayE_bleach1_p051_H010_30.jpg", width=20, height=10, dpi=300)

#bleach 3 , 30 pixels 
abamectin_assayE_bleach1_p053_H07_30 <- easyXpress::viewWell(proc_data_abamectin_30 , proc_img_dir_abamectin, "p053", "H07", boxplot = TRUE)
ggsave(abamectin_assayE_bleach1_p053_H07_30, file="abamectin_assayE_bleach1_p053_H07_30.jpg", width=20, height=10, dpi=300)


#30 pixels emodepside
edge_flaggedD_30_emodepside <- easyXpress::edgeFlag(model_selectedE_30_drugspecificMDHD_emodepside, radius=825, center_x=1024, center_y=1024)
raw_flaggedD_30_emodepside <- easyXpress::setFlags(edge_flaggedD_30_emodepside, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayD_30_emodepside <- easyXpress::process(raw_flaggedD_30_emodepside , Metadata_Plate, Metadata_Well)

save(processed_assayD_30_emodepside, file = "/data/validation_30pixels/processed_assayD_30_emodepside.RData")

proc_data_emodepside_30 <- processed_assayD_30_emodepside[[2]]

edge_flaggedD_50_emodepside <- easyXpress::edgeFlag(model_selectedE_50_drugspecificMDHD_emodepside, radius=825, center_x=1024, center_y=1024)
raw_flaggedD_50_emodepside <- easyXpress::setFlags(edge_flaggedD_50_emodepside, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayD_50_emodepside <- easyXpress::process(raw_flaggedD_50_emodepside , Metadata_Plate, Metadata_Well)

save(processed_assayD_50_emodepside, file = "/data/validation_30pixels/processed_assayD_50_emodepside.RData")

proc_data_emodepside_50 <- processed_assayD_50_emodepside[[2]]

#highest dose concentration
#bleach 1, 50 pixels 
emodepside_assayE_bleach1_p025_H03_50 <- easyXpress::viewWell(proc_data_emodepside_50, proc_img_dir_abamectin, "p025", "H03", boxplot = TRUE)
ggsave(emodepside_assayE_bleach1_p025_H03_50, file="emodepside_assayE_bleach1_p025_H03_50.jpg", width=20, height=10, dpi=300)

#bleach 2, 50 pixels 
emodepside_assayE_bleach2_p027_H10_50 <- easyXpress::viewWell(proc_data_emodepside_50, proc_img_dir_abamectin, "p027", "H10", boxplot = TRUE)
ggsave(emodepside_assayE_bleach2_p027_H10_50, file="emodepside_assayE_bleach2_p027_H10_50.jpg", width=20, height=10, dpi=300)

#bleach 3, 50 pixels 
emodepside_assayE_bleach3_p029_H01_50 <- easyXpress::viewWell(proc_data_emodepside_50, proc_img_dir_abamectin, "p029", "H01", boxplot = TRUE)
ggsave(emodepside_assayE_bleach3_p029_H01_50, file="emodepside_assayE_bleach3_p029_H01_50.jpg", width=20, height=10, dpi=300)

#highest dose concentration
#bleach 1, 30 pixels 
emodepside_assayE_bleach1_p025_H03_30 <- easyXpress::viewWell(proc_data_emodepside_30, proc_img_dir_abamectin, "p025", "H03", boxplot = TRUE)
ggsave(emodepside_assayE_bleach1_p025_H03_30, file="emodepside_assayE_bleach1_p025_H03_30.jpg", width=20, height=10, dpi=300)

#bleach 2, 30 pixels 
emodepside_assayE_bleach2_p027_H10_30 <- easyXpress::viewWell(proc_data_emodepside_30, proc_img_dir_abamectin, "p027", "H10", boxplot = TRUE)
ggsave(emodepside_assayE_bleach2_p027_H10_30, file="emodepside_assayE_bleach2_p027_H10_30.jpg", width=20, height=10, dpi=300)

#bleach 3, 30 pixels 
emodepside_assayE_bleach3_p029_H01_30 <- easyXpress::viewWell(proc_data_emodepside_30, proc_img_dir_abamectin, "p029", "H01", boxplot = TRUE)
ggsave(emodepside_assayE_bleach3_p029_H01_30, file="emodepside_assayE_bleach3_p029_H01_30.jpg", width=20, height=10, dpi=300)


#Levamisole 
#30 pixels Levamisole
edge_flaggedE_30_Levamisole <- easyXpress::edgeFlag(model_selectedE_30_drugspecificMDHD_levamisole, radius=825, center_x=1024, center_y=1024)
raw_flaggedE_30_Levamisole <- easyXpress::setFlags(edge_flaggedE_30_Levamisole, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayE_30_Levamisole <- easyXpress::process(raw_flaggedE_30_Levamisole , Metadata_Plate, Metadata_Well)

save(processed_assayE_30_Levamisole, file = "/data/validation_30pixels/processed_assayE_30_Levamisole.RData")

proc_data_Levamisole_30 <- processed_assayE_30_Levamisole[[2]]

#highest dose concentration

#bleach 1, 30 pixels 
levamisole_assayE_bleach1_p014_H03_30 <- easyXpress::viewWell(proc_data_Levamisole_30 , proc_img_dir_abamectin, "p014", "H03", boxplot = TRUE)
ggsave(levamisole_assayE_bleach1_p014_H03_30, file="levamisole_assayE_bleach1_p014_H03_30.jpg", width=20, height=10, dpi=300)

#bleach 2, 30 pixels 
levamisole_assayE_bleach2_p015_H04_30 <- easyXpress::viewWell(proc_data_Levamisole_30 , proc_img_dir_abamectin, "p015", "H04", boxplot = TRUE)
ggsave(levamisole_assayE_bleach2_p015_H04_30 , file="levamisole_assayE_bleach2_p015_H04_30.jpg", width=20, height=10, dpi=300)

#bleach 3, 30 pixels 
levamisole_assayE_bleach3_p017_H07_30 <- easyXpress::viewWell(proc_data_Levamisole_30 , proc_img_dir_abamectin, "p017", "H07", boxplot = TRUE)
ggsave(levamisole_assayE_bleach3_p017_H07_30, file="levamisole_assayE_bleach3_p017_H07_30.jpg", width=20, height=10, dpi=300)



#Levamisole 
#50 pixels Levamisole
edge_flaggedE_50_Levamisole <- easyXpress::edgeFlag(model_selectedE_50_drugspecificMDHD_levamisole, radius=825, center_x=1024, center_y=1024)
raw_flaggedE_50_Levamisole <- easyXpress::setFlags(edge_flaggedE_50_Levamisole, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayE_50_Levamisole <- easyXpress::process(raw_flaggedE_50_Levamisole , Metadata_Plate, Metadata_Well)

save(processed_assayE_50_Levamisole, file = "/data/validation_30pixels/processed_assayE_50_Levamisole.RData")

proc_data_Levamisole_50 <- processed_assayE_50_Levamisole[[2]]

#highest dose concentration

#bleach 1, 50 pixels 
levamisole_assayE_bleach1_p014_H03_50 <- easyXpress::viewWell(proc_data_Levamisole_50 , proc_img_dir_abamectin, "p014", "H03", boxplot = TRUE)
ggsave(levamisole_assayE_bleach1_p014_H03_50, file="levamisole_assayE_bleach1_p014_H03_50.jpg", width=20, height=10, dpi=500)

#bleach 2, 50 pixels 
levamisole_assayE_bleach2_p015_H04_50 <- easyXpress::viewWell(proc_data_Levamisole_50 , proc_img_dir_abamectin, "p015", "H04", boxplot = TRUE)
ggsave(levamisole_assayE_bleach2_p015_H04_50 , file="levamisole_assayE_bleach2_p015_H04_50.jpg", width=20, height=10, dpi=500)

#bleach 3, 50 pixels 
levamisole_assayE_bleach3_p017_H07_50 <- easyXpress::viewWell(proc_data_Levamisole_50 , proc_img_dir_abamectin, "p017", "H07", boxplot = TRUE)
ggsave(levamisole_assayE_bleach3_p017_H07_50, file="levamisole_assayE_bleach3_p017_H07_50.jpg", width=20, height=10, dpi=500)


#############################
#                           #
#          Assay F          #
#                           #
#############################

# Define experimental directory and file name
datafileF <- "20210319_assayF_Analysis-20220711.RData"

rawF_50<- easyXpress::readXpress(filedir = dirsF, rdafile = datafileF, design = TRUE) %>%
  dplyr::filter(Worm_Length > 50) #remove objects smaller than 50 pixels = 165 um 

rawF_30<- easyXpress::readXpress(filedir = dirsF, rdafile = datafileF, design = TRUE) %>%
  dplyr::filter(Worm_Length > 30) #remove objects smaller than 30 pixels = 100 um 

# Assay F Drugs
# cry5b (remove MDHD)
# eprinomectin (keep MDHD)
# levamisole (keep MDHD)
# monepantel_LY33414916 (remove MDHD)
# monepantel_LY3348298 (keep MDHD)
# piperazine (remove MDHD)
# selamectin (keep MDHD)

#30 pixels 
model_selectedF_30_drugspecificMDHD_cry5b <- easyXpress::modelSelection(rawF_30)%>% 
  dplyr::filter(drug == 'cry5b' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedF_30_drugspecificMDHD_eprinomectin <- easyXpress::modelSelection(rawF_30)%>% 
  dplyr::filter(drug == 'eprinomectin')

model_selectedF_30_drugspecificMDHD_levamisole <- easyXpress::modelSelection(rawF_30)%>% 
  dplyr::filter(drug == 'levamisole')

model_selectedF_30_drugspecificMDHDc_monepantel_LY33414916 <- easyXpress::modelSelection(rawF_30)%>% 
  dplyr::filter(drug == 'monepantel_LY33414916')

model_selectedF_30_drugspecificMDHD_monepantel_LY3348298 <- easyXpress::modelSelection(rawF_30)%>% 
  dplyr::filter(drug == 'monepantel_LY3348298')

model_selectedF_30_drugspecificMDHD_piperazine <- easyXpress::modelSelection(rawF_30)%>% 
  dplyr::filter(drug == 'piperazine' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedF_30_drugspecificMDHD_selamectin <- easyXpress::modelSelection(rawF_30)%>% 
  dplyr::filter(drug == 'selamectin')

model_selectedF_30_drugspecificMDHD  <- model_selectedF_30_drugspecificMDHD_cry5b %>% 
  full_join(model_selectedF_30_drugspecificMDHD_eprinomectin , by = NULL) %>% 
  full_join(model_selectedF_30_drugspecificMDHD_levamisole, by = NULL) %>% 
  full_join(model_selectedF_30_drugspecificMDHDc_monepantel_LY33414916, by = NULL) %>% 
  full_join(model_selectedF_30_drugspecificMDHD_monepantel_LY3348298 , by = NULL) %>% 
  full_join(model_selectedF_30_drugspecificMDHD_piperazine, by = NULL) %>% 
  full_join(model_selectedF_30_drugspecificMDHD_selamectin, by = NULL)

edge_flaggedF_30_drugspecific <- easyXpress::edgeFlag(model_selectedF_30_drugspecificMDHD , radius=825, center_x=1024, center_y=1024)
raw_flaggedF_30_drugspecific <- easyXpress::setFlags(edge_flaggedF_30_drugspecific, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayF_30_drugspecificMDHD <- easyXpress::process(raw_flaggedF_30_drugspecific, Metadata_Plate, Metadata_Well)

save(processed_assayF_30_drugspecificMDHD, file = "/data/processed_assayF_30_drugspecificMDHD.RData")


#50 pixels
model_selectedF_50_drugspecificMDHD_cry5b <- easyXpress::modelSelection(rawF_50)%>% 
  dplyr::filter(drug == 'cry5b' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedF_50_drugspecificMDHD_eprinomectin <- easyXpress::modelSelection(rawF_50)%>% 
  dplyr::filter(drug == 'eprinomectin')

model_selectedF_50_drugspecificMDHD_levamisole <- easyXpress::modelSelection(rawF_50)%>% 
  dplyr::filter(drug == 'levamisole')

model_selectedF_50_drugspecificMDHDc_monepantel_LY33414916 <- easyXpress::modelSelection(rawF_50)%>% 
  dplyr::filter(drug == 'monepantel_LY33414916')

model_selectedF_50_drugspecificMDHD_monepantel_LY3348298 <- easyXpress::modelSelection(rawF_50)%>% 
  dplyr::filter(drug == 'monepantel_LY3348298')

model_selectedF_50_drugspecificMDHD_piperazine <- easyXpress::modelSelection(rawF_50)%>% 
  dplyr::filter(drug == 'piperazine' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedF_50_drugspecificMDHD_selamectin <- easyXpress::modelSelection(rawF_50)%>% 
  dplyr::filter(drug == 'selamectin')

model_selectedF_50_drugspecificMDHD  <- model_selectedF_50_drugspecificMDHD_cry5b %>% 
  full_join(model_selectedF_50_drugspecificMDHD_eprinomectin , by = NULL) %>% 
  full_join(model_selectedF_50_drugspecificMDHD_levamisole, by = NULL) %>% 
  full_join(model_selectedF_50_drugspecificMDHDc_monepantel_LY33414916, by = NULL) %>% 
  full_join(model_selectedF_50_drugspecificMDHD_monepantel_LY3348298 , by = NULL) %>% 
  full_join(model_selectedF_50_drugspecificMDHD_piperazine, by = NULL) %>% 
  full_join(model_selectedF_50_drugspecificMDHD_selamectin, by = NULL)

edge_flaggedF_50_drugspecific <- easyXpress::edgeFlag(model_selectedF_50_drugspecificMDHD , radius=825, center_x=1050, center_y=1050)
raw_flaggedF_50_drugspecific <- easyXpress::setFlags(edge_flaggedF_50_drugspecific, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayF_50_drugspecificMDHD <- easyXpress::process(raw_flaggedF_50_drugspecific, Metadata_Plate, Metadata_Well)

save(processed_assayF_50_drugspecificMDHD, file = "/data/processed_assayF_50_drugspecificMDHD.RData")


#50 pixels eprinomectin
edge_flaggedF_50_eprinomectin <- easyXpress::edgeFlag(model_selectedF_50_drugspecificMDHD_eprinomectin, radius=825, center_x=1024, center_y=1024)
raw_flaggedF_50_eprinomectin <- easyXpress::setFlags(edge_flaggedF_50_eprinomectin, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayF_50_eprinomectin <- easyXpress::process(raw_flaggedF_50_eprinomectin, Metadata_Plate, Metadata_Well)

save(processed_assayF_50_eprinomectin, file = "/data/validation_30pixels/processed_assayF_50_eprinomectin.RData")

proc_data_eprinomectin_50 <- processed_assayF_50_eprinomectin[[2]]

proc_img_dir_F <- '/projects/cellprofiler-nf/projects/20210319_assayF/Analysis-20220711/processed_images'


#highest dose concentration

#bleach 1, 50 pixels 
eprinomectin_assayF_bleach1_p043_H03_50 <- easyXpress::viewWell(proc_data_eprinomectin_50, proc_img_dir_F, "p043", "H03", boxplot = TRUE)
ggsave(eprinomectin_assayF_bleach1_p043_H03_50, file="eprinomectin_assayF_bleach1_p043_H03_50.jpg", width=20, height=10, dpi=500)

#bleach 2, 50 pixels 
eprinomectin_assayF_bleach2_p045_H04_50 <- easyXpress::viewWell(proc_data_eprinomectin_50, proc_img_dir_F, "p045", "H04", boxplot = TRUE)
ggsave(eprinomectin_assayF_bleach2_p045_H04_50, file="eprinomectin_assayF_bleach2_p045_H04_50.jpg", width=20, height=10, dpi=500)

#bleach 3, 50 pixels 
eprinomectin_assayF_bleach2_p047_H07_50 <- easyXpress::viewWell(proc_data_eprinomectin_50, proc_img_dir_F, "p047", "H07", boxplot = TRUE)
ggsave(eprinomectin_assayF_bleach2_p047_H07_50 , file="eprinomectin_assayF_bleach2_p047_H07_50 .jpg", width=20, height=10, dpi=500)


#30 pixels eprinomectin
edge_flaggedF_30_eprinomectin <- easyXpress::edgeFlag(model_selectedF_30_drugspecificMDHD_eprinomectin, radius=825, center_x=1024, center_y=1024)
raw_flaggedF_30_eprinomectin <- easyXpress::setFlags(edge_flaggedF_30_eprinomectin, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayF_30_eprinomectin <- easyXpress::process(raw_flaggedF_30_eprinomectin, Metadata_Plate, Metadata_Well)

save(processed_assayF_30_eprinomectin, file = "/data/validation_30pixels/processed_assayF_30_eprinomectin.RData")

proc_data_eprinomectin_30 <- processed_assayF_30_eprinomectin[[2]]

proc_img_dir_F <- '/projects/20210319_assayF/Analysis-20220711/processed_images'


#highest dose concentration

#bleach 1, 30 pixels 
eprinomectin_assayF_bleach1_p043_H03_30 <- easyXpress::viewWell(proc_data_eprinomectin_30, proc_img_dir_F, "p043", "H03", boxplot = TRUE)
ggsave(eprinomectin_assayF_bleach1_p043_H03_30, file="eprinomectin_assayF_bleach1_p043_H03_30.jpg", width=20, height=10, dpi=300)

#bleach 2, 30 pixels 
eprinomectin_assayF_bleach2_p045_H04_30 <- easyXpress::viewWell(proc_data_eprinomectin_30, proc_img_dir_F, "p045", "H04", boxplot = TRUE)
ggsave(eprinomectin_assayF_bleach2_p045_H04_30, file="eprinomectin_assayF_bleach2_p045_H04_30.jpg", width=20, height=10, dpi=300)

#bleach 3, 30 pixels 
eprinomectin_assayF_bleach2_p047_H07_30 <- easyXpress::viewWell(proc_data_eprinomectin_30, proc_img_dir_F, "p047", "H07", boxplot = TRUE)
ggsave(eprinomectin_assayF_bleach2_p047_H07_30 , file="eprinomectin_assayF_bleach2_p047_H07_30 .jpg", width=20, height=10, dpi=300)



#50 pixels monepantel_LY3348298
edge_flaggedF_50_monepantel_LY3348298 <- easyXpress::edgeFlag(model_selectedF_50_drugspecificMDHD_monepantel_LY3348298, radius=825, center_x=1024, center_y=1024)
raw_flaggedF_50_monepantel_LY3348298 <- easyXpress::setFlags(edge_flaggedF_50_monepantel_LY3348298, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayF_50_monepantel_LY3348298 <- easyXpress::process(raw_flaggedF_50_monepantel_LY3348298, Metadata_Plate, Metadata_Well)

save(processed_assayF_50_monepantel_LY3348298, file = "/data/validation_30pixels/processed_assayF_50_monepantel_LY3348298.RData")

proc_data_monepantel_LY3348298_50 <- processed_assayF_50_monepantel_LY3348298[[2]]

proc_img_dir_F <- '/projects/cellprofiler-nf/projects/20210319_assayF/Analysis-20220711/processed_images'


#highest dose concentration

#bleach 1, 50 pixels 
monepantel_LY3348298_assayF_bleach1_p013_H09_50 <- easyXpress::viewWell(proc_data_monepantel_LY3348298_50, proc_img_dir_F, "p013", "H09", boxplot = TRUE)
ggsave(monepantel_LY3348298_assayF_bleach1_p013_H09_50 , file="monepantel_LY3348298_assayF_bleach1_p013_H09_50 .jpg", width=20, height=10, dpi=500)

#bleach 2, 50 pixels 
monepantel_LY3348298_assayF_bleach2_p015_H04_50 <- easyXpress::viewWell(proc_data_monepantel_LY3348298_50, proc_img_dir_F, "p015", "H04", boxplot = TRUE)
ggsave(monepantel_LY3348298_assayF_bleach2_p015_H04_50, file="monepantel_LY3348298_assayF_bleach2_p015_H04_50.jpg", width=20, height=10, dpi=500)

#bleach 3, 50 pixels 
monepantel_LY3348298_assayF_bleach2_p017_H01_50 <- easyXpress::viewWell(proc_data_monepantel_LY3348298_50, proc_img_dir_F, "p017", "H01", boxplot = TRUE)
ggsave(monepantel_LY3348298_assayF_bleach2_p017_H01_50 , file="monepantel_LY3348298_assayF_bleach2_p017_H01_50.jpg", width=20, height=10, dpi=500)





#30 pixels monepantel_LY3348298
edge_flaggedF_30_monepantel_LY3348298 <- easyXpress::edgeFlag(model_selectedF_30_drugspecificMDHD_monepantel_LY3348298, radius=825, center_x=1024, center_y=1024)
raw_flaggedF_30_monepantel_LY3348298 <- easyXpress::setFlags(edge_flaggedF_30_monepantel_LY3348298, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayF_30_monepantel_LY3348298 <- easyXpress::process(raw_flaggedF_30_monepantel_LY3348298, Metadata_Plate, Metadata_Well)

save(processed_assayF_30_monepantel_LY3348298, file = "/data/validation_30pixels/processed_assayF_30_monepantel_LY3348298.RData")

proc_data_monepantel_LY3348298_30 <- processed_assayF_30_monepantel_LY3348298[[2]]

proc_img_dir_F <- '/projects/20210319_assayF/Analysis-20220711/processed_images'


#highest dose concentration

#bleach 1, 30 pixels 
monepantel_LY3348298_assayF_bleach1_p013_H09_30 <- easyXpress::viewWell(proc_data_monepantel_LY3348298_30, proc_img_dir_F, "p013", "H09", boxplot = TRUE)
ggsave(monepantel_LY3348298_assayF_bleach1_p013_H09_30 , file="monepantel_LY3348298_assayF_bleach1_p013_H09_30 .jpg", width=20, height=10, dpi=300)

#bleach 2, 30 pixels 
monepantel_LY3348298_assayF_bleach2_p015_H04_30 <- easyXpress::viewWell(proc_data_monepantel_LY3348298_30, proc_img_dir_F, "p015", "H04", boxplot = TRUE)
ggsave(monepantel_LY3348298_assayF_bleach2_p015_H04_30, file="monepantel_LY3348298_assayF_bleach2_p015_H04_30.jpg", width=20, height=10, dpi=300)

#bleach 3, 30 pixels 
monepantel_LY3348298_assayF_bleach2_p017_H01_30 <- easyXpress::viewWell(proc_data_monepantel_LY3348298_30, proc_img_dir_F, "p017", "H01", boxplot = TRUE)
ggsave(monepantel_LY3348298_assayF_bleach2_p017_H01_30 , file="monepantel_LY3348298_assayF_bleach2_p017_H01_30.jpg", width=20, height=10, dpi=300)





#############################
#                           #
#          Assay G          #
#                           #
#############################

datafileG <- "20210326_assayG_Analysis-20220709.RData"

rawG_50<- easyXpress::readXpress(filedir = dirsG, rdafile = datafileG, design = TRUE) %>%
  dplyr::filter(Worm_Length > 50) #remove objects smaller than 50 pixels = 165 um 

rawG_30 <- easyXpress::readXpress(filedir = dirsG, rdafile = datafileG, design = TRUE) %>%
  dplyr::filter(Worm_Length > 30) #remove objects smaller than 30 pixels = 100 um 

# Assay G drugs 
# abamectin
# albendazole
# doramectin
# emodepside
# fenbendazole
# ivermectin
# mebendazole
# milbemycin
# morantel
# moxidectin
# thiabendazole

# 30 pixels 
model_selectedG_30_drugspecificMDHD_abamectin <- easyXpress::modelSelection(rawG_30)%>% 
  dplyr::filter(drug == 'abamectin')

model_selectedG_30_drugspecificMDHD_albendazole <- easyXpress::modelSelection(rawG_30)%>% 
  dplyr::filter(drug == 'albendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedG_30_drugspecificMDHD_doramectin <- easyXpress::modelSelection(rawG_30)%>% 
  dplyr::filter(drug == 'doramectin')

model_selectedG_30_drugspecificMDHD_emodepside <- easyXpress::modelSelection(rawG_30)%>% 
  dplyr::filter(drug == 'emodepside')

model_selectedG_30_drugspecificMDHD_fenbendazole <- easyXpress::modelSelection(rawG_30)%>% 
  dplyr::filter(drug == 'fenbendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedG_30_drugspecificMDHD_ivermectin <- easyXpress::modelSelection(rawG_30)%>% 
  dplyr::filter(drug == 'ivermectin')

model_selectedG_30_drugspecificMDHD_mebendazole <- easyXpress::modelSelection(rawG_30)%>% 
  dplyr::filter(drug == 'mebendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedG_30_drugspecificMDHD_milbemycin <- easyXpress::modelSelection(rawG_30)%>% 
  dplyr::filter(drug == 'milbemycin')

model_selectedG_30_drugspecificMDHD_morantel <- easyXpress::modelSelection(rawG_30)%>% 
  dplyr::filter(drug == 'morantel' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedG_30_drugspecificMDHD_moxidectin <- easyXpress::modelSelection(rawG_30)%>% 
  dplyr::filter(drug == 'moxidectin')

model_selectedG_30_drugspecificMDHD_thiabendazole <- easyXpress::modelSelection(rawG_30)%>% 
  dplyr::filter(drug == 'thiabendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedG_30_drugspecificMDHD <- model_selectedG_30_drugspecificMDHD_abamectin %>% 
  full_join(model_selectedG_30_drugspecificMDHD_albendazole, by = NULL) %>% 
  full_join(model_selectedG_30_drugspecificMDHD_doramectin, by = NULL) %>% 
  full_join(model_selectedG_30_drugspecificMDHD_emodepside, by = NULL) %>% 
  full_join(model_selectedG_30_drugspecificMDHD_fenbendazole , by = NULL) %>% 
  full_join(model_selectedG_30_drugspecificMDHD_ivermectin, by = NULL) %>% 
  full_join(model_selectedG_30_drugspecificMDHD_mebendazole, by = NULL) %>% 
  full_join(model_selectedG_30_drugspecificMDHD_milbemycin, by = NULL) %>% 
  full_join(model_selectedG_30_drugspecificMDHD_morantel, by = NULL) %>% 
  full_join(model_selectedG_30_drugspecificMDHD_moxidectin, by = NULL) %>% 
  full_join(model_selectedG_30_drugspecificMDHD_thiabendazole, by = NULL) 

edge_flaggedG_30_drug_specific <- easyXpress::edgeFlag(model_selectedG_30_drugspecificMDHD, radius=825, center_x=1024, center_y=1024)
raw_flaggedG_30_drugspecific <- easyXpress::setFlags(edge_flaggedG_30_drug_specific, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayG_30_drugspecificMDHD <- easyXpress::process(raw_flaggedG_30_drugspecific , Metadata_Plate, Metadata_Well)

save(processed_assayG_30_drugspecificMDHD, file = "/data/processed_assayG_30_drugspecificMDHD.RData")


# 50 pixels 
model_selectedG_50_drugspecificMDHD_abamectin <- easyXpress::modelSelection(rawG_50)%>% 
  dplyr::filter(drug == 'abamectin')

model_selectedG_50_drugspecificMDHD_albendazole <- easyXpress::modelSelection(rawG_50)%>% 
  dplyr::filter(drug == 'albendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedG_50_drugspecificMDHD_doramectin <- easyXpress::modelSelection(rawG_50)%>% 
  dplyr::filter(drug == 'doramectin')

model_selectedG_50_drugspecificMDHD_emodepside <- easyXpress::modelSelection(rawG_50)%>% 
  dplyr::filter(drug == 'emodepside')

model_selectedG_50_drugspecificMDHD_fenbendazole <- easyXpress::modelSelection(rawG_50)%>% 
  dplyr::filter(drug == 'fenbendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedG_50_drugspecificMDHD_ivermectin <- easyXpress::modelSelection(rawG_50)%>% 
  dplyr::filter(drug == 'ivermectin')

model_selectedG_50_drugspecificMDHD_mebendazole <- easyXpress::modelSelection(rawG_50)%>% 
  dplyr::filter(drug == 'mebendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedG_50_drugspecificMDHD_milbemycin <- easyXpress::modelSelection(rawG_50)%>% 
  dplyr::filter(drug == 'milbemycin')

model_selectedG_50_drugspecificMDHD_morantel <- easyXpress::modelSelection(rawG_50)%>% 
  dplyr::filter(drug == 'morantel' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedG_50_drugspecificMDHD_moxidectin <- easyXpress::modelSelection(rawG_50)%>% 
  dplyr::filter(drug == 'moxidectin')

model_selectedG_50_drugspecificMDHD_thiabendazole <- easyXpress::modelSelection(rawG_50)%>% 
  dplyr::filter(drug == 'thiabendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedG_50_drugspecificMDHD <- model_selectedG_50_drugspecificMDHD_abamectin %>% 
  full_join(model_selectedG_50_drugspecificMDHD_albendazole, by = NULL) %>% 
  full_join(model_selectedG_50_drugspecificMDHD_doramectin, by = NULL) %>% 
  full_join(model_selectedG_50_drugspecificMDHD_emodepside, by = NULL) %>% 
  full_join(model_selectedG_50_drugspecificMDHD_fenbendazole , by = NULL) %>% 
  full_join(model_selectedG_50_drugspecificMDHD_ivermectin, by = NULL) %>% 
  full_join(model_selectedG_50_drugspecificMDHD_mebendazole, by = NULL) %>% 
  full_join(model_selectedG_50_drugspecificMDHD_milbemycin, by = NULL) %>% 
  full_join(model_selectedG_50_drugspecificMDHD_morantel, by = NULL) %>% 
  full_join(model_selectedG_50_drugspecificMDHD_moxidectin, by = NULL) %>% 
  full_join(model_selectedG_50_drugspecificMDHD_thiabendazole, by = NULL) 

edge_flaggedG_50_drug_specific <- easyXpress::edgeFlag(model_selectedG_50_drugspecificMDHD, radius=825, center_x=1050, center_y=1050)
raw_flaggedG_50_drugspecific <- easyXpress::setFlags(edge_flaggedG_50_drug_specific, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayG_50_drugspecificMDHD <- easyXpress::process(raw_flaggedG_50_drugspecific, Metadata_Plate, Metadata_Well)

save(processed_assayG_50_drugspecificMDHD, file = "/data/processed_assayG_50_drugspecificMDHD.RData")


#30 pixels Doramectin
edge_flaggedG_30_doramectin <- easyXpress::edgeFlag(model_selectedG_30_drugspecificMDHD_doramectin, radius=825, center_x=1024, center_y=1024)
raw_flaggedG_30_doramectin <- easyXpress::setFlags(edge_flaggedG_30_doramectin, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayG_30_doramectin <- easyXpress::process(raw_flaggedG_30_doramectin, Metadata_Plate, Metadata_Well)

save(processed_assayG_30_doramectin, file = "/data/validation_30pixels/processed_assayG_30_doramectin.RData")

proc_data_doramectin_30 <- processed_assayG_30_doramectin[[2]]

proc_img_dir_G <- '/projects/20210326_assayG/Analysis-20220709/processed_images'


#highest dose concentration

#bleach 1, 30 pixels 
doramectin_assayG_bleach1_p008_H09_30 <- easyXpress::viewWell(proc_data_doramectin_30, proc_img_dir_G, "p008", "H09", boxplot = TRUE)
ggsave(doramectin_assayG_bleach1_p008_H09_30, file="doramectin_assayG_bleach1_p008_H09_30.jpg", width=20, height=10, dpi=500)

#bleach 2, 30 pixels 
doramectin_assayG_bleach2_p009_H10_30 <- easyXpress::viewWell(proc_data_doramectin_30, proc_img_dir_G, "p009", "H10", boxplot = TRUE)
ggsave(doramectin_assayG_bleach2_p009_H10_30 , file="doramectin_assayG_bleach2_p009_H10_30.jpg", width=20, height=10, dpi=500)

#bleach 3, 30 pixels 
doramectin_assayG_bleach3_p011_H07_30 <- easyXpress::viewWell(proc_data_doramectin_30, proc_img_dir_G, "p011", "H07", boxplot = TRUE)
ggsave(doramectin_assayG_bleach3_p011_H07_30  , file="doramectin_assayG_bleach3_p011_H07_30.jpg", width=20, height=10, dpi=500)



#50 pixels Doramectin
edge_flaggedG_50_doramectin <- easyXpress::edgeFlag(model_selectedG_50_drugspecificMDHD_doramectin, radius=825, center_x=1024, center_y=1024)
raw_flaggedG_50_doramectin <- easyXpress::setFlags(edge_flaggedG_50_doramectin, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayG_50_doramectin <- easyXpress::process(raw_flaggedG_50_doramectin, Metadata_Plate, Metadata_Well)

save(processed_assayG_50_doramectin, file = "/data/validation_30pixels/processed_assayG_50_doramectin.RData")

proc_data_doramectin_50 <- processed_assayG_50_doramectin[[2]]

proc_img_dir_G <- '/projects/20210326_assayG/Analysis-20220709/processed_images'


#highest dose concentration

#bleach 1, 50 pixels 
doramectin_assayG_bleach1_p008_H09_50 <- easyXpress::viewWell(proc_data_doramectin_50, proc_img_dir_G, "p008", "H09", boxplot = TRUE)
ggsave(doramectin_assayG_bleach1_p008_H09_50, file="doramectin_assayG_bleach1_p008_H09_50.jpg", width=20, height=10, dpi=500)

#bleach 2, 50 pixels 
doramectin_assayG_bleach2_p009_H10_50 <- easyXpress::viewWell(proc_data_doramectin_50, proc_img_dir_G, "p009", "H10", boxplot = TRUE)
ggsave(doramectin_assayG_bleach2_p009_H10_50 , file="doramectin_assayG_bleach2_p009_H10_50.jpg", width=20, height=10, dpi=500)

#bleach 3, 50 pixels 
doramectin_assayG_bleach3_p011_H07_50 <- easyXpress::viewWell(proc_data_doramectin_50, proc_img_dir_G, "p011", "H07", boxplot = TRUE)
ggsave(doramectin_assayG_bleach3_p011_H07_50  , file="doramectin_assayG_bleach3_p011_H07_50.jpg", width=20, height=10, dpi=500)



#30 pixels Ivermectin
edge_flaggedG_30_ivermectin <- easyXpress::edgeFlag(model_selectedG_30_drugspecificMDHD_ivermectin, radius=825, center_x=1024, center_y=1024)
raw_flaggedG_30_ivermectin <- easyXpress::setFlags(edge_flaggedG_30_ivermectin, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayG_30_ivermectin <- easyXpress::process(raw_flaggedG_30_ivermectin, Metadata_Plate, Metadata_Well)

save(processed_assayG_30_ivermectin, file = "/data/validation_30pixels/processed_assayG_30_ivermectin.RData")

proc_data_ivermectin_30 <- processed_assayG_30_ivermectin[[2]]

proc_img_dir_G <- '/projects/20210326_assayG/Analysis-20220709/processed_images'


#highest dose concentration

#bleach 1, 30 pixels 
ivermectin_assayG_bleach1_p001_H09_30 <- easyXpress::viewWell(proc_data_ivermectin_30, proc_img_dir_G, "p001", "H09", boxplot = TRUE)
ggsave(ivermectin_assayG_bleach1_p001_H09_30, file="ivermectin_assayG_bleach1_p001_H09_30.jpg", width=20, height=10, dpi=500)

#bleach 2, 30 pixels 
ivermectin_assayG_bleach2_p003_H10_30 <- easyXpress::viewWell(proc_data_ivermectin_30, proc_img_dir_G, "p003", "H10", boxplot = TRUE)
ggsave(ivermectin_assayG_bleach2_p003_H10_30, file="ivermectin_assayG_bleach2_p003_H10_30.jpg", width=20, height=10, dpi=500)

#bleach 3, 30 pixels 
ivermectin_assayG_bleach3_p006_H07_30 <- easyXpress::viewWell(proc_data_ivermectin_30, proc_img_dir_G, "p006", "H07", boxplot = TRUE)
ggsave(ivermectin_assayG_bleach3_p006_H07_30, file="ivermectin_assayG_bleach3_p006_H07_30.jpg", width=20, height=10, dpi=500)



#50 pixels Ivermectin
edge_flaggedG_50_ivermectin <- easyXpress::edgeFlag(model_selectedG_50_drugspecificMDHD_ivermectin, radius=825, center_x=1024, center_y=1024)
raw_flaggedG_50_ivermectin <- easyXpress::setFlags(edge_flaggedG_50_ivermectin, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayG_50_ivermectin <- easyXpress::process(raw_flaggedG_50_ivermectin, Metadata_Plate, Metadata_Well)

save(processed_assayG_50_ivermectin, file = "/data/validation_30pixels/processed_assayG_50_ivermectin.RData")

proc_data_ivermectin_50 <- processed_assayG_50_ivermectin[[2]]

proc_img_dir_G <- '/projects/20210326_assayG/Analysis-20220709/processed_images'


#highest dose concentration

#bleach 1, 50 pixels 
ivermectin_assayG_bleach1_p001_H09_50 <- easyXpress::viewWell(proc_data_ivermectin_50, proc_img_dir_G, "p001", "H09", boxplot = TRUE)
ggsave(ivermectin_assayG_bleach1_p001_H09_50, file="ivermectin_assayG_bleach1_p001_H09_50.jpg", width=20, height=10, dpi=500)

#bleach 2, 50 pixels 
ivermectin_assayG_bleach2_p003_H10_50 <- easyXpress::viewWell(proc_data_ivermectin_50, proc_img_dir_G, "p003", "H10", boxplot = TRUE)
ggsave(ivermectin_assayG_bleach2_p003_H10_50, file="ivermectin_assayG_bleach2_p003_H10_50.jpg", width=20, height=10, dpi=500)

#bleach 3, 50 pixels 
ivermectin_assayG_bleach3_p006_H07_50 <- easyXpress::viewWell(proc_data_ivermectin_50, proc_img_dir_G, "p006", "H07", boxplot = TRUE)
ggsave(ivermectin_assayG_bleach3_p006_H07_50, file="ivermectin_assayG_bleach3_p006_H07_50.jpg", width=20, height=10, dpi=500)



#50 pixels milbemycin
edge_flaggedG_50_milbemycin <- easyXpress::edgeFlag(model_selectedG_50_drugspecificMDHD_milbemycin, radius=825, center_x=1024, center_y=1024)
raw_flaggedG_50_milbemycin <- easyXpress::setFlags(edge_flaggedG_50_milbemycin, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayG_50_milbemycin <- easyXpress::process(raw_flaggedG_50_milbemycin, Metadata_Plate, Metadata_Well)

save(processed_assayG_50_milbemycin, file = "/data/validation_30pixels/processed_assayG_50_milbemycin.RData")

proc_data_milbemycin_50 <- processed_assayG_50_milbemycin[[2]]

proc_img_dir_G <- '/projects/20210326_assayG/Analysis-20220709/processed_images'


#highest dose concentration

#bleach 1, 50 pixels 
milbemycin_assayG_bleach1_p020_H03_50 <- easyXpress::viewWell(proc_data_milbemycin_50, proc_img_dir_G, "p020", "H03", boxplot = TRUE)
ggsave(milbemycin_assayG_bleach1_p020_H03_50, file="milbemycin_assayG_bleach1_p020_H03_50.jpg", width=20, height=10, dpi=500)

#bleach 2, 50 pixels 
milbemycin_assayG_bleach2_p021_H04_50 <- easyXpress::viewWell(proc_data_milbemycin_50, proc_img_dir_G, "p021", "H04", boxplot = TRUE)
ggsave(milbemycin_assayG_bleach2_p021_H04_50 , file="milbemycin_assayG_bleach2_p021_H04_50.jpg", width=20, height=10, dpi=500)

#bleach 3, 50 pixels 
milbemycin_assayG_bleach3_p024_H07_50 <- easyXpress::viewWell(proc_data_milbemycin_50, proc_img_dir_G, "p024", "H07", boxplot = TRUE)
ggsave(milbemycin_assayG_bleach3_p024_H07_50, file="milbemycin_assayG_bleach3_p024_H07_50.jpg", width=20, height=10, dpi=500)


#30 pixels milbemycin
edge_flaggedG_30_milbemycin <- easyXpress::edgeFlag(model_selectedG_30_drugspecificMDHD_milbemycin, radius=825, center_x=1024, center_y=1024)
raw_flaggedG_30_milbemycin <- easyXpress::setFlags(edge_flaggedG_30_milbemycin, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayG_30_milbemycin <- easyXpress::process(raw_flaggedG_30_milbemycin, Metadata_Plate, Metadata_Well)

save(processed_assayG_30_milbemycin, file = "/data/validation_30pixels/processed_assayG_30_milbemycin.RData")

proc_data_milbemycin_30 <- processed_assayG_30_milbemycin[[2]]

proc_img_dir_G <- '/projects/20210326_assayG/Analysis-20220709/processed_images'


#highest dose concentration

#bleach 1, 30 pixels 
milbemycin_assayG_bleach1_p020_H03_30 <- easyXpress::viewWell(proc_data_milbemycin_30, proc_img_dir_G, "p020", "H03", boxplot = TRUE)
ggsave(milbemycin_assayG_bleach1_p020_H03_30, file="milbemycin_assayG_bleach1_p020_H03_30.jpg", width=20, height=10, dpi=300)

#bleach 2, 30 pixels 
milbemycin_assayG_bleach2_p021_H04_30 <- easyXpress::viewWell(proc_data_milbemycin_30, proc_img_dir_G, "p021", "H04", boxplot = TRUE)
ggsave(milbemycin_assayG_bleach2_p021_H04_30 , file="milbemycin_assayG_bleach2_p021_H04_30.jpg", width=20, height=10, dpi=300)

#bleach 3, 30 pixels 
milbemycin_assayG_bleach3_p024_H07_30 <- easyXpress::viewWell(proc_data_milbemycin_30, proc_img_dir_G, "p024", "H07", boxplot = TRUE)
ggsave(milbemycin_assayG_bleach3_p024_H07_30, file="milbemycin_assayG_bleach3_p024_H07_30.jpg", width=20, height=10, dpi=300)


#############################
#                           #
#          Assay H          #
#                           #
#############################

#remove all derquantel, niridazole

#Reading in the data
# Define experimental directory and file name
dirsH <- "/projects/20210402_assayH/Analysis-20220710" #data and design match
datafileH <- "20210402_assayH_Analysis-20220710.RData"

rawH_50 <- easyXpress::readXpress(filedir = dirsH, rdafile = datafileH, design = TRUE) %>%
  dplyr::filter(Worm_Length > 50) #remove objects smaller than 50 pixels = 165 um 

rawH_30 <- easyXpress::readXpress(filedir = dirsH, rdafile = datafileH, design = TRUE) %>%
  dplyr::filter(Worm_Length > 30) #remove objects smaller than 30 pixels = 100 um 

# ASSAY H DRUGS
# albendazole (REMOVE MDHD)
# derquantel (KEEP MDHD) (RMEOVE ENTIRELY FROM ASSAY)
# doramectin (KEEP MDHD)
# fenbendazole (REMOVE MDHD)
# ivermectin (KEEP MDHD)
# mebendazole (REMOVE MDHD)
# morantel (REMOVE MDHD)
# niridazole (REMOVE MDHD) (REMOVE ENTIRELY FROM ASSAY )
# oxamniquine (REMOVE MDHD)
# pyrantel_cytrate (REMOVE MDHD)
# thiabendazole (REMOVE MDHD)

model_selectedH_30_drugspecificMDHD_albendazole <- easyXpress::modelSelection(rawH_30)%>% 
  dplyr::filter(drug == 'albendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedH_30_drugspecificMDHD_niridazole <- easyXpress::modelSelection(rawH_30)%>% 
  dplyr::filter(drug == 'niridazole')

model_selectedH_30_drugspecificMDHD_doramectin <- easyXpress::modelSelection(rawH_30)%>% 
  dplyr::filter(drug == 'doramectin')

model_selectedH_30_drugspecificMDHD_fenbendazole  <- easyXpress::modelSelection(rawH_30)%>% 
  dplyr::filter(drug == 'fenbendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedH_30_drugspecificMDHD_ivermectin <- easyXpress::modelSelection(rawH_30)%>% 
  dplyr::filter(drug == 'ivermectin')

model_selectedH_30_drugspecificMDHD_mebendazole  <- easyXpress::modelSelection(rawH_30)%>% 
  dplyr::filter(drug == 'mebendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedH_30_drugspecificMDHD_morantel <- easyXpress::modelSelection(rawH_30)%>% 
  dplyr::filter(drug == 'morantel')

model_selectedH_30_drugspecificMDHD_oxamniquine  <- easyXpress::modelSelection(rawH_30)%>% 
  dplyr::filter(drug == 'oxamniquine' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedH_30_drugspecificMDHD_pyrantel_cytrate  <- easyXpress::modelSelection(rawH_30)%>% 
  dplyr::filter(drug == 'pyrantel_cytrate' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedH_30_drugspecificMDHD_thiabendazole  <- easyXpress::modelSelection(rawH_30)%>% 
  dplyr::filter(drug == 'thiabendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedH_30_assayspecificMDHD <- model_selectedH_30_drugspecificMDHD_albendazole %>% 
  full_join(model_selectedH_30_drugspecificMDHD_doramectin, by = NULL) %>% 
  full_join(model_selectedH_30_drugspecificMDHD_fenbendazole, by = NULL) %>% 
  full_join(model_selectedH_30_drugspecificMDHD_ivermectin, by = NULL) %>% 
  full_join(model_selectedH_30_drugspecificMDHD_mebendazole, by = NULL) %>% 
  full_join(model_selectedH_30_drugspecificMDHD_morantel, by = NULL) %>% 
  full_join(model_selectedH_30_drugspecificMDHD_oxamniquine, by = NULL) %>% 
  full_join(model_selectedH_30_drugspecificMDHD_pyrantel_cytrate, by = NULL) %>% 
  full_join(model_selectedH_30_drugspecificMDHD_thiabendazole, by = NULL) 
  

edge_flaggedH_30_drugspecific <- easyXpress::edgeFlag(model_selectedH_30_assayspecificMDHD, radius=825, center_x=1024, center_y=1024)
raw_flaggedH_30_drugspecific <- easyXpress::setFlags(edge_flaggedH_30_drugspecific, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayH_30_drugspecificMDHD <- easyXpress::process(raw_flaggedH_30_drugspecific, Metadata_Plate, Metadata_Well)

save(processed_assayH_30_drugspecificMDHD, file = "/data/processed_assayH_30_drugspecificMDHD.RData")

edge_flaggedH_30_drugspecific_niridazole <- easyXpress::edgeFlag(model_selectedH_30_drugspecificMDHD_niridazole, radius=825, center_x=1024, center_y=1024)
raw_flaggedH_30_drugspecific_niridazole <- easyXpress::setFlags(edge_flaggedH_30_drugspecific_niridazole, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayH_30_drugspecificMDHD_niridazole <- easyXpress::process(raw_flaggedH_30_drugspecific_niridazole, Metadata_Plate, Metadata_Well)

save(processed_assayH_30_drugspecificMDHD_niridazole, file = "/data/processed_assayH_30_drugspecificMDHD_niridazole.RData")



#50 pixels
model_selectedH_50_drugspecificMDHD_albendazole <- easyXpress::modelSelection(rawH_50)%>% 
  dplyr::filter(drug == 'albendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedH_50_drugspecificMDHD_doramectin <- easyXpress::modelSelection(rawH_50)%>% 
  dplyr::filter(drug == 'doramectin')

model_selectedH_50_drugspecificMDHD_fenbendazole  <- easyXpress::modelSelection(rawH_50)%>% 
  dplyr::filter(drug == 'fenbendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedH_50_drugspecificMDHD_ivermectin <- easyXpress::modelSelection(rawH_50)%>% 
  dplyr::filter(drug == 'ivermectin')

model_selectedH_50_drugspecificMDHD_mebendazole  <- easyXpress::modelSelection(rawH_50)%>% 
  dplyr::filter(drug == 'mebendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedH_50_drugspecificMDHD_morantel <- easyXpress::modelSelection(rawH_50)%>% 
  dplyr::filter(drug == 'morantel')

model_selectedH_50_drugspecificMDHD_oxamniquine  <- easyXpress::modelSelection(rawH_50)%>% 
  dplyr::filter(drug == 'oxamniquine' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedH_50_drugspecificMDHD_pyrantel_cytrate  <- easyXpress::modelSelection(rawH_50)%>% 
  dplyr::filter(drug == 'pyrantel_cytrate' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedH_50_drugspecificMDHD_thiabendazole  <- easyXpress::modelSelection(rawH_50)%>% 
  dplyr::filter(drug == 'thiabendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedH_50_assayspecificMDHD <- model_selectedH_50_drugspecificMDHD_albendazole %>% 
  full_join(model_selectedH_50_drugspecificMDHD_doramectin, by = NULL) %>% 
  full_join(model_selectedH_50_drugspecificMDHD_fenbendazole, by = NULL) %>% 
  full_join(model_selectedH_50_drugspecificMDHD_ivermectin, by = NULL) %>% 
  full_join(model_selectedH_50_drugspecificMDHD_mebendazole, by = NULL) %>% 
  full_join(model_selectedH_50_drugspecificMDHD_morantel, by = NULL) %>% 
  full_join(model_selectedH_50_drugspecificMDHD_oxamniquine, by = NULL) %>% 
  full_join(model_selectedH_50_drugspecificMDHD_pyrantel_cytrate, by = NULL) %>% 
  full_join(model_selectedH_50_drugspecificMDHD_thiabendazole, by = NULL) 


edge_flaggedH_50_drugspecific <- easyXpress::edgeFlag(model_selectedH_50_assayspecificMDHD, radius=825, center_x=1050, center_y=1050)
raw_flaggedH_50_drugspecific <- easyXpress::setFlags(edge_flaggedH_50_drugspecific, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayH_50_drugspecificMDHD <- easyXpress::process(raw_flaggedH_50_drugspecific, Metadata_Plate, Metadata_Well)

save(processed_assayH_50_drugspecificMDHD, file = "/data/processed_assayH_50_drugspecificMDHD.RData")


#############################
#                           #
#          Assay I          # 
#                           #
#############################

#Reading in the data
# Define experimental directory and file name
datafileI <- "20210409_assayI_Analysis-20220709.RData"

rawI_50 <- easyXpress::readXpress(filedir = dirsI, rdafile = datafileI, design = TRUE) %>%
  dplyr::filter(Worm_Length > 50)

rawI_30 <- easyXpress::readXpress(filedir = dirsI, rdafile = datafileI, design = TRUE) %>%
  dplyr::filter(Worm_Length > 30)

# ASSAY I DRUGS
# dec
# doramectin
# eprinomectin
# ivermectin
# monepantel_LY33414916
# monepantel_LY3348298
# piperazine
# pyrantel_cytrate
# selamectin

#30 pixels
model_selectedI_30_drugspecificMDHD_dec  <- easyXpress::modelSelection(rawI_30)%>% 
  dplyr::filter(drug == 'dec' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedI_30_drugspecificMDHD_doramectin  <- easyXpress::modelSelection(rawI_30)%>% 
  dplyr::filter(drug == 'doramectin')

# model_selectedI_assayspecific_eprinomectin  <- easyXpress::modelSelection(rawI_30)%>%
#   dplyr::filter(drug == 'eprinomectin') #assay needs to be removed outlier 

model_selectedI_30_drugspecificMDHD_ivermectin  <- easyXpress::modelSelection(rawI_30)%>% 
  dplyr::filter(drug == 'ivermectin')

model_selectedI_30_drugspecificMDHD_monepantel_LY33414916  <- easyXpress::modelSelection(rawI_30)%>% 
  dplyr::filter(drug == 'monepantel_LY33414916')

model_selectedI_30_drugspecificMDHD_monepantel_LY3348298  <- easyXpress::modelSelection(rawI_30)%>% 
  dplyr::filter(drug == 'monepantel_LY3348298')

# model_selectedI_assayspecific_piperazine  <- easyXpress::modelSelection(rawI_30)%>% 
#   dplyr::filter(drug == 'piperazine' & model_select != 'MDHD_NonOverlappingWorms.model.outputs') #removed, assay outlier 

model_selectedI_30_drugspecificMDHD_pyrantel_cytrate <- easyXpress::modelSelection(rawI_30)%>% 
  dplyr::filter(drug == 'pyrantel_cytrate' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedI_30_drugspecificMDHD_selamectin  <- easyXpress::modelSelection(rawI_30)%>% 
  dplyr::filter(drug == 'selamectin')

model_selectedI_30_assayspecificMDHD <- model_selectedI_30_drugspecificMDHD_dec %>% 
  full_join(model_selectedI_30_drugspecificMDHD_doramectin, by = NULL) %>% 
  # full_join(model_selectedI_assayspecific_eprinomectin, by = NULL) %>% 
  full_join(model_selectedI_30_drugspecificMDHD_ivermectin, by = NULL) %>% 
  full_join(model_selectedI_30_drugspecificMDHD_monepantel_LY33414916, by = NULL) %>% 
  full_join(model_selectedI_30_drugspecificMDHD_monepantel_LY3348298, by = NULL) %>% 
  # full_join(model_selectedI_assayspecific_piperazine, by = NULL) %>% 
  full_join(model_selectedI_30_drugspecificMDHD_pyrantel_cytrate, by = NULL) %>% 
  full_join(model_selectedI_30_drugspecificMDHD_selamectin, by = NULL) 

edge_flaggedI_30_drugspecific <- easyXpress::edgeFlag(model_selectedI_30_assayspecificMDHD, radius=825, center_x=1024, center_y=1024)
raw_flaggedI_30_drugspecific <- easyXpress::setFlags(edge_flaggedI_30_drugspecific, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayI_30_drugspecificMDHD <- easyXpress::process(raw_flaggedI_30_drugspecific, Metadata_Plate, Metadata_Well)

save(processed_assayI_30_drugspecificMDHD, file = "/data/processed_assayI_30_drugspecificMDHD.RData")


#50 pixels 
model_selectedI_50_drugspecificMDHD_dec  <- easyXpress::modelSelection(rawI_50)%>% 
  dplyr::filter(drug == 'dec' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedI_50_drugspecificMDHD_doramectin  <- easyXpress::modelSelection(rawI_50)%>% 
  dplyr::filter(drug == 'doramectin')

# model_selectedI_assayspecific_eprinomectin  <- easyXpress::modelSelection(rawI_50)%>% 
#   dplyr::filter(drug == 'eprinomectin')

model_selectedI_50_drugspecificMDHD_ivermectin  <- easyXpress::modelSelection(rawI_50)%>% 
  dplyr::filter(drug == 'ivermectin')

model_selectedI_50_drugspecificMDHD_monepantel_LY33414916  <- easyXpress::modelSelection(rawI_50)%>% 
  dplyr::filter(drug == 'monepantel_LY33414916')

model_selectedI_50_drugspecificMDHD_monepantel_LY3348298  <- easyXpress::modelSelection(rawI_50)%>% 
  dplyr::filter(drug == 'monepantel_LY3348298')

# model_selectedI_assayspecific_piperazine  <- easyXpress::modelSelection(rawI_50)%>% 
#   dplyr::filter(drug == 'piperazine' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedI_50_drugspecificMDHD_pyrantel_cytrate <- easyXpress::modelSelection(rawI_50)%>% 
  dplyr::filter(drug == 'pyrantel_cytrate' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedI_50_drugspecificMDHD_selamectin  <- easyXpress::modelSelection(rawI_50)%>% 
  dplyr::filter(drug == 'selamectin')

model_selectedI_50_assayspecificMDHD <- model_selectedI_50_drugspecificMDHD_dec %>% 
  full_join(model_selectedI_50_drugspecificMDHD_doramectin, by = NULL) %>% 
  # full_join(model_selectedI_assayspecific_eprinomectin, by = NULL) %>% 
  full_join(model_selectedI_50_drugspecificMDHD_ivermectin, by = NULL) %>% 
  full_join(model_selectedI_50_drugspecificMDHD_monepantel_LY33414916, by = NULL) %>% 
  full_join(model_selectedI_50_drugspecificMDHD_monepantel_LY3348298, by = NULL) %>% 
  # full_join(model_selectedI_assayspecific_piperazine, by = NULL) %>% 
  full_join(model_selectedI_50_drugspecificMDHD_pyrantel_cytrate, by = NULL) %>% 
  full_join(model_selectedI_50_drugspecificMDHD_selamectin, by = NULL) 

edge_flaggedI_50_drugspecific <- easyXpress::edgeFlag(model_selectedI_50_assayspecificMDHD, radius=825, center_x=1050, center_y=1050)
raw_flaggedI_50_drugspecific <- easyXpress::setFlags(edge_flaggedI_50_drugspecific, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayI_50_drugspecificMDHD <- easyXpress::process(raw_flaggedI_50_drugspecific, Metadata_Plate, Metadata_Well)

save(processed_assayI_50_drugspecificMDHD, file = "/data/processed_assayI_50_drugspecificMDHD.RData")


#############################
#                           #
#          Assay J          #
#                           #
#############################

datafileJ <- "20210430_assayJ_Analysis-20220708.RData"

rawJ_50 <- easyXpress::readXpress(filedir = dirsJ, rdafile = datafileJ, design = TRUE) %>%
  dplyr::filter(Worm_Length > 50)

rawJ_30 <- easyXpress::readXpress(filedir = dirsJ, rdafile = datafileJ, design = TRUE) %>%
  dplyr::filter(Worm_Length > 30)

# ASSAY J DRUGS
# abamectin (KEEP MDHD)
# albendazole (REMOVE MDHD)
# doramectin (KEEP MDHD)
# emodepside (KEEP MDHD)
# eprinomectin (KEEP MDHD)
# fenbendazole (REMOVE MDHD)
# ivermectin (KEEP MDHD)
# levamisole (KEEP MDHD)
# milbemycin (KEEP MDHD)
# piperazine (REMOVE MDHD)
# selamectin (KEEP MDHD)

#30 pixels 
model_selectedJ_30_drugspecificMDHD_abamectin  <- easyXpress::modelSelection(rawJ_30)%>% 
  dplyr::filter(drug == 'abamectin')

model_selectedJ_30_drugspecificMDHD_albendazole <- easyXpress::modelSelection(rawJ_30)%>% 
  dplyr::filter(drug == 'albendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedJ_30_drugspecificMDHD_doramectin  <- easyXpress::modelSelection(rawJ_30)%>% 
  dplyr::filter(drug == 'doramectin')

model_selectedJ_30_drugspecificMDHD_emodepside  <- easyXpress::modelSelection(rawJ_30)%>% 
  dplyr::filter(drug == 'emodepside')

model_selectedJ_30_drugspecificMDHD_eprinomectin  <- easyXpress::modelSelection(rawJ_30)%>% 
  dplyr::filter(drug == 'eprinomectin')

model_selectedJ_30_drugspecificMDHD_fenbendazole <- easyXpress::modelSelection(rawJ_30)%>% 
  dplyr::filter(drug == 'fenbendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedJ_30_drugspecificMDHD_ivermectin  <- easyXpress::modelSelection(rawJ_30)%>% 
  dplyr::filter(drug == 'ivermectin')

model_selectedJ_30_drugspecificMDHD_levamisole  <- easyXpress::modelSelection(rawJ_30)%>% 
  dplyr::filter(drug == 'levamisole')

model_selectedJ_30_drugspecificMDHD_milbemycin  <- easyXpress::modelSelection(rawJ_30)%>% 
  dplyr::filter(drug == 'milbemycin')

model_selectedJ_30_drugspecificMDHD_piperazine <- easyXpress::modelSelection(rawJ_30)%>% 
  dplyr::filter(drug == 'piperazine' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedJ_30_drugspecificMDHD_selamectin  <- easyXpress::modelSelection(rawJ_30)%>% 
  dplyr::filter(drug == 'selamectin')

model_selectedJ_30_assayspecificMDHD <- model_selectedJ_30_drugspecificMDHD_abamectin  %>% 
  full_join(model_selectedJ_30_drugspecificMDHD_doramectin , by = NULL) %>% 
  full_join(model_selectedJ_30_drugspecificMDHD_albendazole , by = NULL) %>% 
  full_join(model_selectedJ_30_drugspecificMDHD_emodepside , by = NULL) %>% 
  full_join(model_selectedJ_30_drugspecificMDHD_eprinomectin, by = NULL) %>% 
  full_join(model_selectedJ_30_drugspecificMDHD_fenbendazole, by = NULL) %>% 
  full_join(model_selectedJ_30_drugspecificMDHD_ivermectin, by = NULL) %>% 
  full_join(model_selectedJ_30_drugspecificMDHD_levamisole, by = NULL) %>% 
  full_join(model_selectedJ_30_drugspecificMDHD_milbemycin, by = NULL) %>% 
  full_join(model_selectedJ_30_drugspecificMDHD_piperazine, by = NULL) %>% 
  full_join(model_selectedJ_30_drugspecificMDHD_selamectin, by = NULL)

edge_flaggedJ_30_drugspecific <- easyXpress::edgeFlag(model_selectedJ_30_assayspecificMDHD, radius=825, center_x=1024, center_y=1024)
raw_flaggedJ_30_drugpecific <- easyXpress::setFlags(edge_flaggedJ_30_drugspecific, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayJ_30_drugspecificMDHD <- easyXpress::process(raw_flaggedJ_30_drugpecific, Metadata_Plate, Metadata_Well)

save(processed_assayJ_30_drugspecificMDHD, file = "/data/processed_assayJ_30_drugspecificMDHD.RData")


#50 pixels
model_selectedJ_50_drugspecificMDHD_abamectin  <- easyXpress::modelSelection(rawJ_50)%>% 
  dplyr::filter(drug == 'abamectin')

model_selectedJ_50_drugspecificMDHD_albendazole <- easyXpress::modelSelection(rawJ_50)%>% 
  dplyr::filter(drug == 'albendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedJ_50_drugspecificMDHD_doramectin  <- easyXpress::modelSelection(rawJ_50)%>% 
  dplyr::filter(drug == 'doramectin')

model_selectedJ_50_drugspecificMDHD_emodepside  <- easyXpress::modelSelection(rawJ_50)%>% 
  dplyr::filter(drug == 'emodepside')

model_selectedJ_50_drugspecificMDHD_eprinomectin  <- easyXpress::modelSelection(rawJ_50)%>% 
  dplyr::filter(drug == 'eprinomectin')

model_selectedJ_50_drugspecificMDHD_fenbendazole <- easyXpress::modelSelection(rawJ_50)%>% 
  dplyr::filter(drug == 'fenbendazole' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedJ_50_drugspecificMDHD_ivermectin  <- easyXpress::modelSelection(rawJ_50)%>% 
  dplyr::filter(drug == 'ivermectin')

model_selectedJ_50_drugspecificMDHD_levamisole  <- easyXpress::modelSelection(rawJ_50)%>% 
  dplyr::filter(drug == 'levamisole')

model_selectedJ_50_drugspecificMDHD_milbemycin  <- easyXpress::modelSelection(rawJ_50)%>% 
  dplyr::filter(drug == 'milbemycin')

model_selectedJ_50_drugspecificMDHD_piperazine <- easyXpress::modelSelection(rawJ_50)%>% 
  dplyr::filter(drug == 'piperazine' & model_select != 'MDHD_NonOverlappingWorms.model.outputs')

model_selectedJ_50_drugspecificMDHD_selamectin  <- easyXpress::modelSelection(rawJ_50)%>% 
  dplyr::filter(drug == 'selamectin')

model_selectedJ_50_assayspecificMDHD <- model_selectedJ_50_drugspecificMDHD_abamectin  %>% 
  full_join(model_selectedJ_50_drugspecificMDHD_doramectin , by = NULL) %>% 
  full_join(model_selectedJ_50_drugspecificMDHD_albendazole , by = NULL) %>% 
  full_join(model_selectedJ_50_drugspecificMDHD_emodepside , by = NULL) %>% 
  full_join(model_selectedJ_50_drugspecificMDHD_eprinomectin, by = NULL) %>% 
  full_join(model_selectedJ_50_drugspecificMDHD_fenbendazole, by = NULL) %>% 
  full_join(model_selectedJ_50_drugspecificMDHD_ivermectin, by = NULL) %>% 
  full_join(model_selectedJ_50_drugspecificMDHD_levamisole, by = NULL) %>% 
  full_join(model_selectedJ_50_drugspecificMDHD_milbemycin, by = NULL) %>% 
  full_join(model_selectedJ_50_drugspecificMDHD_piperazine, by = NULL) %>% 
  full_join(model_selectedJ_50_drugspecificMDHD_selamectin, by = NULL)

edge_flaggedJ_50_drugspecific <- easyXpress::edgeFlag(model_selectedJ_50_assayspecificMDHD, radius=825, center_x=1050, center_y=1050)
raw_flaggedJ_50_drugpecific <- easyXpress::setFlags(edge_flaggedJ_50_drugspecific, cluster_flag = TRUE, well_edge_flag = TRUE)
processed_assayJ_50_drugspecificMDHD <- easyXpress::process(raw_flaggedJ_50_drugpecific, Metadata_Plate, Metadata_Well)

save(processed_assayJ_50_drugspecificMDHD, file = "/data/processed_assayJ_50_drugspecificMDHD.RData")



#combine all assays into one .Rdata file. 50 PIXELS
save(processed_assayA_drugMDHD_50, processed_assayB_50_drugspecificMDHD, processed_assayC_50_drugspecificMDHD,
     processed_assayD_50_drugspecificMDHD, processed_assayE_50_drugspecificMDHD, processed_assayF_50_drugspecificMDHD,
     processed_assayG_50_drugspecificMDHD, processed_assayH_50_drugspecificMDHD, processed_assayI_50_drugspecificMDHD, 
     processed_assayJ_50_drugspecificMDHD, file = "/data/all_processed_assays_drugspecificMDHD_50.RData")


#combine all assays into one .Rdata file. 30 Pixels
save(processed_assayA_drugMDHD_30, processed_assayB_30_drugspecificMDHD, processed_assayC_30_drugspecificMDHD,
     processed_assayD_30_drugspecificMDHD, processed_assayE_30_drugspecificMDHD, processed_assayF_30_drugspecificMDHD,
     processed_assayG_30_drugspecificMDHD, processed_assayH_30_drugspecificMDHD, processed_assayI_30_drugspecificMDHD,
     processed_assayJ_30_drugspecificMDHD, file = "/data/all_processed_assays_drugspecificMDHD_30.RData")

