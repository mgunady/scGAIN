library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(rafalib)

library(gridExtra)
library(Rtsne)

library(splatter)
library(scater)
library(edgeR)
library("Matrix")
library(tictoc)

library(scImpute)
library(SCRABBLE)
library(DrImpute)
library(VIPER)
library(SAVER)

source("run_tools.R")

dropout_indexs = c(5)
seed_values = c(20000)

for(drop_index in dropout_indexs){
  
  seed_value <- seed_values[1]

  p <- plot_tsnes(drop_index,
                  seed_value,
                  50,100#, c("gain")
		)
  
  ggsave(filename=paste0("Figures_tsne_",drop_index,"_scater_logn_2.png"),
         plot = p,
         width = 24,
         height = 3*3+2)

  
  #p <- plot_comparison_tsne(drop_index, seed_value, 50,100)
  
  #ggsave(filename=paste0("Figures_tsne_",drop_index,".png"), plot = p, width = 24, height = 3)


}