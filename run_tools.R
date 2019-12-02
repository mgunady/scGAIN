
#library(vsn)
#library(Seurat)
#library(scater)
#library(edgeR)
#library(gridExtra)
#library(R.matlab)
#library(cowplot)
#library(biomaRt)
#library(data.table)
#library(lattice)
#library(VennDiagram)
#library(Rtsne)
#library(DT)
#library(ggpubr)
#library(ggsignif)
#library(scatterplot3d)
#library(ggfortify)
#library(refGenome)
#library(pheatmap)
#library(dendsort)
#library(entropy)
#library(splatter)
#library(mcriPalettes)
#library(plotly)
#library(factoextra)
#library(cluster)
#library(NbClust)
#library(fpc)
#library(class)
#library(SC3)


######## Data Generation ##############
# the function generating the simulation data using bioconductor package Splatter
generate_simulation_splatter <- function(dropout_index, seed_value, nGenes = 800, nCells = 20000){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  #    seed_value: the random seed
  #        nGenes: the number of genes in the simulation data. The default is 800
  
  # Set up the parameters
  params = newSplatParams()
  
  params = setParams(params, list(batchCells = nCells,
                                  nGenes = nGenes,
                                  group.prob = c(0.20, 0.35, 0.45),
                                  de.prob = c(0.045, 0.045, 0.045),
                                  de.facLoc = 0.1,
                                  de.facScale = 0.4)
  )
  
  # Set up the vector of dropout.mid
  dropout_mid = c(4, 5, 5.5, 3.2, 2.3)
  
  # determine if it is a good parameter
  if(dropout_index > length(dropout_mid)){
    
    stop(
      paste0('The dropout_index shold not be greater than ', 
             length(dropout_mid), 
             ' . Please input a proper one.\n')
    )
    
  }
  print(paste0("dropout.mid= ", dropout_mid[dropout_index]))
  
  # Generate the simulation data using Splatter package
  sim = splatSimulateGroups(params,
                            dropout.type = "batch",#"experiment",
                            dropout.shape = -1,
                            dropout.mid = dropout_mid[dropout_index],
                            seed = seed_value)
  
  # genereate the cpm levels of the true simulation data
  data_true = cpm(sim@assays$data$TrueCounts)
  
  data_dropout = data_true
  
  # generate the dropout data based on the counts in sim
  data_dropout[counts(sim) == 0] = 0
  
  # calculate the dropout rate
  percentage_zeros = round(nnzero(data_dropout == 0, na.counted = NA)/
                             (dim(data_dropout)[1]*dim(data_dropout)[2])*100)
  
  
  # generate the bulk RNAseq data
  data_bulk = data.frame(val = rowMeans(data_true))
  
  # define the data list for the simulation data
  # indcluding: data_true: true data
  #          data_dropout: dropout data
  #             data_bluk: bulk data
  #      percentage_zeros: dropout rate
  #                 group: the group label
  
  data = list()
  
  data$data_bulk = data_bulk
  
  data$data_dropout = data_dropout
  
  data$data_true = data_true
  
  data$percentage_zeros = percentage_zeros
  
  data$group = colData(sim)@listData$Group
  
  return(data)
}

# generate the simulation data and save the data
generate_save_data <- function(dropout_index, seed_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # generate the simulation data
  data_simulation = generate_simulation_splatter(dropout_index, seed_value)
  
  # generate the folder saving the simulation data
  dir.create(file.path('simulation_data'), showWarnings = FALSE)
  
  # save the data as RDS format
  saveRDS(data_simulation, 
          file = paste0('simulation_data/simulation_data_dropout_index_',
                        dropout_index, 
                        '_seed_', 
                        seed_value,
                        '.rds')
  )

  # write the data
  write.table(data_simulation$data_dropout,
              paste0('simulation_data/simulation_data_dropout_index_',
                        dropout_index, 
                        '_seed_', 
                        seed_value,
                        '.csv'),
              sep=',',
              row.names = F,
              col.names = F
  )

  
}

######## DrImpute ##############
run_drimpute <- function(dropout_index, seed_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # load the raw data
  data <- readRDS(file = paste0('simulation_data/simulation_data_dropout_index_',
                                dropout_index,
                                '_seed_',
                                seed_value,
                                '.rds')
  )
  
  # build the folder saving the imputed data using DrImpute
  path <- "imputation_drimpute_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using DrImpute
  data_dropout <- as.matrix(data$data_dropout)[,1:1000]
  
  exdata <- DrImpute(data_dropout)
  
  # write the data
  write.table(exdata,
              paste0(path, "drimpute_",dropout_index,"_",seed_value,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}

######## scImpute ##############
run_scimpute <- function(dropout_index, seed_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  
  # load the data
  data = readRDS(file = paste0('simulation_data/simulation_data_dropout_index_',
                               dropout_index,
                               '_seed_',
                               seed_value,
                               '.rds')
  )
  
  # build the folder saving the imputed data using scimpute method
  path = "imputation_scimpute_data/"
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using scImpute
  data_dropout = data$data_dropout[,1:1000]
  
  write.table(data_dropout,paste0(path, 
                                  "dropout_scimpute_",
                                  dropout_index,"_",
                                  seed_value,".csv"),
              sep=',',
              row.names = TRUE,
              col.names = TRUE
  )
  
  file.remove(paste0(path, "scimpute_", dropout_index, "_", seed_value,"_*"))
  
  # run scImpute
  scimpute(
    paste0(path, "dropout_scimpute_",dropout_index,"_",seed_value,".csv"),
    infile = "csv",           
    outfile = "csv",         
    out_dir = paste0(path, "scimpute_", dropout_index, "_", seed_value,"_"),
    drop_thre = 0.5,         
    Kcluster = 2,
    ncores = 2)       
  # 
  # clean the data
  data_dropout = read.table( file = paste0(path, "scimpute_",
                                           dropout_index, "_",
                                           seed_value,
                                           "_scimpute_count.csv") ,
                             header = TRUE, sep=",")
  
  data_dropout$X = NULL
  
  # save the data
  write.table(data_dropout,
              paste0(path, "data_imputation_scimpute_",dropout_index,"_",seed_value,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}

######## MAGIC ##############
# here we use python script to run magic and we post the python script here.
# -----------------------------------------------------------------------------
# import sys
# import os
# import magic 
# import pandas as pd
# 
# cwd = os.getcwd()
# 
# if not os.path.exists(cwd+"/imputation_magic_data"):
#   os.makedirs(cwd+"/imputation_magic_data")
#
# X =pd.read_csv("simulation_data/simulation_data_dropout_index_"+str(dropout_value)+"_seed_"+str(seed_value)+".txt",sep = ' ', header=None)
#
# magic_operator = magic.MAGIC()
# X_magic = magic_operator.fit_transform(X.T)
#
# out_magic = X_magic.T
# out_magic.to_csv(cwd+"/imputation_magic_data/magic_"+str(dropout_value)+"_"+str(seed_value)+".csv", sep = '\t', header= None)
# -----------------------------------------------------------------------------


######## VIPER ##############
run_viper <- function(dropout_index, seed_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # load the raw data
  data <- readRDS(file = paste0('simulation_data/simulation_data_dropout_index_',
                                dropout_index,
                                '_seed_',
                                seed_value,
                                '.rds')
  )
  
  # build the folder saving the imputed data using DrImpute
  path <- "imputation_viper_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using VIPER
  data_dropout <- as.matrix(data$data_dropout)[,1:1000]
  
  exdata <- VIPER(data_dropout, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1, 
                  report = FALSE, outdir = NULL, prefix = NULL)
  
  # write the data
  write.table(exdata$imputed,
              paste0(path, "viper_",dropout_index,"_",seed_value,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}

######## SCRABBLE ##############
run_scrabble <- function(dropout_index, seed_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  
  # load the data
  data = readRDS(file = paste0('simulation_data/simulation_data_dropout_index_',
                               dropout_index,
                               '_seed_',
                               seed_value,
                               '.rds')
  )
  
  
  
  path = "imputation_scrabble_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using DrImpute
  data1 = list()
  
  data1[[1]] = data$data_dropout[,1:1000]
  
  data1[[2]] = data$data_bulk[,1:1000]
  
  # set up the parameters
  parameter = c(1, 1e-06, 1e-04)
  
  # run scrabble
  result = scrabble(data1,
                    parameter = parameter, 
                    nIter = 60,
                    error_out_threshold = 1e-7, 
                    nIter_inner = 100,
                    error_inner_threshold = 1e-5)
  
  # write the data
  write.table(result,
              paste0(path, "scrabble_",dropout_index,"_",seed_value,".csv"),
              sep=',',
              row.names = F,
              col.names = F)
  
}


# Plot tsne function
plot_comparison_tsne <- function(dropout_index,
                                 seed_value, 
                                 initial_dims_value,
                                 perplexity_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  # initial_dims_value: the initial dimensions used in the tsne
  # perplexity_value: the perplexity used in the tnse
  
  data_simulation = readRDS(file = paste0('simulation_data/simulation_data_dropout_index_',
                                          dropout_index,
                                          '_seed_',
                                          seed_value,
                                          '.rds')
  )
  
  # true data
  data_true = data_simulation$data_true
  
  # raw data (dropout data)
  data_dropout = data_simulation$data_dropout
  
  # load imputed data from Drimpute
  data_drimpute = read.table( file = paste0("imputation_drimpute_data/drimpute_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  # load imputed data from scImpute
  data_scimpute = read.table( file = paste0("imputation_scimpute_data/data_imputation_scimpute_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  # load the magic results 
  data_magic = read.csv(paste0("imputation_magic_data/magic_",
                         dropout_index,"_",
                         seed_value,".csv"),
                  header = FALSE,
                  sep = ",")
  
  data_magic$V1 = NULL
  
  data_magic = as.matrix(data_magic)[,1:1000]
  
  # load imputed data from VIPER
  data_viper = read.table( file = paste0("imputation_viper_data/viper_",
                                         dropout_index, "_",
                                         seed_value,
                                         ".csv") ,
                           header = FALSE, sep=","
  )
  
  
  # load imputed data from scrabble
  data_scrabble = read.table( file = paste0("imputation_scrabble_data/scrabble_",
                                            "2", "_",
                                            "1000",
                                            ".csv") ,
                              header = FALSE, sep=","
  )

  # load imputed data from gain
  data_gain = read.table( file = paste0("imputation_gain_data/gain_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )[,1:1000]

  print("Done loading data")
  
  # define the plot list
  pl = list()
  set.seed(1) # Set a seed if you want reproducible results
  
  # calculate the tsne of true data
  true_tsne = Rtsne(t(as.matrix(data_true)[,1:1000]), 
                    initial_dims = initial_dims_value, 
                    perplexity = perplexity_value)
  
  
  pl[[1]] = plot_pca_singlecell(true_tsne$Y,data_simulation$group[1:1000], "True")
  print("Done Truth")
  
  set.seed(1)
  
  # calculate the tsne of dropout data
  dropout_tsne = Rtsne(t(as.matrix(cpm(data_dropout))[,1:1000]), 
                       initial_dims = initial_dims_value, 
                       perplexity = perplexity_value)
  
  pl[[2]] = plot_pca_singlecell(dropout_tsne$Y,data_simulation$group[1:1000], "Dropout")
  print("Done dropout")
  
  set.seed(1)
  
  # calculate the tsne of drimpute imputed data
  drimpute_tsne = Rtsne(t(as.matrix(data_drimpute)[,1:1000]), 
                        initial_dims = initial_dims_value, 
                        perplexity = perplexity_value)
  
  pl[[3]] = plot_pca_singlecell(drimpute_tsne$Y,data_simulation$group[1:1000], "DrImpute")
  print("Done drImpute")
  
  set.seed(1)
  
  # calculate the tsne of scimpute imputed data
  scimpute_tsne = Rtsne(t(as.matrix(data_scimpute)[,1:1000]), 
                        initial_dims = initial_dims_value, 
                        perplexity = perplexity_value)
  
  pl[[4]] = plot_pca_singlecell(scimpute_tsne$Y,data_simulation$group[1:1000], "scImpute")
  print("Done scImpute")

  set.seed(1)
  
  # calculate the tsne of magic imputed data
  magic_tsne = Rtsne(t(as.matrix(data_magic)[,1:1000]), 
                     initial_dims = initial_dims_value, 
                     perplexity = perplexity_value)
  
  pl[[5]] = plot_pca_singlecell(magic_tsne$Y,data_simulation$group[1:1000], "Magic")
  print("Done Magic")

  set.seed(1)
  
  # calculate the tsne of viper imputed data
  viper_tsne = Rtsne(t(as.matrix(data_viper)[,1:1000]), 
                        initial_dims = initial_dims_value, 
                        perplexity = perplexity_value)
  
  pl[[6]] = plot_pca_singlecell(viper_tsne$Y,data_simulation$group[1:1000], "Viper")
  print("Done Viper")

  set.seed(1)
  
  # calculate the tsne of scrabble imputed data
  scrabble_tsne = Rtsne(t(as.matrix(data_scrabble)[,1:1000]), 
                        initial_dims = initial_dims_value, 
                        perplexity = perplexity_value)
  
  #pl[[7]] = plot_pca_singlecell(scrabble_tsne$Y,data_simulation$group[1:1000], "Scrabble")
  #print("Done Scrabble")

  set.seed(1)
  
  # calculate the tsne of Gain imputed data
  print(dim(data_gain))
  gain_tsne = Rtsne(t(as.matrix(data_gain)[,1:1000]), 
                        initial_dims = initial_dims_value, 
                        perplexity = perplexity_value
			)
  
  pl[[7]] = plot_pca_singlecell(gain_tsne$Y,data_simulation$group[1:1000], "Gain")
  print("Done Gain")


  main = grid.arrange(grobs = pl,ncol = 7)
  
  return(main)
  
}

plot_pca_singlecell <- function(Y, groups, name) {
  tsne_plot <- data.frame(x = Y[,1], y = Y[,2], col = groups)
  ntypes = length(unique(groups))
  cols <- palette(brewer.pal(ntypes, "Accent"))[as.fumeric(as.character(groups))] # Use Set3 pallette
  cols <- palette(brewer.pal(ntypes, "Accent"))[as.fumeric(as.character(groups))] # Use Set3 pallette
  #print(unique(as.character(groups)))
  #print(unique(cols))

  p = ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col)) + #, shape=cat)) + facet_wrap(~cat, scale="free") +
  scale_color_manual(values=unique(cols)) +theme_bw() + ggtitle(name) +
  xlab("Dimension 1") + ylab("Dimension 2") + 
  #geom_label(aes(x, y, label=lab),
  #           data=data.frame(x=20, y=37, lab=paste0("R=", format(round(r, 2), nsmall = 2))),alpha=0.7, size=6) +
  theme(legend.position = "top", legend.title = element_blank(),  plot.title = element_text(hjust = 0.5))
  return(p)
}

####################################

# Plot tsne function
plot_tsnes <- function(dropout_index,
                                 seed_value, 
                                 initial_dims_value,
                                 perplexity_value, columns=c("drimpute", "scimpute", "viper", "magic", "gain")){
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  # initial_dims_value: the initial dimensions used in the tsne
  # perplexity_value: the perplexity used in the tnse
  
  data_simulation = readRDS(file = paste0('simulation_data/simulation_data_dropout_index_',
                                          dropout_index,
                                          '_seed_',
                                          seed_value,
                                          '.rds')
  )
  
  # true data
  data_true = data_simulation$data_true
  
  # raw data (dropout data)
  data_dropout = data_simulation$data_dropout
  
if("drimpute" %in% columns) {
  # load imputed data from Drimpute
  data_drimpute = read.table( file = paste0("imputation_drimpute_data/drimpute_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
}
 
if("scimpute" %in% columns) {
  # load imputed data from scImpute
  data_scimpute = read.table( file = paste0("imputation_scimpute_data/data_imputation_scimpute_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
}

if("magic" %in% columns) {  
  # load the magic results 
  data_magic = read.csv(paste0("imputation_magic_data/magic_",
                         dropout_index,"_",
                         seed_value,"_.csv"),
                  header = FALSE,
                  sep = ",")
  
  #print(data_magic[1:6, 1:6])
  data_magic$V1 = NULL
  
  data_magic = as.matrix(data_magic)
}
 
if("viper" %in% columns) { 
  # load imputed data from VIPER
  data_viper = read.table( file = paste0("imputation_viper_data/viper_",
                                         dropout_index, "_",
                                         seed_value,
                                         ".csv") ,
                           header = FALSE, sep=","
  )
}

if("deepimpute" %in% columns) {
  # load imputed data from DeepImpute
  data_deep = read.table( file = paste0("imputation_deepimpute_data/deepimpute_",
                                         dropout_index, "_",
                                         seed_value,
                                         ".csv") ,
                           header = FALSE, sep=","
  )
  data_deep = (data_deep)
  print(dim(data_deep))
}

if("dca" %in% columns) {
  # load imputed data from DCA
  data_dca = read.table( file = paste0("imputation_dca_data/dca_",
                                         dropout_index, "_",
                                         seed_value,
                                         ".csv") ,
                           header = FALSE, sep=","
  )
  data_dca = (data_dca)
  print(dim(data_dca))
  print(data_dca[1:6, 1:6])
}

if("gain" %in% columns) {  
  # load imputed data from gain
  data_gain = read.table( file = paste0("imputation_gain_data/gain_",
                                            dropout_index, "_",
                                            seed_value,
                                            "_10000_reproduced.csv") ,
                              header = FALSE, sep=","
  )
  data_truth_logn = read.table( file = paste0("simulation_data/simulation_data_dropout_index_",
                                            dropout_index, "_seed_",
                                            seed_value,
                                            "_logn_true.csv") ,
                              header = FALSE, sep=","
  )
}

  print("Done loading data")


  # define the plot lis20
  pl = list()
for(i in 1:21) {
pl[[i]] = ggplot() + theme_void()
}

#######################

  pl[[1]] = plot_correlation(as.matrix(data_true), as.matrix(data_true), data_dropout, "Truth")
  print("Done Truth")

  pl[[2]] = plot_correlation(as.matrix(data_dropout), as.matrix(data_true), data_dropout, "Dropout")
  print("Done Dropout")

if("drimpute" %in% columns) {
  pl[[3]] = plot_correlation(myNorm(log(as.matrix(data_drimpute)+1)), myNorm(log(as.matrix(data_true)+1))[,1:1000], data_dropout[,1:1000], "DrImpute")
  print("Done DrImpute")
}
if("scimpute" %in% columns) {
  pl[[4]] = plot_correlation(myNorm(log(as.matrix(data_scimpute)+1)), myNorm(log(as.matrix(data_true)+1))[,1:1000], data_dropout[,1:1000], "scImpute")
  print("Done scImpute")
}
if("magic" %in% columns) {
  pl[[5]] = plot_correlation(myNorm(log(as.matrix(data_magic)+1)), myNorm(log(as.matrix(data_true)+1)), data_dropout, "Magic")
  print("Done Magic")
}
if("viper" %in% columns) {
  pl[[6]] = plot_correlation(myNorm(log(as.matrix(data_viper)+1)), myNorm(log(as.matrix(data_true)+1))[,1:1000], data_dropout[,1:1000], "Viper")
  print("Done Viper")
}

  #pl[[5]] = plot_correlation(myNorm(log(as.matrix(data_deep)+1)), myNorm(log(as.matrix(data_true)+1)), data_dropout, "DeepImpute")
  #print("Done DeepImpute")

  #pl[[6]] = plot_correlation(myNorm(log(as.matrix(data_dca)+1)), myNorm(log(as.matrix(data_true)+1)), data_dropout, "DCA")
  #print("Done DCA")
if("gain" %in% columns) {
  pl[[7]] = plot_correlation(as.matrix(data_gain), as.matrix(data_truth_logn), data_dropout, "scGain")
  print("Done Gain")
}

#######################
  
  #set.seed(1) # Set a seed if you want reproducible results

  pl[[8]] = plot_mytsne(as.matrix(data_true), data_simulation$group)
  print("Done Truth")

  pl[[9]] = plot_mytsne(as.matrix(data_dropout), data_simulation$group)
  print("Done Dropout")
if("drimpute" %in% columns) {
  pl[[10]] = plot_mytsne(as.matrix(data_drimpute), data_simulation$group[1:1000])
  print("Done DrImpute")
}
if("scimpute" %in% columns) {
  pl[[11]] = plot_mytsne(as.matrix(data_scimpute), data_simulation$group[1:1000])
  print("Done scImpute")
}
if("magic" %in% columns) {
  pl[[12]] = plot_mytsne(log(as.matrix(data_magic)+1), data_simulation$group)#, 5000)
  print("Done Magic")
}
if("viper" %in% columns) {
  pl[[13]] = plot_mytsne(as.matrix(data_viper), data_simulation$group[1:1000])
  print("Done Viper")
}
  #pl[[12]] = plot_mytsne(as.matrix(data_deep), data_simulation$group)#, 5000)
  #print("Done DeepImpute")

  #pl[[13]] = plot_mytsne(as.matrix(data_dca), data_simulation$group)#, 5000)
  #print("Done DCA")

if("gain" %in% columns) {
  pl[[14]] = plot_mytsne(as.matrix(data_gain), data_simulation$group)#, 5000)
  print("Done Gain")
}

#######################

  genesVars = rowVars(data_true)
  names(genesVars) = paste0("Gene", seq(dim(data_true)[1]))
  genesVars = genesVars[order(genesVars)]
  genes_f = names(genesVars)[1:100]
  ncells = 1000
  print(genes_f)

  pl[[15]] = plot_mytsne(as.matrix(data_true), data_simulation$group, ncells, genes_f)
  print("Done Truth")

  pl[[16]] = plot_mytsne(as.matrix(data_dropout), data_simulation$group, ncells, genes_f)
  print("Done Dropout")
if("drimpute" %in% columns) {
  pl[[17]] = plot_mytsne(as.matrix(data_drimpute), data_simulation$group[1:1000], ncells, genes_f)
  print("Done DrImpute")
}
if("scimpute" %in% columns) {
  pl[[18]] = plot_mytsne(as.matrix(data_scimpute), data_simulation$group[1:1000], ncells, genes_f)
  print("Done scImpute")
}
if("magic" %in% columns) {
  pl[[19]] = plot_mytsne(as.matrix(data_magic), data_simulation$group, ncells, genes_f)#, 5000)
  print("Done Magic")
}
if("viper" %in% columns) {
  pl[[20]] = plot_mytsne(as.matrix(data_viper), data_simulation$group[1:1000], ncells, genes_f)
  print("Done Viper")
}

  #pl[[12]] = plot_mytsne(as.matrix(data_deep), data_simulation$group, ncells, genes_f)#, 5000)
  #print("Done DeepImpute")

  #pl[[13]] = plot_mytsne(as.matrix(data_dca), data_simulation$group, ncells, genes_f)#, 5000)
  #print("Done DCA")

if("gain" %in% columns) {
  pl[[21]] = plot_mytsne(as.matrix(data_gain), data_simulation$group, ncells, genes_f)#, 5000)
  print("Done Gain")
}

#######################

  main = grid.arrange(grobs = pl,ncol = 7)
  
  return(main)

}

plot_mytsne <- function(data, labels, ncells=1000, genes_f=NA) {
  rownames(data) = paste0("Gene", seq(dim(data)[1]))
  colnames(data) = paste0("Cell", seq(dim(data)[2]))
  genesVars = rowVars(data)
  if(is.na(genes_f)) {
    genes_f = order(genesVars, decreasing = T)[1:100]
  }

  set.seed(1000)
  sceset_imp <- SingleCellExperiment(assays = list(logcounts = data[genes_f,]), colData=data.frame(Group=labels))
  p=plotTSNE(sceset_imp[,1:ncells], colour_by = "Group")
  p = p+theme(legend.position = "none")
  return(p)
}

# Input: genesxcells
myNorm <- function(data) {
  X = data
  X_maxs = apply(X, 2, max)
  print(length(X_maxs))
  X_n = sweep(X, 2, X_maxs, "/")
  return(X_n)
}

plot_correlation <- function(data, data_truth, data_dropout, method, ncells=100) {
  mask = data_dropout
  mask[data_dropout>0] = 1
  mask = 1-mask
  #print(dim(mask))
  #print(dim(data))
  imp_v = as.vector((data*mask)[,1:ncells])
  dropout_truth_v = as.vector((data_truth*mask)[,1:ncells])
  cor = cor.test(imp_v, dropout_truth_v)
  print(cor)
  mse = mean((imp_v-dropout_truth_v)**2)
  print(mse)
  p=data.frame(truth=dropout_truth_v, imp=imp_v) %>% ggplot(aes(truth, imp)) + geom_point(alpha=0.3, shape=16, size=1) + geom_abline(slope=1, color="blue") + theme_bw() + ggtitle(paste0(method, " R=", format(round(cor$estimate, 2), nsmall = 2), " MSE=", format(round(mse, 4), nsmall = 4) ))
  return(p)
}
