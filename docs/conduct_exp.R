rm(list = ls())
library(fda)
library(fdapace)
library(fda.usc)
library(lme4)
library(locpol)
library(rpart)
library(mgcv)
library(refund)
library(matrixcalc)
library(pracma)
library(randomForest)

source("exp_main_tree.R")
source("dat_gen.R")
source("FAM.R")
source("SIM.R")
source("FSIM.R")
source("FAME.R")
source("TREE.R")
source("TREEBOOST.R")
source("utils.R")

source("fdapace/R/predict.FPCA.R")
source("fdapace/R/HandleNumericsAndNAN.R")
source("fdapace/R/GetCEScores.R")
source("fdapace/R/GetINScores.R")
source("fdapace/R/Lwls1D.R")
source("fdapace/R/CreateFolds.R")
source("fdapace/R/CVLwls1D.R")
source("fdapace/R/RcppExports.R")
source("fdapace/R/GCVLwls1D1.R")
source("fdapace/R/Minb.R")

args <- commandArgs(trailingOnly = TRUE) 

make_up <- FALSE

x_type <- c("ferraty") 
nknots <-3; 
g_func_nos <- 1:5
ds <- 1:4
SNRs <- c(5,20)
seeds <- 1:100

# parameters for the competitors 
niter_FPPR <- 15; niters_FAME <- 1:15
nknots_FPPR <- 3;  nknots_FAME <- 3
nbasises_FGAM <- 15

thetas_FRF <- c(1,5,10,15,20,25,30) 

methods <- c("FLM1","FLM2","FAM","FAME","FPPR","FGAM","FRF","RFGroove")


if(length(args)== 0){ # test locally
  #methods <- c("FAME")
  x_type = "mattern"
  #g_func_nos <- 1
  #SNR <- c(5,20)
  thetas_FRF <- 1
  ds <- c(1)
  test <- TRUE
  make_up = FALSE
  RNGkind(sample.kind = "default")
  exp_id <- 1  # row index in the conduct_sheet 
  case_id <- 1 #  0 means only run competitors; 1 means running TFBoost 
  #niter <- 2
  niters_FAME <- 1
  niter_FPPR <- 2; 
  m <- 50
  if(case_id == 0){
    niter <- 1000  # doesn't matter
    ncol_num <- 4
  }
  if(case_id == 1){
    niter <- 1000  # need to be changed 
    ncol_num <- 2
  }
  if(case_id == 2){
    niter <- 1000  #doesn't matter 
    ncol_num <- 8
  }
}else{
  test <- FALSE
  exp_id <- as.numeric(args[1]) 
  case_id <- as.numeric(args[2]) 
  if(case_id == 0){
    niter <- 100  # doesn't matter
    ncol_num <- 4
  }
  if(case_id == 1){
    niter <- 1000  # need to be changed 
    ncol_num <- 2
  }
  if(case_id == 2){
    niter <- 1000  # doesn't matter 
    ncol_num <- 8
  }
}


shrinkage <- 0.05
num_indexs <- c(1,2,3)

all_settings <- expand.grid(nknots, "A", num_indexs)
all_settings <- rbind(all_settings, expand.grid(nknots, "B", num_indexs[1]))
colnames(all_settings) <- c("nknots","type", "num_index")


nmulti <- 5; num_dir <- 200;


#only train competitors 
if(case_id == 0){
  ds <- 1
}


all_exps <- expand.grid(ds, g_func_nos,  seeds, SNRs)
colnames(all_exps) <- c( "d","g_func_no","seed","SNR" )



if(nrow(all_exps)%%ncol_num!=0){
  
  if(nrow(all_exps)%%ncol_num!=0){
    conduct_sheet <- matrix(c(1:nrow(all_exps), rep(NA, ncol_num - nrow(all_exps)%%ncol_num)), ncol = ncol_num)
  }else{
    conduct_sheet <- matrix(c(1:nrow(all_exps)), ncol = ncol_num)
  }
}else{
  conduct_sheet <- matrix(1:nrow(all_exps), ncol = ncol_num)
}

print(dim(conduct_sheet))

if(dir.exists("Results") == FALSE){
  dir.create("Results")
}

if(dir.exists("Results_tmp") == FALSE){
  dir.create("Results_tmp")
}


conduct.exp <- function(exp_id = 1, conduct_sheet){
  
  for(i in 1:ncol(conduct_sheet)) {
    
    print(c(i,"th trial"))
    if(!is.na(conduct_sheet[exp_id,i])){
      
      g_func_no <- all_exps[conduct_sheet[exp_id, i],]$g_func_no
      seed <- all_exps[conduct_sheet[exp_id, i],]$seed
      d <-  all_exps[conduct_sheet[exp_id, i],]$d
      SNR <- all_exps[conduct_sheet[exp_id, i],]$SNR
      
      precision <- 6
   
      control.tree.list <- list()
      for(i in 1:nrow(all_settings)){
        control.tree.list[[i]] <- TFBoost.control(nknot = all_settings$nknots[i], shrinkage = shrinkage, loss = "l2", 
                                                     tree_control = TREE.control(nscreen = 30, nmulti = nmulti, num_dir = num_dir,
                                                                                 num_index = all_settings$num_index[i], tree_type = all_settings$type[i], d = d), trace = TRUE, precision = precision)
      }
      
      if(test){
        dir_name <- paste("Results_tmp/Results_", case_id, "_x_", x_type, "_g_", g_func_no,  "_SNR_", SNR, "_nknot", nknots,  sep = "")
      }else{
        dir_name <- paste("Results/Results_", case_id,"_x_", x_type, "_g_",g_func_no, "_SNR_", SNR, "_nknot", nknots,  sep = "")
      }
      
      if(dir.exists(dir_name) == FALSE){
        dir.create(dir_name)
      }
      
      
      start_time <- Sys.time()
      dat2save <- tryCatch(do.exp(seed, g_func_no, SNR, x_type,  niter = niter, niter_FPPR = niter_FPPR, niters_FAME = niters_FAME, 
                                  control.tree.list = control.tree.list,
                                  nknots_FPPR = nknots_FPPR, nknots_FAME = nknots_FAME, 
                                  nbasises_FGAM = nbasises_FGAM, thetas_FRF = thetas_FRF, case_id = case_id, methods = methods))
      end_time <- Sys.time()
      
      dat2save$time <- end_time - start_time
      print(dat2save$time)
      save(file = paste(dir_name, "/g_func_no_", g_func_no, "_case_id_",case_id, "_d_", d,  "_seed_", seed, ".RData", sep = ""), dat2save)
      
    }
  }
}

conduct.exp(exp_id = exp_id, conduct_sheet) 