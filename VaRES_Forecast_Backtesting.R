################################################
######  GDFM-CHF: VaR and ES application  ###### 
####         Tables and Figures             #### 
################################################
# The codes used to estimate the covariance are in Matlab, 
# this R code is only to perform the backtesting exercices

library(GAS)
library(rugarch)
library(esback)
library(esreg)
library(optimx)
library(xtable)
library(ggplot2)
library(stringr)
library(modelconf)
library(tidyverse)
source("aux_functions/scoring_functions.R")
source("aux_functions/Function_VaR_VQR.R")
source("aux_functions/esr_backtest_modified.R")

pMCS = 0.10
p = 0.05
pMCS = 0.25
#########################################
###           Loading Data            ###
#########################################
### The following *.txt files come from the Covariance_EstimationCodes' outcome
retornos = apply(read.table("OoSdata.txt"),1,mean)
mu = apply(read.table("OoSmu.txt"),1,mean)
rp_GDFM_Boot = read.table("rp_GDFM-CHF_Boot.txt")[-1,]
rp_GDFM = read.table("rp_GDFM-CHF.txt")[-1,]
rp_ABC = read.table("rp_ABC.txt")[-1,]
rp_RM2006 = read.table("rp_RM2006.txt")
rp_DCC = read.table("rp_DCCc.txt")

#########################################
###     Forecasting Risk Measures     ###
#########################################
n = length(retornos)
a1 = 0.010; a2 = 0.025; a5 = 0.050
VaR1m = VaR2m = VaR5m = ES1m = ES2m = ES5m = matrix(0,ncol = 5, nrow = n)
seednumber = 135 
for(l in 1:n){
  set.seed(seednumber+l)
  boot_RM2006m     = mu[l] + sample(rp_RM2006[,l], 5*length(rp_RM2006[,l]), replace = TRUE)
  boot_DCCm        = mu[l] + sample(rp_DCC[,l], 5*length(rp_DCC[,l]), replace = TRUE)
  boot_ABCm        = mu[l] + sample(rp_ABC[,l], 5*length(rp_ABC[,l]), replace = TRUE)
  boot_GDFMm       = mu[l] + sample(rp_GDFM[,l], 5*length(rp_GDFM[,l]), replace = TRUE)
  boot_GDFM_Boot_m = mu[l] + sample(rp_GDFM_Boot[,l], 5*length(rp_GDFM_Boot[,l]), replace = TRUE)
  
  VaR1m[l,] = c(quantile(boot_RM2006m, prob = a1),
                quantile(boot_DCCm, prob = a1),
                quantile(boot_ABCm, prob = a1), 
                quantile(boot_GDFMm, prob = a1), 
                quantile(boot_GDFM_Boot_m, prob = a1))
  VaR2m[l,] = c(quantile(boot_RM2006m, prob = a2),
                quantile(boot_DCCm, prob = a2),
                quantile(boot_ABCm, prob = a2), 
                quantile(boot_GDFMm, prob = a2), 
                quantile(boot_GDFM_Boot_m, prob = a2))
  VaR5m[l,] = c(quantile(boot_RM2006m, prob = a5),
                quantile(boot_DCCm, prob = a5),
                quantile(boot_ABCm, prob = a5), 
                quantile(boot_GDFMm, prob = a5), 
                quantile(boot_GDFM_Boot_m, prob = a5))
  
  ES1m[l,] = c(mean(boot_RM2006m[boot_RM2006m < VaR1m[l,1]]), 
               mean(boot_DCCm[boot_DCCm < VaR1m[l,2]]),
               mean(boot_ABCm[boot_ABCm < VaR1m[l,3]]), 
               mean(boot_GDFMm[boot_GDFMm < VaR1m[l,4]]),
               mean(boot_GDFM_Boot_m[boot_GDFM_Boot_m < VaR1m[l,5]]))
  ES2m[l,] = c(mean(boot_RM2006m[boot_RM2006m < VaR2m[l,1]]), 
               mean(boot_DCCm[boot_DCCm < VaR2m[l,2]]),
               mean(boot_ABCm[boot_ABCm < VaR2m[l,3]]), 
               mean(boot_GDFMm[boot_GDFMm < VaR2m[l,4]]),
               mean(boot_GDFM_Boot_m[boot_GDFM_Boot_m < VaR2m[l,5]]))
  ES5m[l,] = c(mean(boot_RM2006m[boot_RM2006m < VaR5m[l,1]]), 
               mean(boot_DCCm[boot_DCCm < VaR5m[l,2]]),
               mean(boot_ABCm[boot_ABCm < VaR5m[l,3]]), 
               mean(boot_GDFMm[boot_GDFMm < VaR5m[l,4]]),
               mean(boot_GDFM_Boot_m[boot_GDFM_Boot_m < VaR5m[l,5]]))
}


#########################################
###    Backtesting Risk Measures      ###
#########################################
n = ncol(VaR1m)
BackVaRES1 = BackVaRES2 = BackVaRES5 = matrix(0,ncol = 14,nrow = n) 
colnames(BackVaRES1) = colnames(BackVaRES2) = colnames(BackVaRES5) = c("Hits", "UC", "CC", "DQ", "VQ", "MFE", "NZ", "ESR_1", "ESR_2" ,"ESR_3", "QL", "FZG", "NZ", "AL")
row.names(BackVaRES1) = row.names(BackVaRES2) = row.names(BackVaRES5) = c("RM2006", "DCC", "ABC", "GDFM", "GDFM-Boot")

for (i in 1:n){
  set.seed(i+seednumber)
  BackT1 = BacktestVaR(retornos, VaR1m[,i], alpha = a1, Lags = 4)
  BackT2 = BacktestVaR(retornos, VaR2m[,i], alpha = a2, Lags = 4)
  BackT5 = BacktestVaR(retornos, VaR5m[,i], alpha = a5, Lags = 4)
  BackVaRES1[i,] = c(mean(retornos < VaR1m[,i])*100, 
                     BackT1$LRuc[2], BackT1$LRcc[2],BackT1$DQ$pvalue, VaR_VQR(retornos, VaR1m[,i], a1),
                     er_backtest(retornos, VaR1m[,i], ES1m[,i], sigma[,i])$pvalue_onesided_standardized,
                     cc_backtest(retornos, VaR1m[,i], ES1m[,i],  alpha  = a1)$pvalue_twosided_simple, 
                     suppressWarnings(esr_backtest_modified(retornos, VaR1m[,i], ES1m[,i],alpha  = a1, B = 0, version = 1)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos, VaR1m[,i], ES1m[,i],alpha  = a1, B = 0, version = 2)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos, VaR1m[,i], ES1m[,i],alpha  = a1, B = 0, version = 3)$pvalue_onesided_asymptotic),
                     mean(QL(VaR1m[,i],retornos, alpha = a1)),
                     mean(FZG(VaR1m[,i], ES1m[,i], retornos, alpha = a1)),
                     mean(NZ(VaR1m[,i], ES1m[,i], retornos, alpha = a1)),
                     mean(AL(VaR1m[,i], ES1m[,i], retornos, alpha = a1)))
  
  BackVaRES2[i,] = c(mean(retornos < VaR2m[,i])*100, 
                     BackT2$LRuc[2], BackT2$LRcc[2],BackT2$DQ$pvalue, VaR_VQR(retornos, VaR2m[,i], a2),
                     er_backtest(retornos, VaR2m[,i], ES2m[,i], sigma[,i])$pvalue_onesided_standardized,
                     cc_backtest(retornos, VaR2m[,i], ES2m[,i],  alpha  = a2)$pvalue_twosided_simple, 
                     suppressWarnings(esr_backtest_modified(retornos, VaR2m[,i], ES2m[,i],alpha  = a2, B = 0, version = 1)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos, VaR2m[,i], ES2m[,i],alpha  = a2, B = 0, version = 2)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos, VaR2m[,i], ES2m[,i],alpha  = a2, B = 0, version = 3)$pvalue_onesided_asymptotic),
                     mean(QL(VaR2m[,i],retornos, alpha = a2)),
                     mean(FZG(VaR2m[,i], ES2m[,i], retornos, alpha = a2)),
                     mean(NZ(VaR2m[,i], ES2m[,i], retornos, alpha = a2)),
                     mean(AL(VaR2m[,i], ES2m[,i], retornos, alpha = a2)))
  
  
  BackVaRES5[i,] = c(mean(retornos < VaR5m[,i])*100, 
                     BackT5$LRuc[2], BackT5$LRcc[2],BackT5$DQ$pvalue, VaR_VQR(retornos, VaR5m[,i], a5),
                     er_backtest(retornos, VaR5m[,i], ES5m[,i], sigma[,i])$pvalue_onesided_standardized,
                     cc_backtest(retornos, VaR5m[,i], ES5m[,i],  alpha  = a5)$pvalue_twosided_simple, 
                     suppressWarnings(esr_backtest_modified(retornos, VaR5m[,i], ES5m[,i],alpha  = a5, B = 0, version = 1)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos, VaR5m[,i], ES5m[,i],alpha  = a5, B = 0, version = 2)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos, VaR5m[,i], ES5m[,i],alpha  = a5, B = 0, version = 3)$pvalue_onesided_asymptotic),
                     mean(QL(VaR5m[,i],retornos, alpha = a5)),
                     mean(FZG(VaR5m[,i], ES5m[,i], retornos, alpha = a5)),
                     mean(NZ(VaR5m[,i], ES5m[,i], retornos, alpha = a5)),
                     mean(AL(VaR5m[,i], ES5m[,i], retornos, alpha = a5)))
}

VaRES_Table = rbind(BackVaRES1,BackVaRES2,BackVaRES5)
row.names(VaRES_Table) = rep(row.names(BackVaRES1),3)


## MCS 
MQL = QL(VaR1m,retornos, alpha = a1)
MCS_MQL1 = estMCS(MQL, test="t.range", B=5000, l=12)
MFZG = FZG(VaR1m,ES1m,retornos, alpha = a1)
MCS_MFZG1 = estMCS(MFZG, test="t.range", B=5000, l=12)
MNZ = NZ(VaR1m,ES1m,retornos, alpha = a1)
MCS_MNZ1 = estMCS(MNZ, test="t.range", B=5000, l=12)
MAL = AL(VaR1m,ES1m,retornos, alpha = a1)
MCS_MAL1 = estMCS(MAL, test="t.range", B=5000, l=12)

MQL = QL(VaR2m,retornos, alpha = a2)
MCS_MQL2 = estMCS(MQL, test="t.range", B=5000, l=12)
MFZG = FZG(VaR2m,ES2m,retornos, alpha = a2)
MCS_MFZG2 = estMCS(MFZG, test="t.range", B=5000, l=12)
MNZ = NZ(VaR2m,ES2m,retornos, alpha = a2)
MCS_MNZ2 = estMCS(MNZ, test="t.range", B=5000, l=12)
MAL = AL(VaR2m,ES2m,retornos, alpha = a2)
MCS_MAL2 = estMCS(MAL, test="t.range", B=5000, l=12)

MQL = QL(VaR5m,retornos, alpha = a5)
MCS_MQL5 = estMCS(MQL, test="t.range", B=5000, l=12)
MFZG = FZG(VaR5m,ES5m,retornos, alpha = a5)
MCS_MFZG5 = estMCS(MFZG, test="t.range", B=5000, l=12)
MNZ = NZ(VaR5m,ES5m,retornos, alpha = a5)
MCS_MNZ5 = estMCS(MNZ, test="t.range", B=5000, l=12)
MAL = AL(VaR5m,ES5m,retornos, alpha = a5)
MCS_MAL5 = estMCS(MAL, test="t.range", B=5000, l=12)

MCS_MQL = c(MCS_MQL1[,"MCS p-val"],MCS_MQL2[,"MCS p-val"],MCS_MQL5[,"MCS p-val"])
MCS_MFZG = c(MCS_MFZG1[,"MCS p-val"],MCS_MFZG2[,"MCS p-val"],MCS_MFZG5[,"MCS p-val"])
MCS_MNZ = c(MCS_MNZ1[,"MCS p-val"],MCS_MNZ2[,"MCS p-val"],MCS_MNZ5[,"MCS p-val"])
MCS_MAL = c(MCS_MAL1[,"MCS p-val"],MCS_MAL2[,"MCS p-val"],MCS_MAL5[,"MCS p-val"])


#########################################
###     Tables into LaTex format     ###
#########################################

### Adding cell colour
VaRES_Table = VaRES_Table %>% data.frame() %>% mutate(
  UC = ifelse(UC>p,paste0('\\cellcolor{gray!25}',format(round(UC,3),nsmall = 3)), format(round(UC,3),nsmall = 3)),
  CC = ifelse(CC>p,paste0('\\cellcolor{gray!25}',format(round(CC,3),nsmall = 3)),format(round(CC,3),nsmall = 3)),
  DQ = ifelse(DQ>p,paste0('\\cellcolor{gray!25}',format(round(DQ,3),nsmall = 3)),format(round(DQ,3),nsmall = 3)),
  VQ = ifelse(VQ>p,paste0('\\cellcolor{gray!25}',format(round(VQ,3),nsmall = 3)),format(round(VQ,3),nsmall = 3)),
  MFE = ifelse(MFE>p,paste0('\\cellcolor{gray!25}',format(round(MFE,3),nsmall = 3)),format(round(MFE,3),nsmall = 3)),
  NZ = ifelse(NZ>p,paste0('\\cellcolor{gray!25}',format(round(NZ,3),nsmall = 3)),format(round(NZ,3),nsmall = 3)),
  ESR_1 = ifelse(ESR_1>p,paste0('\\cellcolor{gray!25}',format(round(ESR_1,3),nsmall = 3)),format(round(ESR_1,3),nsmall = 3)),
  ESR_2 = ifelse(ESR_2>p,paste0('\\cellcolor{gray!25}',format(round(ESR_2,3),nsmall = 3)),format(round(ESR_2,3),nsmall = 3)),
  ESR_3 = ifelse(ESR_3>p,paste0('\\cellcolor{gray!25}',format(round(ESR_3,3),nsmall = 3)),format(round(ESR_3,3),nsmall = 3)),
  Hits = round(Hits,1),
  QL = ifelse(MCS_MQL>pMCS,paste0('\\cellcolor{gray!25}',format(round(QL,3),nsmall = 3)),format(round(QL,3),nsmall = 3)),
  FZG = ifelse(MCS_MFZG>pMCS,paste0('\\cellcolor{gray!25}',format(round(FZG,3),nsmall = 3)),format(round(FZG,3),nsmall = 3)),
  NZ.1 = ifelse(MCS_MNZ>pMCS,paste0('\\cellcolor{gray!25}',format(round(NZ.1,3),nsmall = 3)),format(round(NZ.1,3),nsmall = 3)),
  AL = ifelse(MCS_MAL>pMCS,paste0('\\cellcolor{gray!25}',format(round(AL,3),nsmall = 3)),format(round(AL,3),nsmall = 3)))


row.names(VaRES_Table) = c("RM20061", "DCC1", "ABC1", "GDFM1", "GDFM-Boot1", "RM20062", "DCC2", "ABC2", "GDFM2","GDFM-Boot2","RM20065", "DCC5", "ABC5", "GDFM5","GDFM-Boot5")
Caption = "One-step-ahead VaR and ES backtesting at 1\\% (top panel), 2.5\\% (middle panel) and 5\\% (bottom panel) risk levels. Out-of-Sample period from December 29, 2014 to July 1, 2020 (1386 obs)."
print(xtable(VaRES_Table, caption = Caption, align = "r|rrrrr|rrrrr|rrrr"), file = "VaRES_Table_m_full.tex", compress = FALSE)
VaRES_Table

######################################################
####   Backtesting Risk Measures-m before COVID    ###
######################################################
## Observation 1309 corresponds to 11 March 2020

for (i in 1:n){
  set.seed(i+seednumber)
  BackT1 = BacktestVaR(retornos[1:1309], VaR1m[1:1309,i], alpha = a1, Lags = 4)
  BackT2 = BacktestVaR(retornos[1:1309], VaR2m[1:1309,i], alpha = a2, Lags = 4)
  BackT5 = BacktestVaR(retornos[1:1309], VaR5m[1:1309,i], alpha = a5, Lags = 4)
  
  BackVaRES1[i,] = c(mean(retornos[1:1309] < VaR1m[1:1309,i])*100, 
                     BackT1$LRuc[2], BackT1$LRcc[2],BackT1$DQ$pvalue, VaR_VQR(retornos[1:1309], VaR1m[1:1309,i], a1),
                     er_backtest(retornos[1:1309], VaR1m[1:1309,i], ES1m[1:1309,i], sigma[1:1309,i])$pvalue_onesided_standardized,
                     cc_backtest(retornos[1:1309], VaR1m[1:1309,i], ES1m[1:1309,i],  alpha  = a1)$pvalue_twosided_simple, 
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR1m[1:1309,i], ES1m[1:1309,i],alpha  = a1, B = 0, version = 1)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR1m[1:1309,i], ES1m[1:1309,i],alpha  = a1, B = 0, version = 2)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR1m[1:1309,i], ES1m[1:1309,i],alpha  = a1, B = 0, version = 3)$pvalue_onesided_asymptotic),
                     mean(QL(VaR1m[1:1309,i],retornos[1:1309], alpha = a1)),
                     mean(FZG(VaR1m[1:1309,i], ES1m[1:1309,i], retornos[1:1309], alpha = a1)),
                     mean(NZ(VaR1m[1:1309,i], ES1m[1:1309,i], retornos[1:1309], alpha = a1)),
                     mean(AL(VaR1m[1:1309,i], ES1m[1:1309,i], retornos[1:1309], alpha = a1)))
  
  BackVaRES2[i,] = c(mean(retornos[1:1309] < VaR2m[1:1309,i])*100, 
                     BackT2$LRuc[2], BackT2$LRcc[2],BackT2$DQ$pvalue, VaR_VQR(retornos[1:1309], VaR2m[1:1309,i], a2),
                     er_backtest(retornos[1:1309], VaR2m[1:1309,i], ES2m[1:1309,i], sigma[1:1309,i])$pvalue_onesided_standardized,
                     cc_backtest(retornos[1:1309], VaR2m[1:1309,i], ES2m[1:1309,i],  alpha  = a2)$pvalue_twosided_simple, 
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR2m[1:1309,i], ES2m[1:1309,i],alpha  = a2, B = 0, version = 1)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR2m[1:1309,i], ES2m[1:1309,i],alpha  = a2, B = 0, version = 2)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR2m[1:1309,i], ES2m[1:1309,i],alpha  = a2, B = 0, version = 3)$pvalue_onesided_asymptotic),
                     mean(QL(VaR2m[1:1309,i],retornos[1:1309], alpha = a2)),
                     mean(FZG(VaR2m[1:1309,i], ES2m[1:1309,i], retornos[1:1309], alpha = a2)),
                     mean(NZ(VaR2m[1:1309,i], ES2m[1:1309,i], retornos[1:1309], alpha = a2)),
                     mean(AL(VaR2m[1:1309,i], ES2m[1:1309,i], retornos[1:1309], alpha = a2)))
  
  
  BackVaRES5[i,] = c(mean(retornos[1:1309] < VaR5m[1:1309,i])*100, 
                     BackT5$LRuc[2], BackT5$LRcc[2],BackT5$DQ$pvalue, VaR_VQR(retornos[1:1309], VaR5m[1:1309,i], a5),
                     er_backtest(retornos[1:1309], VaR5m[1:1309,i], ES5m[1:1309,i], sigma[1:1309,i])$pvalue_onesided_standardized,
                     cc_backtest(retornos[1:1309], VaR5m[1:1309,i], ES5m[1:1309,i],  alpha  = a5)$pvalue_twosided_simple, 
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR5m[1:1309,i], ES5m[1:1309,i],alpha  = a5, B = 0, version = 1)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR5m[1:1309,i], ES5m[1:1309,i],alpha  = a5, B = 0, version = 2)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR5m[1:1309,i], ES5m[1:1309,i],alpha  = a5, B = 0, version = 3)$pvalue_onesided_asymptotic),
                     mean(QL(VaR5m[1:1309,i],retornos[1:1309], alpha = a5)),
                     mean(FZG(VaR5m[1:1309,i], ES5m[1:1309,i], retornos[1:1309], alpha = a5)),
                     mean(NZ(VaR5m[1:1309,i], ES5m[1:1309,i], retornos[1:1309], alpha = a5)),
                     mean(AL(VaR5m[1:1309,i], ES5m[1:1309,i], retornos[1:1309], alpha = a5)))
}

VaRES_Table = rbind(BackVaRES1,BackVaRES2,BackVaRES5)
row.names(VaRES_Table) = rep(row.names(BackVaRES1),3)


## MCS 

MQL = QL(VaR1m[1:1309,],retornos[1:1309], alpha = a1)
MCS_MQL1 = estMCS(MQL, test="t.range", B=5000, l=12)
MFZG = FZG(VaR1m[1:1309,],ES1m[1:1309,],retornos[1:1309], alpha = a1)
MCS_MFZG1 = estMCS(MFZG, test="t.range", B=5000, l=12)
MNZ = NZ(VaR1m[1:1309,],ES1m[1:1309,],retornos[1:1309], alpha = a1)
MCS_MNZ1 = estMCS(MNZ, test="t.range", B=5000, l=12)
MAL = AL(VaR1m[1:1309,],ES1m[1:1309,],retornos[1:1309], alpha = a1)
MCS_MAL1 = estMCS(MAL, test="t.range", B=5000, l=12)

MQL = QL(VaR2m[1:1309,],retornos[1:1309], alpha = a2)
MCS_MQL2 = estMCS(MQL, test="t.range", B=5000, l=12)
MFZG = FZG(VaR2m[1:1309,],ES2m[1:1309,],retornos[1:1309], alpha = a2)
MCS_MFZG2 = estMCS(MFZG, test="t.range", B=5000, l=12)
MNZ = NZ(VaR2m[1:1309,],ES2m[1:1309,],retornos[1:1309], alpha = a2)
MCS_MNZ2 = estMCS(MNZ, test="t.range", B=5000, l=12)
MAL = AL(VaR2m[1:1309,],ES2m[1:1309,],retornos[1:1309], alpha = a2)
MCS_MAL2 = estMCS(MAL, test="t.range", B=5000, l=12)


MQL = QL(VaR5m[1:1309,],retornos[1:1309], alpha = a5)
MCS_MQL5 = estMCS(MQL, test="t.range", B=5000, l=12)
MFZG = FZG(VaR5m[1:1309,],ES5m[1:1309,],retornos[1:1309], alpha = a5)
MCS_MFZG5 = estMCS(MFZG, test="t.range", B=5000, l=12)
MNZ = NZ(VaR5m[1:1309,],ES5m[1:1309,],retornos[1:1309], alpha = a5)
MCS_MNZ5 = estMCS(MNZ, test="t.range", B=5000, l=12)
MAL = AL(VaR5m[1:1309,],ES5m[1:1309,],retornos[1:1309], alpha = a5)
MCS_MAL5 = estMCS(MAL, test="t.range", B=5000, l=12)

MCS_MQL = c(MCS_MQL1[,"MCS p-val"],MCS_MQL2[,"MCS p-val"],MCS_MQL5[,"MCS p-val"])
MCS_MFZG = c(MCS_MFZG1[,"MCS p-val"],MCS_MFZG2[,"MCS p-val"],MCS_MFZG5[,"MCS p-val"])
MCS_MNZ = c(MCS_MNZ1[,"MCS p-val"],MCS_MNZ2[,"MCS p-val"],MCS_MNZ5[,"MCS p-val"])
MCS_MAL = c(MCS_MAL1[,"MCS p-val"],MCS_MAL2[,"MCS p-val"],MCS_MAL5[,"MCS p-val"])


#########################################
###     Tables into LaTex format     ###
#########################################

VaRES_Table = VaRES_Table %>% data.frame() %>% mutate(
  UC = ifelse(UC>p,paste0('\\cellcolor{gray!25}',format(round(UC,3),nsmall = 3)), format(round(UC,3),nsmall = 3)),
  CC = ifelse(CC>p,paste0('\\cellcolor{gray!25}',format(round(CC,3),nsmall = 3)),format(round(CC,3),nsmall = 3)),
  DQ = ifelse(DQ>p,paste0('\\cellcolor{gray!25}',format(round(DQ,3),nsmall = 3)),format(round(DQ,3),nsmall = 3)),
  VQ = ifelse(VQ>p,paste0('\\cellcolor{gray!25}',format(round(VQ,3),nsmall = 3)),format(round(VQ,3),nsmall = 3)),
  MFE = ifelse(MFE>p,paste0('\\cellcolor{gray!25}',format(round(MFE,3),nsmall = 3)),format(round(MFE,3),nsmall = 3)),
  NZ = ifelse(NZ>p,paste0('\\cellcolor{gray!25}',format(round(NZ,3),nsmall = 3)),format(round(NZ,3),nsmall = 3)),
  ESR_1 = ifelse(ESR_1>p,paste0('\\cellcolor{gray!25}',format(round(ESR_1,3),nsmall = 3)),format(round(ESR_1,3),nsmall = 3)),
  ESR_2 = ifelse(ESR_2>p,paste0('\\cellcolor{gray!25}',format(round(ESR_2,3),nsmall = 3)),format(round(ESR_2,3),nsmall = 3)),
  ESR_3 = ifelse(ESR_3>p,paste0('\\cellcolor{gray!25}',format(round(ESR_3,3),nsmall = 3)),format(round(ESR_3,3),nsmall = 3)),
  Hits = round(Hits,1),
  QL = ifelse(MCS_MQL>pMCS,paste0('\\cellcolor{gray!25}',format(round(QL,3),nsmall = 3)),format(round(QL,3),nsmall = 3)),
  FZG = ifelse(MCS_MFZG>pMCS,paste0('\\cellcolor{gray!25}',format(round(FZG,3),nsmall = 3)),format(round(FZG,3),nsmall = 3)),
  NZ.1 = ifelse(MCS_MNZ>pMCS,paste0('\\cellcolor{gray!25}',format(round(NZ.1,3),nsmall = 3)),format(round(NZ.1,3),nsmall = 3)),
  AL = ifelse(MCS_MAL>pMCS,paste0('\\cellcolor{gray!25}',format(round(AL,3),nsmall = 3)),format(round(AL,3),nsmall = 3)))


row.names(VaRES_Table) = c("RM20061", "DCC1", "ABC1", "GDFM1", "GDFM-Boot1", "RM20062", "DCC2", "ABC2", "GDFM2","GDFM-Boot2","RM20065", "DCC5", "ABC5", "GDFM5","GDFM-Boot5")
Caption = "One-step-ahead VaR and ES backtesting at 1\\% (top panel), 2.5\\% (middle panel) and 5\\% (bottom panel) risk levels. Out-of-Sample period from December 29, 2014 to March 11, 2020."
print(xtable(VaRES_Table, caption = Caption, align = "r|rrrrr|rrrrr|rrrr"), file = "VaRES_Table_NoCOVID_m_full.tex", compress = FALSE)






##################################################
####   VaR Figure FULL PERIOD
##################################################


Dates = read.table("data/Prices_GDFM_CHF_VaR_APP2_DATES.txt")[-1,]
Dates = lubridate::ymd(Dates[-c(1:750)]) #Only OoS
which(Dates == "2020-03-12") #1310

VaRm = c(c(VaR1m[,1],VaR1m[,2], VaR1m[,3], VaR1m[,4], VaR1m[,5]),
         c(VaR2m[,1],VaR2m[,2], VaR2m[,3], VaR2m[,4], VaR2m[,5]),
         c(VaR5m[,1],VaR5m[,2], VaR5m[,3], VaR5m[,4], VaR5m[,5]))

risklevel = c(rep("1%",ncol(VaR1m)*nrow(VaR1m)),rep("2.5%",ncol(VaR2m)*nrow(VaR2m)),rep("5%",ncol(VaR5m)*nrow(VaR5m)))
models = rep(c(rep("RM2006", length(VaR1m[,1])), 
               rep("DCCc", length(VaR1m[,2])),
               rep("ABC", length(VaR1m[,3])),
               rep("GDFM-CHF", length(VaR1m[,4])),
               rep("GDFM-CHF-Boot", length(VaR1m[,5]))),3)
dates = rep(Dates, 15)
OoS = rep(retornos,15)

figureVaR = data.frame(dates,OoS,VaRm, risklevel,models)

figureVaR$models = factor(figureVaR$models, levels = c("RM2006", "DCCc", "ABC", "GDFM-CHF", "GDFM-CHF-Boot"))



ggplot(figureVaR) + geom_line(aes(x = dates, y = OoS), colour="gray29") + 
  geom_vline(xintercept = as.Date("2020-03-12"), linetype = "dotted") + 
  geom_line(aes(x = dates, y = VaRm), colour = "red2", linetype = "dashed") + facet_grid(risklevel~models, scale = "free") +
  ylab(" ") + xlab(" ") + theme_bw() + 
  theme(legend.position = "none")



##################################################
####   VaR Figure 2018-2019 only
##################################################


Dates = read.table("data/Prices_GDFM_CHF_VaR_APP2_DATES.txt")[-1,]
Dates = lubridate::ymd(Dates[-c(1:750)]) #Only OoS
which(Dates == "2020-03-12") #1310



VaRm = c(c(VaR1m[,1],VaR1m[,2], VaR1m[,3], VaR1m[,4], VaR1m[,5]),
         c(VaR2m[,1],VaR2m[,2], VaR2m[,3], VaR2m[,4], VaR2m[,5]),
         c(VaR5m[,1],VaR5m[,2], VaR5m[,3], VaR5m[,4], VaR5m[,5]))

risklevel = c(rep("1%",ncol(VaR1m)*nrow(VaR1m)),rep("2.5%",ncol(VaR2m)*nrow(VaR2m)),rep("5%",ncol(VaR5m)*nrow(VaR5m)))
models = rep(c(rep("RM2006", length(VaR1m[,1])), 
               rep("DCCc", length(VaR1m[,2])),
               rep("ABC", length(VaR1m[,3])),
               rep("GDFM-CHF", length(VaR1m[,4])),
               rep("GDFM-CHF-Boot", length(VaR1m[,5]))),3)
dates = rep(Dates, 15)
OoS = rep(retornos,15)

figureVaR = data.frame(dates,OoS,VaRm, risklevel,models) %>% filter(dates>'2017-12-31', dates < '2019-07-01')

figureVaR$models = factor(figureVaR$models, levels = c("RM2006", "DCCc", "ABC", "GDFM-CHF","GDFM-CHF-Boot"))

ggplot(figureVaR) + geom_line(aes(x = dates, y = OoS), colour="gray29") + 
  geom_vline(xintercept = as.Date("2020-03-12"), linetype = "dotted") + 
  geom_line(aes(x = dates, y = VaRm), colour = "red2", linetype = "dashed") + 
  scale_x_date(breaks = as.Date(c("2018-02-01","2018-10-01","2019-05-01")),
               labels = c("2018-02", "2018-10", "2019-05")) +
  facet_grid(risklevel~models, scale = "free") +
  ylab(" ") + xlab(" ") + theme_bw() + 
  theme(legend.position = "none", axis.text.x=element_text(angle=360))








