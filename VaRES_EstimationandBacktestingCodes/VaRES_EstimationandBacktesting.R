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
source("Optimizations.R")
source("Function_VaR_VQR.R")
source("esr_backtest_modified.R")

pMCS = 0.10
#########################################
###           Loading Data            ###
#########################################
### The following *.txt files come from the Covariance_EstimationCodes' outcome

retornos = apply(read.table("OoSdata.txt"),1,mean)
mu = apply(read.table("OoSmu.txt"),1,mean)

e_RM2006 = read.table("e_RM2006.txt")
e_DCC = read.table("e_DCCc.txt")
e_ABC = read.table("e_ABC.txt")[-1,]
e_GDFM = read.table("e_GDFM-CHF.txt")[-1,]

s_RM2006 = read.table("s_RM2006.txt")
s_DCC = read.table("s_DCCc.txt")
s_ABC = read.table("s_ABC.txt")
s_GDFM = read.table("s_GDFM-CHF.txt")


sigma = cbind(as.numeric(s_RM2006[1,]),as.numeric(s_DCC[1,]),as.numeric(s_ABC[1,]),as.numeric(s_GDFM[1,]))

#########################################
###     Estimating Risk Measures      ###
#########################################
n = length(retornos)
a1 = 0.010
a2 = 0.025
a5 = 0.050

a_1 = 0.990
a_2 = 0.975
a_5 = 0.950
quant1 = c(1-a_1, 1-0.75*a_1 - 0.25, 1-0.5*a_1-0.5, 1-0.25*a_1-0.75)
quant2 = c(1-a_2, 1-0.75*a_2 - 0.25, 1-0.5*a_2-0.5, 1-0.25*a_2-0.75)
quant5 = c(1-a_5, 1-0.75*a_5 - 0.25, 1-0.5*a_5-0.5, 1-0.25*a_5-0.75)

VaR1 = VaR2 = VaR5 = ES1 = ES2 = ES5 = ES1ap = ES2ap = ES5ap = matrix(0,ncol = 4, nrow = n)

for(l in 1:n){
  set.seed(7656+l)
  boot_RM2006 = mu[l] + sample(e_RM2006[,l], 3*length(e_RM2006[,l]), replace = TRUE)*s_RM2006[1,l]
  boot_DCC    = mu[l] + sample(e_DCC[,l], 3*length(e_DCC[,l]), replace = TRUE)*s_DCC[1,l]
  boot_ABC    = mu[l] + sample(e_ABC[,l], 3*length(e_ABC[,l]), replace = TRUE)*s_ABC[1,l]
  boot_GDFM   = mu[l] + sample(e_GDFM[,l], 3*length(e_GDFM[,l]), replace = TRUE)*s_GDFM[1,l]
  
  
  VaR1[l,] = c(quantile(boot_RM2006, prob = a1),quantile(boot_DCC, prob = a1),quantile(boot_ABC, prob = a1), quantile(boot_GDFM, prob = a1))
  VaR2[l,] = c(quantile(boot_RM2006, prob = a2),quantile(boot_DCC, prob = a2),quantile(boot_ABC, prob = a2), quantile(boot_GDFM, prob = a2))
  VaR5[l,] = c(quantile(boot_RM2006, prob = a5),quantile(boot_DCC, prob = a5),quantile(boot_ABC, prob = a5), quantile(boot_GDFM, prob = a5))

  ## First Way To estimate the ES:
  ES1[l,] = c(mean(boot_RM2006[boot_RM2006 < VaR1[l,1]]), 
              mean(boot_DCC[boot_DCC < VaR1[l,2]]),
              mean(boot_ABC[boot_ABC < VaR1[l,3]]),
              mean(boot_GDFM[boot_GDFM < VaR1[l,4]]))
  
  ES2[l,] = c(mean(boot_RM2006[boot_RM2006 < VaR2[l,1]]), 
              mean(boot_DCC[boot_DCC < VaR2[l,2]]),
              mean(boot_ABC[boot_ABC < VaR2[l,3]]),
              mean(boot_GDFM[boot_GDFM < VaR2[l,4]]))
  
  ES5[l,] = c(mean(boot_RM2006[boot_RM2006 < VaR5[l,1]]), 
              mean(boot_DCC[boot_DCC < VaR5[l,2]]),
              mean(boot_ABC[boot_ABC < VaR5[l,3]]),
              mean(boot_GDFM[boot_GDFM < VaR5[l,4]]))
}

#########################################
###    Backtesting Risk Measures      ###
#########################################
n = ncol(VaR1)
a1 = 0.010
a2 = 0.025
a5 = 0.050

BackVaRES1 = BackVaRES2 = BackVaRES5 = matrix(0,ncol = 6 + 8,nrow = n) 
colnames(BackVaRES1) = colnames(BackVaRES2) = colnames(BackVaRES5) = c("Hits", "UC", "CC", "DQ", "VQ", "MFE", "NZ", "ESR_1", "ESR_2" ,"ESR_3", "QL", "FZG", "NZ", "AL")
row.names(BackVaRES1) = row.names(BackVaRES2) = row.names(BackVaRES5) = c("RM2006", "DCC", "ABC", "GDFM")


for (i in 1:n){
  set.seed(i+7153)
  BackT1 = BacktestVaR(retornos, VaR1[,i], alpha = a1, Lags = 4)
  BackT2 = BacktestVaR(retornos, VaR2[,i], alpha = a2, Lags = 4)
  BackT5 = BacktestVaR(retornos, VaR5[,i], alpha = a5, Lags = 4)
  
  
  BackVaRES1[i,] = c(mean(retornos < VaR1[,i])*100, 
                     BackT1$LRuc[2], BackT1$LRcc[2],BackT1$DQ$pvalue, VaR_VQR(retornos, VaR1[,i], a1),
                     er_backtest(retornos, VaR1[,i], ES1[,i], sigma[,i])$pvalue_onesided_standardized,
                     cc_backtest(retornos, VaR1[,i], ES1[,i],  alpha  = a1)$pvalue_twosided_simple, 
                     suppressWarnings(esr_backtest_modified(retornos, VaR1[,i], ES1[,i],alpha  = a1, B = 0, version = 1)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos, VaR1[,i], ES1[,i],alpha  = a1, B = 0, version = 2)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos, VaR1[,i], ES1[,i],alpha  = a1, B = 0, version = 3)$pvalue_onesided_asymptotic),
                     mean(QL(VaR1[,i],retornos, alpha = a1)),
                     mean(FZG(VaR1[,i], ES1[,i], retornos, alpha = a1)),
                     mean(NZ(VaR1[,i], ES1[,i], retornos, alpha = a1)),
                     mean(AL(VaR1[,i], ES1[,i], retornos, alpha = a1)))
  
  BackVaRES2[i,] = c(mean(retornos < VaR2[,i])*100, 
                     BackT2$LRuc[2], BackT2$LRcc[2],BackT2$DQ$pvalue, VaR_VQR(retornos, VaR2[,i], a2),
                     er_backtest(retornos, VaR2[,i], ES2[,i], sigma[,i])$pvalue_onesided_standardized,
                     cc_backtest(retornos, VaR2[,i], ES2[,i],  alpha  = a2)$pvalue_twosided_simple, 
                     suppressWarnings(esr_backtest_modified(retornos, VaR2[,i], ES2[,i],alpha  = a2, B = 0, version = 1)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos, VaR2[,i], ES2[,i],alpha  = a2, B = 0, version = 2)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos, VaR2[,i], ES2[,i],alpha  = a2, B = 0, version = 3)$pvalue_onesided_asymptotic),
                     mean(QL(VaR2[,i],retornos, alpha = a2)),
                     mean(FZG(VaR2[,i], ES2[,i], retornos, alpha = a2)),
                     mean(NZ(VaR2[,i], ES2[,i], retornos, alpha = a2)),
                     mean(AL(VaR2[,i], ES2[,i], retornos, alpha = a2)))
  
  
  BackVaRES5[i,] = c(mean(retornos < VaR5[,i])*100, 
                     BackT5$LRuc[2], BackT5$LRcc[2],BackT5$DQ$pvalue, VaR_VQR(retornos, VaR5[,i], a5),
                     er_backtest(retornos, VaR5[,i], ES5[,i], sigma[,i])$pvalue_onesided_standardized,
                     cc_backtest(retornos, VaR5[,i], ES5[,i],  alpha  = a5)$pvalue_twosided_simple, 
                     suppressWarnings(esr_backtest_modified(retornos, VaR5[,i], ES5[,i],alpha  = a5, B = 0, version = 1)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos, VaR5[,i], ES5[,i],alpha  = a5, B = 0, version = 2)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos, VaR5[,i], ES5[,i],alpha  = a5, B = 0, version = 3)$pvalue_onesided_asymptotic),
                     mean(QL(VaR5[,i],retornos, alpha = a5)),
                     mean(FZG(VaR5[,i], ES5[,i], retornos, alpha = a5)),
                     mean(NZ(VaR5[,i], ES5[,i], retornos, alpha = a5)),
                     mean(AL(VaR5[,i], ES5[,i], retornos, alpha = a5)))


}

VaRES_Table = rbind(BackVaRES1,BackVaRES2,BackVaRES5)
row.names(VaRES_Table) = rep(row.names(BackVaRES1),3)

p = 0.05

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
                       QL = format(round(QL,3),nsmall = 3),
                       FZG = format(round(FZG,3),nsmall = 3),
                       NZ.1 = format(round(NZ.1,3),nsmall = 3),
                       AL = format(round(AL,3),nsmall = 3))





#MCS (not informative)



pMCS = 0.1

MQL = QL(VaR1,retornos, alpha = a1)
estMCS.quick(MQL, test="t.range", B=5000, l=12, alpha = pMCS)
MFZG = FZG(VaR1,ES1,retornos, alpha = a1)
estMCS.quick(MFZG, test="t.range", B=5000, l=12, alpha = pMCS)
MNZ = NZ(VaR1,ES1,retornos, alpha = a1)
estMCS.quick(MNZ, test="t.range", B=5000, l=12, alpha = pMCS)
MAL = AL(VaR1,ES1,retornos, alpha = a1)
estMCS.quick(MAL, test="t.range", B=5000, l=12, alpha = pMCS)

MQL = QL(VaR2,retornos, alpha = a2)
estMCS.quick(MQL, test="t.range", B=5000, l=12, alpha = pMCS)
MFZG = FZG(VaR2,ES2,retornos, alpha = a2)
estMCS.quick(MFZG, test="t.range", B=5000, l=12, alpha = pMCS)
MNZ = NZ(VaR2,ES2,retornos, alpha = a2)
estMCS.quick(MNZ, test="t.range", B=5000, l=12, alpha = pMCS)
MAL = AL(VaR2,ES2,retornos, alpha = a2)
estMCS.quick(MAL, test="t.range", B=5000, l=12, alpha = pMCS)


MQL = QL(VaR5,retornos, alpha = a5)
estMCS.quick(MQL, test="t.range", B=5000, l=12, alpha = pMCS)
MFZG = FZG(VaR5,ES5,retornos, alpha = a5)
estMCS.quick(MFZG, test="t.range", B=5000, l=12, alpha = pMCS)
MNZ = NZ(VaR5,ES5,retornos, alpha = a5)
estMCS.quick(MNZ, test="t.range", B=5000, l=12, alpha = pMCS)
MAL = AL(VaR5,ES5,retornos, alpha = a5)
estMCS.quick(MAL, test="t.range", B=5000, l=12, alpha = pMCS)











                       

row.names(VaRES_Table) = c("RM20061", "DCC1", "ABC1", "GDFM1","RM20062", "DCC2", "ABC2", "GDFM2","RM20065", "DCC5", "ABC5", "GDFM5")

Caption = "One-step-ahead VaR and ES backtesting at 1\\% (top panel), 2.5\\% (middle panel) and 5\\% (bottom panel) risk levels. Out-of-Sample period from September 18, 2019 to July 1, 2020 (1386 obs)."
print(xtable(VaRES_Table, caption = Caption, align = "r|rrrrr|rrrrr|rrrr"), file = "VaRES_Table.tex", compress = FALSE)






##################################################
####   VaR Figure FULL PERIOD
##################################################


Dates = read.table("Prices_GDFM_CHF_VaR_APP2_DATES.txt")[-1,]
Dates = lubridate::ymd(Dates[-c(1:750)]) #Only OoS
which(Dates == "2020-03-12") #1310



VaR = c(c(VaR1[,1],VaR1[,2], VaR1[,3], VaR1[,4]),
        c(VaR2[,1],VaR2[,2], VaR2[,3], VaR2[,4]),
        c(VaR5[,1],VaR5[,2], VaR5[,3], VaR5[,4]))
risklevel = c(rep("1%",ncol(VaR1)*nrow(VaR1)),rep("2.5%",ncol(VaR2)*nrow(VaR2)),rep("5%",ncol(VaR5)*nrow(VaR5)))
models = rep(c(rep("RM2006", length(VaR1[,1])), 
         rep("DCCc", length(VaR1[,2])),
         rep("ABC", length(VaR1[,3])),
         rep("GDFM-CHF", length(VaR1[,4]))),3)
dates = rep(Dates, 12)
OoS = rep(retornos,12)

figureVaR = data.frame(dates,OoS,VaR, risklevel,models)

figureVaR$models = factor(figureVaR$models, levels = c("RM2006", "DCCc", "ABC", "GDFM-CHF"))


          
ggplot(figureVaR) + geom_line(aes(x = dates, y = OoS), colour="gray29") + 
  geom_vline(xintercept = as.Date("2020-03-12"), linetype = "dotted") + 
  geom_line(aes(x = dates, y = VaR), colour = "red2", linetype = "dashed") + facet_grid(risklevel~models, scale = "free") +
  ylab("Porfolio returns") + xlab(" ") + theme_bw() + 
  theme(legend.position = "none")



##################################################
####   VaR Figure 2018-2019 only
##################################################


Dates = read.table("Prices_GDFM_CHF_VaR_APP2_DATES.txt")[-1,]
Dates = lubridate::ymd(Dates[-c(1:750)]) #Only OoS
which(Dates == "2020-03-12") #1310



VaR = c(c(VaR1[,1],VaR1[,2], VaR1[,3], VaR1[,4]),
        c(VaR2[,1],VaR2[,2], VaR2[,3], VaR2[,4]),
        c(VaR5[,1],VaR5[,2], VaR5[,3], VaR5[,4]))
risklevel = c(rep("1%",ncol(VaR1)*nrow(VaR1)),rep("2.5%",ncol(VaR2)*nrow(VaR2)),rep("5%",ncol(VaR5)*nrow(VaR5)))
models = rep(c(rep("RM2006", length(VaR1[,1])), 
               rep("DCCc", length(VaR1[,2])),
               rep("ABC", length(VaR1[,3])),
               rep("GDFM-CHF", length(VaR1[,4]))),3)
dates = rep(Dates, 12)
OoS = rep(retornos,12)

figureVaR = data.frame(dates,OoS,VaR, risklevel,models) %>% filter(dates>'2017-12-31', dates < '2019-07-01')

figureVaR$models = factor(figureVaR$models, levels = c("RM2006", "DCCc", "ABC", "GDFM-CHF"))

ggplot(figureVaR) + geom_line(aes(x = dates, y = OoS), colour="gray29") + 
  geom_vline(xintercept = as.Date("2020-03-12"), linetype = "dotted") + 
  geom_line(aes(x = dates, y = VaR), colour = "red2", linetype = "dashed") + facet_grid(risklevel~models, scale = "free") +
  ylab("Porfolio returns") + xlab(" ") + theme_bw() + 
  theme(legend.position = "none")


  
##################################################
####   Tables before march 12
##################################################




for (i in 1:n){
  set.seed(i+123)
  BackT1 = BacktestVaR(retornos[1:1309], VaR1[1:1309,i], alpha = a1, Lags = 4)
  BackT2 = BacktestVaR(retornos[1:1309], VaR2[1:1309,i], alpha = a2, Lags = 4)
  BackT5 = BacktestVaR(retornos[1:1309], VaR5[1:1309,i], alpha = a5, Lags = 4)
  
  
  
  BackVaRES1[i,] = c(mean(retornos[1:1309] < VaR1[1:1309,i])*100, 
                     BackT1$LRuc[2], BackT1$LRcc[2],BackT1$DQ$pvalue, VaR_VQR(retornos[1:1309], VaR1[1:1309,i], a1),
                     er_backtest(retornos[1:1309], VaR1[1:1309,i], ES1[1:1309,i], sigma[1:1309,i])$pvalue_onesided_standardized,
                     cc_backtest(retornos[1:1309], VaR1[1:1309,i], ES1[1:1309,i],  alpha  = a1)$pvalue_twosided_simple, 
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR1[1:1309,i], ES1[1:1309,i],alpha  = a1, B = 0, version = 1)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR1[1:1309,i], ES1[1:1309,i],alpha  = a1, B = 0, version = 2)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR1[1:1309,i], ES1[1:1309,i],alpha  = a1, B = 0, version = 3)$pvalue_onesided_asymptotic),
                     mean(QL(VaR1[1:1309,i],retornos[1:1309], alpha = a1)),
                     mean(FZG(VaR1[1:1309,i], ES1[1:1309,i], retornos[1:1309], alpha = a1)),
                     mean(NZ(VaR1[1:1309,i], ES1[1:1309,i], retornos[1:1309], alpha = a1)),
                     mean(AL(VaR1[1:1309,i], ES1[1:1309,i], retornos[1:1309], alpha = a1)))
  
  BackVaRES2[i,] = c(mean(retornos[1:1309] < VaR2[1:1309,i])*100, 
                     BackT2$LRuc[2], BackT2$LRcc[2],BackT2$DQ$pvalue, VaR_VQR(retornos[1:1309], VaR2[1:1309,i], a2),
                     er_backtest(retornos[1:1309], VaR2[1:1309,i], ES2[1:1309,i], sigma[1:1309,i])$pvalue_onesided_standardized,
                     cc_backtest(retornos[1:1309], VaR2[1:1309,i], ES2[1:1309,i],  alpha  = a2)$pvalue_twosided_simple, 
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR2[1:1309,i], ES2[1:1309,i],alpha  = a2, B = 0, version = 1)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR2[1:1309,i], ES2[1:1309,i],alpha  = a2, B = 0, version = 2)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR2[1:1309,i], ES2[1:1309,i],alpha  = a2, B = 0, version = 3)$pvalue_onesided_asymptotic),
                     mean(QL(VaR2[1:1309,i],retornos[1:1309], alpha = a2)),
                     mean(FZG(VaR2[1:1309,i], ES2[1:1309,i], retornos[1:1309], alpha = a2)),
                     mean(NZ(VaR2[1:1309,i], ES2[1:1309,i], retornos[1:1309], alpha = a2)),
                     mean(AL(VaR2[1:1309,i], ES2[1:1309,i], retornos[1:1309], alpha = a2)))
  
  
  BackVaRES5[i,] = c(mean(retornos[1:1309] < VaR5[1:1309,i])*100, 
                     BackT5$LRuc[2], BackT5$LRcc[2],BackT5$DQ$pvalue, VaR_VQR(retornos[1:1309], VaR5[1:1309,i], a5),
                     er_backtest(retornos[1:1309], VaR5[1:1309,i], ES5[1:1309,i], sigma[1:1309,i])$pvalue_onesided_standardized,
                     cc_backtest(retornos[1:1309], VaR5[1:1309,i], ES5[1:1309,i],  alpha  = a5)$pvalue_twosided_simple, 
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR5[1:1309,i], ES5[1:1309,i],alpha  = a5, B = 0, version = 1)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR5[1:1309,i], ES5[1:1309,i],alpha  = a5, B = 0, version = 2)$pvalue_twosided_asymptotic),
                     suppressWarnings(esr_backtest_modified(retornos[1:1309], VaR5[1:1309,i], ES5[1:1309,i],alpha  = a5, B = 0, version = 3)$pvalue_onesided_asymptotic),
                     mean(QL(VaR5[1:1309,i],retornos[1:1309], alpha = a5)),
                     mean(FZG(VaR5[1:1309,i], ES5[1:1309,i], retornos[1:1309], alpha = a5)),
                     mean(NZ(VaR5[1:1309,i], ES5[1:1309,i], retornos[1:1309], alpha = a5)),
                     mean(AL(VaR5[1:1309,i], ES5[1:1309,i], retornos[1:1309], alpha = a5)))
  
  
}

VaRES_Table = rbind(BackVaRES1,BackVaRES2,BackVaRES5)
row.names(VaRES_Table) = rep(row.names(BackVaRES1),3)

p = 0.05

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
  QL = format(round(QL,3),nsmall = 3),
  FZG = format(round(FZG,3),nsmall = 3),
  NZ.1 = format(round(NZ.1,3),nsmall = 3),
  AL = format(round(AL,3),nsmall = 3))



row.names(VaRES_Table) = c("RM20061", "DCC1", "ABC1", "GDFM1","RM20062", "DCC2", "ABC2", "GDFM2","RM20065", "DCC5", "ABC5", "GDFM5")

Caption = "One-step-ahead VaR and ES backtesting at 1\\% (top panel), 2.5\\% (middle panel) and 5\\% (bottom panel) risk levels. Out-of-Sample period from September 18, 2019 to July 1, 2020 (1386 obs)."
print(xtable(VaRES_Table, caption = Caption, align = "r|rrrrr|rrrrr|rrrr"), file = "VaRES_Table_NoCOVID.tex", compress = FALSE)





