###   Real Data    Runchao Jiang

rm(list=ls())


### read in data and combine
all <- read.table("all.dat", header=F, sep="", 
                  col.names=c("pidnum","age","wtkg","hemo","homo","drugs","karnof","oprior","z30","zprior","preanti",
                              "race","gender","str2","strat","symptom","treat","offtrt","cd40","cd420","cd496","r",
                              "cd80","cd820"))

events <- read.table("events.dat", header=F, sep="", 
                     col.names=c("pidnum","delta","Ttilde","A"))

data0 <- merge(all, events, by="pidnum", all=T)
# dim(data0)
# class(data0)
# sapply(data0, class)

# A = 0: ZDV
# A = 1: ZDV + ddi
# A = 2: ZDV + ZAL
# A = 3: ddi



### remove incomplete data and treatments of no interest
keep <- (!is.na(rowSums(data0))) & (data0$A %in% c(1,2))
data1 <- data0[keep,]
data1$A[data1$A==2] <- 0 # recoding the treatments  
# dim(data1)
# 1-mean(data1$delta) # this is the proportion of people censored





### set the approporiate class for each variable
for (i in setdiff(names(data1), 
                  c("age","wtkg","karnof","preanti","cd40","cd420","cd496","cd80", "cd820", "delta","Ttilde","A")))
  data1[[i]] <- as.factor(data1[[i]])
# sapply(data1, class)

data1$fA <- as.factor(data1$A)
# sapply(data1, class)





### re-order the data by increasing Ttilde
index <- order(data1$Ttilde, decreasing=F)  
data <- data1[index,]





### save the raw data
save(data, file="RealData-RawData.RData")
# rm(list=setdiff(ls(), "DataList0"))











rm(list=ls())
library(survival)

load("RealData-RawData.RData")


### Build DataList with selected variables
DataList <- list()
DataList$Ttilde <- data$Ttilde
DataList$delta <- data$delta
DataList$A <- data$A
DataList$X <- model.matrix( ~ karnof+cd40+age, data)[,-1]
Neta <- ncol(DataList$X)+1




### save the DataList
save(DataList, Neta, file="RealData-SelectVar.RData")










rm(list=ls())
library(compiler)
enableJIT(3)
library(survival)
library(rgenoud)

load("RealData-SelectVar.RData")
source("RealData-Functions-v3.R")





Frd <- function(t0, Datalist, Neta)
{
  # t0 = 400, 600, 800, 1000
  ps <- Fps(DataList)
  reg <- Freg(DataList)
  PC <- FPC(DataList)
  prep <- Fprep(DataList, t0, ps, reg, PC)
  
  result <- matrix(NA, nrow=2, ncol=10, 
                   dimnames=list(c("I","A"),
                                 c("int", "Karnof", "CD40", "Age", "Stilde","se",
                                   "S1", "diff1.se", "S0", "diff0.se")))
  
  result[1, 1:(Neta+1)] <- Genetic.IPWE(t0, DataList, ps, T, Neta)
  result[1, "se"] <- IPWE.se(result[1, 1:Neta], t0, DataList, ps)
  
  result[2, 1:(Neta+1)] <- Genetic.AIPWE(t0, DataList, prep, T, Neta)
  result[2, "se"] <- AIPWE.se(result[2, 1:Neta], t0, DataList, ps, reg, PC)
  
  result[1, c("S1","S0")] <- result[2, c("S1","S0")] <- Fsub(DataList, t0)
  
  result[1, c("diff1.se", "diff0.se")] <- IPWE.diff.se(result[1, 1:Neta], t0, DataList, ps)[c("diff.se.surv.1", "diff.se.surv.0")]
  result[2, c("diff1.se", "diff0.se")] <- AIPWE.diff.se(result[2, 1:Neta], t0, DataList, ps, reg, PC)[c("diff.se.surv.1", "diff.se.surv.0")]
  
  return(result)
}


pre0 <- matrix(NA, nrow=8, ncol=10, 
               dimnames=list(paste0(rep(c(400, 600, 800, 1000), each=2), ",", rep(c("I","A"),4)),
                             c("int", "Karnof", "CD40", "Age", "Stilde","se",
                               "S1", "diff1.se", "S0", "diff0.se")))

# The below code takes a really long time to run. Use "ReadData-Analysis-Result.Rdata" if you can
pre0[1:2, ] <- Frd(t0=400, DataList, Neta)
pre0[3:4, ] <- Frd(t0=600, DataList, Neta)
pre0[5:6, ] <- Frd(t0=800, DataList, Neta)
pre0[7:8, ] <- Frd(t0=1000, DataList, Neta)


save(pre0, file="RealData-Analysis-Result.Rdata")





############################################
###   make tables


rm(list=ls())
library(xtable)
load("RealData-Analysis-Result.RData")




pre2 <- matrix(NA, nrow=8, ncol=4,
               dimnames=list(NULL, c("T1.Lower", "T1.Upper", "T0.Lower", "T0.Upper")))
pre2[,1] <- pre0[,"Stilde"]-pre0[,"S1"]-1.96*pre0[,"diff1.se"]
pre2[,2] <- pre0[,"Stilde"]-pre0[,"S1"]+1.96*pre0[,"diff1.se"]
pre2[,3] <- pre0[,"Stilde"]-pre0[,"S0"]-1.96*pre0[,"diff0.se"]
pre2[,4] <- pre0[,"Stilde"]-pre0[,"S0"]+1.96*pre0[,"diff0.se"]



pre1 <- data.frame(
  "t" = c("$400$", "", "$600$", "", "$800$", "", "$1000$", ""),
  "Method" = rep(c("I", "A"), 4),
  "Intercept" = pre0[, "int"],
  "Karnof" = pre0[, "Karnof"],
  "CD40" = pre0[, "CD40"],
  "Age" = pre0[, "Age"],
  "$\\tilde{S}(t; \\tilde{\\eta})$" = paste0("$", formatC(pre0[, "Stilde"], digits=3, format="f"), "\\;(", formatC(pre0[, "se"], digits=3, format="f"), ")$"),
  
  "Trt 1" = paste0("$(", formatC(pre2[, "T1.Lower"], digits=3, format="f"), ",", formatC(pre2[, "T1.Upper"], digits=3, format="f"), ")$"),
  "Trt 0" = paste0("$(", formatC(pre2[, "T0.Lower"], digits=3, format="f"), ",", formatC(pre2[, "T0.Upper"], digits=3, format="f"), ")$"),
  check.names = FALSE, stringsAsFactors = FALSE
)



print.xtable(xtable(pre1, digits=3, 
                    caption = "Real Data",
                    label = "tab:realdata"),
             hline.after = c(-1, 0, 8),
             sanitize.colnames.function=function(x) x, 
             sanitize.text.function=function(x) x, 
             include.rownames=FALSE)




