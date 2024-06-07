# CD Ayasse

# Testing non-Rasch and Rasch IRT Models under different 
#   underlying data conditions/structures
#     - Gaded Response Model (GRM; Samejima)
#     - Rasch: Partial Credit Model (PCM)
#     - Graded Partial Credit Model (GPCM)

# Varying:
#     - Item loadings/slopes
#         - Variability - constrained in Rasch, not in GRM
#         - Strength
#     - Sample size
#     - Number of items
#     - Number of response categories
# Also:
#     - Presence and strength of item local dependence (LD)
#     - Unidimensional vs not

rm(list=ls())
gc() 

options(rstudio.help.showDataPreview = FALSE) #temp while R fixes a bug


# Save folder tree name for data, figures, tables: ----------------------------- ####

# main.folder.name 
# code.name 
# tables.name 
# figures.name 
# object.name 
# data.name 
# genfntn.name 





# Load packages & custom functions: -------------------------------------------- ####

library(tidyverse)
library(tidyr) #data wrangling
library(dplyr) #data wrangling
library(mirt) #for fitting models
library(faux) #correlated variables
library(doParallel) #for parallel
library(parallel) #for random number stream
library(foreach) #for parallel
# library(doRNG) #to make fully reproducible parallel loops

# General-use custom functions:
source(paste0(genfntn.name,"CDA_R2_Function_dump_df_mat_to_file.R")) #custom function for printing table to .doc file
source(paste0(genfntn.name,"CDA_IRT_ScoringStats_Fntns.R")) #custom function for bifactor/IRT scoring statistics  
source(paste0(genfntn.name,"CDA_OrdItemsDat_Fntn.R")) #custom function for creating dataset with ordinal items
source(paste0(genfntn.name,"CDA_ConvertUnitDifftime_Fntn.R")) #custom function to convert unit of difftime objects (since R sometimes doesn't automatically)
source(paste0(genfntn.name,"CDA_FileNameCreate_Fntn.R")) #custom fntn filename str, GENERAL 

# Custom functions that are specific to this project or script:
source(paste0(code.name,"IRTSims_GRMvsRasch_ParallelCondsOutside_Fntn.R")) #function for parallelization for THIS script/sim


# Create a new class to use for outputting 2 results in parallel loops:
resultsAndData <- function(resultOut=NULL, resultData=NULL)
{
  me <- list(
    resultOut = resultOut,
    resultData = resultData
  )
  
  # Set the name for the class
  class(me) <- append(class(me),"resultsAndData")
  return(me)
}





# Version for all conditions etc, set up for specific batches: ----------------- ####
# --- Set up conditions, parameters, etc.: ------------------------------------- ####


# WHICH TO RUN?:
run.currently <- c("Uni","2Dim") # "2Dim" # "Uni" # "LD" # 
sub.conds <- TRUE # FALSE # 
if (sub.conds==TRUE) {
  cur.sample.run <- 500 # c(25, 100, 250) #c(25,100,250,500) #
  cur.avg.gf.run <- c(1,2) # 1 #
}


# Parallel properties:
no_cores <- detectCores() -2 # *(3/4) # *(5/6) # *(1/2) # *(1/4) # *(1/3) # *(1/4) # *(2/3) # use cores minus some of available
no_cores <- round(no_cores)
message(paste0("Using ",no_cores," cores out of ",detectCores()," total cores."))


# Set some values: --- #

#   To help with run time, optionally decide which models to run:
models.run <- c("GRM","PCM","GPCM") #c("GRM","PCM") #c("GRM","PCM","RSM")

#     Decide whether to increase the number of iterations from default settings:
num.cycles <- NULL # 5000 # 2000 # if NULL will remain default; if set as a numeric value, 
#                   will run that many iterations to fit IRT models
if (!is.null(num.cycles)) {
  n.cyc.str <- paste0("NCyc",num.cycles,"_")
} else { n.cyc.str <- "" }

#   Only save datasets up to a certain number of replications?
#     (To conserve memory/computational resources)
max.repl.dat.rbind <- 100 # num.repls

#   Randomization, parallelization, looping:
seed.start <- 07062024+1
add.seed.new <- 1717
num.repls <- 500 # 100 # 250 # 10 # 
num.out.repls <- 1 # 2 # 

#   Sample and items:
sample.size.vec <-  c(25,100,250,500) # c(25, 100, 250) # 
number.items.vec <- c(20, 12, 5) #c(5,12,20) #c(5,10,15,20)
num.resp.cat.vec <- c(3, 7) #c(3,5,7)

#   Item parameters and intercepts:
#     Item intercepts depend on the number of response categories:
int.3 <- c(-1,1); int.5 <- c(-2,-0.75,0.75,2); int.7 <- c(-3,-2,-0.75,0.75,2,3)
intercepts <- list(int.3, int.5, int.7)
names(intercepts) <- c("3","5","7")
#     Vary the average General Factor (GF) loading and GF loading variability 
#       (since Rasch assumes constant):
avg.gf.loads <- c(1, 2, 3) # c(1, 3) # 
vary.gf.loads <- c(0, 0.125, 0.5, 0.75) # c(0, 0.125, 0.25, 0.5, 0.75) 


# Set up local dependence and unidimensionality conditions so that we test
#   each in the standard/unchanged/expected condition of the other:

#   Local Dependence between Items?:
ld.load.vec <- c("*0","*1.0","*2") # c("*0","*0.5","*1.0","*2")
ld.items <- c(2, 5)

#   Unidimensional?
uni.cond.vec <- c("Unidimensional","Two.Fac.Half.Items") #,"Unrelated")

#   Local Dependence and Unidimensionality, as we will run them:
uni.ld.df <- data.frame(
  "LD.Load"=c(ld.load.vec, rep(ld.load.vec[1], times=(length(uni.cond.vec)-1))),
  "Uni.Condition"=c(rep(uni.cond.vec[1], times=length(ld.load.vec)), uni.cond.vec[-1]))



# Number of conditions:
num.conds.all <- nrow(uni.ld.df) * length(sample.size.vec) * 
  length(number.items.vec) * length(num.resp.cat.vec) * 
  length(avg.gf.loads) * length(vary.gf.loads)
msg.all.conds <- paste0("Total number of condition-combinations is: ",num.conds.all)
message(msg.all.conds)



# Run Separately: Unidimensional + No LD; Two-Dimensional (No LD); With LD (Unid.)
if (length(run.currently)==1) {
  if (run.currently == "Uni") {
    uni.ld.df.run <- uni.ld.df[which(uni.ld.df$LD.Load=="*0" & 
                                       uni.ld.df$Uni.Condition=="Unidimensional"),]
  } else if (run.currently == "2Dim") {
    uni.ld.df.run <- uni.ld.df[which(uni.ld.df$LD.Load=="*0" & 
                                       uni.ld.df$Uni.Condition=="Two.Fac.Half.Items"),]
  } else if (run.currently == "LD") {
    uni.ld.df.run <- uni.ld.df[which(uni.ld.df$LD.Load!="*0" & 
                                       uni.ld.df$Uni.Condition=="Unidimensional"),]
  }
} else if (length(run.currently)==2) {
  if ("Uni" %in% run.currently & "2Dim" %in% run.currently) {
    uni.ld.df.run <- uni.ld.df[which(
      uni.ld.df$LD.Load=="*0" & (uni.ld.df$Uni.Condition=="Unidimensional" |
                                   uni.ld.df$Uni.Condition=="Two.Fac.Half.Items")),]
  } else if ("Uni" %in% run.currently & "LD" %in% run.currently) {
    uni.ld.df.run <- uni.ld.df[which(
      (uni.ld.df$LD.Load=="*0" & uni.ld.df$Uni.Condition=="Unidimensional") |
        (uni.ld.df$LD.Load!="*0" & uni.ld.df$Uni.Condition=="Unidimensional")),]
  } else if ("2Dim" %in% run.currently & "LD" %in% run.currently) {
    uni.ld.df.run <- uni.ld.df[which(
      (uni.ld.df$LD.Load=="*0" & uni.ld.df$Uni.Condition=="Two.Fac.Half.Items") |
        (uni.ld.df$LD.Load!="*0" & uni.ld.df$Uni.Condition=="Unidimensional")),]
  } else {
    warning("run.currently conditions (",paste0(run.currently, collapse=","),
            ") not all recognized")
  } 
} else if (length(run.currently)==3) {
  if ("Uni" %in% run.currently & "2Dim" %in% run.currently & "LD" %in% run.currently) {
    uni.ld.df.run <- uni.ld.df[which(
      (uni.ld.df$LD.Load=="*0" & (uni.ld.df$Uni.Condition=="Unidimensional" |
                                    uni.ld.df$Uni.Condition=="Two.Fac.Half.Items")) |
        (uni.ld.df$LD.Load!="*0" & uni.ld.df$Uni.Condition=="Unidimensional")),]
  } else {
    warning("run.currently conditions (",paste0(run.currently, collapse=", "),
            ") not all recognized")
  }
}
#   Check that run.currently matches uni.ld.df.run:
message(paste0("run.currently: ",paste0(run.currently,collapse=", "),
               "\nConditions running:"))
print(uni.ld.df.run)

#   Updated: Number of conditions:
num.conds.all <- nrow(uni.ld.df.run) * length(sample.size.vec) * 
  length(number.items.vec) * length(num.resp.cat.vec) * 
  length(avg.gf.loads) * length(vary.gf.loads)
msg.all.conds <- paste0("Total number of condition-combinations is: ",num.conds.all)
message(msg.all.conds)

message.num.cores <- paste0("(Cores: ",no_cores," / ",
                            parallel::detectCores(),"; Replications: ",num.repls,
                            " inner x ",num.out.repls," outer)")
message.longertime.conds <- paste0(
  "Run levels for big-time-commitment variables: \n   Num.Resp.Cat: ",
  paste0(num.resp.cat.vec, collapse=", "),"\n   Num.Items: ",
  paste0(number.items.vec, collapse=", "),"\n  (Sample.Size: ",
  paste0(sample.size.vec, collapse=", "),")")
message(paste0(msg.all.conds, "\n", message.num.cores,
               " \n", message.longertime.conds))


if (sub.conds==FALSE) {
  
  number.items.vec.run <- number.items.vec
  num.resp.cat.vec.run <- num.resp.cat.vec
  sample.size.vec.run <- sample.size.vec
  avg.gf.loads.run <- avg.gf.loads
  vary.gf.loads.run <- vary.gf.loads
  
  run.seq.cur <- ""
}

message(paste0(msg.all.conds, "\n", message.num.cores,
               " \n", message.longertime.conds))





# ------ Choose Subset of Conditions: ------------------------------------------ ####


if (sub.conds==TRUE) {
  
  
  # View all original conditions (besides LD & Dimensionality):
  number.items.vec
  vary.gf.loads
  num.resp.cat.vec
  sample.size.vec
  avg.gf.loads
  
  
  
  # Narrow in on specific conditions to run, from the condition variables:
  
  #   Pick conditions (besides LD & Dimensionality):
  sample.size.vec.run <- sample.size.vec[which(sample.size.vec %in% cur.sample.run)] # sample.size.vec[c(2:3)] # sample.size.vec[1] # sample.size.vec[4] # sample.size.vec # 
  avg.gf.loads.run <- avg.gf.loads[which(avg.gf.loads %in% cur.avg.gf.run)] # avg.gf.loads # avg.gf.loads[3] # avg.gf.loads[c(1:2)] # 
  num.resp.cat.vec.run <- num.resp.cat.vec # !!!NOTE: resp.cat=7 [2] produces BIG increase in run time # num.resp.cat.vec[1] # num.resp.cat.vec[2] # num.resp.cat.vec[c(1,2)]
  number.items.vec.run <- number.items.vec # !!!NOTE: n.items=20 [1] produces BIG increase in run time # number.items.vec[3] # number.items.vec[1] # number.items.vec[c(2:3)] # number.items.vec[3] # number.items.vec[2] # number.items.vec # order: 20,12,5
  vary.gf.loads.run <- vary.gf.loads # 
  #                   ^#vary.gf.loads[c(1,2)] # vary.gf.loads[c(4:length(vary.gf.loads))] # 
  #                   ^#vary.gf.loads[-3] # vary.gf.loads # vary.gf.loads[(length(vary.gf.loads)-1)] # 
  #                   ^#vary.gf.loads[c(1,2,length(vary.gf.loads))] # vary.gf.loads[c(3:(length(vary.gf.loads)-1))] # 
  
  #   View current/chosen conditions (besides LD & Dimensionality):
  number.items.vec.run
  vary.gf.loads.run
  num.resp.cat.vec.run
  sample.size.vec.run
  avg.gf.loads.run
  
  
  
  # Updated: Number of conditions:
  num.conds.all <- nrow(uni.ld.df.run) * length(sample.size.vec.run) * 
    length(number.items.vec.run) * length(num.resp.cat.vec.run) * 
    length(avg.gf.loads.run) * length(vary.gf.loads.run)
  msg.all.conds <- paste0("Total number of condition-combinations is: ",num.conds.all)
  message(msg.all.conds)
  
  message.num.cores <- paste0("(Cores: ",no_cores," / ",
                              parallel::detectCores(),"; Replications: ",num.repls,
                              " inner x ",num.out.repls," outer)")
  message.longertime.conds <- paste0(
    "Run levels for big-time-commitment variables: \n   Num.Resp.Cat: ",
    paste0(num.resp.cat.vec.run, collapse=", "),"\n   Num.Items: ",
    paste0(number.items.vec.run, collapse=", "),"\n  (Sample.Size: ",
    paste0(sample.size.vec.run, collapse=", "),")")
  message(paste0(msg.all.conds, "\n", message.num.cores,
                 " \n", message.longertime.conds))
  
  
  
  # Planned Sequence for Conditions:
  
  run.seq.cur <- NULL
  if (length(avg.gf.loads.run)==1 & length(sample.size.vec.run)==1) {
    if (avg.gf.loads.run==1 & sample.size.vec.run==250) {
      run.seq.cur <- "_Seq1"
    } else if (avg.gf.loads.run==1 & sample.size.vec.run==100) {
      run.seq.cur <- "_Seq2"
    } else if (avg.gf.loads.run==1 & sample.size.vec.run==25) {
      run.seq.cur <- "_Seq3"
    } else if (avg.gf.loads.run==3 & sample.size.vec.run==250) {
      run.seq.cur <- "_Seq4"
    } else if (avg.gf.loads.run==3 & sample.size.vec.run==100) {
      run.seq.cur <- "_Seq5"
    } else if (avg.gf.loads.run==3 & sample.size.vec.run==25) {
      run.seq.cur <- "_Seq6"
    } else if (avg.gf.loads.run==2 & sample.size.vec.run==250) {
      run.seq.cur <- "_Seq7"
    } else if (avg.gf.loads.run==2 & sample.size.vec.run==100) {
      run.seq.cur <- "_Seq8"
    } else if (avg.gf.loads.run==2 & sample.size.vec.run==25) {
      run.seq.cur <- "_Seq9"
    }
  } else if (length(avg.gf.loads.run)==1 & length(sample.size.vec.run)>=3) {
    if (avg.gf.loads.run==1) {
      run.seq.cur <- "_Seq1,2,3"
    } else if (avg.gf.loads.run==3) {
      run.seq.cur <- "_Seq4,5,6"
    } else if (avg.gf.loads.run==2) {
      run.seq.cur <- "_Seq7,8,9"
    }
    if (length(sample.size.vec.run)==4) {
      run.seq.cur <- paste0(run.seq.cur, ",plN500")
    } 
  } else if (length(avg.gf.loads.run)==2 & length(sample.size.vec.run)>=3) {
    if (1 %in% avg.gf.loads.run & 2 %in% avg.gf.loads.run) {
      run.seq.cur <- "_Seq1t3,7t9"
    } else if (1 %in% avg.gf.loads.run & 3 %in% avg.gf.loads.run) {
      run.seq.cur <- "_Seq1t3,4t6"
    } else if (2 %in% avg.gf.loads.run & 3 %in% avg.gf.loads.run) {
      run.seq.cur <- "_Seq4t6,7t9"
    }
    if (length(sample.size.vec.run)==4) {
      run.seq.cur <- paste0(run.seq.cur, ",plN500")
    }
  } else if (length(avg.gf.loads.run)==1 & length(sample.size.vec.run)==1 &
             sample.size.vec.run==500) {
    if (avg.gf.loads.run==1) {
      run.seq.cur <- "_Seq10"
    } else if (avg.gf.loads.run==3) {
      run.seq.cur <- "_Seq11"
    } else if (avg.gf.loads.run==2) {
      run.seq.cur <- "_Seq12"
    }
  } else if (length(avg.gf.loads.run)==2 & length(sample.size.vec.run)==1 &
             sample.size.vec.run==500) {
    if (1 %in% avg.gf.loads.run & 2 %in% avg.gf.loads.run) {
      run.seq.cur <- "_Seq10,12"
    } else if (1 %in% avg.gf.loads.run & 3 %in% avg.gf.loads.run) {
      run.seq.cur <- "_Seq10,11"
    } else if (2 %in% avg.gf.loads.run & 3 %in% avg.gf.loads.run) {
      run.seq.cur <- "_Seq11,12"
    }
  } else {
    run.seq.cur <- "_SeqOth"
    warning("Current chosen conditions are not perfectly matching an existing numbered sequence!")
  }
  # run.seq.cur <- "_Seq4t9"
  message(paste0("Running: ",paste0(run.currently, collapse=", ")))
  message(paste0("Sequence: ",run.seq.cur))
  
  # Unidimensional:
  #   Run 1: avg.gf==1; sample.size==250 [all others run] -> 24 conditions
  #   Run 2: avg.gf==1; sample.size==100 [all others run] -> 24 conditions
  #   Run 3: avg.gf==1; sample.size==25 [all others run] -> 24 conditions
  #   Run 4: avg.gf==3; sample.size==250 [all others run] -> 24 conditions
  #   Run 5: avg.gf==3; sample.size==100 [all others run] -> 24 conditions
  #   Run 6: avg.gf==3; sample.size==25 [all others run] -> 24 conditions
  #   Run 7: avg.gf==2; sample.size==250 [all others run] -> 24 conditions
  #   Run 8: avg.gf==2; sample.size==100 [all others run] -> 24 conditions
  #   Run 9: avg.gf==2; sample.size==25 [all others run] -> 24 conditions
  
  # 2-Dimensional:
  #   Run 1: avg.gf==1; sample.size==250 [all others run] -> 24 conditions
  #   Run 2: avg.gf==1; sample.size==100 [all others run] -> 24 conditions
  #   Run 3: avg.gf==1; sample.size==25 [all others run] -> 24 conditions
  #   Run 4: avg.gf==3; sample.size==250 [all others run] -> 24 conditions
  #   Run 5: avg.gf==3; sample.size==100 [all others run] -> 24 conditions
  #   Run 6: avg.gf==3; sample.size==25 [all others run] -> 24 conditions
  #   Run 7: avg.gf==2; sample.size==250 [all others run] -> 24 conditions
  #   Run 8: avg.gf==2; sample.size==100 [all others run] -> 24 conditions
  #   Run 9: avg.gf==2; sample.size==25 [all others run] -> 24 conditions
  
  # LD = 1 or LD = 2:
  #   Run 1: avg.gf==1; sample.size==250 [all others run] -> 48 conditions
  #   Run 2: avg.gf==1; sample.size==100 [all others run] -> 48 conditions
  #   Run 3: avg.gf==1; sample.size==25 [all others run] -> 48 conditions
  #   Run 4: avg.gf==3; sample.size==250 [all others run] -> 48 conditions
  #   Run 5: avg.gf==3; sample.size==100 [all others run] -> 48 conditions
  #   Run 6: avg.gf==3; sample.size==25 [all others run] -> 48 conditions
  #   Run 7: avg.gf==2; sample.size==250 [all others run] -> 48 conditions
  #   Run 8: avg.gf==2; sample.size==100 [all others run] -> 48 conditions
  #   Run 9: avg.gf==2; sample.size==25 [all others run] -> 48 conditions
  
} else if (sub.conds==FALSE) {
  run.seq.cur <- ""
}


#   Check that run.currently matches uni.ld.df.run:
message(paste0("run.currently: ",paste0(run.currently,collapse=", "),
               "\nConditions running:"))
print(uni.ld.df.run)





# --- PARALLEL: Loop through conditions & replications: ------------------------ ####


RNGkind("L'Ecuyer-CMRG")
set.seed(seed=seed.start, kind="L'Ecuyer-CMRG")
cur.new.seed <- .Random.seed
# cur.new.seed <- seed.start 
temp.out.ct <- 0


# Set up output:
col.out <- c("Uni.Condition","LD.Load","Num.Items","Num.Resp.Cat","Sample.Size",
             "Avg.GF","Vary.GF",
             "LD.Load.Val","Model","Repl","Out.Repl","Cur.Seed",
             "Converged","RMSEA","RMSEA_5","RMSEA_95","TLI","CFI",
             "Num.Item.Pairs","Adj.Alpha.Level","Count.P.LD.G2","Count.P.LD.X2",
             "Count.AdjP.LD.G2","Count.AdjP.LD.X2","Count.Gt5LE10.LD.zG2",
             "Count.Gt5LE10.LD.zX2","Count.Gt10.LD.zG2","Count.Gt10.LD.zX2",
             "Items2.5.LD.zX2","Items2.5.LD.zG2","Items2.5.LD.P.X2","Items2.5.LD.P.G2",
             "Items2.5.LD.PltAdjAlpha.X2","Items2.5.LD.PltAdjAlpha.G2",
             "Omega","ECV.Overall","All.Items.Corr.Cat",
             "All.Items.Same.Cat","Num.Diff.Obs.Cat","Min.Obs.Cat","Max.Obs.Cat",
             "Item","Gen.Fac.Load.Val","Loads.On.Factor","Obs.Num.Resp.Cat",
             "S_X2","RMSEA.S_X2","p.S_X2",
             "a",paste0("b",1:(max(num.resp.cat.vec.run)-1)),"g","u","c") #adding the extra coef cols for pcm & rsm, jic
output.all <- as.data.frame(matrix(NA, nrow=0, ncol=length(col.out)))
colnames(output.all) <- col.out


# Secondary output, time per replication/condition:
col.time <- c("Uni.Condition","LD.Load","Num.Items","Num.Resp.Cat",
              "Sample.Size","Avg.GF","Vary.GF","LD.Load.Val","Num.Repl",
              "Out.Repl","No.Cores","Seed.Start.Cur","Time.AllRepl.Val",
              "Time.AllRepl.Unit","Time.PerRepl.Val","Time.PerRepl.Unit")
time.all <- as.data.frame(matrix(
  NA, nrow=(num.conds.all*num.out.repls), ncol=length(col.time)))
colnames(time.all) <- col.time
counter.time <- 0


# Create a list for the dataframes of each condition, rbinded over replications:
dat.list.all <- vector(mode="list", length=num.conds.all)
dat.temp <- NULL



# Loop through conditions & replications:
start.time.all <- Sys.time()
for (out.repl in 1:num.out.repls) {
  cond.count <- 0
  start.time.outrepl <- Sys.time()
  
  for (uni.ld.num in 1:nrow(uni.ld.df.run)) {
    ld.load <- uni.ld.df.run[uni.ld.num,"LD.Load"]
    uni.cond <- uni.ld.df.run[uni.ld.num,"Uni.Condition"]
    for (number.items in number.items.vec.run) {
      for (num.resp.cat in num.resp.cat.vec.run) {
        for (sample.size in sample.size.vec.run) {
          for (avg.gf in avg.gf.loads.run) {
            for (vary.gf in vary.gf.loads.run) {
              
              if (cond.count > 0) {
                cur.new.seed <- .Random.seed + add.seed.new
              }
              start.time.cond <- Sys.time()
              cond.count <- cond.count + 1
              counter.time <- counter.time + 1
              
              
              # Get numbers:
              number.items.num <- which(number.items.vec.run == number.items) 
              num.resp.cat.num <- which(num.resp.cat.vec.run == num.resp.cat)
              sample.size.num <- which(sample.size.vec.run == sample.size) 
              avg.gf.num <- which(avg.gf.loads.run == avg.gf)
              vary.gf.num <- which(vary.gf.loads.run == vary.gf)
              cond.num.all <- nrow(uni.ld.df.run) * length(number.items.vec.run) * 
                length(num.resp.cat.vec.run) * length(sample.size.vec.run) * 
                length(avg.gf.loads.run) * length(vary.gf.loads.run)
              
              # Put conditions info together into string (for message):
              conditions.str <- paste0("LD loading: ",ld.load," times avg.gf load;",
                                       " Unidimensional condition: ",uni.cond,
                                       " (looped together: ",uni.ld.num," / ",
                                       nrow(uni.ld.df.run),")",
                                       "\nNumber items: ",number.items,
                                       " (",number.items.num," / ",
                                       length(number.items.vec.run),")",
                                       "\nNumber response categories: ",num.resp.cat,
                                       " (",num.resp.cat.num," / ",
                                       length(num.resp.cat.vec.run),")",
                                       "\nSample size: ",sample.size,
                                       " (",sample.size.num," / ",
                                       length(sample.size.vec.run),")",
                                       "\nAverage GF loadings: ",avg.gf,
                                       " (",avg.gf.num," / ",length(avg.gf.loads.run),")",
                                       "\nVary GF loadings: ",vary.gf,
                                       " (",vary.gf.num," / ",length(vary.gf.loads.run),")")
              
              # Get time differences:
              cur.time.diff <- start.time.cond-start.time.all
              cur.time.outrepl.diff <- start.time.cond-start.time.outrepl
              # Estimate total times:
              est.tottime.outrepl <- cur.time.outrepl.diff / ((cond.count-1)/num.conds.all) #estimated total time for THIS out.repl
              if (out.repl == 1) { #If currently on first out.repl, estimate total time based on just this out.repl info
                est.tottime.all <- est.tottime.outrepl*num.out.repls
              } else if (out.repl > 1) { #If currently on 2nd or later out.repl, estimate total time considering prior out.repls
                time.prevrepls.diff <- start.time.outrepl-start.time.all #how long PRIOR out.repls (ALL) took
                if (((cond.count-1)/num.conds.all) >= (1/3)) { #If at least 1/3rd of the way through the current out.repl...
                  est.single.outrepl.time <- (time.prevrepls.diff + est.tottime.outrepl) / out.repl #...Estimate single out.repl time based on BOTH previous out.repls and current out.repl (so far)
                } else { #If less than 1/3rd of the way through the current out.repl...
                  est.single.outrepl.time <- (time.prevrepls.diff) / (out.repl-1) #...Estimate single out.repl time based on JUST the previous out.repls
                }
                est.tottime.all <- est.single.outrepl.time*num.out.repls #Estimate single out.repl time based on JUST the previous out.repls
              }
              # Estimate the time things will finish:
              est.endtime.all <- start.time.all + est.tottime.all
              est.endtime.outrepl <- start.time.outrepl + est.tottime.outrepl
              
              # Convert total times to different unit, if needed:
              est.tottime.all <- convert.unit.difftime(est.tottime.all)
              est.tottime.outrepl <- convert.unit.difftime(est.tottime.outrepl)
              
              # Print messages with status info:
              message(paste0("\n------------------------------\n",
                             "------------------------------\n\n",
                             "NOW STARTING:\n Running: ",paste0(run.currently, collapse=", "),
                             "\n Sequence, if relevant: ",run.seq.cur,
                             "\n Starting Seed: ",seed.start,
                             "\n Latest new seed: ",paste0(cur.new.seed,collapse=", "),
                             "\n",conditions.str,"   (STARTING Overall Condition: ",
                             cond.count," / ",num.conds.all,"; have COMPLETED: ",
                             (cond.count-1)," / ",num.conds.all,", or ",
                             round((cond.count-1)/num.conds.all,digits=2),")"))
              message(paste0("   Current out.repl: ",out.repl," / ",num.out.repls))
              message(paste0("\nTOTAL time so far is: ",
                             round(cur.time.diff,digits=2)," ",units(cur.time.diff),
                             " (Total Start date & time: ",start.time.all,")"))
              message(paste0("   Projected TOTAL time (ALL repls, ALL conds): ",
                             round(est.tottime.all,digits=2)," ",units(est.tottime.all)))
              message(paste0("   Projected TOTAL end time (ALL repls, ALL conds):\n      * ",
                             est.endtime.all," * "))
              message(paste0("   Time left: ",round(est.endtime.all-Sys.time(),digits=2),
                             " ",units(est.endtime.all-Sys.time())))
              message(paste0(" (Number of cores being used: ",no_cores," / ",detectCores(),")\n"))
              if (num.out.repls > 1) {
                message(paste0("Time so far on this Out.Repl is: ",
                               round(cur.time.outrepl.diff,digits=2)," ",
                               units(cur.time.outrepl.diff),
                               " (current out.repl start date & time: ",start.time.outrepl,")"))
                message(paste0("   Projected time, this out.repl: ",round(est.tottime.outrepl,digits=2),
                               " ",units(est.tottime.outrepl)))
                message(paste0("   Projected end time, this out.repl:\n      * ",est.endtime.outrepl," * "))
                message(paste0("   Time left: ",round(est.endtime.outrepl-Sys.time(),digits=2),
                               " ",units(est.endtime.outrepl-Sys.time())))
                message(paste0(" (Number of cores being used: ",no_cores," / ",detectCores(),")\n"))
              }
              
              
              
              
              # Get loadings and item intercepts: --- #
              
              #   GF loadings:
              min.cur <- avg.gf - (avg.gf*vary.gf)
              max.cur <- avg.gf + (avg.gf*vary.gf)
              gen.fac.loads <- seq(from=min.cur, to=max.cur, length.out=number.items)
              
              #   Item intercepts:
              item.ints <- intercepts[[as.character(num.resp.cat)]]
              
              # Local Dependence (if relevent):
              if (ld.load != "*0") {
                # Get the loading value for the local dependence:
                ld.load.val <- eval(rlang::parse_expr(paste0(avg.gf, ld.load)))
                # Create a vector of loadings and put that value in for the correct items:
                ld.fac.loads <- rep(0, length.out=number.items)
                ld.fac.loads[ld.items] <- ld.load.val
              }
              
              # Create item parameters df:
              if (ld.load == "*0") {
                item.params <- data.frame("Gen.Fac"=gen.fac.loads)
              } else if (ld.load != "*0") {
                item.params <- data.frame("Gen.Fac"=gen.fac.loads,
                                          "LD.Fac"=ld.fac.loads)
              }
              for (item.int.num in 1:length(item.ints)) {
                item.int.cur <- item.ints[item.int.num]
                new.int.name <- paste0("Int.",item.int.num)
                item.params[[new.int.name]] <- rep(item.ints[item.int.num], times=number.items)
              }
              
              
              
              
              # Parallel Loop: --- #
              
              
              # set up parallelization of the replications:
              cl <- makeCluster(no_cores)
              clusterSetRNGStream(cl, iseed=cur.new.seed)
              registerDoParallel(cl)
              getDoParWorkers()
              getDoParName()
              
              # now loop through replications:
              out.temp <- foreach(repl = 1:num.repls, #foreach works with the parallelization process
                                  # .combine='rbind',
                                  .errorhandling='pass',
                                  .packages=c('foreach','doParallel','parallel',
                                              'dplyr','tidyr','mirt','faux',
                                              'stats')
                                  # ,.verbose=TRUE
              )  %dopar% {
                
                grmvrasch.parallel.condsout.fntn(
                  repl = repl,
                  col.out = col.out,
                  models.run = models.run, 
                  uni.cond = uni.cond,
                  ld.load = ld.load,
                  sample.size = sample.size,
                  number.items = number.items, 
                  num.resp.cat = num.resp.cat, 
                  item.params = item.params,
                  num.cycles = num.cycles)
                
              } #end of looping through replications
              
              # Rbind Results:
              results.temp <- foreach(repl = 1:num.repls,
                                      .combine = 'rbind'
              ) %do% {
                out.temp[[repl]]$resultOut
              } #end of looping thru repls to save results
              
              # Rbind Data:
              if (num.repls > max.repl.dat.rbind) {
                max.dat.repl <- max.repl.dat.rbind
              } else { max.dat.repl <- num.repls }
              dat.temp <- foreach(repl = 1:max.dat.repl,
                                  .combine = 'rbind'
              ) %do% {
                out.temp[[repl]]$resultData
              } #end of looping thru repls to save datasets
              
              # Stop the parallelization:
              stopCluster(cl)
              
              # Add some condition information (that wasn't relevant for the parallel fntn):
              #   Results output:
              results.temp$Out.Repl <- out.repl
              results.temp$Avg.GF <- avg.gf
              results.temp$Vary.GF <- vary.gf
              if (ld.load == "*0") {
                results.temp$LD.Load.Val <- 0
              } else if (ld.load != "*0") {
                results.temp$LD.Load.Val <- ld.load.val
              }
              #   Data output:
              dat.temp$Out.Repl <- out.repl
              dat.temp$Avg.GF <- avg.gf
              dat.temp$Vary.GF <- vary.gf
              if (ld.load == "*0") {
                dat.temp$LD.Load.Val <- 0
              } else if (ld.load != "*0") {
                dat.temp$LD.Load.Val <- ld.load.val
              }
              
              # Rbind the current output in the overall output:
              output.all <- rbind(output.all, results.temp)
              
              # Put the rbinded datasets into the data list:
              dat.list.all[[cond.count]] <- dat.temp
              
              
              
              # Message of Time for This Condition:
              end.time.cond <- Sys.time()
              total.time.cond <- end.time.cond - start.time.cond
              time.per.repl.cur <- total.time.cond / num.repls
              #   Convert time per replication to a different unit, if needed:
              time.per.repl.cur <- convert.unit.difftime(time.per.repl.cur)
              #   Put message together and print to console:
              message.str <- paste0(
                "Total time for this one condition (",num.repls," replications): ",
                round(total.time.cond, digits=2)," ",units(total.time.cond),
                "\n  (Started: ",start.time.cond,"\n   Ended: ",end.time.cond,
                "\n   Time per 1 replication: ",round(time.per.repl.cur,digits=2),
                " ",units(time.per.repl.cur),")")
              if (!is.null(num.cycles)) {
                message.str <- paste0(
                  "Note: Running ",num.cycles," iterations to fit IRT model(s), NOT default.\n",
                  message.str)
              }
              message(message.str)
              
              
              # Add time info to the time output df:
              time.all[counter.time,"Uni.Condition"] <- uni.cond
              time.all[counter.time,"LD.Load"] <- ld.load
              time.all[counter.time,"Num.Items"] <- number.items
              time.all[counter.time,"Num.Resp.Cat"] <- num.resp.cat
              time.all[counter.time,"Sample.Size"] <- sample.size
              time.all[counter.time,"Avg.GF"] <- avg.gf
              time.all[counter.time,"Vary.GF"] <- vary.gf
              if (ld.load == "*0") { time.all[counter.time,"LD.Load.Val"] <- 0
              } else if (ld.load != "*0") {
                time.all[counter.time,"LD.Load.Val"] <- ld.load.val }
              time.all[counter.time,"Num.Repl"] <- num.repls
              time.all[counter.time,"Out.Repl"] <- out.repl
              time.all[counter.time,"Seed.Start.Cur"] <- paste0(cur.new.seed, collapse=",")
              time.all[counter.time,"No.Cores"] <- no_cores
              time.all[counter.time,"Time.AllRepl.Val"] <- total.time.cond
              time.all[counter.time,"Time.AllRepl.Unit"] <- units(total.time.cond)
              time.all[counter.time,"Time.PerRepl.Val"] <- time.per.repl.cur
              time.all[counter.time,"Time.PerRepl.Unit"] <- units(time.per.repl.cur)
              
              
              
            } #end loop through Vary GF loadings
            
            # Get the temporary count for the output:
            # if (length(unique(vary.gf.loads.run)) > 1) {
            #   temp.out.ct <- which(vary.gf.loads.run == vary.gf)
            #   temp.out.ct.str <- toString(temp.out.ct)
            # } else 
            if (length(unique(number.items.vec.run)) > 1) {
              temp.out.ct <- which(number.items.vec.run == number.items) 
              temp.out.ct.str <- toString(temp.out.ct)
            } else if (length(unique(sample.size.vec.run)) > 1) {
              temp.out.ct <- which(sample.size.vec.run == sample.size) 
              temp.out.ct.str <- toString(temp.out.ct)
            } else if (length(unique(avg.gf.loads.run)) > 1) {
              temp.out.ct <- which(avg.gf.loads.run == avg.gf)
              temp.out.ct.str <- toString(temp.out.ct)
            } else if (length(unique(num.resp.cat.vec.run)) > 1) {
              temp.out.ct <- which(num.resp.cat.vec.run == num.resp.cat)
              temp.out.ct.str <- toString(temp.out.ct)
            } else {
              temp.out.ct <- temp.out.ct + 1
              temp.out.ct.str <- toString(temp.out.ct)
              if (temp.out.ct >= 5) {
                temp.out.ct <- 0
              }
            }
            
            # Get filename info (use function):
            cond.vec <- c("Uni.Condition","LD.Load","Model","Num.Items",
                          "Num.Resp.Cat","Sample.Size","Avg.GF","Vary.GF")
            nicknames.vec <- c("DimC","LD","Md","NItem","NRsp","N",
                               "AvgGF","VrGF")
            print.min.max.vec <- FALSE # c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
            filename.part.str <- sims.filename.fntn(
              dataframe=output.all, 
              vec.of.cond.vars=cond.vec,
              nicknames.cond.vars.vec=nicknames.vec,
              print.min.max=print.min.max.vec )
            filename.part.str <- gsub(filename.part.str, pattern="*", 
                                      replacement="", fixed=T)
            
            # Put together the total name of the output as a string:
            temp.out.name.str <- paste0(
              object.name,"TEMP",temp.out.ct.str,"GRMvRasch",
              paste0(run.currently,collapse=","),run.seq.cur,
              "_",filename.part.str,"_MxRepls",
              max(output.all$Repl)*max(output.all$Out.Repl),"_Sd",seed.start,
              "_",n.cyc.str,Sys.Date(),".rds")
            
            # Save the temporary output (save all output so far):
            saveRDS(output.all, file=temp.out.name.str)
            
            # Ensure temp.out.name.str doesn't get held over:
            temp.out.name.str <- NULL
            
            # Save temp datasets list: 
            saveRDS(object=dat.list.all, paste0(
              data.name,"TEMP",temp.out.ct.str,"GRMvRaschDat",
              paste0(run.currently, collapse=","),run.seq.cur,
              "_",filename.part.str,"_DatRpls",max.repl.dat.rbind,"_Sd",seed.start,
              "_",Sys.Date(),".rds"))
            
            #   Re-null the temporary dataframe:
            dat.temp <- NULL
            
            
            
          } #end loop through Avg GF loadings
        } #end loop through sample sizes
      } #end loop through number of response categories
    } #end loop through number of items
  } #end loop through LD loadings & unidimensionality conditions
  
} #end loop through out.repl

end.time.all <- Sys.time()
total.time.all <- end.time.all - start.time.all

time.per.condcombo.repl <- total.time.all / (cond.num.all*num.repls*num.out.repls)
time.per.condcombo.repl <- convert.unit.difftime(time.per.condcombo.repl)
if (min(output.all$Num.Resp.Cat) != max(output.all$Num.Resp.Cat)) { 
  n.rspcat.str.notes <- paste0(min(output.all$Num.Resp.Cat),", ",max(output.all$Num.Resp.Cat))
} else { n.rspcat.str.notes <- min(output.all$Num.Resp.Cat) }
notes.message <- paste0("Total time to run ALL conditions: ",
                        round(total.time.all, digits=2)," ",
                        units(total.time.all)," (number of cores used: ",
                        no_cores,")\nNumber of total condition combinations ",
                        "was: ",cond.num.all,"\nNumber of total replications ",
                        "per condition was: ",(num.repls*num.out.repls),
                        "\n(That is ",time.per.condcombo.repl," ",
                        units(time.per.condcombo.repl)," per single ",
                        "replication per single condition combination.)",
                        "\nNumber of response categories conditions included: ",
                        n.rspcat.str.notes,".\nStarting seed: ",
                        seed.start,"\nDate finished running: ",Sys.Date())
if (!is.null(num.cycles)) {
  notes.message <- paste0(
    "Note: Running ",num.cycles," iterations to fit IRT model(s), NOT default.\n",
    notes.message)
}
output.all$Notes <- NA
output.all[1,"Notes"] <- notes.message


# Get filename info (use function):
vec.of.cond.vars <- c("Uni.Condition","LD.Load","Num.Items",
                      "Num.Resp.Cat","Sample.Size","Avg.GF","Vary.GF")
nicknames.cond.vars.vec <- c("DimC","LD","NItm","NRsp","N",
                             "AvgGF","VrGF")
print.min.max.vec <- FALSE # c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE)
filename.part.str <- sims.filename.fntn(
  dataframe = output.all, 
  vec.of.cond.vars = vec.of.cond.vars,
  nicknames.cond.vars.vec = nicknames.cond.vars.vec,
  print.min.max = print.min.max.vec)
filename.part.str <- gsub(filename.part.str, pattern="*", replacement="", fixed=T)
#   Add model info:
mod.info <- paste0(length(unique(output.all$Model)),"Mod")
filename.part.str <- paste0(mod.info,"_",filename.part.str)

# Put together the total name of the output as a string:
out.name.all.str <- paste0(object.name,"OutGRMvRasch",
                           paste0(run.currently,collapse=","),run.seq.cur,"_",
                           filename.part.str,"_MxRepls",
                           max(output.all$Repl)*max(output.all$Out.Repl),
                           "_Sd",seed.start,"_",n.cyc.str,Sys.Date(),".rds")

# Save ALL the output:
saveRDS(output.all, file=out.name.all.str)

# Save the datasets list:
saveRDS(object=dat.list.all, paste0(
  data.name,"GRMvRasch",paste0(run.currently,collapse=","),run.seq.cur,"_",
  filename.part.str,"_DatRpls",max.repl.dat.rbind,"_Sd",seed.start,"_",
  Sys.Date(),".rds"))


# Save the TIME output:
saveRDS(time.all, paste0(
  object.name,"TimeGRMvRasch",paste0(run.currently,collapse=","),run.seq.cur,
  "_",filename.part.str,"_MxRpls",max(output.all$Repl)*max(output.all$Out.Repl),
  "_Sd",seed.start,"_",n.cyc.str,Sys.Date(),".rds"))


# Print message to console of how long total it took to run, as well as
#   the number of cores used, the total condition combinations, and the 
#   replications per condition:
message(notes.message)



# Check actual replications:
temp.out.sum <- output.all %>% 
  dplyr::group_by(Uni.Condition, LD.Load, Num.Items, Num.Resp.Cat, Sample.Size, 
                  Avg.GF, Vary.GF) %>%
  summarise(Max.Repl=max(Repl, na.rm=T),
            Max.Out.Repl=max(Out.Repl, na.rm=T))
temp.out.sum <- as.data.frame(temp.out.sum)
temp.out.sum <- dplyr::mutate(temp.out.sum, Max.All.Repl = Max.Repl * Max.Out.Repl)
min.all.repl <- min(temp.out.sum$Max.All.Repl, na.rm=T)
max.all.repl <- max(temp.out.sum$Max.All.Repl, na.rm=T)

#   Add info to notes:
output.all[1,"Notes"] <- paste0(output.all[1,"Notes"],"\n",
                                "Replications range from: ",min.all.repl,
                                " to ",max.all.repl)

#   Re-save the output:
#     Put together the string info using replications:
if (min.all.repl != max.all.repl) {
  repl.str.temp <- paste0(min.all.repl,",",max.all.repl)
} else { repl.str.temp <- min.all.repl }
#     Put together the entire name of the output file:
out.name.all.str <- paste0(object.name,"OutGRMvRasch",
                           paste0(run.currently,collapse=","),run.seq.cur,"_",
                           filename.part.str,"_Repls",repl.str.temp,
                           "_Sd",seed.start,"_",n.cyc.str,Sys.Date(),".rds")
#     Save ALL the output, using the replication info:
saveRDS(output.all, file=out.name.all.str)




