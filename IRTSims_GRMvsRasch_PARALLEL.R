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
source(paste0(genfntn.name,"CDA_IRT_ScoringStats_Fntns.R")) #custom function for bifactor/IRT scoring statistics  
source(paste0(genfntn.name,"CDA_OrdItemsDat_Fntn.R")) #custom function for creating dataset with ordinal items

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


# Parallel properties:
no_cores <- detectCores() -2 # use cores minus some of available
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

#   How many datasets to save
#     (If needed - To conserve memory/computational resources)
max.repl.dat.rbind <- num.repls

#   Randomization, parallelization, looping:
seed.start <- 07062024+1
add.seed.new <- 1717
num.repls <- 500 
num.out.repls <- 1 

#   Sample and items:
sample.size.vec.run <-  c(25,100,250,500) 
number.items.vec.run <- c(20, 12, 5) 
num.resp.cat.vec.run <- c(3, 7) 

#   Item parameters and intercepts:
#     Item intercepts depend on the number of response categories:
int.3 <- c(-1,1); int.5 <- c(-2,-0.75,0.75,2); int.7 <- c(-3,-2,-0.75,0.75,2,3)
intercepts <- list(int.3, int.5, int.7)
names(intercepts) <- c("3","5","7")
#     Vary the average General Factor (GF) loading and GF loading variability:
avg.gf.loads.run <- c(1, 2, 3) 
vary.gf.loads.run <- c(0, 0.125, 0.25, 0.5, 0.75) 


# Set up local dependence and unidimensionality conditions so that we test
#   each in the standard/unchanged/expected condition of the other:

#   Local Dependence between Items?:
ld.load.vec <- c("*0","*1.0","*2") # c("*0","*0.5","*1.0","*2")
ld.items <- c(2, 5)

#   Unidimensional?
uni.cond.vec <- c("Unidimensional","Two.Fac.Half.Items") 

#   Local Dependence and Unidimensionality, as we will run them:
uni.ld.df.run <- data.frame(
  "LD.Load"=c(ld.load.vec, rep(ld.load.vec[1], times=(length(uni.cond.vec)-1))),
  "Uni.Condition"=c(rep(uni.cond.vec[1], times=length(ld.load.vec)), uni.cond.vec[-1]))



# Number of conditions:
num.conds.all <- nrow(uni.ld.df.run) * length(sample.size.vec.run) * 
  length(number.items.vec.run) * length(num.resp.cat.vec.run) * 
  length(avg.gf.loads.run) * length(vary.gf.loads.run)
msg.all.conds <- paste0("Total number of condition-combinations is: ",num.conds.all)
message(msg.all.conds)




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




