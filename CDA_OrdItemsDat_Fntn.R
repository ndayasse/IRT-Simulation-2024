# C.D. Ayasse
# 4 Jan 2023
# Functions for creation of ordinal item data

# Some updates:
#   14 Nov 2023: changed object name "latent.vars" to "params.latent.vars"
#                 (intended to reduce confusion)
#   30 Nov 2023: enable custom item names to be supplied
#   08 Dec 2023: allow an option where the inputted single latent var can be 
#                 the threshold used instead of producing random uniform


# Will create a dataset of ordinal items:
#   - Can vary the number of either latent variables or the number of response 
#       categories for the items.
#   - Can provide a dataframe with latent variables (latent.vars.df), OR can 
#       provide information for the function to compute latent variables within.
#   - MUST provide EITHER parameters object or BOTH params.latent.vars AND 
#       params.intercepts objects. Note that if parameters object is supplied, 
#       will default to using that over the other two options.
# 
#   "parameters" object is assumed to be set up such that each item has its own row,
#       the first X columns have loadings for each of the X latent variables,
#       and the last Y columns have the Y params.intercepts (Y is one less than the 
#       number of response categories)
#   "params.latent.vars" object is assumed to have one row per each item with columns 
#       that are all the latent variables to be used.
#   "params.intercepts" object is assumed to have one row per each item with columns
#       that contain each of the params.intercepts (one less than the number of response categories).
# 
#   "latent.vars.df" default is NULL, in which case function will compute the latent vars.
#       If provided, would be a dataframe structured as one row per subject with 
#       each column being a latent variable.
#   "latent.var.threshold" default is FALSE, in which case function will create 
#       a relationship between item(s) and inputted or within-function
#       simulated latent variables, but will also simulate a random uniform
#       variable to use as the threshold.
#       If TRUE, will skip the step of simulating a random uniform variable and
#       and instead use the latent variable directly as the threshold.
#   "latent.var.thresh.name" only used if @latent.var.threshold is TRUE;
#       default is NULL, in which case the function will use the left-most 
#       latent variable column as the threshold; can provide the name of the
#       desired column as a string and it will use that variable in the latent
#       variable dataframe as the threshold.
#   "ind.categ.latent.vars" is the index or indices (column number) for latent variables that should be categorical
#   "num.categ.latent.vars" is the number of category options for the categorical vars:
#       if one number, that will be applied across all categ vars; 
#       otherwise must provide a vector of numbers for each, i.e.,
#       length(num.categ.latent.vars)==length(ind.categ.latent.vars)
#       If "num.categ.latent.vars" is NOT provided but "ind.categ.latent.vars" is,
#         will assume binary/2 categories.
#   "prob.categ.latent.vars" provides the probability weights for categories, if desired.
#       If not provided, probabilities will all be assumed to be evenly split across categories.
#       If no num.categ.latent.vars is provided, then must just provide a vector of 2 numbers (that add up to 1).
#       Otherwise, if num.categ.latent.vars is provided, 
#         MUST be a list of length==length(ind.categ.latent.vars),
#         And each element within the list MUST be a vector of 
#         probabilities (that add up to 1) such that their length matches the 
#         corresponding number of categories provided in num.categ.latent.vars
#         for that latent variable.
#   "corr.latent.val" default is NULL, in which case latent variables are NOT correlated
#       Can be single number, vars\*vars matrix, vars\*vars vector; passed 
#           directlyr to faux::rnorm_multi ("r" in that function)
# 
#   "subj.id.name" default is NULL, only relevant if "latent.vars.df" is supplied.
#       If supplied, is the name of the column in "latent.vars.df" dataframe
#       that gives the subject ID. If not provided, will assume all columns
#       are latent variables and no subject ID was provided, and will create one.
#   "item.name.vec" default is NULL; can be supplied as a vector of strings
#       that will be used to name the simulated items -- length MUST match
#       number.items

ord.items.dat <- function(N, number.items, number.resp.cat,
                          parameters=NULL, params.latent.vars=NULL, params.intercepts=NULL,
                          ind.categ.latent.vars=NULL, num.categ.latent.vars=NULL,
                          prob.categ.latent.vars=NULL, corr.latent.val=NULL,
                          latent.var.mean=0, latent.var.sd=1, latent.vars.df=NULL,
                          latent.var.threshold=FALSE, latent.var.thresh.name=NULL,
                          subj.id.name=NULL, item.name.vec=NULL) {
  
  library(faux)
  
  
  # re-name the columns of the parameters object(s) to make using it easier:
  
  if (!is.null(parameters)) {
    
    if (is.null(corr.latent.val) | is.null(ind.categ.latent.vars)) { #if there are NO latent var correlations OR no categ latent vars
      
      int.name.vec <- paste0("Int.",1:(number.resp.cat-1))
      latent.var.num <- ncol(parameters) - length(int.name.vec)
      
      orig.latent.var.names <- colnames(parameters)[1:latent.var.num]
      latent.var.name.vec <- paste0("Latent.Var.",1:latent.var.num)
      colnames(parameters) <- c(latent.var.name.vec, int.name.vec)
      
    } else if (!is.null(corr.latent.val) & !is.null(ind.categ.latent.vars)) { #if there ARE latent var correlations AND there are categorical vars:
      
      # create the names of the params.intercepts, which also gives the number:
      int.name.vec <- paste0("Int.",1:(number.resp.cat-1))
      
      # get the number of latent variables:
      latent.var.num <- ncol(parameters) - length(int.name.vec)
      # vector of the indices of all the latent variables within "parameters"
      #       (since parameters has the latent vars first, params.intercepts last):
      ind.all.latent.vars <- c(1:latent.var.num) 
      # vector of the indices of the continuous latent variables:
      ind.cont.latent.vars <- ind.all.latent.vars[!(ind.all.latent.vars %in% 
                                                      ind.categ.latent.vars)]
      # vector of the indices of the intercepts:
      ind.ints <- c((latent.var.num+1):ncol(parameters))
      # re-order columns of parameters so that continuous latent vars are first,
      #       categorical latent vars are next, and intercepts are last:
      parameters <- parameters[ , c(ind.cont.latent.vars, ind.categ.latent.vars,
                                    ind.ints)]
      
      # save the original column names (in the NEW order):
      orig.latent.var.names <- colnames(parameters)[1:latent.var.num]
      
      # re-name parameters columns (temp), for easier working with in the function:
      latent.var.name.vec <- paste0("Latent.Var.",1:latent.var.num)
      colnames(parameters) <- c(latent.var.name.vec, int.name.vec)
      
    }
    
  } else if (!is.null(params.latent.vars)) { #if parameters are supplied as params.latent.vars and params.intercepts separately
    
    if (is.null(corr.latent.val) | is.null(ind.categ.latent.vars)) { #if there are NO latent var correlations OR not categorical latent vars
      
      orig.latent.var.names <- colnames(params.latent.vars)
      latent.var.name.vec <- paste0("Latent.Var.",1:ncol(params.latent.vars))
      colnames(params.latent.vars) <- c(latent.var.name.vec)
      latent.var.num <- ncol(params.latent.vars)
      
      if (!is.null(params.intercepts)) {
        if (ncol(params.intercepts)!=(number.resp.cat-1)) {
          warning(paste0("The provided number of response categories does not match the",
                         " number of columns in the intercepts object (ncol should be number.resp.cat-1)"))
        }
        
        # re-name params.intercepts as well, for easier use later:
        int.name.vec <- paste0("Int.",1:ncol(params.intercepts))
        colnames(params.intercepts) <- c(int.name.vec)
        
      } else if (is.null(params.intercepts)) {
        stop("Must provide EITHER parameters object OR BOTH params.latent.vars AND params.intercepts object")
      }
      
      parameters <- cbind(params.latent.vars, params.intercepts)
      
    } else if (!is.null(corr.latent.val) & !is.null(ind.categ.latent.vars)) { #if there ARE corrs b/w latent vars AND there are categ latent vars
      
      # get the number of latent variables:
      latent.var.num <- ncol(params.latent.vars)
      # vector of the indices of all the latent variables:
      ind.all.latent.vars <- c(1:ncol(params.latent.vars)) 
      # vector of the indices of the continuous latent variables:
      ind.cont.latent.vars <- ind.all.latent.vars[!(ind.all.latent.vars %in% 
                                                      ind.categ.latent.vars)]
      # re-order columns of params.latent.vars so that continuous latent vars are first,
      #       categorical latent vars are last:
      params.latent.vars <- params.latent.vars[ , c(ind.cont.latent.vars, ind.categ.latent.vars)]
      
      # save the original column names (in the NEW order):
      orig.latent.var.names <- colnames(params.latent.vars)
      
      # re-name params.latent.vars columns (temp), for easier working with in the function:
      latent.var.name.vec <- paste0("Latent.Var.",1:latent.var.num)
      colnames(params.latent.vars) <- c(latent.var.name.vec)
      
      if (!is.null(params.intercepts)) {
        if (ncol(params.intercepts)!=(number.resp.cat-1)) {
          warning(paste0("The provided number of response categories does not match the",
                         " number of columns in the params.intercepts object (ncol should be number.resp.cat-1)"))
        }
        
        # re-name intercepts as well, for easier use later:
        int.name.vec <- paste0("Int.",1:ncol(params.intercepts))
        colnames(params.intercepts) <- c(int.name.vec)
        
      } else if (is.null(params.intercepts)) {
        stop("Must provide EITHER parameters object OR BOTH params.latent.vars AND params.intercepts object")
      }
      
      parameters <- cbind(params.latent.vars, params.intercepts)
      
    } #end of whether there are/aren't latent var correlations or categ latent vars
    
  } #end of supplied parameters vs supplied params.latent.vars and params.intercepts
  # print(parameters)
  
  
  
  # Create a blank matrix to put item values into later: --- #
  if (is.null(latent.vars.df)) { #if did NOT supply a dataframe of latent variables
    
    temp.items <- as.data.frame(matrix(data=NA, nrow=N, ncol=number.items))
    colnames(temp.items) <- paste0("Item_",1:ncol(temp.items))
    #   Create a df starting with subject ID's:
    dat.ord.all <- data.frame(
      'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')),
                      each = 1, length.out = N),
      stringsAsFactors=F)
    #   Put those together:
    dat.ord.all <- cbind(dat.ord.all, temp.items)
    
  } else if (!is.null(latent.vars.df)) { #if DID supply a dataframe of latent variables
    
    if (is.null(subj.id.name)) {
      
      temp.items <- as.data.frame(matrix(data=NA, nrow=N, ncol=number.items))
      colnames(temp.items) <- paste0("Item_",1:ncol(temp.items))
      #   Create a df starting with subject ID's:
      dat.ord.all <- data.frame(
        'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')),
                        each = 1, length.out = N),
        stringsAsFactors=F)
      #   Put those together:
      dat.ord.all <- cbind(dat.ord.all, temp.items)
      
    } else if (!is.null(subj.id.name)) {
      
      dat.ord.all <- as.data.frame(matrix(data=NA, nrow=N, ncol=number.items))
      colnames(dat.ord.all) <- paste0("Item_",1:ncol(dat.ord.all))
      
    }
    
  }
  
  
  
  # Start computing the latent variables: --------- #
  
  if (is.null(latent.vars.df)) { #only if did NOT supply dataframe of latent variables
    
    temp.latent <- as.data.frame(matrix(NA, nrow=N, ncol=latent.var.num))
    colnames(temp.latent) <- latent.var.name.vec
    
    count.categ.latent.var <- 0 #initialize counting the categorical latent variables
    count.cont.latent.var <- 0 #initialize counting the continuous latent variables, for when correlating them
    
    if (is.null(ind.categ.latent.vars)) { #if no indices provided for categorical latent vars (assume all are continuous)
      
      if (is.null(corr.latent.val)) { #if no latent variables are meant to be correlated
        for (cur.latent.var.name in latent.var.name.vec) {
          temp.latent[[cur.latent.var.name]] <- rnorm(n=N, mean=latent.var.mean, 
                                                      sd=latent.var.sd)
        }
      } else if (!is.null(corr.latent.val)) {
        # print(paste0("Corr value is ",corr.latent.val))
        num.cont.latent.vars <- length(latent.var.name.vec)
        cont.latent.vars <- as.data.frame(
          faux::rnorm_multi(n=N, vars=num.cont.latent.vars,
                            mu=latent.var.mean, sd=latent.var.sd,
                            r=corr.latent.val,
                            empirical=F,
                            varnames=latent.var.name.vec))
        temp.latent[,c(1:num.cont.latent.vars)] <- cont.latent.vars
        # print(cor.test(temp.latent[,1], temp.latent[,2]))
      }
      
    } else if (!is.null(ind.categ.latent.vars)) { #if indices were provided for categorical latent variables
      
      # Provide warnings only once, up here (outside of loop through categ latent vars):
      #   If provided num.categ.latent.vars but is NOT length of 1 NOR does its 
      #     length equal the length of ind.categ.latent.vars.
      if (!is.null(num.categ.latent.vars)) {
        if (length(num.categ.latent.vars)!=1) {
          if (length(num.categ.latent.vars)!=length(ind.categ.latent.vars)) { 
            warning(paste0("If num.categ.latent.vars is provided, must be EITHER\nlength of 1 ",
                           "(and applied across all categorical variables)\nOR ",
                           "must have a length equal to the length of ind.categ.latent.vars.\n",
                           " The supplied num.categ.latent.vars does not meet these criteria;\n",
                           "function will revert to using the first element of num.categ.latent.vars",
                           " across all categorical variables."))
          }
        }
      }
      #   If user has NOT provided number of categories (so assuming binary) 
      #     but HAS provided probabilities, but the length of the 
      #     provided probabilities is not correct.
      if (is.null(num.categ.latent.vars)) {
        if (!is.null(prob.categ.latent.vars)) { 
          if (length(prob.categ.latent.vars)!=2) { 
            warning(paste0("If number of categories for the latent var is NOT provided, the assumption is binary.\n",
                           "User can either provide no probabilities (assume an even split),\nor can provide ",
                           "a vector of exactly two probabilities that add up to 1.\nYou provided the incorrect ",
                           "length of probabilities (either <2 or >2).\n",
                           "Function will default to assuming an even split until fixed."))
          }
        }
      }
      #   If the number of categories WAS provided, AND the user HAS provided
      #     probabilities, BUT the length of the provided probabilities list
      #     does NOT match the length of ind.categ.latent.vars.
      if (!is.null(num.categ.latent.vars)) { 
        if (!is.null(prob.categ.latent.vars)) { 
          if (length(prob.categ.latent.vars)!=length(ind.categ.latent.vars)) { 
            warning(paste0("If providing probabilities, the list length must",
                           " equal the length of ind.categ.latent.vars,\n",
                           "unless no numbers of categories are provided and",
                           " all variables are assumed binary.\nIn which case",
                           " a single vector of length 2 may be provided.\n",
                           "Since the length of the list is not appropriate,",
                           " probabilities will be assumed even across categories."))
          }
        }
      }
      
      if (is.null(corr.latent.val)) { #if no latent variables are meant to be correlated
        
        for (cur.latent.var.num in 1:length(latent.var.name.vec)) { #loop through all latent vars (continuous and categ)
          
          # save the name of the current latent variable:
          cur.latent.var.name <- latent.var.name.vec[cur.latent.var.num] 
          
          if (cur.latent.var.num %in% ind.categ.latent.vars) { #if current latent var is marked categorical
            
            # Increase the counter, to use for accessing num.categ.latent.vars and prob.categ.latent.vars:
            count.categ.latent.var <- count.categ.latent.var + 1
            
            if (is.null(num.categ.latent.vars)) { #if number of categories NOT provided, assume binary
              
              if (is.null(prob.categ.latent.vars)) { #if user has NOT provided the probabilities for categories
                
                temp.latent[[cur.latent.var.name]] <- sample(x=c(0,1), size=N, replace=T)
                
              } else if (!is.null(prob.categ.latent.vars)) { #if user HAS provided probabilities (still assuming binary)
                
                if (length(prob.categ.latent.vars)!=2) { #if the length is *NOT* correct/INcorrect
                  temp.latent[[cur.latent.var.name]] <- sample(x=c(0,1), size=N, replace=T) #assume even probabilities
                } else { #if the length IS Correct
                  temp.latent[[cur.latent.var.name]] <- sample(x=c(0,1), size=N, 
                                                               replace=T, prob=prob.categ.latent.vars)
                }
                
              } #end of if user HAS provided probabilities
              
            } else if (!is.null(num.categ.latent.vars)) { #if number of categories WAS provided
              
              # Get the number of categories for THIS latent var:
              if (length(num.categ.latent.vars)==1) {
                cur.num.categ <- num.categ.latent.vars
              } else if (length(num.categ.latent.vars)!=length(ind.categ.latent.vars)) { #if provided num.categ.latent.vars but is not length of 1 NOR length==index vec
                cur.num.categ <- num.categ.latent.vars[1]
              } else {
                cur.num.categ <- num.categ.latent.vars[count.categ.latent.var]
              }
              
              if (is.null(prob.categ.latent.vars)) { #if user has NOT provided the probabilities for categories
                
                sample.vec <- 0:(cur.num.categ-1)
                temp.latent[[cur.latent.var.name]] <- sample(x=sample.vec, size=N, replace=T)
                
              } else if (!is.null(prob.categ.latent.vars)) { #if user HAS provided probabilities
                
                if (length(prob.categ.latent.vars)==length(ind.categ.latent.vars)) { #must match the length of ind.categ.latent.vars
                  
                  # save the current vector of probabilities:
                  cur.prob.vec <- prob.categ.latent.vars[[count.categ.latent.var]]
                  # create the vector of category possibilities:
                  sample.vec <- 0:(cur.num.categ-1)
                  # create the latent variable:
                  temp.latent[[cur.latent.var.name]] <- sample(x=sample.vec, size=N, 
                                                               replace=T, prob=cur.prob.vec)
                  
                } else { #if length of prob does NOT match the length of ind.categ.latent.vars
                  
                  sample.vec <- 0:(cur.num.categ-1)
                  temp.latent[[cur.latent.var.name]] <- sample(x=sample.vec, size=N, replace=T)
                  
                } #end of if length prob matches/does not match length of ind.categ.latent.vars
                
              }
              
            } #end if number of categories was not/was provided
            
          } else { #if current latent var is NOT categorical
            
            for (cur.latent.var.name in latent.var.name.vec) {
              temp.latent[[cur.latent.var.name]] <- rnorm(n=N, mean=latent.var.mean, 
                                                          sd=latent.var.sd)
            }
            
          } #end of if the current latent variable is/is not a categorical one
          
        } #end of loop through latent variables
        
      } else if (!is.null(corr.latent.val)) { #if cont latent vars ARE correlated
        
        # save the number of continuous latent vars:
        num.cont.latent.vars <- length(ind.cont.latent.vars)
        
        # create that many latent variables that are correlated:
        cont.latent.vars <- as.data.frame(
          faux::rnorm_multi(n=N, vars=num.cont.latent.vars,
                            mu=latent.var.mean, sd=latent.var.sd,
                            r=corr.latent.val,
                            empirical=F,
                            varnames=latent.var.name.vec))
        #   add those onto the temp.latent dataframe:
        temp.latent[,c(1:num.cont.latent.vars)] <- cont.latent.vars
        
        for (cur.latent.var.num in (num.cont.latent.vars+1):length(latent.var.name.vec)) { #loop through all CATEG latent vars
          
          # save the name of the current latent variable:
          cur.latent.var.name <- latent.var.name.vec[cur.latent.var.num] 
          
          # Increase the counter, to use for accessing num.categ.latent.vars and prob.categ.latent.vars:
          count.categ.latent.var <- count.categ.latent.var + 1
          
          if (is.null(num.categ.latent.vars)) { #if number of categories NOT provided, assume binary
            
            if (is.null(prob.categ.latent.vars)) { #if user has NOT provided the probabilities for categories
              
              temp.latent[[cur.latent.var.name]] <- sample(x=c(0,1), size=N, replace=T)
              
            } else if (!is.null(prob.categ.latent.vars)) { #if user HAS provided probabilities (still assuming binary)
              
              if (length(prob.categ.latent.vars)!=2) { #if the length is *NOT* correct/INcorrect
                temp.latent[[cur.latent.var.name]] <- sample(x=c(0,1), size=N, replace=T) #assume even probabilities
              } else { #if the length IS Correct
                temp.latent[[cur.latent.var.name]] <- sample(x=c(0,1), size=N, 
                                                             replace=T, prob=prob.categ.latent.vars)
              }
              
            } #end of if user HAS provided probabilities
            
          } else if (!is.null(num.categ.latent.vars)) { #if number of categories WAS provided
            
            # Get the number of categories for THIS latent var:
            if (length(num.categ.latent.vars)==1) {
              cur.num.categ <- num.categ.latent.vars
            } else if (length(num.categ.latent.vars)!=length(ind.categ.latent.vars)) { #if provided num.categ.latent.vars but is not length of 1 NOR length==index vec
              cur.num.categ <- num.categ.latent.vars[1]
            } else {
              cur.num.categ <- num.categ.latent.vars[count.categ.latent.var]
            }
            
            if (is.null(prob.categ.latent.vars)) { #if user has NOT provided the probabilities for categories
              
              sample.vec <- 0:(cur.num.categ-1)
              temp.latent[[cur.latent.var.name]] <- sample(x=sample.vec, size=N, replace=T)
              
            } else if (!is.null(prob.categ.latent.vars)) { #if user HAS provided probabilities
              
              if (length(prob.categ.latent.vars)==length(ind.categ.latent.vars)) { #must match the length of ind.categ.latent.vars
                
                # save the current vector of probabilities:
                cur.prob.vec <- prob.categ.latent.vars[[count.categ.latent.var]]
                # create the vector of category possibilities:
                sample.vec <- 0:(cur.num.categ-1)
                # create the latent variable:
                temp.latent[[cur.latent.var.name]] <- sample(x=sample.vec, size=N, 
                                                             replace=T, prob=cur.prob.vec)
                
              } else { #if length of prob does NOT match the length of ind.categ.latent.vars
                
                sample.vec <- 0:(cur.num.categ-1)
                temp.latent[[cur.latent.var.name]] <- sample(x=sample.vec, size=N, replace=T)
                
              } #end of if length prob matches/does not match length of ind.categ.latent.vars
              
            }
            
          } #end if number of categories was not/was provided
          
        } #end of loop through latent variables
        
      } #end of if cont latent vars areN'T/ARE correlated
      
    } #end of if indices for categ latent vars were provided
    
    
    # If DID supply dataframe of latent variables:
  } else if (!is.null(latent.vars.df)) {
    
    orig.latent.var.names <- colnames(latent.vars.df)
    if (!is.null(subj.id.name)) {
      full.orig.latent.var.names <- orig.latent.var.names
      orig.latent.var.names <- 
        orig.latent.var.names[which(orig.latent.var.names!=subj.id.name)] 
    }
    latent.var.name.vec <- paste0("Latent.Var.",1:latent.var.num)
    if (!is.null(subj.id.name)) {
      col.num.subj <- which(full.orig.latent.var.names==subj.id.name)
      col.not.subj.nums <- which(full.orig.latent.var.names!=subj.id.name)
      if (col.num.subj!=1) {
        latent.vars.df <- latent.vars.df[,c(col.num.subj, col.not.subj.nums)]
      }
      colnames(latent.vars.df) <- c(subj.id.name, latent.var.name.vec)
    } else if (is.null(subj.id.name)) {
      colnames(latent.vars.df) <- latent.var.name.vec
    }
    
    temp.latent <- latent.vars.df
    
  } #end of if DID supply dataframe of latent vars
  
  
  
  # ACTUAL PROCESS OF SIMULATING THE ORDINAL ITEM(S): ------ #
  
  
  # Put together the overall df (subject ID, empty ordinal items, latent variables):
  dat.ord.all <- cbind(dat.ord.all, temp.latent)
  
  
  # Figure out what is being used for the threshold (a latent var or from runif):
  #   If using runif for threshold(s):
  if (!latent.var.threshold) { 
    
    # Use runif to simulate thresholds per person per item:
    thresholds <- as.data.frame(matrix(data=NA, nrow=N, ncol=number.items))
    colnames(thresholds) <- paste0("Threshold",1:number.items)
    for (thresh in 1:number.items) {
      cur.threshold.name <- paste0("Threshold",thresh)
      thresholds[[cur.threshold.name]] <- runif(n=N, min=0, max=1)
    }
    
    #   If using a latent variable, no column name provided:
  } else if (latent.var.threshold & is.null(latent.var.thresh.name)) {
    
    first.latent.var.name <- latent.var.name.vec[1]
    
    # Put the first latent variable in as the threshold for each item:
    thresholds <- as.data.frame(matrix(data=NA, nrow=N, ncol=number.items))
    colnames(thresholds) <- paste0("Threshold",1:number.items)
    for (thresh in 1:number.items) {
      cur.threshold.name <- paste0("Threshold",thresh)
      thresholds[[cur.threshold.name]] <- dat.ord.all[,first.latent.var.name]
    }
    
    #   If using a latent variable, with column name provided:
  } else if (latent.var.threshold & !is.null(latent.var.thresh.name)) {
    
    # Figure out what the correct latent variable is currently being called:
    latent.var.thresh.index <- which(orig.latent.var.names == latent.var.thresh.name)
    cur.latent.var.name <- latent.var.name.vec[latent.var.thresh.index]
    
    # Put THAT latent variable in as the threshold for each item:
    thresholds <- as.data.frame(matrix(data=NA, nrow=N, ncol=number.items))
    colnames(thresholds) <- paste0("Threshold",1:number.items)
    for (thresh in 1:number.items) {
      cur.threshold.name <- paste0("Threshold",thresh)
      thresholds[[cur.threshold.name]] <- dat.ord.all[,cur.latent.var.name]
    }
    
  }
  
  
  # item probabilities:
  item.probs <- as.data.frame(matrix(data=NA, nrow=N, ncol=number.items*(number.resp.cat-1)))
  prob.names <- paste0("Prob",1:(number.resp.cat-1))
  colnames(item.probs) <- paste0(rep(paste0("Item",1:number.items),
                                     each=(number.resp.cat-1)),
                                 "_",prob.names)
  
  # item prob differences:
  item.pdiffs <- as.data.frame(matrix(data=NA, nrow=N*number.items, 
                                      ncol=(number.resp.cat+1)))
  pdiff.names <- paste0("PDiff",1:number.resp.cat)
  colnames(item.pdiffs) <- c("Item",pdiff.names)
  item.pdiffs$Item <- rep(paste0("Item",1:number.items), 
                          length.out=nrow(item.pdiffs))
  
  
  # loop through items, put info in:
  for (item in 1:number.items) {
    
    # item probabilities:
    prob.names.vec <- paste0("Item",item,"_","Prob",1:(number.resp.cat-1))
    
    #   calculate the part where we multiply each latent var by its parameter and then add those cols together:
    store.new.vec <- matrix(0, nrow=N, ncol=1) #initialize to 0 so we can add onto it
    for (cur.latent.var.num in 1:latent.var.num) {
      
      cur.latent.var.name <- latent.var.name.vec[cur.latent.var.num]
      constant.param <- parameters[item,cur.latent.var.name]
      col.data <- dat.ord.all[,cur.latent.var.name]
      
      # calculate this current new vector:
      cur.new.vec <- constant.param*col.data
      
      # add this vector onto the previous ones:
      store.new.vec <- store.new.vec + cur.new.vec
    } #end of loop through latent vars
    
    #   calculate item probs:
    for (cur.prob.num in 1:length(prob.names.vec)) {
      cur.int.name <- paste0("Int.",cur.prob.num)
      item.probs[[prob.names.vec[cur.prob.num]]] <- 
        1/(1+exp(-1*(parameters[item,cur.int.name] + store.new.vec)))
    } #end of calc item probs
    
    
    
    # item probability differences:
    cur.item.name <- paste0("Item",item)
    pdiff.names.vec <- paste0("PDiff",1:number.resp.cat)
    
    for (cur.pdiff.num in 1:length(pdiff.names.vec)) {
      cur.pdiff.name <- pdiff.names.vec[cur.pdiff.num]
      
      if (cur.pdiff.num==1) { #if on the first one
        item.pdiffs[which(item.pdiffs[,"Item"]==cur.item.name),cur.pdiff.name] <- 
          1 - item.probs[[prob.names.vec[cur.pdiff.num]]]
        
      } else if (cur.pdiff.num==length(pdiff.names.vec)) { #if on the last one
        item.pdiffs[which(item.pdiffs[,"Item"]==cur.item.name),cur.pdiff.name] <- 
          item.probs[[prob.names.vec[cur.pdiff.num-1]]]
        
      } else { #if on one of the middle ones
        item.pdiffs[which(item.pdiffs[,"Item"]==cur.item.name),cur.pdiff.name] <- 
          item.probs[[prob.names.vec[cur.pdiff.num-1]]] - item.probs[[prob.names.vec[cur.pdiff.num]]]
      }
    } #end of loop through pdiff
    
    
    # create categorical variables:
    cur.item.name.dat <- paste0("Item_",item)
    cur.threshold.name <- paste0("Threshold",item)
    dat.ord.all[[cur.item.name.dat]] <- 1 #initialize the new item as "1" since 
    #       that is what we start with and then add on to based 
    #       on threshold comparison
    
    for (cur.pdiff.num in 1:length(pdiff.names.vec)) {
      cur.pdiff.name <- pdiff.names.vec[cur.pdiff.num]
      pdiff.nums.touse <- 1:cur.pdiff.num
      pdiff.names.touse <- paste0("PDiff",pdiff.nums.touse)
      # print(cur.pdiff.name)
      # print(pdiff.names.touse)
      
      # Get the part of the item.pdiffs dataset that is just THIS CURRENT item:
      item.pdiffs.curitem <- subset(item.pdiffs, Item==cur.item.name)
      
      # Now get just the PDiff column(s) we want:
      item.pdiffs.curitem.curcols <- dplyr::select(item.pdiffs.curitem, 
                                                   c(all_of(pdiff.names.touse)))
      
      # Now add those columns together to get the new column:
      cur.new.col <- rowSums(item.pdiffs.curitem.curcols)
      
      # Now calculate compared to thresholds, to add together later:
      if (cur.pdiff.num==1) { #if on the first one, do >=
        cur.thresh.compare <- thresholds[[cur.threshold.name]] >= cur.new.col
      } else { #if on any other, do just >
        cur.thresh.compare <- thresholds[[cur.threshold.name]] > cur.new.col
      }
      
      # Now add that info onto the current item:
      dat.ord.all[[cur.item.name.dat]] <- dat.ord.all[[cur.item.name.dat]] + cur.thresh.compare
      
    } #end of loop through cdiff
    
  } #end of loop through items
  
  # # Checking:
  # temp.list <- list(dat.ord.all, latent.var.num, latent.var.name.vec)
  # names(temp.list) <- c("Data","Number of Latent Variables","Latent Variable Names")
  # return (temp.list)
  
  
  # If item.name.vec was supplied, change the item names before returning the object:
  if (!is.null(item.name.vec)) {
    if (length(item.name.vec) != number.items) {
      warning("Length of item.name.vec MUST match number.items; ignoring item.name.vec")
    } else {
      colnames(dat.ord.all)[which(colnames(dat.ord.all) %in% 
                                    paste0("Item_",1:number.items))] <- item.name.vec
    }
  } 
  
  
  # create a version that is just the items for fitting models:
  if (is.null(subj.id.name)) {
    dat.ord.items <- dplyr::select(dat.ord.all, -c(USUBJID))
  } else if (!is.null(subj.id.name)) {
    dat.ord.items <- dplyr::select(dat.ord.all, -c(all_of(subj.id.name)))
  }
  dat.ord.items <- dplyr::select(dat.ord.items, -c(all_of(latent.var.name.vec)))

  # print(colnames(dat.ord.all))
  # print(colnames(dat.ord.items))


  # Put the original latent variable names back in the dataset to return:
  colnames(dat.ord.all)[(number.items+2):ncol(dat.ord.all)] <- orig.latent.var.names


  # return the data objects: the ordinal items, with and without the latent variables:
  data.list.out <- list(dat.ord.all, dat.ord.items)
  names(data.list.out) <- c("Data.Items.Latent","Data.Items.Only")
  return(data.list.out)
  
} #end function
