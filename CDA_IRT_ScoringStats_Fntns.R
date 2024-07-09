# CD Ayasse
# 9 Dec 2022
# Functions for calculation of Omega, OmegaH, ECV, and Ratio for IRT Models



# Omega: ----------------------------------------------------------------------- ####

omega.fntn <- function(irt.model) {
  
  # save the summary of the model:
  sum.mod <- summary(irt.model, verbose=F)
  
  # create an empty list to save things into:
  var.list <- vector(mode="list", length=ncol(sum.mod$rotF))
  
  # loop through and save for each factor:
  for (i in 1:ncol(sum.mod$rotF)) {
    var.list[i] <- sum(sum.mod$rotF[,i])^2
    names(var.list)[[i]] <- colnames(sum.mod$rotF)[i]
  }
  
  # save the error:
  sum.h2 <- sum(1-sum.mod$h2)
  
  # put the equation together as a string:
  equation.str <- paste0("(", paste0(var.list, collapse=" + "), ")", "/",
                         "(", paste0(var.list, collapse=" + "), " + ",sum.h2,")")
  
  # evaluate the string/calculate omega:
  omega <- eval(parse(text=equation.str))
  
  # return omega:
  return(omega)
  
}

omega.fntn.alt <- function(sum.mod.rotF, sum.mod.h2) {
  
  # create an empty list to save things into:
  var.list <- vector(mode="list", length=ncol(sum.mod.rotF))
  
  # loop through and save for each factor:
  for (i in 1:ncol(sum.mod.rotF)) {
    var.list[i] <- sum(sum.mod.rotF[,i])^2
    names(var.list)[[i]] <- colnames(sum.mod.rotF)[i]
  }
  
  # save the error:
  sum.h2 <- sum(1-sum.mod.h2)
  
  # put the equation together as a string:
  equation.str <- paste0("(", paste0(var.list, collapse=" + "), ")", "/",
                         "(", paste0(var.list, collapse=" + "), " + ",sum.h2,")")
  
  # evaluate the string/calculate omega:
  omega <- eval(parse(text=equation.str))
  
  # return omega:
  return(omega)
  
}


# OmegaH: ---------------------------------------------------------------------- ####

omega.h.fntn <- function(irt.model) {
  
  # save a summary of the model:
  sum.mod <- summary(irt.model, verbose=F)
  
  # create an empty list to save things into:
  factor.omegah.list <- vector(mode="list", length=ncol(sum.mod$rotF))
  
  for (factor in 1:ncol(sum.mod$rotF)) {
    
    # create an empty list to save things into:
    var.list <- vector(mode="list", length=ncol(sum.mod$rotF))
    
    # loop through and save those:
    for (i in 1:ncol(sum.mod$rotF)) {
      var.list[i] <- sum(sum.mod$rotF[,i])^2
      names(var.list)[[i]] <- colnames(sum.mod$rotF)[i]
    }
    
    # calculate error:
    sum.h2 <- sum(1-sum.mod$h2)
    
    # put the equation string together:
    equation.str <- paste0("(", var.list[[factor]], ")", "/",
                           "(", paste0(var.list, collapse=" + "), " + ",sum.h2,")")
    
    # evaluate equation/calculate omega H:
    omega.h <- eval(parse(text=equation.str))
    
    # save in the list:
    factor.omegah.list[[factor]] <- omega.h
    names(factor.omegah.list)[[factor]] <- colnames(sum.mod$rotF)[factor]
    
  }
  
  # return omega H list:
  return(factor.omegah.list)
  
}


omega.h.fntn.alt <- function(sum.mod.rotF, sum.mod.h2) {
  
  # create an empty list to save things into:
  factor.omegah.list <- vector(mode="list", length=ncol(sum.mod.rotF))
  
  for (factor in 1:ncol(sum.mod.rotF)) {
    
    # create an empty list to save things into:
    var.list <- vector(mode="list", length=ncol(sum.mod.rotF))
    
    # loop through and save those:
    for (i in 1:ncol(sum.mod.rotF)) {
      var.list[i] <- sum(sum.mod.rotF[,i])^2
      names(var.list)[[i]] <- colnames(sum.mod.rotF)[i]
    }
    
    # calculate error:
    sum.h2 <- sum(1-sum.mod.h2)
    
    # put the equation string together:
    equation.str <- paste0("(", var.list[[factor]], ")", "/",
                           "(", paste0(var.list, collapse=" + "), " + ",sum.h2,")")
    
    # evaluate equation/calculate omega H:
    omega.h <- eval(parse(text=equation.str))
    
    # save in the list:
    factor.omegah.list[[factor]] <- omega.h
    names(factor.omegah.list)[[factor]] <- colnames(sum.mod.rotF)[factor]
    
  }
  
  # return omega H list:
  return(factor.omegah.list)
  
}


# ECV: ------------------------------------------------------------------------- ####

ecv.fntn <- function(irt.model) {
  
  # save summary of the model:
  sum.mod <- summary(irt.model, verbose=F)
  
  # create an empty list to store things in:
  var.list <- vector(mode="list", length=ncol(sum.mod$rotF))
  
  # loop through and calculate:
  for (i in 1:ncol(sum.mod$rotF)) {
    var.list[i] <- sum(sum.mod$rotF[,i]^2)
    names(var.list)[[i]] <- colnames(sum.mod$rotF)[i]
  }
  
  # put together the equation as a string:
  equation.str <- paste0("(", var.list[[1]], ")", "/",
                         "(", paste0(var.list, collapse=" + "), ")")
  
  # evaluate the equation/calculate ECV, overall:
  ecv.overall <- eval(parse(text=equation.str))
  
  # calculate ECV for each item:
  ecv.i <- (sum.mod$rotF^2)[ ,1] / rowSums(sum.mod$rotF^2)
  
  # put those together:
  ecv.all <- c(ecv.overall, ecv.i)
  names(ecv.all)[1] <- "Overall"
  
  
  # return ECV:
  return(ecv.all)
  
}

ecv.fntn.alt <- function(sum.mod.rotF) {
  
  # create an empty list to store things in:
  var.list <- vector(mode="list", length=ncol(sum.mod.rotF))
  
  # loop through and calculate:
  for (i in 1:ncol(sum.mod.rotF)) {
    var.list[i] <- sum(sum.mod.rotF[,i]^2)
    names(var.list)[[i]] <- colnames(sum.mod.rotF)[i]
  }
  
  # put together the equation as a string:
  equation.str <- paste0("(", var.list[[1]], ")", "/",
                         "(", paste0(var.list, collapse=" + "), ")")
  
  # evaluate the equation/calculate ECV, overall:
  ecv.overall <- eval(parse(text=equation.str))
  
  # calculate ECV for each item:
  ecv.i <- (sum.mod.rotF^2)[ ,1] / rowSums(sum.mod.rotF^2)
  
  # put those together:
  ecv.all <- c(ecv.overall, ecv.i)
  names(ecv.all)[1] <- "Overall"
  
  
  # return ECV:
  return(ecv.all)
  
}


# OmegaH/omega ratio: ---------------------------------------------------------- ####

omega.ratio.fntn <- function(omega.h.gen=NULL, omega.gen=NULL, irt.model=NULL) {
  # Note: omega.h.gen is the omega H for the general factor
  # (omega.gen is one number, for the model)
  
  
  if (nargs()==1 & typeof(omega.h.gen)=="S4") {
    irt.model <- omega.h.gen
  }
  if (missing(irt.model) & (missing(omega.h.gen)|missing(omega.gen))) {
    stop("Must supply EITHER the irt.model, or must supply BOTH omega.h.gen AND omega.gen")
  }

  # If inputted the model:
  if (!missing(irt.model)) { #if irt.model was supplied...
    if (typeof(irt.model)=="S4") { #...and is the correct type
      # save a summary of the model:
      sum.mod <- summary(irt.model, verbose=F)

      # Calc Omega H from the model (just general/first factor):
      #     create an empty list to save things into:
      var.list <- vector(mode="list", length=ncol(sum.mod$rotF))
      #     loop through and save those:
      for (i in 1:ncol(sum.mod$rotF)) {
        var.list[i] <- sum(sum.mod$rotF[,i])^2
        names(var.list)[[i]] <- colnames(sum.mod$rotF)[i]
      }
      #     calculate error:
      sum.h2 <- sum(1-sum.mod$h2)
      #     put the equation string together:
      equation.str <- paste0("(", var.list[[1]], ")", "/", #FIRST/GENERAL FACTOR ONLY
                             "(", paste0(var.list, collapse=" + "), " + ",sum.h2,")")
      #     evaluate equation/calculate omega.gen H:
      omega.h.gen <- eval(parse(text=equation.str))
      
      # Calc Omega from the model (just general/first factor):
      #     create an empty list to save things into:
      var.list <- vector(mode="list", length=ncol(sum.mod$rotF))
      #     loop through and save for each factor:
      for (i in 1:ncol(sum.mod$rotF)) {
        var.list[i] <- sum(sum.mod$rotF[,i])^2
        names(var.list)[[i]] <- colnames(sum.mod$rotF)[i]
      }
      #     save the error:
      sum.h2 <- sum(1-sum.mod$h2)
      #     put the equation together as a string:
      equation.str <- paste0("(", paste0(var.list, collapse=" + "), ")", "/",
                             "(", paste0(var.list, collapse=" + "), " + ",sum.h2,")")
      #     evaluate the string/calculate omega.gen:
      omega.gen <- eval(parse(text=equation.str))
      
    } #end of if the model is the correct type
  } #end of if provided the full model
  
  gen.fac.explained.var <- omega.h.gen / omega.gen
  res.fac.explained.var <- 1 - (omega.h.gen / omega.gen)
  
  out.list <- list(gen.fac.explained.var, res.fac.explained.var)
  names(out.list) <- c("Proportion of the total variance explained by the general factor",
                       "Proportion of the total variance explained by the residual/specific factors")
  
  return(out.list)
  
}




