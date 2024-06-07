# CD Ayasse

# Function to enable parallelization of Simulations for GRM vs Rasch


grmvrasch.parallel.condsout.fntn <- 
  function(col.out, repl, uni.cond=NULL, ld.load=NULL, 
           sample.size, number.items, num.resp.cat, item.params,
           models.run=c("GRM","PCM","GPCM"), #c("GRM","PCM","RSM"), 
           alpha=0.05, num.cycles=NULL) {
    
    
    
    # RNGkind("L'Ecuyer-CMRG") #INCLUDING THIS IN THE FUNCTION CAUSES PROBLEMS FOR DORNG! (works for the older parallel workaround)
    #   save current state of RNG:
    cur.seed = get(".Random.seed", .GlobalEnv)
    
    
    # Create output for THIS replication:
    output.temp <- as.data.frame(matrix(NA, nrow=length(models.run)*number.items, 
                                        ncol=length(col.out)))
    colnames(output.temp) <- col.out
    #   ^nrow=3*number.items because testing 3 models plus all item info/stats
    
    # Put basic info in output:
    if (!is.null(uni.cond)) {
      output.temp[,"Uni.Condition"] <- uni.cond
    }
    if (!is.null(ld.load)) {
      output.temp[,"LD.Load"] <- ld.load
    }
    output.temp[,"Num.Items"] <- number.items
    output.temp[,"Num.Resp.Cat"] <- num.resp.cat
    output.temp[,"Sample.Size"] <- sample.size
    output.temp[,"Repl"] <- repl
    output.temp[,"Cur.Seed"] <- paste0(cur.seed, collapse=",")
    
    
    # Get basic info that may not be provided:
    if (is.null(ld.load)) {
      ld.load <- "*0"
    }
    if (is.null(uni.cond)) {
      uni.cond <- "Unidimensional"
    }
    
    
    
    
    # --- Simulate starting df, latent variable, and ordinal items: ---------------- ####
    
    
    
    # Create starting dataframe: --- #
    
    #   Create a df with subject ID's:
    start.dat <- data.frame(
      'USUBJID' = rep(paste0('Subject_', formatC(1:sample.size, width=4, flag='0')),
                      each=1, length.out=sample.size),
      stringsAsFactors=F)
    
    #   If Unidimensional:
    if (uni.cond == "Unidimensional") {
      #     Create one General Factor:
      # latent.vars.df <- as.data.frame("Gen.Fac"=rnorm(n=sample.size, mean=0, sd=1))
      # colnames(latent.vars.df) <- "Gen.Fac"
      latent.vars.df <- as.data.frame(
        faux::rnorm_multi(n=sample.size, vars=1, mu=0, sd=1, r=0, empirical=F,
                          varnames=c("Gen.Fac")))
    } else if (uni.cond == "Two.Fac.Half.Items") {
      #     Create two Factors:
      latent.vars.df <- as.data.frame(
        faux::rnorm_multi(n=sample.size, vars=2, mu=0, sd=1, r=0, empirical=F, 
                          varnames=c("Gen.Fac","Fac.2")))
    } else if (uni.cond == "Unrelated") {
      #     Create one factor per item:
      latent.vars.df <- as.data.frame(
        faux::rnorm_multi(n=sample.size, vars=number.items, mu=0, sd=1, r=0, empirical=F, 
                          varnames=paste0("Fac.",1:number.items)))
    }
    
    #   ALSO ADD the LD/Nuisance factor IF RELEVANT:
    if (!is.null(ld.load) & ld.load!="*0") {
      latent.vars.df$LD.Fac <- rnorm(n=sample.size, mean=0, sd=1)
    }
    
    #   Add the factor(s) to the dataframe:
    start.dat <- cbind(start.dat, latent.vars.df)
    
    # #     Examine/check the general factor/latent variable briefly:
    # summary(latent.vars.df)
    # hist(latent.vars.df$Gen.Fac)
    
    
    
    
    # Create Ordinal Items: --- #
    
    #     First, shuffle the factor loadings:
    gen.fac.loads.notshuffled <- item.params$Gen.Fac
    gen.fac.loads <- sample(gen.fac.loads.notshuffled, replace=F, 
                            size=length(gen.fac.loads.notshuffled))
    item.params$Gen.Fac <- gen.fac.loads
    
    #     Change dimensionality, if relevant:
    if (uni.cond != "Unidimensional") {
      if (uni.cond == "Two.Fac.Half.Items") {
        fac1.seq <- seq(from=1, to=number.items, by=2)
        fac2.seq <- seq(from=2, to=number.items, by=2)
        item.params$Fac.2 <- 0
        item.params[fac2.seq,"Fac.2"] <- item.params[fac2.seq,"Gen.Fac"]
        item.params[fac2.seq,"Gen.Fac"] <- 0
        item.params <- item.params[,c(which(colnames(item.params)=="Gen.Fac"),
                                      which(colnames(item.params)=="Fac.2"),
                                      which(colnames(item.params)!="Gen.Fac" &
                                              colnames(item.params)!="Fac.2"))]
        # Label which factor each item loads on:
        label.fac.items <- data.frame("Item"=paste0("Item_",1:number.items),
                                      "Loads.On.Factor"=NA)
        label.fac.items[fac1.seq,"Loads.On.Factor"] <- "Gen.Fac"
        label.fac.items[fac2.seq,"Loads.On.Factor"] <- "Fac.2"
      } else if (uni.cond == "Unrelated") {
        temp.fac.params.mat <- matrix(0, nrow=number.items, ncol=number.items)
        diag(temp.fac.params.mat) <- gen.fac.loads
        temp.fac.params <- as.data.frame(temp.fac.params.mat)
        colnames(temp.fac.params) <- paste0("Fac.",1:number.items)
        item.params <- cbind(temp.fac.params, 
                             item.params[,which(grepl(colnames(item.params), 
                                                      pattern="Int.", fixed=T))])
        # Label which factor each item loads on:
        label.fac.items <- data.frame("Item"=paste0("Item_",1:number.items),
                                      "Loads.On.Factor"=paste0("Fac.",1:number.items))
      }
    } else if (uni.cond == "Unidimensional") {
      # Label which factor each item loads on:
      label.fac.items <- data.frame("Item"=paste0("Item_",1:number.items),
                                    "Loads.On.Factor"=rep("Gen.Fac", length.out=number.items))
    }
    
    #     Run function:
    dat.ord.list <- ord.items.dat(N=sample.size, number.items=number.items, 
                                  number.resp.cat=num.resp.cat,
                                  parameters=item.params,
                                  corr.latent.val=NULL, #don't want to specify the latent var correlations
                                  latent.vars.df=start.dat, #supply the latent variables to the function
                                  subj.id.name="USUBJID") #supply the subject ID var name so it doesn't make its own
    #     Save the two outputted datasets separately:
    dat.ord.items <- dat.ord.list$Data.Items.Only
    dat.ord.all <- dat.ord.list$Data.Items.Latent
    
    # #     Examine data/check briefly:
    # xtabs(~Item_1, dat.ord.items, addNA=T)
    # xtabs(~Item_2, dat.ord.items, addNA=T)
    # xtabs(~Item_3, dat.ord.items, addNA=T)
    # xtabs(~Item_4, dat.ord.items, addNA=T)
    # xtabs(~Item_5, dat.ord.items, addNA=T)
    # hist(dat.ord.items$Item_1)
    # hist(dat.ord.items$Item_2)
    # hist(dat.ord.items$Item_5)
    
    
    
    
    
    # --- Run (Unidimensional) Models, Get Fit info, Get Scoring Stats: ------------ ####
    
    
    #   Make sure item responses are a dataframe:
    dat.ord.items <- as.data.frame(dat.ord.items)
    
    
    #   Prep some info for output: 
    #     Get row indices for each model:
    row.indices <- seq(from=1, to=nrow(output.temp), by=number.items)
    if ("GRM" %in% models.run) {
      grm.rows <- c(row.indices[1]:(row.indices[1]+number.items-1))
      if ("PCM" %in% models.run) {
        pcm.rows <- c(row.indices[2]:(row.indices[2]+number.items-1))
        if ("GPCM" %in% models.run) {
          gpcm.rows <- c(row.indices[3]:nrow(output.temp))
        }
        # if ("RSM" %in% models.run) {
        #   rsm.rows <- c(row.indices[3]:nrow(output.temp))
        # }
      }
    } else {
      if ("PCM" %in% models.run) {
        pcm.rows <- c(row.indices[1]:(row.indices[1]+number.items-1))
        if ("GPCM" %in% models.run) {
          gpcm.rows <- c(row.indices[2]:nrow(output.temp))
        }
        # if ("RSM" %in% models.run) {
        #   rsm.rows <- c(row.indices[2]:nrow(output.temp))
        # }
      } else {
        if ("GPCM" %in% models.run) {
          gpcm.rows <- c(row.indices[1]:nrow(output.temp))
        }
        # if ("RSM" %in% models.run) {
        #   rsm.rows <- c(row.indices[1]:nrow(output.temp))
        # }
      }
    }
    #     Save the names of threshold columns:
    b.name.vec <- paste0("b",1:(num.resp.cat-1))
    #     For adjusted LD P-values: Put total number of item pairs in output:
    output.temp <- dplyr::mutate(output.temp, Num.Item.Pairs = 
                                   (((Num.Items*Num.Items)-Num.Items)/2) )
    output.temp <- dplyr::mutate(output.temp, Adj.Alpha.Level = alpha / Num.Item.Pairs )
    if (length(unique(output.temp$Num.Item.Pairs)) == 1) {
      adj.alpha.level <- alpha / unique(output.temp$Num.Item.Pairs)[1]
    }
    
    
    
    #   Run GRM:
    if ("GRM" %in% models.run) {
      if (is.null(num.cycles)) {
        grm.uni.mod <- mirt::mirt(data=dat.ord.items, model=1, itemtype="graded", verbose=F)
      } else {
        grm.uni.mod <- mirt::mirt(data=dat.ord.items, model=1, itemtype="graded", 
                                  verbose=F, technical=list(NCYCLES=num.cycles))
      }
      #     Need tryCatch statement for GRM Fit because it occasionally errors ("Error: could not invert orthogonal complement matrix")
      fit.grm.uni <- tryCatch(
        {mirt::M2(grm.uni.mod, type="C2")},
        error=function(e) { return(NA) }
      )
      if (length(fit.grm.uni) == 1 ) {
        if (is.na(fit.grm.uni)) {
          fit.grm.uni <- vector(mode="list", length=5)
          names(fit.grm.uni) <- c("RMSEA","RMSEA_5","RMSEA_95","TLI","CFI")
          fit.grm.uni[1:5] <- NA
        }
      }
      #     Get both G2 and X2 statistics (Chen & Thissen, 1997) and as approximately standardized:
      res.g2.grm.uni <- mirt::residuals(grm.uni.mod, type="LDG2", df.p=T, 
                                        verbose=F, approx.z=T)
      res.x2.grm.uni <- mirt::residuals(grm.uni.mod, type="LD", df.p=T, 
                                        verbose=F, approx.z=T)
      #         Save just the approximately standardized G2 and X2 stats separately:
      res.g2.grm.uni.ld <- res.g2.grm.uni$LD
      res.g2.grm.uni.ld[upper.tri(res.g2.grm.uni.ld, diag=T)] <- NA
      res.x2.grm.uni.ld <- res.x2.grm.uni$LD
      res.x2.grm.uni.ld[upper.tri(res.x2.grm.uni.ld, diag=T)] <- NA
      #         Save just info for item2-item5 pairing:
      #           X2 stat for item2-item5:
      res.x2.grm.uni.ld.df <- as.data.frame(res.x2.grm.uni.ld)
      res.x2.grm.uni.ld.df <- tibble::rownames_to_column(res.x2.grm.uni.ld.df,
                                                         var="Item")
      x2.ld.grm.items.2.5 <- res.x2.grm.uni.ld.df[
        which(res.x2.grm.uni.ld.df$Item=="Item_5"), "Item_2"]
      #           G2 stat for item2-item5:
      res.g2.grm.uni.ld.df <- as.data.frame(res.g2.grm.uni.ld)
      res.g2.grm.uni.ld.df <- tibble::rownames_to_column(res.g2.grm.uni.ld.df,
                                                         var="Item")
      g2.ld.grm.items.2.5 <- res.g2.grm.uni.ld.df[
        which(res.g2.grm.uni.ld.df$Item=="Item_5"), "Item_2"]
      #           P-value for X2 stat for item2-item5:
      res.x2.grm.uni.p <- res.x2.grm.uni$df.p
      res.x2.grm.uni.p[lower.tri(res.x2.grm.uni.p, diag=T)] <- NA
      res.x2.grm.uni.p.df <- as.data.frame(res.x2.grm.uni.p)
      res.x2.grm.uni.p.df <- tibble::rownames_to_column(res.x2.grm.uni.p.df,
                                                        var="Item")
      p.x2.ld.grm.items.2.5 <- res.x2.grm.uni.p.df[
        which(res.x2.grm.uni.p.df$Item=="Item_2"), "Item_5"]
      #           P-value for G2 stat for item2-item5:
      res.g2.grm.uni.p <- res.g2.grm.uni$df.p
      res.g2.grm.uni.p[lower.tri(res.g2.grm.uni.p, diag=T)] <- NA
      res.g2.grm.uni.p.df <- as.data.frame(res.g2.grm.uni.p)
      res.g2.grm.uni.p.df <- tibble::rownames_to_column(res.g2.grm.uni.p.df,
                                                        var="Item")
      p.g2.ld.grm.items.2.5 <- res.g2.grm.uni.p.df[
        which(res.g2.grm.uni.p.df$Item=="Item_2"), "Item_5"]
      #     Save coefs:
      coef.grm.uni <- coef(grm.uni.mod, IRTpars=T, simplify=T)
      coef.grm.uni.df <- as.data.frame(coef.grm.uni$items)
      coef.grm.uni.df <- tibble::rownames_to_column(coef.grm.uni.df, var="Item")
      #     Save scoring statistics and item fit info:
      omega.grm.uni <- omega.fntn(grm.uni.mod)
      ecv.grm.uni <- ecv.fntn(grm.uni.mod)
      ecv.grm.uni <- ecv.grm.uni[["Overall"]]
      itemfit.grm.uni <- mirt::itemfit(grm.uni.mod)
      
      #     GRM: Put model results into output:
      output.temp[grm.rows,"Model"] <- "GRM"
      output.temp[grm.rows,"Converged"] <- ifelse(grm.uni.mod@OptimInfo$converged==T,
                                                  yes=1, no=0)
      output.temp[grm.rows,"RMSEA"] <- fit.grm.uni$RMSEA
      output.temp[grm.rows,"RMSEA_5"] <- fit.grm.uni$RMSEA_5
      output.temp[grm.rows,"RMSEA_95"] <- fit.grm.uni$RMSEA_95
      output.temp[grm.rows,"TLI"] <- fit.grm.uni$TLI
      output.temp[grm.rows,"CFI"] <- fit.grm.uni$CFI
      output.temp[grm.rows,"Count.P.LD.G2"] <- sum(res.g2.grm.uni.p.df < alpha, na.rm=T)
      output.temp[grm.rows,"Count.P.LD.X2"] <- sum(res.x2.grm.uni.p.df < alpha, na.rm=T)
      output.temp[grm.rows,"Count.AdjP.LD.G2"] <- paste0(sum(res.g2.grm.uni.p.df < adj.alpha.level, na.rm=T),collapse=",")
      output.temp[grm.rows,"Count.AdjP.LD.X2"] <- paste0(sum(res.x2.grm.uni.p.df < adj.alpha.level, na.rm=T),collapse=",") 
      # output.temp[grm.rows,"Count.P.LD.G2"] <- sum(res.g2.grm.uni$df.p < alpha, na.rm=T)
      # output.temp[grm.rows,"Count.P.LD.X2"] <- sum(res.x2.grm.uni$df.p < alpha, na.rm=T)
      # output.temp[grm.rows,"Count.AdjP.LD.G2"] <- paste0(sum(res.g2.grm.uni$df.p < adj.alpha.level, na.rm=T),collapse=",")
      # output.temp[grm.rows,"Count.AdjP.LD.X2"] <- paste0(sum(res.x2.grm.uni$df.p < adj.alpha.level, na.rm=T),collapse=",")
      output.temp[grm.rows,"Count.Gt5LE10.LD.zG2"] <- sum((res.g2.grm.uni.ld > 5 & 
                                                             res.g2.grm.uni.ld <= 10), na.rm=T)
      output.temp[grm.rows,"Count.Gt5LE10.LD.zX2"] <- sum((res.x2.grm.uni.ld > 5 & 
                                                             res.x2.grm.uni.ld <= 10), na.rm=T)
      output.temp[grm.rows,"Count.Gt10.LD.zG2"] <- sum(res.g2.grm.uni.ld > 10, na.rm=T)
      output.temp[grm.rows,"Count.Gt10.LD.zX2"] <- sum(res.x2.grm.uni.ld > 10, na.rm=T)
      output.temp[grm.rows,"Items2.5.LD.zX2"] <- x2.ld.grm.items.2.5
      output.temp[grm.rows,"Items2.5.LD.zG2"] <- g2.ld.grm.items.2.5
      output.temp[grm.rows,"Items2.5.LD.P.X2"] <- p.x2.ld.grm.items.2.5
      output.temp[grm.rows,"Items2.5.LD.P.G2"] <- p.g2.ld.grm.items.2.5
      output.temp[grm.rows,"Items2.5.LD.PltAdjAlpha.X2"] <- ifelse(p.x2.ld.grm.items.2.5 < adj.alpha.level,
                                                                   yes=1, no=0)
      output.temp[grm.rows,"Items2.5.LD.PltAdjAlpha.G2"] <- ifelse(p.g2.ld.grm.items.2.5 < adj.alpha.level,
                                                                   yes=1, no=0)
      output.temp[grm.rows,"Omega"] <- omega.grm.uni
      output.temp[grm.rows,"ECV.Overall"] <- ecv.grm.uni
      
      #     GRM: Put item info results into output:
      output.temp[grm.rows,"Item"] <- itemfit.grm.uni$item
      output.temp[grm.rows,"Gen.Fac.Load.Val"] <- gen.fac.loads
      output.temp[grm.rows,"Loads.On.Factor"] <- label.fac.items$Loads.On.Factor
      output.temp[grm.rows,"S_X2"] <- itemfit.grm.uni$S_X2
      output.temp[grm.rows,"RMSEA.S_X2"] <- itemfit.grm.uni$RMSEA.S_X2
      output.temp[grm.rows,"p.S_X2"] <- itemfit.grm.uni$p.S_X2
      if (all(output.temp[grm.rows,"Item"] == coef.grm.uni.df$Item)) {
        output.temp[grm.rows,"a"] <- coef.grm.uni.df$a
        for (b.name in b.name.vec) {
          if (b.name %in% colnames(coef.grm.uni.df)) {
            output.temp[grm.rows,b.name] <- coef.grm.uni.df[[b.name]]
          }
        } #end of loop through b.name.vec
      } else {
        for (item.name in unique(coef.grm.uni.df$Item)) {
          out.item.row <- which(output.temp$Item==item.name & output.temp$Model=="GRM")
          coef.item.row <- which(coef.grm.uni.df$Item==item.name)
          output.temp[out.item.row,"a"] <- coef.grm.uni.df[coef.item.row,"a"] 
          for (b.name in b.name.vec) {
            if (b.name %in% colnames(coef.grm.uni.df)) {
              output.temp[out.item.row,b.name] <- coef.grm.uni.df[coef.item.row,b.name]
            }
          } #end of loop through b.name.vec
        }
      }
    } #end of if "GRM" is in models.run
    
    
    
    #   Run PCM (Rasch):
    if ("PCM" %in% models.run) {
      if (is.null(num.cycles)) {
        pcm.uni.mod <- mirt::mirt(data=dat.ord.items, model=1, itemtype="Rasch", verbose=F)
      } else {
        pcm.uni.mod <- mirt::mirt(data=dat.ord.items, model=1, itemtype="Rasch", 
                                  verbose=F, technical=list(NCYCLES=num.cycles))
      }
      fit.pcm.uni <- mirt::M2(pcm.uni.mod, type="C2")
      #     Get both G2 and X2 statistics (Chen & Thissen, 1997) and as approximately standardized:
      res.g2.pcm.uni <- mirt::residuals(pcm.uni.mod, type="LDG2", df.p=T, 
                                        verbose=F, approx.z=T)
      res.x2.pcm.uni <- mirt::residuals(pcm.uni.mod, type="LD", df.p=T, 
                                        verbose=F, approx.z=T)
      #         Save just the approximately standardized G2 and X2 stats separately:
      res.g2.pcm.uni.ld <- res.g2.pcm.uni$LD
      res.g2.pcm.uni.ld[upper.tri(res.g2.pcm.uni.ld, diag=T)] <- NA
      res.x2.pcm.uni.ld <- res.x2.pcm.uni$LD
      res.x2.pcm.uni.ld[upper.tri(res.x2.pcm.uni.ld, diag=T)] <- NA
      #         Save just info for item2-item5 pairing:
      #           X2 stat for item2-item5:
      res.x2.pcm.uni.ld.df <- as.data.frame(res.x2.pcm.uni.ld)
      res.x2.pcm.uni.ld.df <- tibble::rownames_to_column(res.x2.pcm.uni.ld.df,
                                                         var="Item")
      x2.ld.pcm.items.2.5 <- res.x2.pcm.uni.ld.df[
        which(res.x2.pcm.uni.ld.df$Item=="Item_5"), "Item_2"]
      #           G2 stat for item2-item5:
      res.g2.pcm.uni.ld.df <- as.data.frame(res.g2.pcm.uni.ld)
      res.g2.pcm.uni.ld.df <- tibble::rownames_to_column(res.g2.pcm.uni.ld.df,
                                                         var="Item")
      g2.ld.pcm.items.2.5 <- res.g2.pcm.uni.ld.df[
        which(res.g2.pcm.uni.ld.df$Item=="Item_5"), "Item_2"]
      #           P-value for X2 stat for item2-item5:
      res.x2.pcm.uni.p <- res.x2.pcm.uni$df.p
      res.x2.pcm.uni.p[lower.tri(res.x2.pcm.uni.p, diag=T)] <- NA
      res.x2.pcm.uni.p.df <- as.data.frame(res.x2.pcm.uni.p)
      res.x2.pcm.uni.p.df <- tibble::rownames_to_column(res.x2.pcm.uni.p.df,
                                                        var="Item")
      p.x2.ld.pcm.items.2.5 <- res.x2.pcm.uni.p.df[
        which(res.x2.pcm.uni.p.df$Item=="Item_2"), "Item_5"]
      #           P-value for G2 stat for item2-item5:
      res.g2.pcm.uni.p <- res.g2.pcm.uni$df.p
      res.g2.pcm.uni.p[lower.tri(res.g2.pcm.uni.p, diag=T)] <- NA
      res.g2.pcm.uni.p.df <- as.data.frame(res.g2.pcm.uni.p)
      res.g2.pcm.uni.p.df <- tibble::rownames_to_column(res.g2.pcm.uni.p.df,
                                                        var="Item")
      p.g2.ld.pcm.items.2.5 <- res.g2.pcm.uni.p.df[
        which(res.g2.pcm.uni.p.df$Item=="Item_2"), "Item_5"]
      #     Save coefs:
      coef.pcm.uni <- coef(pcm.uni.mod, IRTpars=T, simplify=T)
      coef.pcm.uni.df <- as.data.frame(coef.pcm.uni$items)
      coef.pcm.uni.df <- tibble::rownames_to_column(coef.pcm.uni.df, var="Item")
      #     Save scoring statistics and item fit info:
      omega.pcm.uni <- omega.fntn(pcm.uni.mod)
      ecv.pcm.uni <- ecv.fntn(pcm.uni.mod)
      ecv.pcm.uni <- ecv.pcm.uni[["Overall"]]
      itemfit.pcm.uni <- mirt::itemfit(pcm.uni.mod)
      
      #     Rasch.PCM: Put model results into output:
      output.temp[pcm.rows,"Model"] <- "Rasch.PCM"
      output.temp[pcm.rows,"Converged"] <- ifelse(pcm.uni.mod@OptimInfo$converged==T,
                                                  yes=1, no=0)
      output.temp[pcm.rows,"RMSEA"] <- fit.pcm.uni$RMSEA
      output.temp[pcm.rows,"RMSEA_5"] <- fit.pcm.uni$RMSEA_5
      output.temp[pcm.rows,"RMSEA_95"] <- fit.pcm.uni$RMSEA_95
      output.temp[pcm.rows,"TLI"] <- fit.pcm.uni$TLI
      output.temp[pcm.rows,"CFI"] <- fit.pcm.uni$CFI
      output.temp[pcm.rows,"Count.P.LD.G2"] <- sum(res.g2.pcm.uni.p.df < alpha, na.rm=T)
      output.temp[pcm.rows,"Count.P.LD.X2"] <- sum(res.x2.pcm.uni.p.df < alpha, na.rm=T)
      output.temp[pcm.rows,"Count.AdjP.LD.G2"] <- paste0(sum(res.g2.pcm.uni.p.df < adj.alpha.level, na.rm=T),collapse=",")
      output.temp[pcm.rows,"Count.AdjP.LD.X2"] <- paste0(sum(res.x2.pcm.uni.p.df < adj.alpha.level, na.rm=T),collapse=",")
      # output.temp[pcm.rows,"Count.P.LD.G2"] <- sum(res.g2.pcm.uni$df.p < alpha, na.rm=T)
      # output.temp[pcm.rows,"Count.P.LD.X2"] <- sum(res.x2.pcm.uni$df.p < alpha, na.rm=T)
      # output.temp[pcm.rows,"Count.AdjP.LD.G2"] <- paste0(sum(res.g2.pcm.uni$df.p < adj.alpha.level, na.rm=T),collapse=",")
      # output.temp[pcm.rows,"Count.AdjP.LD.X2"] <- paste0(sum(res.x2.pcm.uni$df.p < adj.alpha.level, na.rm=T),collapse=",")
      output.temp[pcm.rows,"Count.Gt5LE10.LD.zG2"] <- sum((res.g2.pcm.uni.ld > 5 & 
                                                             res.g2.pcm.uni.ld <= 10), na.rm=T)
      output.temp[pcm.rows,"Count.Gt5LE10.LD.zX2"] <- sum((res.x2.pcm.uni.ld > 5 & 
                                                             res.x2.pcm.uni.ld <= 10), na.rm=T)
      output.temp[pcm.rows,"Count.Gt10.LD.zG2"] <- sum(res.g2.pcm.uni.ld > 10, na.rm=T)
      output.temp[pcm.rows,"Count.Gt10.LD.zX2"] <- sum(res.x2.pcm.uni.ld > 10, na.rm=T)
      output.temp[pcm.rows,"Items2.5.LD.zX2"] <- x2.ld.pcm.items.2.5
      output.temp[pcm.rows,"Items2.5.LD.zG2"] <- g2.ld.pcm.items.2.5
      output.temp[pcm.rows,"Items2.5.LD.P.X2"] <- p.x2.ld.pcm.items.2.5
      output.temp[pcm.rows,"Items2.5.LD.P.G2"] <- p.g2.ld.pcm.items.2.5
      output.temp[pcm.rows,"Items2.5.LD.PltAdjAlpha.X2"] <- ifelse(p.x2.ld.pcm.items.2.5 < adj.alpha.level,
                                                                   yes=1, no=0)
      output.temp[pcm.rows,"Items2.5.LD.PltAdjAlpha.G2"] <- ifelse(p.g2.ld.pcm.items.2.5 < adj.alpha.level,
                                                                   yes=1, no=0)
      output.temp[pcm.rows,"Omega"] <- omega.pcm.uni
      output.temp[pcm.rows,"ECV.Overall"] <- ecv.pcm.uni
      
      #     Rasch.PCM: Put item info results into output:
      output.temp[pcm.rows,"Item"] <- itemfit.pcm.uni$item
      output.temp[pcm.rows,"Gen.Fac.Load.Val"] <- gen.fac.loads
      output.temp[pcm.rows,"Loads.On.Factor"] <- label.fac.items$Loads.On.Factor
      output.temp[pcm.rows,"S_X2"] <- itemfit.pcm.uni$S_X2
      output.temp[pcm.rows,"RMSEA.S_X2"] <- itemfit.pcm.uni$RMSEA.S_X2
      output.temp[pcm.rows,"p.S_X2"] <- itemfit.pcm.uni$p.S_X2
      if (all(output.temp[pcm.rows,"Item"] == coef.pcm.uni.df$Item)) {
        output.temp[pcm.rows,"a"] <- coef.pcm.uni.df$a
        for (b.name in b.name.vec) {
          if (b.name %in% colnames(coef.pcm.uni.df)) {
            output.temp[pcm.rows,b.name] <- coef.pcm.uni.df[[b.name]]
          }
        } #end of loop through b.name.vec
        if ("b" %in% colnames(coef.pcm.uni.df)) { #IF not all categories endorsed, 
          #   when num.resp.cat==3, PCM will move the single threshold to "b" col, 
          #   with no num after (whereas GRM just moves the thresholds down one) -- fix it same as GRM!
          output.temp[pcm.rows,"b1"] <- coef.pcm.uni.df[["b"]]
        }
        if ("g" %in% colnames(coef.pcm.uni.df)) { #IF not all categories endorsed, 
          #   PCM MAY have a "g" and a "u" column (when just one threshold, I think)
          output.temp[pcm.rows,"g"] <- coef.pcm.uni.df[["g"]]
        }
        if ("u" %in% colnames(coef.pcm.uni.df)) { #IF not all categories endorsed, 
          #   PCM MAY have a "g" and a "u" column (when just one threshold, I think)
          output.temp[pcm.rows,"u"] <- coef.pcm.uni.df[["u"]]
        }
      } else {
        for (item.name in unique(coef.pcm.uni.df$Item)) {
          out.item.row <- which(output.temp$Item==item.name & output.temp$Model=="Rasch.PCM")
          coef.item.row <- which(coef.pcm.uni.df$Item==item.name)
          output.temp[out.item.row,"a"] <- coef.pcm.uni.df[coef.item.row,"a"] 
          for (b.name in b.name.vec) {
            if (b.name %in% colnames(coef.pcm.uni.df)) {
              output.temp[out.item.row,b.name] <- coef.pcm.uni.df[coef.item.row,b.name]
            }
          } #end of loop through b.name.vec
          if ("b" %in% colnames(coef.pcm.uni.df)) { #IF not all categories endorsed, 
            #   when num.resp.cat==3, PCM will move the single threshold to "b" col, 
            #   with no num after (whereas GRM just moves the thresholds down one) -- fix it same as GRM!
            output.temp[out.item.row,"b1"] <- coef.pcm.uni.df[coef.item.row,"b"]
          }
          if ("g" %in% colnames(coef.pcm.uni.df)) { #IF not all categories endorsed, 
            #   PCM MAY have a "g" and a "u" column (when just one threshold)
            output.temp[out.item.row,"g"] <- coef.pcm.uni.df[coef.item.row,"g"]
          }
          if ("u" %in% colnames(coef.pcm.uni.df)) { #IF not all categories endorsed, 
            #   PCM MAY have a "g" and a "u" column (when just one threshold)
            output.temp[out.item.row,"u"] <- coef.pcm.uni.df[coef.item.row,"u"]
          }
        }
      }
    } #end of if "PCM" is in models.run
    
    
    
    #   Run GPCM:
    if ("GPCM" %in% models.run) {
      if (is.null(num.cycles)) {
        gpcm.uni.mod <- mirt::mirt(data=dat.ord.items, model=1, itemtype="graded", verbose=F)
      } else {
        gpcm.uni.mod <- mirt::mirt(data=dat.ord.items, model=1, itemtype="graded", 
                                   verbose=F, technical=list(NCYCLES=num.cycles))
      }
      #     tryCatch statement in case occasionally errors:
      fit.gpcm.uni <- tryCatch(
        {mirt::M2(gpcm.uni.mod, type="C2")},
        error=function(e) { return(NA) }
      )
      if (length(fit.gpcm.uni) == 1 ) {
        if (is.na(fit.gpcm.uni)) {
          fit.gpcm.uni <- vector(mode="list", length=5)
          names(fit.gpcm.uni) <- c("RMSEA","RMSEA_5","RMSEA_95","TLI","CFI")
          fit.gpcm.uni[1:5] <- NA
        }
      }
      #     Get both G2 and X2 statistics (Chen & Thissen, 1997) and as approximately standardized:
      res.g2.gpcm.uni <- mirt::residuals(gpcm.uni.mod, type="LDG2", df.p=T, 
                                         verbose=F, approx.z=T)
      res.x2.gpcm.uni <- mirt::residuals(gpcm.uni.mod, type="LD", df.p=T, 
                                         verbose=F, approx.z=T)
      #         Save just the approximately standardized G2 and X2 stats separately:
      res.g2.gpcm.uni.ld <- res.g2.gpcm.uni$LD
      res.g2.gpcm.uni.ld[upper.tri(res.g2.gpcm.uni.ld, diag=T)] <- NA
      res.x2.gpcm.uni.ld <- res.x2.gpcm.uni$LD
      res.x2.gpcm.uni.ld[upper.tri(res.x2.gpcm.uni.ld, diag=T)] <- NA
      #         Save just info for item2-item5 pairing:
      #           X2 stat for item2-item5:
      res.x2.gpcm.uni.ld.df <- as.data.frame(res.x2.gpcm.uni.ld)
      res.x2.gpcm.uni.ld.df <- tibble::rownames_to_column(res.x2.gpcm.uni.ld.df,
                                                          var="Item")
      x2.ld.gpcm.items.2.5 <- res.x2.gpcm.uni.ld.df[
        which(res.x2.gpcm.uni.ld.df$Item=="Item_5"), "Item_2"]
      #           G2 stat for item2-item5:
      res.g2.gpcm.uni.ld.df <- as.data.frame(res.g2.gpcm.uni.ld)
      res.g2.gpcm.uni.ld.df <- tibble::rownames_to_column(res.g2.gpcm.uni.ld.df,
                                                          var="Item")
      g2.ld.gpcm.items.2.5 <- res.g2.gpcm.uni.ld.df[
        which(res.g2.gpcm.uni.ld.df$Item=="Item_5"), "Item_2"]
      #           P-value for X2 stat for item2-item5:
      res.x2.gpcm.uni.p <- res.x2.gpcm.uni$df.p
      res.x2.gpcm.uni.p[lower.tri(res.x2.gpcm.uni.p, diag=T)] <- NA
      res.x2.gpcm.uni.p.df <- as.data.frame(res.x2.gpcm.uni.p)
      res.x2.gpcm.uni.p.df <- tibble::rownames_to_column(res.x2.gpcm.uni.p.df,
                                                         var="Item")
      p.x2.ld.gpcm.items.2.5 <- res.x2.gpcm.uni.p.df[
        which(res.x2.gpcm.uni.p.df$Item=="Item_2"), "Item_5"]
      #           P-value for G2 stat for item2-item5:
      res.g2.gpcm.uni.p <- res.g2.gpcm.uni$df.p
      res.g2.gpcm.uni.p[lower.tri(res.g2.gpcm.uni.p, diag=T)] <- NA
      res.g2.gpcm.uni.p.df <- as.data.frame(res.g2.gpcm.uni.p)
      res.g2.gpcm.uni.p.df <- tibble::rownames_to_column(res.g2.gpcm.uni.p.df,
                                                         var="Item")
      p.g2.ld.gpcm.items.2.5 <- res.g2.gpcm.uni.p.df[
        which(res.g2.gpcm.uni.p.df$Item=="Item_2"), "Item_5"]
      #     Save coefs:
      coef.gpcm.uni <- coef(gpcm.uni.mod, IRTpars=T, simplify=T)
      coef.gpcm.uni.df <- as.data.frame(coef.gpcm.uni$items)
      coef.gpcm.uni.df <- tibble::rownames_to_column(coef.gpcm.uni.df, var="Item")
      #     Save scoring statistics and item fit info:
      omega.gpcm.uni <- omega.fntn(gpcm.uni.mod)
      ecv.gpcm.uni <- ecv.fntn(gpcm.uni.mod)
      ecv.gpcm.uni <- ecv.gpcm.uni[["Overall"]]
      itemfit.gpcm.uni <- mirt::itemfit(gpcm.uni.mod)
      
      #     GPCM: Put model results into output:
      output.temp[gpcm.rows,"Model"] <- "GPCM"
      output.temp[gpcm.rows,"Converged"] <- ifelse(gpcm.uni.mod@OptimInfo$converged==T,
                                                   yes=1, no=0)
      output.temp[gpcm.rows,"RMSEA"] <- fit.gpcm.uni$RMSEA
      output.temp[gpcm.rows,"RMSEA_5"] <- fit.gpcm.uni$RMSEA_5
      output.temp[gpcm.rows,"RMSEA_95"] <- fit.gpcm.uni$RMSEA_95
      output.temp[gpcm.rows,"TLI"] <- fit.gpcm.uni$TLI
      output.temp[gpcm.rows,"CFI"] <- fit.gpcm.uni$CFI
      output.temp[gpcm.rows,"Count.P.LD.G2"] <- sum(res.g2.gpcm.uni.p.df < alpha, na.rm=T)
      output.temp[gpcm.rows,"Count.P.LD.X2"] <- sum(res.x2.gpcm.uni.p.df < alpha, na.rm=T)
      output.temp[gpcm.rows,"Count.AdjP.LD.G2"] <- paste0(sum(res.g2.gpcm.uni.p.df < adj.alpha.level, na.rm=T),collapse=",")
      output.temp[gpcm.rows,"Count.AdjP.LD.X2"] <- paste0(sum(res.x2.gpcm.uni.p.df < adj.alpha.level, na.rm=T),collapse=",")
      # output.temp[gpcm.rows,"Count.P.LD.G2"] <- sum(res.g2.gpcm.uni$df.p < alpha, na.rm=T)
      # output.temp[gpcm.rows,"Count.P.LD.X2"] <- sum(res.x2.gpcm.uni$df.p < alpha, na.rm=T)
      # output.temp[gpcm.rows,"Count.AdjP.LD.G2"] <- paste0(sum(res.g2.gpcm.uni$df.p < adj.alpha.level, na.rm=T),collapse=",")
      # output.temp[gpcm.rows,"Count.AdjP.LD.X2"] <- paste0(sum(res.x2.gpcm.uni$df.p < adj.alpha.level, na.rm=T),collapse=",")
      output.temp[gpcm.rows,"Count.Gt5LE10.LD.zG2"] <- sum((res.g2.gpcm.uni.ld > 5 & 
                                                              res.g2.gpcm.uni.ld <= 10), na.rm=T)
      output.temp[gpcm.rows,"Count.Gt5LE10.LD.zX2"] <- sum((res.x2.gpcm.uni.ld > 5 & 
                                                              res.x2.gpcm.uni.ld <= 10), na.rm=T)
      output.temp[gpcm.rows,"Count.Gt10.LD.zG2"] <- sum(res.g2.gpcm.uni.ld > 10, na.rm=T)
      output.temp[gpcm.rows,"Count.Gt10.LD.zX2"] <- sum(res.x2.gpcm.uni.ld > 10, na.rm=T)
      output.temp[gpcm.rows,"Items2.5.LD.zX2"] <- x2.ld.gpcm.items.2.5
      output.temp[gpcm.rows,"Items2.5.LD.zG2"] <- g2.ld.gpcm.items.2.5
      output.temp[gpcm.rows,"Items2.5.LD.P.X2"] <- p.x2.ld.gpcm.items.2.5
      output.temp[gpcm.rows,"Items2.5.LD.P.G2"] <- p.g2.ld.gpcm.items.2.5
      output.temp[gpcm.rows,"Items2.5.LD.PltAdjAlpha.X2"] <- ifelse(p.x2.ld.gpcm.items.2.5 < adj.alpha.level,
                                                                    yes=1, no=0)
      output.temp[gpcm.rows,"Items2.5.LD.PltAdjAlpha.G2"] <- ifelse(p.g2.ld.gpcm.items.2.5 < adj.alpha.level,
                                                                    yes=1, no=0)
      output.temp[gpcm.rows,"Omega"] <- omega.gpcm.uni
      output.temp[gpcm.rows,"ECV.Overall"] <- ecv.gpcm.uni
      
      #     GPCM: Put item info results into output:
      output.temp[gpcm.rows,"Item"] <- itemfit.gpcm.uni$item
      output.temp[gpcm.rows,"Gen.Fac.Load.Val"] <- gen.fac.loads
      output.temp[gpcm.rows,"Loads.On.Factor"] <- label.fac.items$Loads.On.Factor
      output.temp[gpcm.rows,"S_X2"] <- itemfit.gpcm.uni$S_X2
      output.temp[gpcm.rows,"RMSEA.S_X2"] <- itemfit.gpcm.uni$RMSEA.S_X2
      output.temp[gpcm.rows,"p.S_X2"] <- itemfit.gpcm.uni$p.S_X2
      if (all(output.temp[gpcm.rows,"Item"] == coef.gpcm.uni.df$Item)) {
        output.temp[gpcm.rows,"a"] <- coef.gpcm.uni.df$a
        for (b.name in b.name.vec) {
          if (b.name %in% colnames(coef.gpcm.uni.df)) {
            output.temp[gpcm.rows,b.name] <- coef.gpcm.uni.df[[b.name]]
          }
        } #end of loop through b.name.vec
      } else {
        for (item.name in unique(coef.gpcm.uni.df$Item)) {
          out.item.row <- which(output.temp$Item==item.name & output.temp$Model=="GPCM")
          coef.item.row <- which(coef.gpcm.uni.df$Item==item.name)
          output.temp[out.item.row,"a"] <- coef.gpcm.uni.df[coef.item.row,"a"] 
          for (b.name in b.name.vec) {
            if (b.name %in% colnames(coef.gpcm.uni.df)) {
              output.temp[out.item.row,b.name] <- coef.gpcm.uni.df[coef.item.row,b.name]
            }
          } #end of loop through b.name.vec
        }
      }
    } #end of if "GPCM" is in models.run
    
    
    
    #   Check whether any items are missing response categories:
    #       (Doing this after input model & item info because item names now in output)
    
    #     Does every item have the same number of response categories?: 
    #       (Could still have some missing, though unlikely)
    lengths.vec <- vector(mode="numeric", length=number.items)
    for (item.num in 1:number.items) {
      item <- paste0("Item_",item.num)
      lengths.vec[item.num] <- length(xtabs(~dat.ord.items[[item]]))
      output.temp[which(output.temp$Item==item),"Obs.Num.Resp.Cat"] <- 
        length(xtabs(~dat.ord.items[[item]]))
    }
    all.items.same.cat <- ifelse(length(unique(lengths.vec))==1,
                                 yes=1, no=0)
    
    #     Does every item have the *correct* number of response categories?:
    #       (Can only be true if they all have the same number of categories)
    if (all.items.same.cat==1) {
      all.items.correct.cat <- ifelse(unique(lengths.vec)==num.resp.cat,
                                      yes=1, no=0)
    } else {
      all.items.correct.cat <- 0
    }
    
    #     If differing numbers of response categories, how many different ones?:
    num.diff.resp.cat.items <- length(unique(lengths.vec))
    
    #     Min and Max number of response categories:
    min.resp.cat.items <- min(lengths.vec, na.rm=T)
    max.resp.cat.items <- max(lengths.vec, na.rm=T)
    
    #     Put that info into the output:
    output.temp[,"All.Items.Corr.Cat"] <- all.items.correct.cat
    output.temp[,"All.Items.Same.Cat"] <- all.items.same.cat
    output.temp[,"Num.Diff.Obs.Cat"] <- num.diff.resp.cat.items
    output.temp[,"Min.Obs.Cat"] <- min.resp.cat.items
    output.temp[,"Max.Obs.Cat"] <- max.resp.cat.items
    
    
    
    
    # Dataset - Add condition info to dataset:
    dat.ord.all$Repl <- repl
    dat.ord.all$Cur.Seed <- paste0(cur.seed, collapse=",")
    if (!is.null(uni.cond)) {
      dat.ord.all$Uni.Condition <- uni.cond
    }
    if (!is.null(ld.load)) {
      dat.ord.all$LD.Load <- ld.load
    }
    dat.ord.all$Sample.Size <- sample.size
    dat.ord.all$Num.Items <- number.items
    dat.ord.all$Num.Resp.Cat <- num.resp.cat
    
    
    
    
    # Return the results using the custom class:
    result <- resultsAndData()
    result$resultOut <- output.temp
    result$resultData <- dat.ord.all
    return(result)
    
    
    
    
  }
