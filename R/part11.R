# Script contains helper functions for drug response prediction

# label the drugs with proper u
giveDrugLabel3 = function(drid, dtab=drugs, ctab=conctab) {
  vapply(strsplit(drid, "_"), function(x) {
    k <- paste(x[1:2], collapse="_")
    paste0(dtab[k, "name"], " ",
           switch(x[3], "1:5"="c1:5", "4:5"="c4:5",
                  paste0(ctab[k, as.integer(x[3])], " \u00B5","M")))
  }, character(1))
}

# select features from the lpd object to use as response and predictors in the model
featureSelectionForLasso = function(objective, predictors, lpd) {
  
  dimobj = dimnames(objective)
  objectiveName = dimobj[[1]][1]
  
  #design matix
  fx = t(exprs(lpd)[predictors, ])

    # check the objective
  stopifnot(identical(rownames(fx), dimobj[[2]]))
  fy = objective[,,drop=TRUE]
  
  # select predictors
  fx = fx[, (fData(lpd)[predictors, "type"] %in% c("Methylation_Cluster",
                                                   "IGHV","gen", "pretreat")),
          drop=FALSE]
  
  #  make sure that there are no NAs in the table
  stopifnot(all(!is.na(fx)))
  
  # print table with numer of predictors from each group
  # message("Objective: ", objectiveName, "\n")
  prefix = paste(ifelse(dataType %in% unique(fData(lpd)[colnames(fx), "type"]),
                        names(dataType), "_"), collapse="")
  
  n = length(unique(fy))
  family = "gaussian"
  return(list(fx=fx, fy=fy, objectiveName=objectiveName, prefix=prefix,
              family=family))
}

# plotting the predictions
plotPredictions = function(fx, fy, objective, pred, coeffs, lpd, nm, lim,
                           objectiveName) {

  # PREPARE DATA FOR PLOTTING
  
  stopifnot( identical(dim(pred), c(length(fy), 1L)), identical(rownames(pred),
                                                                names(fy)) )
  ordy <- order(fy, pred[, 1])
  # design matrix of selected predictors (unscaled)
  mat  <- t(fx[ordy, names(coeffs), drop=FALSE])
  # human-readable names where available & drugs with concentrations
  nicename = names(coeffs) %>% `names<-`(names(coeffs))
  idx = grepl("D_", nicename)
  nicename[idx] = giveDrugLabel3(nicename[idx])
  nicename[!idx] = fData(lpd)[names(nicename)[!idx], "name"] %>%
    `names<-`(names(nicename)[!idx])
  nicename = ifelse(!is.na(nicename) & (nicename!=""), nicename, names(nicename))
  nicename = gsub("Methylation_Cluster", "Methylation cluster", nicename)
  rownames(mat) = nicename[rownames(mat)]
  # prepare plot title
  title=nm
  
  # CREATE INDIVIDUAL PARTS OF THE FIGURE  
  
  # bar plot
  stopifnot(all(coeffs<lim & coeffs>-lim))
  part1df = data.frame(coeffs,
                       nm=factor(names(coeffs), levels=names(rev(coeffs))))
  part1df$col= ifelse(rownames(part1df)=="IGHV", "I",
                      ifelse(rownames(part1df)=="Methylation_Cluster", "M",
                             ifelse(rownames(part1df)=="Pretreatment", "P","G")))
  
  part1 = ggplot(data=part1df, aes(x=nm, y=coeffs, fill=col)) +
    geom_bar(stat="identity", colour="black", position = "identity",
             width=.66, size=0.2) +theme_bw() +
    geom_hline(yintercept=0, size=0.3) + scale_x_discrete(expand=c(0,0.5)) +
    scale_y_continuous(expand=c(0,0)) + coord_flip(ylim=c(-lim,lim)) +
    theme(panel.grid.major=element_blank(),
          panel.background=element_blank(),
          panel.grid.minor = element_blank(),
          axis.text=element_text(size=8),
          panel.border=element_blank()) +
    xlab("") + ylab("Model Coefficients") +
    geom_vline(xintercept=c(0.5), color="black", size=0.6)+
    scale_fill_manual(c("M", "I", "G", "P"),
                      values=c(M=coldef[["M"]][2],
                               I=coldef[["I"]],
                               G=coldef[["G"]],
                               P=coldef[["P"]]))
  
  
  # heat map
  # mat contains selected predictors with status for each patient
  # (e.g. 0/1 for mutations and IGHV, 0/0.5/1 for M)
  # to assign colors Gosia added 5 to meth and 2 to IGHV values,resuting in
  # 5 LP, 5.5 IP, 6 HP, 2 unmut IGHV or mut, 4 IGHV mut, 3 mut, 7 pre-treatment
  idx = grep("Methylation cluster", rownames(mat))
  mat[idx,] = mat[idx,]+5
  ## ighv
  idx = grep("IGHV", rownames(mat))
  mat[idx,] = (mat[idx,]*2)+2
  ## gene
  rnm = sapply(rownames(mat), function(nm) strsplit(nm," ")[[1]][1])
  idx = rownames(fData(lpd))[fData(lpd)$type=="gen" & rownames(lpd) %in% rnm]
  idx = match(idx, rnm)
  mat[idx,] = mat[idx,]+2
  ## pretreat
  idx = grep("Pretreatment", rownames(mat))
  mat[idx,] = ifelse(mat[idx,]==0,2,7)
  
  mat2 = meltWholeDF(mat)
  mat2$Measure = factor(mat2$Measure, levels=sort(unique(mat2$Measure)))
  mat2$X = factor(mat2$X, levels=colnames(mat))
  mat2$Y = factor(mat2$Y, levels=rev(rownames(mat)))
  
  part2 = ggplot(mat2, aes(x=X, y=Y, fill=Measure)) +
    geom_tile() + theme_bw() +
    scale_fill_manual(name="Mutated",
                      values=c(`2`="gray96", `3`=paste0(coldef["G"], "E5"),
                               `5`=coldef[["M"]][1], `5.5`=coldef[["M"]][2],
                               `6`=coldef[["M"]][3], `7`=coldef[["P"]],
                               `4`=paste0(coldef["I"],"E5")), guide=FALSE) +
    scale_y_discrete(expand=c(0,0)) +
    theme(axis.text.y=element_text(hjust=0, size=14),
          axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          panel.border=element_rect(colour="gainsboro"),
          plot.title=element_text(size=12),
          legend.title=element_text(size=12),
          legend.text=element_text(size=12),
          panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    xlab("patients") + ylab("") + ggtitle(title)
  if(length(levels(mat2$Y)) > 1) {
    part2 = part2 + geom_hline(yintercept=seq(1.5, length(levels(mat2$Y)), 1),
                               colour="gainsboro", size=0.2)
  }
  
  # scatter plot
  mat3 = fy[colnames(mat)]
  mat3 = data.frame(X=factor(names(mat3), levels=names(mat3)), Y=mat3*100)
  Yrange = range(mat3$Y)
  Yhangs = diff(Yrange)*0.05
  Ylims = c(Yrange[1]-Yhangs, Yrange[2]+Yhangs)
  Yran = diff(Yrange)
  Ybreaks = if(Yran<=15) 5 else if(Yran>15 & Yran<=30) 10 else if(Yran>30 & Yran<=40) 15 else if(Yran>40 & Yran<=60) 20 else 40
  
  part4 = ggplot(mat3, aes(x=X, y=Y)) +
    geom_point(shape=21, fill="dimgrey", colour="black", size=1.2) +
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_text(size=8),
          panel.border=element_rect(colour="dimgrey", size=0.1),
          panel.background=element_rect(fill="gray96")) +
    ylab("Viability [%]") +
    scale_y_continuous(limits=c(Yrange[1]-Yhangs, Yrange[2]+Yhangs),
                       breaks=seq(-200,200,Ybreaks))
  
  
  # MERGE PARTS OF THE FIGURE
  # construct the gtable
  wdths = c(4, 0.5, 0.05*ncol(mat), 0.05)
  hghts = c(0.3, 0.25*nrow(mat), 0.1, 0.4, 0.35)
  gt = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))
  ## make grobs
  gg1 = ggplotGrob(part1)
  gg2 = ggplotGrob(part2)
  gg4 = ggplotGrob(part4)
  ## fill in the gtable
  gt = gtable_add_grob(gt, gtable_filter(gg1, "panel"), 2, 1) # add boxplot
  gt = gtable_add_grob(gt, gtable_filter(gg2, "panel"), 2, 3) # add heatmap
  gt = gtable_add_grob(gt, gtable_filter(gg4, "panel"), 4, 3) # add scatterplot
  gt = gtable_add_grob(gt, gtable_filter(gg2, "xlab"), 5, 3)
  gt = gtable_add_grob(gt, gg2$grobs[[whichInGrob(gg2, "title")]], 1, 3)
  gt = gtable_add_grob(gt, gg4$grobs[[whichInGrob(gg4, "axis-l")]], 4, 2)
  gt = gtable_add_grob(gt, gg1$grobs[[whichInGrob(gg1, "axis-b")]], 3, 1)
  gt = gtable_add_grob(gt, gtable_filter(gg1, "xlab-b"), 4, 1)
  
  # now complicated: add the axis-l labels from gg2
  gia = which(gg2$layout$name == "axis-l")
  gga = gg2$grobs[[gia]]
  gax = gga$children[[2]]
  gax$widths = rev(gax$widths)
  gax$grobs = rev(gax$grobs)
  gt = gtable_add_cols(gt, gg2$widths[gg2$layout[gia, ]$l])
  gt = gtable_add_grob(gt, gax, 2, 5)
  
  wdth = convertUnit(gt$widths, "in", val=TRUE)[5]
  # add column on the right of appropriate width
  maxwdth = 2
  gt = gtable_add_cols(gt, unit(maxwdth-wdth, "in"))
  
  return(list(plot=gt, width=sum(wdths)+maxwdth, height=sum(hghts),
              name=make.names(objectiveName)))
}


# main function to fit Lasso and produce barplots to find genetic
# determinants of drug response
doLasso = function(objective, predictors, lpd,suffix="",
                  nm=NA, lim=0.21, ncv=10, nfolds=10, std=FALSE,
                  adaLasso = TRUE) {
  
  #construct design and response matrix
  out = featureSelectionForLasso(objective, predictors, lpd)
  
  # for simplification
  fy = out[["fy"]]
  fx = out[["fx"]]
  family = out[["family"]]
  objectiveName = out[["objectiveName"]]
  prefix = out[["prefix"]]
  fxdim = dim(fx)
  print(sprintf("Prediction for: %s; #samples: %d; #features: %d",
                objectiveName, fxdim[1], fxdim[2]))
  
  # adaptive lasso for a more stable feature selection
  set.seed(19087)
  if(adaLasso){
    if(ncol(fx)>= nrow(fx)) {
      RidgeFit <- cv.glmnet(fx, fy, alpha = 0, standardize = std,
                            family = family, nfolds=10)
      # wRidge <- pmin(1/abs((coef(RidgeFit, s = RidgeFit$lambda.min))), 1e+300)
      wRidge <- 1/abs(coef(RidgeFit, s = RidgeFit$lambda.min))
      wRidge <- wRidge[-1]
      weights <- wRidge
    } else {
      lmFit <- lm(fy ~ fx)
      # wLM <- pmin(1/abs(coefficients(lmFit)[-1]), 1e+300)
      wLM <- 1/abs(coefficients(lmFit)[-1])
      weights <- wLM
    }
    excludedFeatures <- which(weights==Inf)
  } else {
    weights <- rep(1, ncol(fx))
    excludedFeatures <- NULL
  }
  
  #perform repeated cross-validation to find an optimal penalisatio
  #parameter minimizing the cross-validated MSE
  cv.out <- cvr.glmnet(Y=fy, X=fx, family=family,
                       alpha=1, nfolds=nfolds,
                       ncv=ncv, standardize=std,
                       exclude = excludedFeatures,
                       type.measure = "mse", penalty.factor = weights)
  
  # #fit Lasso model for optimal lambda                   
  fit = glmnet(y=fy, x=fx, family=family,
               alpha=1,
               standardize=std,
               exclude = excludedFeatures,
               lambda=cv.out$lambda, penalty.factor = weights)
  
  #get optimal lambda and corresponding predictors with coefficients
  lambda <- cv.out$lambda[which.min(cv.out$cvm)]
  coeffs <- coef(fit, lambda)
  coeffs_all <- coeffs
  coeffs <- as.vector(coeffs) %>%
    `names<-`(rownames(coeffs))   # cast from sparse matrix to ordinary vector
  coeffs <- coeffs[ coeffs!=0 ]
  
  # remove intercept term
  stopifnot(names(coeffs)[1]=="(Intercept)")
  if (length(coeffs) > 1) {
    coeffs <- coeffs[-1]
  } else {
    print("No (0) predictors for given parameters!")
    return(0)
  }
  coeffs <- sort(coeffs)
  
  # Residuals in the model
  pred <- predict(fit, newx = fx, s = lambda, type = "response")
  residuals <- pred[,1]-fy
  
  plot = plotPredictions(fx, fy, objective, pred, coeffs, lpd, nm, lim,
                         objectiveName)
  
  return(list(residuals=residuals, coeffs=coeffs_all, plot=plot))
}

# Make list of predictors for the given lpd
makeListOfPredictors = function(lpd) {
  return(list(
    predictorsM = rownames(fData(lpd))[fData(lpd)$type=="Methylation_Cluster"],
    predictorsG = rownames(fData(lpd))[fData(lpd)$type=="gen"],
    predictorsI = rownames(fData(lpd))[fData(lpd)$type=="IGHV"],
    predictorsP = rownames(fData(lpd))[fData(lpd)$type=="pretreat"]
  ))
}


# Pre-process data: explore what is available & feature selection
prepareLPD = function(lpd, minNumSamplesPerGroup, withMC=TRUE) {
  
  # PRETREATMENT
  # update the expression set by adding row about pretreatment
  pretreated <- t(matrix(ifelse(patmeta[colnames(lpd),
                                        "IC50beforeTreatment"], 0, 1),
                         dimnames=list(colnames(lpd), "Pretreatment"))) 
  fdata_pretreat <- data.frame(name=NA, type="pretreat", id=NA, subtype=NA,
                               row.names="Pretreatment")
  lpd <- ExpressionSet(assayData=rbind(exprs(lpd), pretreated),
                          phenoData=new("AnnotatedDataFrame", data=pData(lpd)),
                          featureData=new("AnnotatedDataFrame",
                                          rbind(fData(lpd), fdata_pretreat)))
  # METHYLATION
  exprs(lpd)[fData(lpd)$type=="Methylation_Cluster",] =
    exprs(lpd)[fData(lpd)$type=="Methylation_Cluster",]/2
  
  # IGHV: changing name from Uppsala to IGHV
  rownames(lpd)[which(fData(lpd)$type=="IGHV")] = "IGHV"
  
  # GENETICS
  # changing name od del13q any to del13q & remove the rest of del13q (mono, single)
  idx = which(rownames(lpd) %in% c("del13q14_bi", "del13q14_mono"))
  lpd = lpd[-idx,]
  rownames(lpd)[which(rownames(lpd)=="del13q14_any")] = "del13q14"
  
  # remove CHROMOTHRYPSIS
  if("Chromothripsis" %in% rownames(lpd))
    lpd = lpd[-which(rownames(lpd)=="Chromothripsis"),]
  
  # SELECT GOOD SAMPLES
  idx = !is.na(exprs(lpd)["IGHV",])
  if(withMC) idx = idx & !is.na(exprs(lpd)["Methylation_Cluster",])
  # cut out the data to have information about Methylation_Cluster and IGHV for all samples
  lpd = lpd[, idx]
  # for the genetics - remove the genes which do not have enough samples
  which2remove = names(
    which(!apply(exprs(lpd)[rownames(lpd)[fData(lpd)$type %in%
                                            c("gen")],], 1, function(cl) {
    if(all(is.na(cl))) return(FALSE)
    if(sum(is.na(cl)) >= 0.1*length(cl)) return(FALSE)
    tmp = table(cl)
    return(length(tmp)==2 & all(tmp>=minNumSamplesPerGroup))
  })))
  lpd = lpd[-match(which2remove, rownames(lpd)),]
  
  # for the ones with NA put 0 instead
  featOther = rownames(lpd)[fData(lpd)$type %in%
                              c("IGHV", "Methylation_Cluster", "gen", "viab",
                                "pretreat")]
  tmp = exprs(lpd)[featOther,]
  tmp[is.na(tmp)] = 0
  exprs(lpd)[featOther,] = tmp
  
  return(lpd)
}


# wrapper functions to do Lasso model fitting, plotting and prediction
makePredictions = function(drs, frq, lpd, predictorList, lim, std=FALSE,
                           adaLasso = TRUE) {
  res = lapply(names(drs), function(typ) {
    setNames(lapply(drs[[typ]], function(dr) {
      if(typ=="1:5") nm <- paste0(drugs[dr, "name"],
                                  " (average of all concentrations)")
        else if(typ=="4:5") nm <- paste0(drugs[dr, "name"], " (average of ",
                                         paste(conctab[dr,4:5]*1000,
                                               collapse = " and "), " nM)")
      # G & I & M & P
      doLasso(exprs(lpd)[grepl(dr, rownames(lpd)) &
                           fData(lpd)$subtype==typ,, drop=FALSE], 
              predictors=with(predictorList,
                              c(predictorsI, predictorsG, predictorsM,
                                predictorsP)), 
              lpd=lpd,
              suffix=paste0("_","th0", "_c",gsub(":","-",typ)),
              nm=nm, lim=lim)
    }), nm=drs[[typ]])
  })
  return(res)
}


# Function to plot the legends
makeLegends = function(legendFor, colors=coldef) {
  
  # select the colors needed
  colors = colors[names(colors) %in% legendFor]
  
  nleg = length(colors)
  wdths = rep(2, length.out=nleg); hghts = c(2)
  gtl = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))
  n=1

  # M
  if("M" %in% names(colors)) {
    Mgg = ggplot(data=data.frame(x=1, y=factor(c("LP","IP","HP"),
                                               levels=c("LP","IP","HP"))),
                 aes(x=x, y=y, fill=y)) + geom_tile() +
      scale_fill_manual(name="Methylation cluster",
                        values=setNames(colors[["M"]], nm=c("LP","IP","HP"))) +
      theme(legend.title=element_text(size=12),
            legend.text=element_text(size=12))
    gtl = gtable_add_grob(gtl, gtable_filter(ggplotGrob(Mgg), "guide-box"), 1, n)
    n = n+1
  }
  
  # I
  if("I" %in% names(colors)) {
    Igg = ggplot(data=data.frame(x=1,
                                 y=factor(c("unmutated","mutated"),
                                          levels=c("unmutated","mutated"))),
                 aes(x=x, y=y, fill=y)) + geom_tile() +
      scale_fill_manual(name="IGHV",
                        values=setNames(c("gray96",
                                          paste0(colors["I"], c("E5"))),
                                        nm=c("unmutated","mutated"))) +
      theme(legend.title=element_text(size=12),
            legend.text=element_text(size=12))
    gtl = gtable_add_grob(gtl, gtable_filter(ggplotGrob(Igg), "guide-box"), 1, n)
    n = n+1
  }
  
  # G
  if("G" %in% names(colors)) {
    Ggg = ggplot(data=data.frame(x=1,
                                 y=factor(c("wild type","mutated"),
                                          levels=c("wild type","mutated"))),
                 aes(x=x, y=y, fill=y)) + geom_tile() +
      scale_fill_manual(name="Gene",
                        values=setNames(c("gray96",
                                          paste0(colors["G"], c("E5"))),
                                        nm=c("wild type","mutated"))) +
      theme(legend.title=element_text(size=12),
            legend.text=element_text(size=12))
    gtl = gtable_add_grob(gtl, gtable_filter(ggplotGrob(Ggg), "guide-box"), 1, n)
    n = n+1
  }
  
  # P
  if("P" %in% names(colors)) {
    Pgg = ggplot(data=data.frame(x=1,
                                 y=factor(c("no","yes"),
                                          levels=c("no","yes"))),
                 aes(x=x, y=y, fill=y)) + geom_tile() +
      scale_fill_manual(name="Pretreatment",
                        values=setNames(c(colors[["P"]], "white"),
                                        nm=c("yes","no"))) +
      theme(legend.title=element_text(size=12),
            legend.text=element_text(size=12))
    gtl = gtable_add_grob(gtl, gtable_filter(ggplotGrob(Pgg), "guide-box"), 1, n)
    n = n+1
  }
  
  return(list(plot=gtl, width=sum(wdths), height=sum(hghts)))
}

