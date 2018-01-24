################################################################################
# VOLCANO PLOTS
################################################################################

#_______________________________________________________________________________
## color by number of sign concentration steps / one dot per drug
#-------------------------------------------------------------------------------

## MAIN FUNCTION
ggvolc2 = function(df, title, Ycut, color=NA, xlab, maxX, maxY, expY, hghBox,
                   axisMarkY, breaksX, arrowLength, Xhang, minConc, fixedHght) {

  # quiets concerns of R CMD check "no visible binding for global variable"
  X=NULL; Y=NULL; labX=NULL; labY=NULL; Label=NULL;
  hjust=NULL; x=NULL; y=NULL; xend=NULL; yend=NULL;
  Color=NULL; .x=NULL; Fill=NULL                 
                   
  # color palette
  pal = setNames(rev(c(brewer.pal(11, "RdYlBu")[1:5], "#000000",
                       brewer.pal(11, "RdBu")[7:11])), nm=1:11)

  # check if dataq.frame have required columns
  stopifnot(all(c("X","Y","Label","Grey") %in% colnames(df)))

  # check if colors are for the palette # if not use the default
  col.default = c(pal[2], pal[9])
  color = if(all(color %in% colors())) color else  col.default

  # x axis
  Xlims = c(-maxX, maxX)

  # y axis
  Ylims = c(-0.05, maxY)

  # direction of the effect
  df$Direction = 0
  df$Direction[df$Y>=Ycut & df$X<0] = -1
  df$Direction[df$Y>=Ycut & df$X>0] = 1

  # add column saying if association is significant
  df$IsSignificant = df$Y>=Ycut & !df$Grey

  # decide what should have label
  df$haveLabel = df$IsSignificant & df$signConc >= minConc

  # positions for labels
  calcY = function(y.org) {
    rng = c(Ycut+0.1,maxY-0.1)
    if(length(y.org)==1) {
      y.org = rng[1]+(rng[2]-rng[1])/2
    } else {
      inc = (rng[2]-rng[1])/length(y.org)
      newY = rng[1]+inc*(1:length(y.org))
      y.org[order(y.org)] = newY
    }
    y.org
  }
  df$labY[df$haveLabel] = calcY(y.org=df$Y[df$haveLabel])
  df$labX = with(df, ifelse(haveLabel,
                            ifelse(Direction==1, Xlims[2]-Xhang,
                                   Xlims[1]+Xhang), NA))
  df$hjust = ifelse(sign(df$labX)==1, 0, 1)
  labNo = length(df$Y[df$haveLabel])
  df$Color = factor(6+(df$signConc*df$Direction), levels=1:11)
  # make grey dots really grey
  df$Color[df$Grey] = 6

  # construct the plot
  if(sum(df$haveLabel)>0) {
    # construct the labels
    df2 = df[df$haveLabel,]
    gg = ggplot() +
      geom_segment(data=df2, aes(x=X, y=Y, xend=labX, yend=labY),
                   colour="darkgrey", alpha=0.7, linetype="longdash", size=0.1) +
      geom_text(data=df2, mapping=aes(x=labX, y=labY, label=Label, hjust=hjust),
                size=2.5) # dotted, size=2.5
    # and arrows
    gg = gg +
      geom_segment(data=data.frame(x=0, y=0, xend=-arrowLength, yend=0),
                   aes(x=x,y=y,xend=xend,yend=yend),
                   arrow=arrow(length = unit(0.25, "cm"), type="open"),
                   colour=color[1], size=0.8) +
      geom_segment(data=data.frame(x=0, y=0, xend=arrowLength, yend=0),
                   aes(x=x,y=y,xend=xend,yend=yend),
                   arrow=arrow(length = unit(0.25, "cm"), type="open"),
                   colour=color[2], size=0.8)
  } else {
    gg = ggplot()
  }
  gg = gg + geom_hline(yintercept=Ycut, colour="dimgrey", linetype="dashed") +
    geom_vline(xintercept=0, colour="dimgrey", size=0.1) +
    geom_point(data=df, aes(x=X, y=Y, fill=Color, color=Color), shape=21,
               size=3) +
    scale_fill_manual(values=col2hex(pal, alpha=0.7)) +
    scale_color_manual(values=pal) + theme_bw() +
    xlab(xlab) + ggtitle(title) +
    scale_y_continuous(
      expression(italic(p)*"-value"), breaks=seq(0,maxY,axisMarkY),
      labels=math_format(expr=10^.x)(-seq(0,maxY,axisMarkY)), limits=Ylims,
      expand = c(expY,0)) +
    scale_x_continuous(limits=Xlims, breaks=breaksX, labels=percentAxisScale) +
    theme(axis.text.y=element_text(size=8), axis.text.x=element_text(size=8),
          axis.title.x=element_text(size=8), axis.title.y=element_text(size=8),
          plot.title= element_text(size=8))

  # construct the gtable
  wdths = c(0.2, 0.4, 3.5*(Xlims[2]-Xlims[1]), 0.2)
  hghts = c(0.3, ifelse(is.na(fixedHght), hghBox*nrow(df2), fixedHght), 0.3, 0.2)
  gt = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))
  ## make grobs
  ggr = ggplotGrob(gg)
  ## fill in the gtable
  gt = gtable_add_grob(gt, gtable_filter(ggr, "panel"), 2, 3)
  gt = gtable_add_grob(gt, ggr$grobs[[whichInGrob(ggr, "axis-l")]], 2, 2)
  gt = gtable_add_grob(gt, ggr$grobs[[whichInGrob(ggr, "axis-b")]], 3, 3)
  gt = gtable_add_grob(gt, ggr$grobs[[whichInGrob(ggr, "xlab-b")]], 4, 3)
  gt = gtable_add_grob(gt, ggr$grobs[[whichInGrob(ggr, "ylab-l")]], 2, 1)
  gt = gtable_add_grob(gt, gtable_filter(ggr, "title"), 1, 3) # title

  # legend # it should adjust depending on minConc
  pal = setNames(pal[-ceiling(length(pal)/2)], nm=c(-5:-1,1:5))

  fakeDF = data.frame(X=1:length(pal), Y=LETTERS[1:length(pal)], Fill=names(pal))

  gl = ggplot() +
    geom_bar(data=fakeDF, aes(x=X, y=Y, fill=Fill),
             stat="identity", position="identity") +
    scale_fill_manual(name="Number of\nsignificant\nconcentrations",
                      values=pal,
                      labels=c(`-5`="", `-4`="", `-3`="", `-2`="", `-1`="",
                               `1`="1 conc.", `2`="2 conc.", `3`="3 conc.",
                               `4`="4 conc.", `5`="5 conc"), drop = FALSE) +
    guides(fill=guide_legend(ncol=2)) +
    theme(legend.text=element_text(size=8),
          legend.title=element_text(size=8, face="bold"),
          legend.title.align=0.5)

  # construct the gtable
  wdthsL = c(2)
  hghtsL = 2
  gtL = gtable(widths=unit(wdthsL, "in"), heights=unit(hghtsL, "in"))
  ## make grobs
  ggl = ggplotGrob(gl)
  ## fill in the gtable
  gtL = gtable_add_grob(gtL, ggl$grobs[[whichInGrob(ggl, "guide-box")]], 1, 1)

  return(list("figure"=list(width=sum(wdths), height=sum(hghts), plot=gt),
              "legend"=list(width=sum(wdthsL), height=sum(hghtsL), plot=gtL)))
}


## RUNNING FUNCTION

run.ggvolcGr2 = function(results, effects, screen, mts, fdr, maxX, maxY,
                         expY, hghBox, axisMarkY, breaksX, arrowLength,
                         Xhang, minConc, dtab=BloodCancerMultiOmics2017::drugs, fixedHght=NA) {

  # select appropriate data to plot
  results = results[[screen]]
  effects = effects[[screen]]

  # merge the results and effects
  reseff = merge(results, effects, by=c("DrugID","TestFac","FacDr"))

  # filter the data to only those lines which will be plotted
  reseff = reseff[reseff$TestFac %in% mts,]

  # mark significant cases
  reseff$fdr = reseff$adj.pval <= fdr
  # find out the p-value threshold
  cut = max(reseff$pval[reseff$fdr])

  # iterate on each mutation
  waste = tapply(1:nrow(reseff), reseff$TestFac, function(idx) {
    # mutation to be plotted
    mt = reseff[idx[1], "TestFac"]
    # select the data to be plotted
    re = reseff[idx,]
    ## for each drug select the association with the biggest difference in effects
    re = do.call(rbind, tapply(1:nrow(re), re$DrugID2, function(i) {
      tmp = re[i,]
      if(all(tmp$fdr==FALSE)) {
        return(cbind(tmp[which.max(abs(tmp$WM)),], signConc=0))
      } else {
        tmp = tmp[tmp$fdr,]
        return(cbind(tmp[which.max(abs(tmp$WM)),], signConc=nrow(tmp)))
      }
    }))

    # create input data frame
    plotDF = with(re, data.frame(X=(-1)*WM, Y=-log10(pval),
                                 Label=dtab[DrugID2,"name"],
                                 Grey= mean.0 <0.1 & mean.1 <0.1 & fdr,
                                 signConc))
    # maxY
    if(is.na(maxY)) maxY=max(ceiling(plotDF$Y))
    # additional paremeters (labels, colors)
    xlab = "Difference of effects"
    # run the plotting function
    ggvolc2(df=plotDF, title=mt, Ycut=-log10(cut), color=NA, xlab=xlab,
            maxX=maxX, maxY=maxY, expY=expY, hghBox=hghBox,
            axisMarkY=axisMarkY, breaksX=breaksX, arrowLength=arrowLength,
            Xhang=Xhang, minConc=minConc, fixedHght=fixedHght)
  })

  # names(waste) = paste(names(waste), screen, sep=".")
  waste
}


################################################################################
# HEAT MAP
################################################################################

ggheat = function(results, effects, dtab=BloodCancerMultiOmics2017::drugs, ctab=BloodCancerMultiOmics2017::conctab) {

  # quiets concerns of R CMD check "no visible binding for global variable"
  diag.conc=NULL; Drug2=NULL; Fill=NULL; x=NULL; fill=NULL

  # COLORS
  # color palette
  pal = c(brewer.pal(11, "PiYG")[2], brewer.pal(11, "RdBu")[10]) # only 2 colors!
  # border of whole plot
  cPanb = "black"
  sPanb = 0.5
  # minor & major grid
  cGrid = c("grey", "dimgrey")

  # apply the FDR threshold
  FDR = 0.1

  plotDF = do.call(rbind, lapply(names(results), function(nm) {
    df = merge(results[[nm]], effects[[nm]], by=c("FacDr", "TestFac", "DrugID"))
    df$plot.pval = ifelse(df$adj.pval <= FDR, df$pval*sign(df$WM), 1)
    df$plot.simple.pval = ifelse(df$adj.pval <= FDR, 1*sign(df$WM), 0)
    df$plot.simple.pval = factor(df$plot.simple.pval, levels=c(-1,1,0))
    df
  }))

  # add column which orders the coulumns of the plot
  plotDF$diag = sapply(plotDF$TestFac, function(tf)
    strsplit(tf, ".", fixed=TRUE)[[1]][[2]])
  plotDF$conc = sapply(plotDF$DrugID, function(drid) {
    splits = strsplit(drid, "-")[[1]]
    colnames(ctab)[which(ctab[splits[1],] == splits[2])]
  })
  plotDF$diag.conc = paste(plotDF$diag, plotDF$conc, sep=".")
  diagOrder = names(colDiagS)[names(colDiagS) %in% unique(plotDF$diag)]
  plotDF$diag.conc = factor(plotDF$diag.conc,
                            levels=as.vector(sapply(diagOrder,
                                                    function(d)
                                                      paste(d, paste0("c", 1:5),
                                                            sep="."))))

  # names of drugs / order them
  plotDF$Drug2 = dtab[plotDF$DrugID2, "name"]
  plotDF$Drug2 = factor(plotDF$Drug2, levels=dtab[rev(mainDrOrd),"name"])

  # data frame for vline
  vltab = table(sapply(unique(plotDF$diag), function(x) grep(x, diagAmt)))
  vldf = data.frame(x=seq(5.5, length(levels(plotDF$diag.conc)), 5), col=1)
  vldf$col[cumsum(vltab[1:(length(vltab)-1)])] = 2 # remove latest
  vldf$col = factor(vldf$col, levels=c(1,2))

  plotDF$Fill = log10div(plotDF$plot.pval) # transform the values to log10 scale
  plotDF$Fill = pmax(pmin(plotDF$Fill, 12), -12) # censore
  # color palette
  pal = rev(c(brewer.pal(11,"PiYG")[1:6], brewer.pal(11,"RdBu")[7:11]))
  pal[6] = "gray95"

  gMain2 = ggplot() + geom_tile(data=plotDF,
                                aes(x=diag.conc, y=Drug2, fill=Fill)) +
    scale_fill_gradientn(expression(italic(p)*"-value"),
                         colours=pal, limits=c(-12, 12),
                         breaks=c(-12,-8,-4,0,4,8,12), labels=exp10div) +
    theme_bw() +
    geom_vline(xintercept=seq(0.5, length(levels(plotDF$diag.conc))+1, 1),
               color="white", size=1.5) +
    geom_hline(yintercept=seq(0.5, length(levels(plotDF$Drug2))+1, 1),
               color="white", size=1.5) +
    geom_vline(data=vldf, aes(xintercept=x, color=col)) +
    scale_color_manual(values=c("1"=cGrid[1], "2"=cGrid[2])) +
    coord_equal() + xlab("") + ylab("") +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    theme(axis.title=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=10),
          axis.ticks.x=element_blank(),
          panel.border=element_rect(colour=cPanb, fill=NA, size=sPanb),
          legend.key=element_rect(color="black"),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10, face="bold")) +
    guides(color=FALSE)

  #################################################

  # concentrations
  concDF = data.frame(x=factor(levels(plotDF$diag.conc)),
                      fill=factor(rep(paste0("c",1:5)), levels=paste0("c",1:5)))
  concpal = setNames(paste0("grey", seq(35, 95, 15)), nm=paste0("c",1:5))
  gConc = ggplot() + geom_tile(data=concDF, aes(x=x, y=1, fill=fill)) +
    scale_fill_manual("Drug\nconcentration", values=concpal) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) + theme_bw() +
    theme(axis.title=element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank(), legend.text=element_text(size=10),
          legend.title=element_text(size=10, face="bold"),
          legend.key=element_rect(color="black"),
          panel.border=element_rect(colour=cPanb, fill=NA, size=sPanb)) +
    geom_vline(data=vldf, aes(xintercept=x, color=col)) +
    scale_color_manual(values=c("1"=cGrid[1], "2"=cGrid[2])) +
    guides(color=FALSE)

  # annotation of cell of origin
  diagOrd = sapply(unique(plotDF$diag), function(x) grep(x, diagAmt))
  diagOrd = diagOrd[diagOrder]
  annoDF = data.frame(x=names(diagAmt)[diagOrd],
                      diag=factor(names(diagOrd), levels=names(diagOrd)))
  annoDF$x = factor(annoDF$x, levels=unique(annoDF$x))
  # vline
  vldf2 = data.frame(x=seq(1.5, length(levels(annoDF$diag)), 1), col=1)
  vldf2$col[cumsum(table(diagOrd)[1:(length(table(diagOrd))-1)])] = 2
  vldf2$col = factor(vldf2$col, levels=c(1,2))
  gAnno = ggplot() + geom_tile(data=annoDF, aes(x=diag, y=1, fill=x)) +
    scale_fill_manual("Disease origin", values=colDiagL) +
    scale_y_continuous(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
    theme_bw() +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10, face="bold"),
          panel.border=element_rect(colour=cPanb, fill=NA, size=sPanb),
          legend.key=element_rect(color="black")) +
    geom_vline(data=vldf2, aes(xintercept=x, color=col)) +
    scale_color_manual(values=c("1"=cGrid[1], "2"=cGrid[2])) +
    guides(color=FALSE) +
    geom_text(data=annoDF, aes(label=diag, y=1, x=diag))


  ## MERGE PLOTS INTO ONE FIGURE

  # prepare grobs
  ggM = ggplotGrob(gMain2)
  ggA = ggplotGrob(gAnno)
  ggC = ggplotGrob(gConc)

  # prepare gtable
  nX = length(levels(plotDF$diag.conc))
  nY = length(levels(plotDF$Drug2))
  unM = 0.15 # unit for main heat map
  unA = 0.25 # unit for y-axis of the annotation - diagnosis
  unC = 0.15 # unit for y-axis of the concentration
  sp = 0.0 # space

  wdths = c(1.5, nX*unM, 2)
  hghts = c(unA, sp, nY*unM, unC)
  gt = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))

  # add pieces to gtable
  gt = gtable_add_grob(gt, ggA$grobs[[whichInGrob(ggA, "panel")]], 1, 2)
  gt = gtable_add_grob(gt, ggM$grobs[[whichInGrob(ggM, "panel")]], 3, 2)
  gt = gtable_add_grob(gt, ggC$grobs[[whichInGrob(ggC, "panel")]], 4, 2)
  gt = gtable_add_grob(gt, ggM$grobs[[whichInGrob(ggM, "axis-l")]], 3, 1)

  # legend cell origin & concentration & effect
  wdthsL = c(2,2,2)
  hghtsL = 2
  gtL = gtable(widths=unit(wdthsL, "in"), heights=unit(hghtsL, "in"))
  gtL = gtable_add_grob(gtL, ggA$grobs[[whichInGrob(ggA, "guide-box")]], 1, 1)
  gtL = gtable_add_grob(gtL, ggM$grobs[[whichInGrob(ggM, "guide-box")]], 1, 2)
  gtL = gtable_add_grob(gtL, ggC$grobs[[whichInGrob(ggC, "guide-box")]], 1, 3)

  return(list("figure"=list(width=sum(wdths), height=sum(hghts), plot=gt),
              "legend"=list(width=sum(wdthsL), height=sum(hghtsL), plot=gtL)))
}


################################################################################
# BEE SWARM PLOTS
################################################################################


beeF <- function(drug, mut, cs, diag, y1, y2, custc,
                 lpd=BloodCancerMultiOmics2017::lpdAll,
                 ctab=BloodCancerMultiOmics2017::conctab,
                 dtab=BloodCancerMultiOmics2017::drugs) {

  col1 <- vector(); col2 <- vector()
  dr <- lpd[ , lpd$Diagnosis %in% diag   ]
  p = t.test( exprs(dr)[drug,] ~ exprs(dr)[mut,], var.equal = TRUE)$p.value

  #clonsize
  af <- fData(BloodCancerMultiOmics2017::mutCOM)[colnames(dr), paste0(mut, "cs")]

  #Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('coral1','blue4'))

  # This adds a cont. color code
  if (cs==T) {
    col2 <- rbPal(100)[as.numeric(cut( af, breaks = 100 ))]
  } else {
    col2 <- "blue4"
  }

  if (custc==T) {
    col1 <- ifelse(exprs(dr)[mut,]==1, col2,"coral1")
  } else {
    col1 <- ifelse(exprs(dr)[mut,]==1, "green","magenta")
  }

  beeswarm( exprs(dr)[drug,] ~ exprs(dr)[mut,],
            method = 'swarm',
            pch = 19, pwcol = col1, cex = 1.6,
            xlab = '', ylab = 'Viability', cex.axis=1.6, cex.lab=1.9,
            ylim=c(y1,y2),
            labels = c("AF=0", "AF>0"),
            main=(paste0( giveDrugLabel(drug, ctab, dtab), " ~ ", mut,
                          "\n (p = ",
                          digits = format.pval(p,
                                               max(1, getOption("digits") - 4)),
                          ")")),
            cex.main=2.0,
            bty="n"
  )

  boxplot(exprs(dr)[drug,] ~ exprs(dr)[mut,],  add = T, names = c("",""),
          col="#0000ff22", axes = 0, outline=FALSE)
}


# BEESWARM FOR PRETREATMENT
beePretreatment = function(lpd, drug, y1, y2, fac, val, name) {

  dr <- lpd[  , exprs(lpd)[fac,] %in% val]
  pretreat <- BloodCancerMultiOmics2017::patmeta[colnames(dr),
    "IC50beforeTreatment"]
  p = t.test( exprs(dr)[drug,] ~ pretreat, var.equal = TRUE)$p.value

  beeswarm( exprs(dr)[drug,] ~ pretreat,
            method = 'swarm',
            pch = 19,  cex = 1.2, pwcol=ifelse(pretreat, "blue4", "coral1"),
            xlab = '', ylab = 'Viability', cex.axis=1.6, cex.lab=1.8,
            labels = c("pre-treated", "not treated"),
            main=paste0(
              name," (p = ", digits = format.pval(
                p, max(1, getOption("digits") - 4)), ")"),  ylim=c(y1,y2),
            cex.main=2,
            bty="n"
  )
  boxplot(exprs(dr)[drug,] ~ pretreat,  add = T, names = c("",""),
          col="#0000ff22", axes = 0, outline=FALSE)
}
