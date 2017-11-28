################################################################################
# Function which:
# - calculates the correlations between drug response profiles
# - produces correlation matrix (either symmetric or 2 combined*)
#
# * for example, for two sets of samples coming from patients with different
#   diagnosis
################################################################################

makeCorrHeatmap = function(mt, mt2=NA, colsc, concNo="one", ctab=conctab,
                           dtab=drugs) {
  
  # cluster the drugs
  mt = mt[apply(mt, 1, sd)!=0, ]
  mt = cor(t(mt))
  hc = hclust(as.dist(1-mt))
  hc = callback(hc, mt)
  ord = rev(hc$label[hc$order])
  
  # prepare data in LF
  mt = meltWholeDF(mt)
  mt$NameX = factor(mt$Y, levels=rev(ord))
  if(concNo=="one") {
    mt$NameY = factor(dtab[mt$X,"name"], levels=dtab[ord,"name"])
  } else if(concNo=="all") {
    mt$NameY = factor(giveDrugLabel(mt$X, ctab, dtab),
                      levels=giveDrugLabel(ord, ctab, dtab))
  }
  
  # take care of the secoond heatmap if present
  if(!is.na(mt2)[1]) {
    # obtain the correlation coefficients
    mt2 = mt2[apply(mt2, 1, sd)!=0, ]
    mt2 = cor(t(mt2))
    
    # order the matrix in the same way as mt is ordered
    mt2 = mt2[levels(mt$NameX), levels(mt$NameX)]
    # put 0 in lower triangle of matrix, include the diagonal
    mt2[upper.tri(mt2, diag=TRUE)] = NA
    
    mt2 = meltWholeDF(mt2)
    # remove NA rows
    mt2 = mt2[!is.na(mt2$Measure),]
    # insert top half of matrix data to plotting data frame
    idx = match(paste(mt2$X, mt2$Y, sep="."), paste(mt$X, mt$Y, sep="."))
    mt[idx, "Measure"] = mt2$Measure
  }
  
  # color the drug names on the axis
  pathColor = pace:::pathColor
  pathColor["Other"] = "black"
  if(concNo=="one") {
    drcol = unname(pathColor[match(dtab[rev(levels(mt$NameX)),"pathway"],
                                   names(pathColor))])
  } else if(concNo=="all") {
    drcol = unname(pathColor[match(dtab[rev(substring(levels(mt$NameX), 1, 5)),
                                        "pathway"], names(pathColor))])
  }
  
  # construct the main plot
  mainPlot = ggplot() +
    geom_tile(data=mt, aes(x=NameX, y=NameY, fill=Measure)) +
    theme_bw() +
    theme(axis.text.y = element_text(size=14, colour=drcol),
          panel.background=element_rect(fill=NA),
          legend.title=element_text(size=14),
          legend.text=element_text(size=14)) +
    scale_fill_gradientn(name="Correlation coefficient", limits=c(-1,1),
                         colours=colsc) +
    xlab("") + ylab("") + coord_fixed()
  
  # construct the dendrogram
  dhc = as.dendrogram(hc)
  ddata = dendro_data(dhc, type = "rectangle")
  ddata$segments = ddata$segments - 0.5
  nXY = nrow(ddata$labels)
  dendPlot = ggplot(segment(ddata)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    theme_dendro() +
    scale_x_continuous(limits=c(0, nXY), expand = c(0,0)) +
    scale_y_reverse()
  
  # construct annotation for drug labels
  annoPlot = ggplot(data.frame(x=1, y=factor(names(pathColor),
                                             levels=names(pathColor))),
                    aes(x, y, fill=y)) +
    geom_tile() +
    scale_fill_manual("Targeted pathway", values=pathColor) +
    theme(legend.title=element_text(size=14), legend.text=element_text(size=14))
  
  # construct gtable
  wdths = c(3, 0.22*nXY, 0.1)
  hghts = c(0.1, 0.22*nXY, 1)
  
  gt = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))
  
  # make grobs
  heat = ggplotGrob(mainPlot)
  denX = ggplotGrob(dendPlot)
  anno = ggplotGrob(annoPlot)
  
  # fill in the gtable
  gt = gtable_add_grob(gt, gtable_filter(heat, "panel"), 2, 2)
  gt = gtable_add_grob(gt, denX$grobs[[whichInGrob(denX, "panel")]], 3, 2)
  gt = gtable_add_grob(gt, heat$grobs[[whichInGrob(heat, "axis-l")]], 2, 1)
  
  # gtable with legends
  wdthsl = c(3, 3)
  hghtsl = c(2.5)
  gtl = gtable(widths=unit(wdthsl, "in"), heights=unit(hghtsl, "in"))
  gtl = gtable_add_grob(gtl, gtable_filter(heat, "guide-box"), 1, 1)
  gtl = gtable_add_grob(gtl, gtable_filter(anno, "guide-box"), 1, 2)
  
  return(list("figure"=list(width=sum(wdths), height=sum(hghts), plot=gt),
              "legend"=list(width=sum(wdthsl), height=sum(hghtsl), plot=gtl)))
}