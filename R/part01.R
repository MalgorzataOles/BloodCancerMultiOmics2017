################################################################################
# Function which plots the drug characteristics
################################################################################

plotPathways = function(dat) {
  
  ordM = sort(table(dat$group), decreasing=TRUE)
  ordS = tapply(dat$"target_category", dat$group, function(pth) {
    sort(table(pth), decreasing=TRUE)
  })
  
  ocur  = ordS[names(ordM)]
  
  # transform the list to LF df
  tmp = do.call(rbind, lapply(names(ocur), function(pathgroup) {
    data.frame(Group=pathgroup,
               Pathway=names(ocur[[pathgroup]]),
               No=as.vector(unname(ocur[[pathgroup]])))
  }))
  tmp$Group = factor(tmp$Group, levels=rev(names(ocur)))
  # order the targets by occurance and leave Other last
  lev = sort(tapply(tmp$No, tmp$Pathway, function(x) sum(x)), decreasing=TRUE)
  lev = c(lev[-which(names(lev)=="Other")], lev["Other"])
  tmp$Pathway = factor(tmp$Pathway, levels=names(lev))
  
  # range to be plotted on the x axis
  widthmax = max(lev)+1
  
  g = ggplot(tmp, aes(x=Pathway, y=No, fill=Group)) +
    geom_bar(width=0.6, stat="identity") +
    theme_bw() + scale_fill_manual(values=pace:::typeColor, name="Drug type") +
    scale_y_continuous(breaks=seq(0,20,2), expand = c(0,0),
                       limits = c(0,widthmax)) +
    xlab("") + ylab("Number of drugs") +
    theme(axis.text.x=element_text(size=12, angle = 45, vjust = 1, hjust=1),
          axis.text.y=element_text(size=12), axis.title.y=element_text(size=12),
          legend.text=element_text(size=12), legend.title=element_text(size=12))
  
  # make gtable
  hghts = c(0.2,0.22*widthmax,1.7)
  wdths = c(1,0.3,0.2,0.26*length(levels(tmp$Pathway)),0.1)
  
  gg = ggplotGrob(g)
  
  gt = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))
  # fill in the gtable
  gt = gtable_add_grob(gt, gtable_filter(gg, "panel"), 2, 4) # panel
  gt = gtable_add_grob(gt, gg$grobs[[whichInGrob(gg, "axis-b")]], 3, 4) # x axis
  gt = gtable_add_grob(gt, gg$grobs[[whichInGrob(gg, "axis-l")]], 2, 3) # y axis
  gt = gtable_add_grob(gt, gg$grobs[[whichInGrob(gg, "ylab-l")]], 2, 2)
  
  # make legend
  wdthsl = c(2)
  hghtsl = c(1.5)
  gtl = gtable(widths=unit(wdthsl, "in"), heights=unit(hghtsl, "in"))
  gtl = gtable_add_grob(gtl, gg$grobs[[whichInGrob(gg, "guide-box")]], 1, 1)
  
  return(list("figure"=list(width=sum(wdths), height=sum(hghts), plot=gt),
              "legend"=list(width=sum(wdthsl), height=sum(hghtsl), plot=gtl)))
}

################################################################################
# Function which plots the patient characteristics as a bar plot
################################################################################
plotPatientStat = function(pats, gap, ptab=patmeta) {
  
  # create plotting data.frame with Diagnosis, Origin and number of cases 
  plotDF = data.frame(table(ptab[pats,"Diagnosis"]))
  colnames(plotDF) = c("Diagnosis","NO")
  plotDF$Diagnosis = as.character(plotDF$Diagnosis)
  plotDF$Origin = 
    names(pace:::diagAmt)[unlist(sapply(plotDF$Diagnosis,
                                        function(x) grep(x, pace:::diagAmt)))]
  
  # set the order of Diagnosis
  ord = smunlist(
    tapply(1:nrow(plotDF), plotDF$Origin,
           function(idx) plotDF$Diagnosis[idx[order(plotDF[idx,"NO"], 
           decreasing=TRUE)]])[names(pace:::diagAmt)])
  plotDF$Diagnosis = factor(plotDF$Diagnosis, levels=ord)
  
  # adjustments for gap
  if(any(plotDF$NO>gap[1] & plotDF$NO<gap[2]))
    stop("Gap is wrongly defined")
  idx = plotDF$NO > gap[2]
  plotDF$NO[idx] = plotDF$NO[idx] - (gap[2]-gap[1])
  
  # round the ceiling to tens (find the latest break point)
  xlimits = c(0, moround(max(plotDF$NO),5)) 
  xbreaks = seq(0, xlimits[2], 10)
  # labels for breaks
  xlabels = ifelse(xbreaks>gap[1], xbreaks + (gap[2]-gap[1]), xbreaks)
  
  g = ggplot() + geom_bar(data=plotDF, aes(x=Diagnosis, y=NO, fill=Origin),
                          stat="identity", colour="black", size=0.1, width=.5) +
    theme_bw() + scale_x_discrete() +
    scale_fill_manual(values=pace:::colDiagL) +
    scale_y_continuous(breaks=xbreaks, labels=xlabels, expand=c(0,0),
                       limits=xlimits) +
    geom_hline(yintercept=c(gap[1]+5,gap[1]+5.5), linetype="dashed", size=0.3) +
    xlab("") + ylab("") +
    theme(axis.text.x=element_text(size=12, angle=45, hjust=1),
          axis.text.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          legend.key.size=unit(0.2,"in"),
          legend.title.align=0.5, legend.text.align=0,
          panel.border=element_rect(color="black", size=0.1))
  
  # construct the gtable
  wdths = c(0.4, 0.4*length(levels(plotDF$Diagnosis)), 0.1) 
  hghts = c(0.2, 0.04*max(xbreaks), 1)
  gt = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))
  ## make grobs
  gg = ggplotGrob(g)
  ## fill in the gtable
  gt = gtable_add_grob(gt, gtable_filter(gg, "panel"), 2, 2)
  gt = gtable_add_grob(gt, gg$grobs[[whichInGrob(gg, "axis-l")]], 2, 1) # y axis
  gt = gtable_add_grob(gt, gg$grobs[[whichInGrob(gg, "axis-b")]], 3, 2) # x axis
  
  # make legend
  wdthsl = c(2)
  hghtsl = c(1.5)
  gtl = gtable(widths=unit(wdthsl, "in"), heights=unit(hghtsl, "in"))
  gtl = gtable_add_grob(gtl, gg$grobs[[whichInGrob(gg, "guide-box")]], 1, 1)
  
  return(list("figure"=list(width=sum(wdths), height=sum(hghts), plot=gt),
              "legend"=list(width=sum(wdthsl), height=sum(hghtsl), plot=gtl)))
  
}

################################################################################
# Function which plots the legends in one row
################################################################################
drawLegends = function(plobj, lng=5, w=2, h=2) { #, alone=FALSE
  
  gt = gtable(widths=unit(rep(w, lng), "in"), heights=unit(h, "in"))
  plotlen = length(plobj)
  
  if(plotlen>lng)
    stop("Number of objects to plot exceeds the number of available slots!")

  for(po in 1:plotlen) {
    gt = gtable_add_grob(gt,
                         plobj[[po]]$grobs[[whichInGrob(plobj[[po]],
                                                        "guide-box")]], 1, po)
  }
  
  grid.draw(gt)
}