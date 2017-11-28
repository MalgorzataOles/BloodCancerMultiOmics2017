defineResponseGroups = function(lpd) {
  
  z_factor <- qnorm(0.05, lower.tail = FALSE)
  
  drugnames <- c( "ibrutinib", "everolimus", "selumetinib")
  ib  <- "D_002_4:5" 
  ev  <- "D_063_4:5" 
  se  <- "D_012_4:5"
  stopifnot(identical(fData(lpd)[c(ib, ev, se), "name"], drugnames))
  
  df  <- exprs(lpd)[c(ib, ev, se), lpd$Diagnosis=="CLL"] %>%
    t %>% data.frame %>% `colnames<-`(drugnames)
  df$PatientID=rownames(df)
  mdf <- melt(df)
  
  # Determine standard deviation using mirror method
  pmdf = mdf[mdf$value >= 1,]
  ssd  <- mean( (pmdf$value - 1) ^ 2 ) ^ 0.5
  
  # Normal density
  dn <- tibble(
    x = seq(min(mdf$value), max(mdf$value), length.out = 100),
    y = dnorm(x, mean = 1, sd = ssd) * 2 * nrow(pmdf) / nrow(mdf) 
  )
  
  # Setting up the threshold
  thresh   <- 1 - z_factor * ssd
  
  # Decision rule
  df <- mutate(df,
               group = ifelse(ibrutinib < thresh, "BTK",
                              ifelse(everolimus < thresh, "mTOR",
                                     ifelse(selumetinib < thresh, "MEK",
                                     "none")))
  )
  
  return(data.frame(df[,c("ibrutinib","everolimus","selumetinib","group")],
                    row.names=df$PatientID))
}