###########################################################
### ADDITIONAL FUNCTIONS
###########################################################


# Wide to long format conversion
meltWholeDF = function(df) {
  data.frame(X=rep(colnames(df), each=nrow(df)),
             Y=rep(rownames(df), times=ncol(df)),
             Measure=as.vector(as.matrix(df)))
}

# Smart unlist - preserves the names of the vectors as they were
smunlist = function(li) {
  setNames(unlist(li),
           nm=unname(rapply(li,
                            function(i)
                              if(!is.null(names(i)))
                                names(i)
                            else
                              rep("", length(i)))))
}

# Function that gives nice drug names out of IDs (D_001-40 or D_001_5 type)
giveDrugLabel = function(drid, ctab, dtab) {
  vapply(strsplit(drid, "-"), function(x) {
    if(length(x)==2) { # ID type: D_001-40
      paste0(dtab[x[1],"name"], " ", x[2], " \u00B5", "M")
    } else if(length(x)==1) { # ID type: D_001_5
      x = unlist(strsplit(drid, "_"))
      k = paste(x[1:2], collapse="_")
      paste0(dtab[k, "name"], " ", ctab[k, as.integer(x[3])], " \u00B5","M")
    }}, character(1))
}

# Round up to the nearest 5
moround = function(x,base) {
    base*ceiling(x/base)
}

# capitalize the first letter
toCaps = function(word) {
  paste0(toupper(substring(word,1,1)), substring(word,2,nchar(word)))
}

# out of IDs likes D_001_1, strip the trailing '_1'
stripConc <- function(x) 
  vapply(strsplit(x, "_"), function(x)
    paste(x[-length(x)], collapse="_"), character(1))

# treshold an array from below and above
deckel <- function(x, lower = -Inf, upper = +Inf)
  ifelse(x<lower, lower, ifelse(x>upper, upper, x))

# log10 scale labels in ggplot2, use: scale_x_log10(labels=scientific_10)
scientific_10 = function(x) {
    x = scientific_format()(x)
    parse(text=ifelse(x=="1e+00", "1   ", gsub("1e", "10^", x)))
}

# labels of volcano x axis scale
percentAxisScale = function(x) {
  x*100
}

# Function which computes log10 and returns with the sign of the input value
log10div = function(x) {
    sign(x)*log10(abs(x))  
}

# Function for axis labels of p-values going in two directions
# (sensitive/resistant)
exp10div = function(x) {
    x = -abs(x)
    x = paste0("10^", x)
    x = gsub("10^0", "1", x, fixed=TRUE)
    parse(text=x)
}

# change color names to hex with alpha level
col2hex = function(cols, alpha=1, names=NA) {
  tmp = col2rgb(cols) 
  max = 255
  tmp = apply(tmp, 2, function(t)
    rgb(red=t[1], green=t[2], blue=t[3], maxColorValue=max, alpha=alpha*max))
  if (all(!is.na(names)) && length(names)==length(tmp))
      tmp = setNames(tmp,nm=names)
  tmp
}

# safe match
safeMatch <- function(x, ...) {
  rv <- match(x, ...)
  if (any(is.na(rv)))
    stop(sprintf("`match` failed to match %s", paste(x[is.na(rv)],
                                                     collapse=", ")))
  rv 
}

# find out the index of appropriate layer in a grob table
whichInGrob = function(grob, layer) {
  match(layer, grob[["layout"]][["name"]])
}