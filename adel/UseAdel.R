#
# User front End and utilities for using adel
#
#
#
#Run adel for several dates and a list of plants. Returns a list of string
#
runAdel <- function(x,p) {
  out <- vector("list",length(x))
  for (t in seq(out)) {
    kinlist <- lapply(p,function(plant) kinLvis(kinL(x[t],plant)))
    desc <- getdesc(kinlist,p)
    out[[t]] <- genString(desc)
  }
  out
}
#
#set Adel from parameters 
#
setAdeluser <- function(devT,geoLeaf,geoAxe,nplants,seed=NULL) {
  setAdel(devT$axeT,devT$dimT,devT$phenT,devT$earT,devT$ssisenT,geoLeaf,geoAxe,nplants,seed)
}
#
                                        #build devT from csv parameter files
#
devTcsv <- function(axeTfile,dimTfile,phenTfile,earTfile=NULL,ssisenTfile=NULL,type=1) {
  reader <- get(ifelse(type==1,"read.csv","read.csv2"))
    #type 1 : "." for decimal, "," for separator
    #type 2:  "," for decimal, ";" for separator
  axeT <- reader(axeTfile)
  
  if (!is.null(earTfile))
    earT <- reader(earTfile)
  else
    earT <- NULL
    
  if (!is.null(ssisenTfile))
    ssisenT <- reader(ssisenTfile)
  else
    ssisenT <- NULL
    
  list(axeT = axeT,
       dimT = reader(dimTfile),
       phenT = reader(phenTfile),
       earT = earT,
       ssisenT = ssisenT)
}
#
readCsv <- function(file,type=1) {
  reader <- get(ifelse(type==1,"read.csv","read.csv2"))
    #type 1 : "." for decimal, "," for separator
    #type 2:  "," for decimal, ";" for separator
  reader(file)
}
#
#
#geoAxe from parameter
#
genGeoAxe <- function(azTM = 75,dazT = 5,incBmM = 2,dincBm = 2,incT = 60,dincT = 5,depMax = 7) {
  list(
       azT = function(a) {
         ifelse(a == 0,
                runif(1) * 360,#plant azimuth
                azTM + (runif(1) - .5) * dazT)
       },
       incT = function(a) {
         ifelse(a == 0,
                incBmM + (runif(1) - .5) * dincBm,
                incT + (runif(1) - .5) * dincT)
       },
       dredT = function(a) {
         #1.5 is an offset to avoid tiller superposed to mainstem
         ifelse(a == 0,
                0,
                1.5 + runif(1) * (depMax-1.5))
       }
       )
}
#
#geoLeaf from parameter (prevoir aussi une boite freeGeomAxe et freegeoLeaf)
#
genGeoLeaf <- function(ntoplim = 4,dazTop = 60,dazBase = 30,topIndex=TRUE) {
   list(
        Azim = function(a,n,ntop) {
          ifelse(ntop <= ntoplim,
                 180 + dazTop * (runif(1) - .5),
                 180 + dazBase * (runif(1) - .5))
               },
        Lindex = function(a,n,ntop) {
          ifelse(topIndex,
                 ntop + 1,
                 n)}
        )
 }
