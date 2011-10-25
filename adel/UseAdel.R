#
# User front End and utilities for using adel
#
#
#
#Run adel for several dates and a list of plants. Returns a list of string
#
runAdel <- function(dates,plants,pars = list('senescence_leaf_shrink' = 0.5,'startLeaf' = -0.4, 'endLeaf' =1.6, 'stemLeaf' = 1.2,'epsillon' = 1e-6)) {
  out <- vector("list",length(dates))
  for (i in seq(out)) {
    kinlist <- lapply(plants,function(plant) kinLvis(kinL(dates[i],plant,pars),pars))
    desc <- getdesc(kinlist,plants,pars)
    chn <- genString(desc,pars)
    out[[i]] <- cbind(date=dates[i],desc)
  }
  out
}
#
#set Adel from parameters 
#
setAdeluser <- function(devT,geoLeaf,geoAxe,nplants,seed=NULL,xy_db=NULL,sr_db=NULL) {
  setAdel(devT$axeT,devT$dimT,devT$phenT,devT$earT,devT$ssisenT,geoLeaf,geoAxe,nplants,seed,xy_db,sr_db)
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
#
# Compute surfaces from lengths in canopy table
#
leafSurface <- function(shape,scL,scW,from=0,to=1) {
  if (is.na(scL) | is.na(scW) | is.null(shape))
    res <- NA
  else if (to <= from | scL == 0 | scW == 0)
    res <- 0
  else {
    shape[,1] <- shape[,1] / max(shape[,1]) * scL
    shape[,2] <- shape[,2] / max(shape[,2]) * scW
    ws <- approxfun(shape[,1],shape[,2],rule=2)
    res <- integrate(ws,from * scL,to *scL)$value
  }
  res
}
#
canL2canS <- function(canT,sr_db,leaf_shrink=NULL) {
  res <- canT
  if (all(canT$LcIndex <= 1))#leaf shape has not been set by setAdel
    {
      print("canL2canS : can't compute Blade surfaces : SR data not connected to setAdel!!!")
      res[,c("Lv","Lsen","Ll")] <- NA
    }
  else {
    rg <- ifelse(is.na(canT$LcType) | canT$LcType == 0, 1,canT$LcType)
    res[,"Ll"] <- sapply(seq(nrow(canT)),function(x) leafSurface(sr_db[[rg[x]]],canT$Ll[x],canT$Lw[x]))
    #
    base <- canT$Ll - canT$Lsen - canT$Lv
    top <- canT$Ll - canT$Lsen
    res[,"Lv"] <- sapply(seq(nrow(canT)),function(x) leafSurface(sr_db[[rg[x]]],canT$Ll[x],canT$Lw[x],base[x] / canT$Ll[x],top[x] / canT$Ll[x]))
    #
    base <- canT$Ll - canT$Lsen
    res[,"Lsen"] <- sapply(seq(nrow(canT)),function(x) leafSurface(sr_db[[rg[x]]],canT$Ll[x],canT$Lw[x]*canT$LsenShrink[x],base[x] / canT$Ll[x],1))
  }
  names(res)[match(c("Ll","Lsen","Lv"),names(res))] <- c("SLl","SLsen","SLv")
  #
  res[,c("Gl","Gv","Gsen")] <- res[,c("Gl","Gv","Gsen")] * pi * res$Gd^2
  res[,c("El","Ev","Esen")] <- res[,c("El","Ev","Esen")] * pi * res$Ed^2
  names(res)[match(c("Gl","Gsen","Gv","El","Ev","Esen"),names(res))] <- c("SGl","SGsen","SGv","SEl","SEv","SEsen")
  #
  res
}
             
                         
    
    
