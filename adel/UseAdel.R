#
# User front End and utilities for using adel
#
#
#
#Run adel for several dates and a list of plants. Returns a list of string
#
runAdel <- function(dates,plants,pars = list('senescence_leaf_shrink' = 0.5,'startLeaf' = -0.4, 'endLeaf' =1.6, 'stemLeaf' = 1.2,'epsillon' = 1e-6)) {
  out <- vector("list",length(dates))
  #deals with python-flatten lists
  if ("axeT" %in% names(plants)) {
    plants <- list(plants)
    names(plants) <- plants[[1]]$refp
  }
  for (i in seq(out)) {
    kinlist <- lapply(plants,function(plant) kinLvis(kinL(dates[i],plant,pars),pars))
    desc <- getdesc(kinlist,plants,pars)
    #chn <- genString(desc,pars)
    if (!is.null(desc))
      out[[i]] <- cbind(date=dates[i],desc)
  }
  out
}
#
#set Adel from parameters 
#
setAdeluser <- function(devT,geoLeaf,geoAxe,nplants,sample='random',seed=NULL,xy_db=NULL,sr_db=NULL) {
  setAdel(devT$axeT,devT$dimT,devT$phenT,devT$earT,devT$ssisenT,geoLeaf,geoAxe,nplants,sample,seed,xy_db,sr_db)
}
#
                                        #build devT from csv parameter files
#
devTcsv <- function(axeTfile,dimTfile,phenTfile,earTfile=NULL,ssisenTfile=NULL,type=1) {
  reader <- get(ifelse(type==1,"read.csv","read.csv2"))
    #type 1 : "." for decimal, "," for separator
    #type 2:  "," for decimal, ";" for separator
  
  #conversion nouvelle nomencalture
  #axeT
  axeT <- reader(axeTfile)
  conv <- c("plant","axe","nf","end","disp","dimIndex","phenIndex","earIndex","emf1","ligf1","senf1","dispf1")
  names(conv) <- c("id_plt","id_axis","N_phytomer","TT_stop_axis","TT_del_axis","id_dim","id_phen","id_ear","TT_em_phytomer1","TT_col_phytomer1","TT_sen_phytomer1","TT_del_phytomer1")
  if (all(conv %in% colnames(axeT)))
    colnames(axeT)[colnames(axeT) %in% conv] <- names(conv)[na.omit(match(colnames(axeT),conv))]
  else if (!all(names(conv) %in% colnames(axeT)))
    stop(paste("axeT : missing data: ",paste(names(conv)[!names(conv) %in% colnames(axeT)],collapse=" ")))
  #dimT
  dimT <- reader(dimTfile)
  conv <- c("index","nrel","Ll","Lw","Gl","Gd","El","Ed")
  names(conv) <- c("id_dim","index_rel_phytomer","L_blade","W_blade","L_sheath","W_sheath","L_internode","W_internode")
  if (all(conv %in% colnames(dimT)))
    colnames(dimT)[colnames(dimT) %in% conv] <- names(conv)[na.omit(match(colnames(dimT),conv))]
  else if (!all(names(conv) %in% colnames(dimT)))
    stop(paste("dimT : missing data: ",paste(names(conv)[!names(conv) %in% colnames(dimT)],collapse=" ")))
  # phenT
  phenT = reader(phenTfile)
  conv <- c("index","nrel","tip","col","ssi","disp")
  names(conv) <- c("id_phen","index_rel_phytomer","dTT_em_phytomer","dTT_col_phytomer","dTT_sen_phytomer","dTT_del_phytomer")
  if (all(conv %in% colnames(phenT)))
    colnames(phenT)[colnames(phenT) %in% conv] <- names(conv)[na.omit(match(colnames(phenT),conv))]
  else if (!all(names(conv) %in% colnames(phenT)))
    stop(paste("phenT : missing data: ",paste(names(conv)[!names(conv) %in% colnames(phenT)],collapse=" ")))
  # earT
  if (!is.null(earTfile)) {
    earT <- reader(earTfile)
    conv <- c("index","em_ear","em_ped","end_gf","l_ped","d_ped","l_ear","Sp_ear","l_ear_awn")
    names(conv) <- c("id_ear","dTT_em_ear","dTT_em_peduncle","TT_z92","L_peduncle","W_peduncle","L_ear","A_ear","L_spike")
    if (all(conv %in% colnames(earT)))
    colnames(earT)[colnames(earT) %in% conv] <- names(conv)[na.omit(match(colnames(earT),conv))]
  else if (!all(names(conv) %in% colnames(earT)))
    stop(paste("earT : missing data: ",paste(names(conv)[!names(conv) %in% colnames(earT)],collapse=" ")))
  }
  else
    earT <- NULL
    
  if (!is.null(ssisenTfile))
    ssisenT <- reader(ssisenTfile)
  else
    ssisenT <- NULL
  list(axeT = axeT,
       dimT = dimT,
       phenT = phenT,
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
canL2canS <- function(canL,sr_db,leaf_shrink=NULL) {
  canS <- canL
  if (all(canS$LcIndex <= 1)) {#leaf shape has not been set by setAdel
    print("canL2canS : can't compute Blade surfaces : SR data not connected to setAdel!!!")
  }
  else {
    rg <- ifelse(is.na(canS$LcType) | canS$LcType == 0, 1,canS$LcType)
    canS$SLl <- sapply(seq(nrow(canS)),function(x) leafSurface(sr_db[[rg[x]]],canS$Ll[x],canS$Lw[x]))
                                        #
    base <- canS$Ll - canS$Lsen - canS$Lv
    top <- canS$Ll - canS$Lsen
    canS$Slv <- sapply(seq(nrow(canS)),function(x) leafSurface(sr_db[[rg[x]]],canS$Ll[x],canS$Lw[x],base[x] / canS$Ll[x],top[x] / canS$Ll[x]))
    #
    base <- canS$Ll - canS$Lsen
    canS$SLsen <- sapply(seq(nrow(canS)),function(x) leafSurface(sr_db[[rg[x]]],canS$Ll[x],canS$Lw[x]*canS$LsenShrink[x],base[x] / canS$Ll[x],1))
  }
  #names(res)[match(c("Ll","Lsen","Lv"),names(res))] <- c("SLl","SLsen","SLv")
  #
  for (w in c("Gl","Gv","Gsen"))
    canS[[paste("S",w,sep="")]] <- canS[[w]] * pi * canS$Gd
  for (w in c("El","Ev","Esen"))
    canS[[paste("S",w,sep="")]] <- canS[[w]] * pi * canS$Ed
  canS
}
             
                         
    
    
