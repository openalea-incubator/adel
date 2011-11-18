#
#
#          R routine for converting a dataframe representation of canopy to String for CanMTG interpreter
#
#
# conventions tissue_ttype :
ttype <- list(green_lamina=1,sen_lamina=2,green_sheath=3,sen_sheath=4,green_internode=5, sen_internode=6,green_peduncle = 7, sen_peduncule=8, green_ear=9, sen_ear=10)
# 1 = Green lamina
# 2 = Senescent Lamina
# 3 = Green sheath
# 4 = Senescent sheath1,3,
# 5 = Green internode
# 6 = senescent Internode
# 7 = green peduncule
# 8 = senescent peduncule
# 9 = green ear
# 10 = senescent ear
#
#
StemElement <- function(po, length, dbase, dtop,prec=3) {
  paste("StemElement(",paste(c(po,round(c(length,dbase,dtop),prec)),collapse=","),")",sep="")
}
#
LeafElement <- function(po, Lf, l, wM, srb, srtop,index,seed,incB,prec=3) {
  paste("LeafElement(",paste(c(po,round(c(Lf,l,wM,srb,srtop),prec),index,round(c(seed,incB),prec)),collapse=","),")",sep="")
  #paste("LeafElement(",paste(c(po,round(c(Lf,l,wM,srb,srtop),prec),index,round(c(seed),prec)),collapse=","),")",sep="") # incB removed untill debugging completed
}
indexedLeafElement <- function(po, Lf, l, wM, srb, srtop,index,seed,incB,tindex,prec=3) {
  paste("LeafElement(",paste(c(po,round(c(Lf,l,wM,srb,srtop),prec),index,round(c(seed,incB),prec),tindex),collapse=","),")",sep="")
  #paste("LeafElement(",paste(c(po,round(c(Lf,l,wM,srb,srtop),prec),index,round(c(seed),prec)),collapse=","),")",sep="") # incB removed untill debugging completed
}
#
pitch <- function(inc) paste("+(",round(inc,2),")",sep="")
roll <- function(ang) paste("/(",round(ang,2),")",sep="")
up <- function(ang) paste("^(",round(ang,2),")",sep="")
#
Internode <- function(lgreen,lsen,diam,po,pos,epsillon) {
  chn <- ""
  if (lgreen > epsillon)
    chn <- paste(chn,StemElement(ttype$green_internode,lgreen,diam,diam),sep="")
  if (lsen > epsillon)
    chn <- paste(chn,StemElement(ttype$sen_internode,lsen,diam,diam))
  chn
}
#
Sheath <- function(lgreen,lsen,diam,po,pos,epsillon) {
  chn <- ""
  if (lgreen > epsillon)
    chn <- paste(chn,StemElement(ttype$green_sheath,lgreen,diam,diam),sep="")
  if (lsen > epsillon)
    chn <- paste(chn,StemElement(ttype$sen_sheath,lsen,diam,diam))
  chn
}
#
Blade <- function(lv, lsen, Lf, wM, lsenshrink,lctype, lcindex, incB, po, pos, epsillon,randindex) {
  chn <- ""
  lgreen = lv - lsen
  if (randindex) {#lcindex has not been set by adel)
    if (lgreen > epsillon)
      chn <- paste(chn,
                   LeafElement(ttype$green_lamina, Lf, lv, wM, 0, lgreen / lv,lctype, lcindex, incB))
    if (lsen > epsillon)
      chn <- paste(chn,
                   LeafElement(ttype$sen_lamina, Lf, lv, wM * lsenshrink, lgreen / lv, 1, lctype, lcindex, incB))
  } else {
    if (lgreen > epsillon)
      chn <- paste(chn,
                   indexedLeafElement(ttype$green_lamina, Lf, lv, wM, 0, lgreen / lv,lctype, lcindex, incB,lcindex))
    if (lsen > epsillon)
      chn <- paste(chn,
                   indexedLeafElement(ttype$sen_lamina, Lf, lv, wM * lsenshrink, lgreen / lv, 1, lctype, lcindex, incB,lcindex))
  }
 
  chn
}
#
Metamer <- function(dat,epsillon,azcum,axil = NULL) {
  
  azm <- (azcum + dat$Laz) %% 360
  
  chn <- "newMetamer"
  if (abs(dat$Einc) > 0)
    chn <- paste(chn,
                 up(dat$Einc))
  
  if (!is.null(axil)) {
    azaxil <- 0
    chn <- paste(chn,
                 "[",
                 roll(azm),
                 "newAxe")
    for (n in axil$numphy) {
      chn <- paste(chn,Metamer(axil[axil$numphy == n,],epsillon,azaxil))
      azaxil = azaxil + axil$Laz[axil$numphy == n]
    }
    chn <- paste(chn,
                 "]")
  }
  
 if (dat$Ev > epsillon)
  chn <- paste(chn,
                 Internode(dat$Ev-dat$Esen,dat$Esen,dat$Ed,dat$Epo,dat$Epos,epsillon))
  if (abs(dat$Ginc) > 0)
    chn <- paste(chn,
                 up(dat$Ginc))
  if (dat$Gv > epsillon)
    chn <- paste(chn,
                 Sheath(dat$Gv-dat$Gsen,dat$Gsen,dat$Gd,dat$Gpo,dat$Gpos,epsillon))
  if (dat$Lv > epsillon) 
    chn <- paste(chn,
                 "[",
                 roll(azm),
                 Blade(dat$Lv, dat$Lsen, dat$Ll,dat$Lw,dat$LsenShrink,dat$LcType,dat$LcIndex,dat$Linc,dat$Lpo,dat$Lpos,epsillon,all(dat$LcIndex <= 1)),
                 "]")
  chn
}
#
genString <- function(can,pars=list("epsillon" = 1e-6)) {
  epsillon = pars$epsillon
  chn <- ""
  for (p in unique(can$plant)) {
    pl <- can[can$plant == p,]
    axes <- unique(pl$axe)
    axe <- pl[pl$axe == axes[1],]
    azcum <- 0
    chn <- paste(chn,
                 "[ newPlant newAxe")
    for (m in axe$numphy) {
      if (m %in% axes)
        chn <- paste(chn,Metamer(axe[axe$numphy == m,],epsillon,azcum,pl[pl$axe == m,]))
      else
        chn <- paste(chn,Metamer(axe[axe$numphy == m,],epsillon,azcum))
      azcum <- azcum + axe$Laz[axe$numphy == m]
    }
   chn <- paste(chn,"]")
  }
  chn
}

