#
#
#          R routine for converting a dataframe representation of canopy to String for CanMTG interpreter
#
#
# conventions tissue_ttype :
ttype <- list(green_lamina=1,sen_lamina=2,green_sheath=3,sen_sheath=4,green_internode=5, sen_internode=6,green_peduncle = 7, sen_peduncule=8, green_ear=9, sen_ear=10, green_awn = 11, sen_awn = 12)
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
# 11 = green awn
# 12 = senescent awn
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
#  if (lgreen > epsillon) # no condition, as internode should always be there to allow bearing tillers
    chn <- paste(chn,StemElement(ttype$green_internode,lgreen,diam,diam),sep="")
  #if (lsen > epsillon)
    chn <- paste(chn,StemElement(ttype$sen_internode,lsen,diam,diam))
  chn
}
#
Peduncle <- function(lgreen,lsen,diam,po,pos,epsillon) {
  chn <- ""
  if (lgreen > epsillon)
    chn <- paste(chn,StemElement(ttype$green_peduncle,lgreen,diam,diam),sep="")
  if (lsen > epsillon)
    chn <- paste(chn,StemElement(ttype$sen_peduncle,lsen,diam,diam))
  chn
}
#
Ear <- function(lgreen,lsen,diam,po,pos,epsillon) {
  chn <- ""
  if (lgreen > epsillon)
    chn <- paste(chn,StemElement(ttype$green_ear,lgreen,diam,diam),sep="")
  if (lsen > epsillon)
    chn <- paste(chn,StemElement(ttype$sen_ear,lsen,diam,diam))
  chn
}
#
Awn <- function(lgreen,lsen,diam,po,pos,epsillon) {
  chn <- ""
  if (lgreen > epsillon)
    chn <- paste(chn,StemElement(ttype$green_awn,lgreen,diam,diam),sep="")
  if (lsen > epsillon)
    chn <- paste(chn,StemElement(ttype$sen_awn,lsen,diam,diam))
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

  #if (dat$Ev > epsillon) #not tested to ensure persistence of stems even with 0 length (for tiller attachment)
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
  
  if (!is.null(axil))
    for (i in seq(along=unique(axil$axe_id))) {
      ax <- axil[axil$axe_id == unique(axil$axe_id)[i],]
      azaxil <- 0
      if (i > 1)
        azaxil = azaxil + 170
      chn <- paste(chn,
                 "[",
                 roll(azm),
                 "newAxe")
      nf <- max(1,nrow(ax) - 3)
      for (n in ax$numphy[1:nf]) {
        chn <- paste(chn,Metamer(ax[ax$numphy == n,],epsillon,azaxil))
        azaxil = azaxil + ax$Laz[ax$numphy == n]
      }
	  if (nf > 3) {
                                        #add Peduncle
      if (ax$Ev[nf + 1] > epsillon)
        chn <- paste(chn,
                 Peduncle(ax$Ev[nf + 1]-ax$Esen[nf + 1],ax$Esen[nf + 1],ax$Ed[nf + 1],ax$Epo[nf + 1],ax$Epos[nf + 1],epsillon))
                                        #add ear
      if (ax$Ev[nf + 2] > epsillon)
        chn <- paste(chn,
                 Ear(ax$Ev[nf + 2]-ax$Esen[nf + 2],ax$Esen[nf + 2],ax$Ed[nf + 2],ax$Epo[nf + 2],ax$Epos[nf + 2],epsillon))
    #add awn
      if (ax$Ev[nf + 3] > epsillon)
        chn <- paste(chn,
                 Awn(ax$Ev[nf + 3]-ax$Esen[nf + 3],ax$Esen[nf + 3],ax$Ed[nf + 3],ax$Epo[nf + 3],ax$Epos[nf + 3],epsillon))
      }
	  chn <- paste(chn,
                 "]")
    } 
  chn
}
#
genString <- function(can,pars=list("epsillon" = 1e-6)) {
  epsillon = pars$epsillon
  chn <- ""
  if (!is.null(can)) {#no plant is emerged
    for (p in unique(can$plant)) {
      pl <- can[can$plant == p,]
      axes <- unique(pl$axe)
      axe <- pl[pl$axe == 0,]#main stem
      if (nrow(axe) > 0) {
        azcum <- 0
        chn <- paste(chn,
                 "[ newPlant newAxe")
                                        # stringify main stem and delegates to Metamer axilary branches, if any
        nf <- max(nrow(axe) - 3,1)
        for (m in axe$numphy[1:nf]) {
          if (m %in% axes)
            chn <- paste(chn,Metamer(axe[axe$numphy == m,],epsillon,azcum,pl[pl$axe == m,]))
          else
            chn <- paste(chn,Metamer(axe[axe$numphy == m,],epsillon,azcum))
          azcum <- azcum + axe$Laz[axe$numphy == m]
        }
		if (nf > 3) {
                                        #add Peduncle
        if (axe$Ev[nf + 1] > epsillon)
          chn <- paste(chn,
                 Peduncle(axe$Ev[nf + 1]-axe$Esen[nf + 1],axe$Esen[nf + 1],axe$Ed[nf + 1],axe$Epo[nf + 1],axe$Epos[nf + 1],epsillon))
                                        #add ear
        if (axe$Ev[nf + 2] > epsillon)
          chn <- paste(chn,
                 Ear(axe$Ev[nf + 2]-axe$Esen[nf + 2],axe$Esen[nf + 2],axe$Ed[nf + 2],axe$Epo[nf + 2],axe$Epos[nf + 2],epsillon))
                                        #add awn
        if (axe$Ev[nf + 3] > epsillon)
          chn <- paste(chn,
                 Awn(axe$Ev[nf + 3]-axe$Esen[nf + 3],axe$Esen[nf + 3],axe$Ed[nf + 3],axe$Epo[nf + 3],axe$Epos[nf + 3],epsillon))
        }
		chn <- paste(chn,"]")
      }
    }
  }
  chn
}

