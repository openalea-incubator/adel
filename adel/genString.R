#
#
#          R routine for converting a dataframe representation of canopy to String for CanMTG interpreter
#
#
StemElement <- function(po, length, dbase, dtop,prec=3) {
  paste("StemElement(",paste(c(po,round(c(length,dbase,dtop),prec)),collapse=","),")",sep="")
}
#
LeafElement <- function(po, Lf, l, wM, srb, srtop,index,seed,incB,prec=3) {
  paste("LeafElement(",paste(c(po,round(c(Lf,l,wM,srb,srtop),prec),index,round(c(seed,incB),prec)),collapse=","),")",sep="")
}
#
pitch <- function(inc) paste("+(",round(inc,2),")",sep="")
roll <- function(ang) paste("/(",round(ang,2),")",sep="")
#
Internode <- function(lgreen,lsen,diam,po,pos,epsillon) {
  chn <- ""
  if (lgreen > epsillon)
    chn <- paste(chn,StemElement(po,lgreen,diam,diam),sep="")
  if (lsen > epsillon)
    chn <- paste(chn,StemElement(pos,lsen,diam,diam))
  chn
}
#
Sheath <- function(lgreen,lsen,diam,po,pos,epsillon) {
  chn <- ""
  if (lgreen > epsillon)
    chn <- paste(chn,StemElement(po,lgreen,diam,diam),sep="")
  if (lsen > epsillon)
    chn <- paste(chn,StemElement(pos,lsen,diam,diam))
  chn
}
#
Blade <- function(lv, lsen, Lf, wM, lctype, lcindex, incB, po, pos, epsillon) {
  chn <- ""
  lgreen = lv - lsen
  if (lgreen > epsillon)
    chn <- paste(chn,
                 LeafElement(po, Lf, lv, wM, 0, lgreen / lv,lctype, lcindex, incB))
  if (lsen > epsillon)
    chn <- paste(chn,
                 LeafElement(pos, Lf, lv, wM / 2, lgreen / lv, 1, lctype, lcindex, incB))
  chn
}
#
Metamer <- function(dat,epsillon) {
  chn <- "newMetamer"
  if (abs(dat$Einc) > 0)
    chn <- paste(chn,
                 pitch(dat$Einc))
  if (dat$Ev > epsillon)
    chn <- paste(chn,
                 Internode(dat$Ev-dat$Esen,dat$Esen,dat$Ed,dat$Epo,dat$Epos,epsillon))
  if (abs(dat$Ginc) > 0)
    chn <- paste(chn,
                 pitch(dat$Ginc))
  if (dat$Gv > epsillon)
    chn <- paste(chn,
                 Sheath(dat$Gv-dat$Gsen,dat$Gsen,dat$Gd,dat$Gpo,dat$Gpos,epsillon))
  if (dat$Lv > epsillon) 
    chn <- paste(chn,
                 "[",
                 roll(dat$Laz),
                 Blade(dat$Lv, dat$Lsen, dat$Ll,dat$Lw,dat$LcType,dat$LcIndex,dat$Linc,dat$Lpo,dat$Lpos,epsillon),
                 "]")
  chn
}
#
genString <- function(can,epsillon = 1e-6) {
  chn <- ""
  for (p in unique(can$plant)) {
    pl <- can[can$plant == p,]
    chn <- paste(chn,
                 "[ newPlant")
    for (a in unique(pl$axe)) {
      axe <- pl[pl$axe == a,]
      if (a == 0)
        azAxe <- 0
      else
        azAxe <- pl$Laz[pl$axe == 0 & pl$numphy == a]
      incB <- axe$Einc[axe$numphy == 1]
      axe$Einc[axe$numphy == 1] <- 0#do not turn again at the base of phytomer 1
      chn <- paste(chn,
                   roll(azAxe),#to make axe emerge from within leaf
                   "[ newAxe",
                   pitch(incB))
      last <- max(which(apply(axe[,c("Lv","Gv","Ev")],1,sum) > epsillon))
      for (n in 1:last)
        chn <- paste(chn,Metamer(axe[axe$numphy == n,],epsillon))
      chn <- paste(chn,"]")
      }
    chn <- paste(chn,"]")
  }
  chn
}

