
############### Simulations ####################################################

ReadOut <- function (folder) {

folder=paste(folder,"/",sep="")
  Nph.tab <<- vector("list",5)
  Nstruct.tab <<- vector("list",5)
  Agreen.tab <<- vector("list",5)
  Nphm2.tab <<- vector("list",5)
  DMrem.tab <<- vector("list",5)
DMstruct.tab <<- vector("list",5)
  PAR.tab <<- vector("list",5)
GAI.tab <<- vector("list",5)
  DegRate.tab <<- vector("list",5)
SynthRate.tab <<- vector("list",5)
DMDegRate.tab <<- vector("list",5)
DMSynthRate.tab <<- vector("list",5)
Photo.tab <<- vector("list",5)
  
  ty <<- 0
  for (type in c("Lamina","Sheath","Internode","Peduncle","Chaff")){
    ty <<- ty+1
    Nstruct.tab[[ty]] <<- read.table(paste(folder,"NstructPred",type,".txt",sep=""),header=FALSE)
    Nph.tab[[ty]] <<- read.table(paste(folder,"NphPred",type,".txt",sep=""),header=FALSE)
    Agreen.tab[[ty]] <<- read.table(paste(folder,"AgreenPred",type,".txt",sep=""),header=FALSE)
    Nphm2.tab[[ty]] <<- read.table(paste(folder,"Nphm2Pred",type,".txt",sep=""),header=FALSE)
    DMrem.tab[[ty]] <<- read.table(paste(folder,"DMremPred",type,".txt",sep=""),header=FALSE)
    DMstruct.tab[[ty]] <<- read.table(paste(folder,"DMstructPred",type,".txt",sep=""),header=FALSE)
    PAR.tab[[ty]] <<- read.table(paste(folder,"PARPred",type,".txt",sep=""),header=FALSE)
    GAI.tab[[ty]] <<- read.table(paste(folder,"GAIPred",type,".txt",sep=""),header=FALSE)
    DegRate.tab[[ty]] <<-  read.table(paste(folder,"NphDegRatePred",type,".txt",sep=""),header=FALSE)
    SynthRate.tab[[ty]] <<-  read.table(paste(folder,"NphSynthRatePred",type,".txt",sep=""),header=FALSE)
    DMDegRate.tab[[ty]] <<- read.table(paste(folder,"DMremDegRatePred",type,".txt",sep=""),header=FALSE)
    DMSynthRate.tab[[ty]] <<- read.table(paste(folder,"DMremSynthRatePred",type,".txt",sep=""),header=FALSE)
    Photo.tab[[ty]] <<- read.table(paste(folder,"PhotoPred",type,".txt",sep=""),header=FALSE)
  }
Nmob.vect <<- read.table(paste(folder,"NmobPred.txt",sep=""),header=FALSE)
Ngrain.tab <<- read.table(paste(folder,"NgrainPred.txt",sep=""),header=FALSE)
MSgrain.tab <<- read.table(paste(folder,"MSgrainPred.txt",sep=""),header=FALSE)
ProdDem.vect <<- read.table(paste(folder,"ProdDemPred.txt",sep=""),header=FALSE)
Absorb.vect <<- read.table(paste(folder,"AbsorbPred.txt",sep=""),header=FALSE)
AreaGreenTot.vect <<- read.table(paste(folder,"AgreenTot.txt",sep=""),header=FALSE)
DMremroot.vect <<- read.table(paste(folder,"DMremRootPred.txt",sep=""),header=FALSE)
NremrootDegRate.vect <<- read.table(paste(folder,"NremDegRatePred.txt",sep=""),header=FALSE)
NremRoot.vect <<- read.table(paste(folder,"NremPredRoot.txt",sep=""),header=FALSE)
NstructRoot.vect <<- read.table(paste(folder,"NstructPredRoot.txt",sep=""),header=FALSE)
DMstructRoot.vect <<- read.table(paste(folder,"DMstructRootPred.txt",sep=""),header=FALSE)
RootSynth.vect <<- read.table(paste(folder,"RootSynthPred.txt",sep=""),header=FALSE)
}

# Simulations:

NAgreenDynamics <- function(){

  par(mfrow=c(2,5))
  type.vect <- c("Lamina","Sheath","Internode","Peduncle-Chaff")
  
  for (ty in 1:4){
    plot(Nph.tab[[ty]][,1],Nph.tab[[ty]][,2],ylim=c(0,0.008),type="l",col=1,main=type.vect[ty])
    if (ty!=4){
      lines(Nph.tab[[ty]][,1],Nph.tab[[ty]][,3],col=2)
      lines(Nph.tab[[ty]][,1],Nph.tab[[ty]][,4],col=3)
      lines(Nph.tab[[ty]][,1],Nph.tab[[ty]][,5],col=4)
      lines(Nstruct.tab[[ty]][,1],Nstruct.tab[[ty]][,2],lty=2,col=1)
      lines(Nstruct.tab[[ty]][,1],Nstruct.tab[[ty]][,3],lty=2,col=2)
      lines(Nstruct.tab[[ty]][,1],Nstruct.tab[[ty]][,4],lty=2,col=3)
      lines(Nstruct.tab[[ty]][,1],Nstruct.tab[[ty]][,5],lty=2,col=4)
    }else{
      lines(Nph.tab[[ty+1]][,1],Nph.tab[[ty+1]][,2],col="pink")
      lines(Nstruct.tab[[ty]][,1],Nstruct.tab[[ty]][,2],lty=2,col=1)
      lines(Nstruct.tab[[ty]][,1],Nstruct.tab[[ty]][,2],lty=2,col="pink")
    }
  }
  
  plot(Ngrain.tab[,1],Ngrain.tab[,3],type="l",col="brown",main="Grain")
  lines(Ngrain.tab[,1],Ngrain.tab[,2],col=1)
  lines(NstructRoot.vect[,1],NstructRoot.vect[,2],lty=2,col="brown")
  
  for (ty in 1:4){
    plot(Agreen.tab[[ty]][,1],Agreen.tab[[ty]][,2],ylim=c(0,0.004),type="l",col=1,main=type.vect[ty])
    if (ty!=4){
      lines(Agreen.tab[[ty]][,1],Agreen.tab[[ty]][,3],col=2)
      lines(Agreen.tab[[ty]][,1],Agreen.tab[[ty]][,4],col=3)
      lines(Agreen.tab[[ty]][,1],Agreen.tab[[ty]][,5],col=4)
    }else{
      lines(Agreen.tab[[ty+1]][,1],Agreen.tab[[ty+1]][,2],col="pink")
    }
    
  }

  plot(Absorb.vect[,1],Absorb.vect[,2],ylim=c(0,0.004),type="l",col=1,main="Absorption racinaire & Nmob")
  lines(Nmob.vect[,1],Nmob.vect[,2],col=3)
  
}


Absorb <- function(){

  par(mfrow=c(1,5))
  plot(Absorb.vect[,1],Absorb.vect[,2],ylim=c(0,0.002),type="l",col=1,main="Absorption")
  plot(Absorb.vect[,1],Absorb.vect[,3],type="l",col=3,ylim=c(0,1), main="N&Ceffect")
  lines(Absorb.vect[,1],Absorb.vect[,4],col=2)
  plot(Nmob.vect[,1],Nmob.vect[,3],col=3,type="l",main="ConcNmob")
plot(RootSynth.vect[,1],RootSynth.vect[,2],type="l",main="RootSynthRate")
plot(ProdDem.vect[,1],ProdDem.vect[,2],ylim=c(0,0.015),type="l",main="Production")

}


DMdynamics <- function(){

  par(mfrow=c(3,4))
  type.vect <- c("Lamina","Sheath","Internode","Peduncle-Chaff")
  
  for (ty in 1:4){
    plot(DMrem.tab[[ty]][,1],DMrem.tab[[ty]][,2],ylim=c(0,0.5),type="l",col=1,main=type.vect[ty])
    lines(DMstruct.tab[[ty]][,1],DMstruct.tab[[ty]][,2],lty=2,col=1)
    if (ty!=4){
      lines(DMrem.tab[[ty]][,1],DMrem.tab[[ty]][,3],col=2)
      lines(DMrem.tab[[ty]][,1],DMrem.tab[[ty]][,4],col=3)
      lines(DMrem.tab[[ty]][,1],DMrem.tab[[ty]][,5],col=4)
      lines(DMstruct.tab[[ty]][,1],DMstruct.tab[[ty]][,3],lty=2,col=2)
      lines(DMstruct.tab[[ty]][,1],DMstruct.tab[[ty]][,4],lty=2,col=3)
      lines(DMstruct.tab[[ty]][,1],DMstruct.tab[[ty]][,5],lty=2,col=4)
    }else{
      lines(DMrem.tab[[ty+1]][,1],DMrem.tab[[ty+1]][,2],col="pink")
      lines(DMstruct.tab[[ty+1]][,1],DMstruct.tab[[ty+1]][,2],lty=2,col="pink")
    }
  }

  for (ty in 1:4){
    plot(PAR.tab[[ty]][,1],PAR.tab[[ty]][,2],ylim=c(0,10000000),type="l",col=1,main=type.vect[ty])
    if (ty!=4){
      lines(PAR.tab[[ty]][,1],PAR.tab[[ty]][,3],col=2)
      lines(PAR.tab[[ty]][,1],PAR.tab[[ty]][,4],col=3)
      lines(PAR.tab[[ty]][,1],PAR.tab[[ty]][,5],col=4)
    }else{
      lines(PAR.tab[[ty+1]][,1],PAR.tab[[ty+1]][,2],col="pink")
    }
  }
  plot(MSgrain.tab[,1],MSgrain.tab[,2],type="l",main="DMgrain")
  plot(ProdDem.vect[,1],ProdDem.vect[,2],ylim=c(0,0.15),type="l",col=1,main="Production")
  plot(ProdDem.vect[,1],ProdDem.vect[,3],ylim=c(0,10),type="l",col=2,main="Demand")
  plot(ProdDem.vect[,1],ProdDem.vect[,4],ylim=c(0,0.05),type="l",col=1,main="QonD")
  
}


