#
#                Setting and dressing of adel plants
#
#
#Utilities foor reconstructing a plant from parameters
#
#extract or replicate a desired number of plant from a canopy table
#
setCanopy <- function(canT, nplants=1, randomize = TRUE, seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)
  out <- vector("list",nplants)
  plantdb <- by(canT,list(canT$plant),function(x) x)
  if (randomize)
    plnb <- ceiling(runif(nplants) * length(plantdb))
  else
    plnb <- rep(seq(length(plantdb)),len=nplants)
  for (p in seq(out)) {
    pT <- plantdb[[plnb[p]]]
    pT$plant <- p
    if (randomize)
      pT$Laz <- (pT$Laz + 360 * runif(1)) %% 360
    out[[p]] <- pT
  }
  do.call("rbind",out)
}
  
#
predictDim <- function(dimT,index,nf) {
  nout <- seq(nf)/nf
  dim <- dimT[dimT$index == index,]
  out <- vector("list",ncol(dimT)-2)
  names(out) <- colnames(dim)[-match(c('index','nrel'),colnames(dim))]
  for (w in names(out))
    out[[w]] <- approx(dim$nrel,dim[,w],nout,rule=2)$y
  data.frame(do.call("cbind",out))
}
#
predictPhen <- function(phenT,index,nf,datesf1) {
  if (!"disp"%in%colnames(phenT))
    print("Missing input for leaf desapearance in phenT (see docAdel.txt")
  nout <- c(0,seq(nf))/nf
  phen <- phenT[phenT$index == index,]
  out <- vector("list",ncol(phenT) -2)
  names(out) <- colnames(phen)[-match(c('index','nrel'),colnames(phen))]
  for (i in 1:4) {
    w <- names(out)[i]
    out[[w]] <- approx(phen$nrel,phen[,w],nout,rule=2)$y + unlist(datesf1)[i] 
  }
  data.frame(cbind(n=c(0,seq(nf)),do.call("cbind",out)))
}
#
#peduncle elongation
#
predictPed <- function(pheno,phyto,index,nf,earT) {
  par <- earT[index,]
  rate <- 0
  if (par$em_ped - par$em_ear > 0)
    rate <- par$l_ear / (par$em_ped - par$em_ear)
  start <- pheno[pheno$n==nf,"col"] + par$em_ped - phyto[nf,"Gl"] / ifelse(rate == 0,1,rate)
  end <- start + par$l_ped / ifelse(rate == 0,1,rate)
  data.frame(startPed=start,endPed=end,senPed=pheno[pheno$n==nf,"col"] + par$end_gf)
}

  
#setAdel performs the dressing (geometry, tiller number ...) of plants from parameters and duplicate them for a given number of outputed plants
#
setAdel <- function(axeT,dimT,phenT,earT,ssisenT,geoLeaf,geoAxe,nplants=1,seed=NULL) {
  #verif inputs et completion default values
  if (!is.null(seed))
    set.seed(seed)
  if (is.null(earT)) {
    iear <- grep("earIndex",colnames(axeT),fixed=TRUE)
    if (length(iear) > 0)
      axeT <- axeT[,-iear]
    axeT <- cbind(axeT,earIndex = 1)
    earT <- data.frame(index = 1,em_ear = 100, em_ped = 200, end_gf = 1000, l_ped = 0, d_ped = 0, l_ear = 0, Sp_ear = 0, l_ear_awn = 0)
  }
  
  if (is.null(ssisenT))
    ssisenT <- data.frame(ndel=1:4,rate1=0.07,dssit1=c(0,1.2,2.5,3),dssit2=c(1.2,2.5,3.7,4))


  if (!"incB"%in%colnames(dimT)) 
    dimT <- cbind(dimT,incB = -999,dincB = 0)

  useAzim <- FALSE
  if (!"pAngle"%in%colnames(dimT)) {
    useAzim <- TRUE
    dimT <- cbind(dimT,pAngle = -999,dpAngle = 0)
  }
  
  out <- vector("list",nplants)
  plantdb <- by(axeT,list(axeT$plant),function(x) x)
  #sampling nplants in the database
  plnb <- ceiling(runif(nplants) * length(unique(axeT$plant)))
  for (p in seq(out)) {
    #axeTable from axeT and geoAxe or dimT if azim/azdev in colums
    pT <- plantdb[[plnb[p]]]
    axeTable <- data.frame(axe = pT$axe,
                           nf = pT$nf,
                           emf1 = pT$emf1,
                           end = pT$end,
                           disp = pT$disp,
                           azT = sapply(pT$axe,geoAxe$azT),
                           incT = sapply(pT$axe,geoAxe$incT),
                           dredT = sapply(pT$axe,geoAxe$dredT)
                           )
    #phytoT from dimT and geoleaf
    nomsdim <- c("Ll","Lw","Gl","Gd","El","Ed","pAngle","dpAngle","incB","dincB")
    #row = phytomer number + 3 (peduncle, ear, awns)
    nfM <- max(pT$nf) + 3
    phytoT <- array(NA,dim=c(nfM,length(nomsdim)+3,nrow(pT)),dimnames=list(seq(nfM),c(nomsdim,"Azim","Lindex","Lseed"),pT$axe))
    for (a in seq(nrow(pT))) {
      nf <- pT$nf[a]
      phytoT[seq(nf),nomsdim,a] <- unlist(predictDim(dimT,pT$dimIndex[a],nf)[,nomsdim])
      phytoT[seq(nf),"incB",a] <- phytoT[seq(nf),"incB",a] + (runif(nf) - .5) * phytoT[seq(nf),"dincB",a]
      if (useAzim) 
        phytoT[seq(nf),"Azim",a] <- sapply(seq(nf),function(n) geoLeaf$Azim(a,n,nf-n))
      else
        phytoT[seq(nf),"Azim",a] <- phytoT[seq(nf),"pAngle",a] + (runif(nf) - .5) * phytoT[seq(nf),"dpAngle",a]
      #phytoT[seq(nf),"Lindex",a] <- sapply(seq(nf),function(n) geoLeaf$Lindex(a,n,nf-n))
      phytoT[seq(nf),"Lseed",a] <- runif(nf)
      phytoT[(nf+1):(nf+3),,a] <- 0
      phytoT[(nf+1):(nf+2),"El",a] <- unlist(earT[pT$earIndex[a],c("l_ped","l_ear")])
      phytoT[nf+3,"El",a] <- unlist(earT[pT$earIndex[a],"l_ear_awn"] - earT[pT$earIndex[a],"l_ear"])
      phytoT[nf+1,"Ed",a] <- unlist(earT[pT$earIndex[a],"d_ped"])
      if (earT[pT$earIndex[a],"l_ear"] > 0)
        phytoT[(nf+2):(nf+3),"Ed",a] <- unlist(earT[pT$earIndex[a],"Sp_ear"]/earT[pT$earIndex[a],"l_ear"])
    }
    #tip, collar and ssi : table input for oppenapprox
    phenoT <- vector("list",nrow(pT))
    names(phenoT) <- pT$axe
    if (!"dispf1"%in%colnames(pT))
      print("Missing input for leaf desapearance in axeT (see docAdel.txt")
    for (a in seq(nrow(pT)))
      phenoT[[a]] <- predictPhen(phenT,pT$phenIndex[a],pT$nf[a],pT[a,c("emf1","ligf1","senf1","dispf1")])
    # date of start and end of peduncle elongation
    pedT <- vector("list",nrow(pT))
    names(pedT) <- pT$axe
    for (a in seq(nrow(pT)))
      pedT[[a]] <- predictPed(phenoT[[a]],phytoT[,,a],pT$earIndex[a],pT$nf[a],earT)
    
    out[[p]] <- list(axeT = axeTable,phytoT=phytoT,pheno=phenoT,pedT = pedT,ssisenT=ssisenT,geoLeaf=geoLeaf)
  }
  out
}

