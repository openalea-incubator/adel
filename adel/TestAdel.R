#
#
#           Script for using/debugging Adel.R from within R
#
#
# Will create .RData in the directory. Don't commit it !
#
#
#
#Load R files for Adel.R
#
AdelRfiles <- c("Adel.R","setAdel.R","UseAdel.R","genString.R")
sapply(AdelRfiles,source)
#
#run the model with csv files
#
#read parameters from csv
manip <- c("Ca0N")
chem <- c("./data/")
pars <- devTcsv(paste(chem,"axeT",manip,".csv",sep=""),
                paste(chem,"dimT",manip,".csv",sep=""),
                paste(chem,"phenT",manip,".csv",sep=""),
                paste(chem,"earT",manip,".csv",sep=""),
                paste(chem,"ssi2sen.csv",sep=""))
#
#avatar
pars <- devTcsv(paste(chem,"axeTSoissonsNormal_10plants.csv",sep=""),
                paste(chem,"dimTSoissons.csv",sep=""),
                paste(chem,"phenTNormalPhyllochron.csv",sep=""),
                paste(chem,"earT",manip,".csv",sep=""),
                paste(chem,"ssi2sen.csv",sep=""))
#Specify geoLeaf & geoAxe
geoLeaf <- genGeoLeaf()
geoAxe <- genGeoAxe()
#generate a list of plant to simulate from parameters
pl <- setAdel(pars$axeT,pars$dimT,pars$phenT,pars$earT,pars$ssisenT,geoLeaf,geoAxe,nplants=1)
#run the model as a whole from plant list to AleaChn
chn <- runAdel(1200,pl)[[1]]
#
#
#Run step by step the model
#
#
#Predict kinetic of extension
#
kinExt <- lapply(pl,function(plant) kinL(1000,plant))
#
# idem for visible parts
kinExtv <- lapply(kinExt,kinLvis)
#
#Dataframe representing the canopy
#
canopy <- getdesc(kinExtv,pl)
#
#string
#
chn <- genString(canopy)
#
#
#Test setAdelArv
#
load("./data/Arvalis.RData")
source("./data/Rfun.R")
source("ArvalisToAdel.R")
