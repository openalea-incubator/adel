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
devT <- devTcsv(paste(chem,"axeT",manip,".csv",sep=""),
                paste(chem,"dimT",manip,".csv",sep=""),
                paste(chem,"phenT",manip,".csv",sep=""),
                paste(chem,"earT",manip,".csv",sep=""),
                paste(chem,"ssi2sen.csv",sep=""))
#
#avatar
#devT <- devTcsv(paste(chem,"axeTSoissonsNormal_10plants.csv",sep=""),
#                paste(chem,"dimTSoissons.csv",sep=""),
#                paste(chem,"phenTNormalPhyllochron.csv",sep=""),
#                paste(chem,"earT",manip,".csv",sep=""),
#                paste(chem,"ssi2sen.csv",sep=""))
#Specify geoLeaf & geoAxe
geoLeaf <- genGeoLeaf()
geoAxe <- genGeoAxe()
#
attach(paste(chem,"So99.RData",sep=""))
xydb <- get("So99")
detach()
attach(paste(chem,"SRSo.RData",sep=""))
srdb <- get("SRSo")
detach()
#
#generate a list of plant to simulate from parameters
pl <- setAdel(devT$axeT,devT$dimT,devT$phenT,devT$earT,devT$ssisenT,geoLeaf,geoAxe,nplants=1,xy_db=xydb,sr_db=srdb)
#run the model as a whole from plant list to AleaChn
canopy <- runAdel(1200,pl)[[1]]
#
Scanopy <- canL2canS(canopy,srdb)
#
chn <- genString(canopy)
#
#
#Run step by step the model
#
#
#Predict kinetic of extension
#
kinExt <- lapply(pl,function(plant) kinL(2000,plant))
#kinExt[[1]][[1]][1,,]
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
#
# simulates how adel param is relaed to HS in case collar app = endleaf
#
phyl <- 110
endLeaf <- 1.6
nF <- 11
#
time <- seq(0,1500) # time since emergence leaf 1
leaves <- matrix(0,ncol=nF,nrow=length(time))
for (i in seq(nF)) {
  xem <- phyl * (i - 1)
  xlig <- xem + endLeaf * phyl
  l <- leaves[,i]
  l[time >= xem] <- (time[time >= xem] - xem) / phyl / endLeaf
  l[time >= xlig] <- 1
  leaves[,i] <- l
}
#
hs <- apply(leaves,1,sum)
plot(time / phyl,hs,xlab='phyllochronic time since emergence 1(adel)/HS_0', ylab='Leaf #',type='l',cex.lab=1.5)
abline(1,1,lwd=2,col=4)
abline(1-endLeaf,1,lwd=2,col=3)
abline(-.5/1.6,1,lwd=2,col=2)
legend(8,4,c('HS','Tips','Collars','HS_lin (Tips - 1.31)'),col=c(1,4,3,2),lty=1,lwd=c(1,2,2,2))
