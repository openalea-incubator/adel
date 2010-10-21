library("rsm")

ReadOut <- function (folder) {
folder=paste(folder,"/",sep="")

DOE <<- read.table(paste(folder,"MATRIX5.txt",sep=""),header=FALSE)
}

#PBNEMA <- coded.data(DOE, Ng~V1, DMg~V2, DMr~V3, Nr~V4, DMc~V5, Nc~V6, DMin.3~V7, Nin.3~V8, DMin.2~V9, Nin.2~V10, DMin.1~V11, Nin.1~V12, DMin~V13, Nin~V14, DMln.3~V15, Nln.3~V16, DMln.2~V17, Nln.2~V18, DMln.1~V19, Nln.1~V20, DMln~V21, Nln~V22, DMp~V23, Np~V24, DMsn.3~V25, Nsn.3~V26, DMsn.2~V27, Nsn.2~V28, DMsn.1~V29, Nsn.1~V30, DMsn~V31, Nsn~V32, d1~V33, d2~V34, d3~V35, d4~V36, d5~V37, d6~V38, d7~V39, Ngrain~V40, Nlamina~V41, Nstem~V42, Nchaff~V43, Nmobile~V44, Nculm~V45, Photoarea~V46, Totalphoto~V47, DMculm~V48, DMlamina~V49, DMchaff~V50, DMstem~V51, DMgrain~V52)


#EFFECTSNgrain <- rsm(Ngrain ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = PB.NEMA)
#sink("EFFECTS.Ngrain.txt")
#summary(EFFECTS.Ngrain)
#sink()


#EFFECTSNgrain <- rsm(V40 ~ FO(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, V21, V22, V23, V24, V25, V26, V27, V28, V29, V30, V31, V32), data = DOE)


cod <- function(){	
PBNEMA<-coded.data(DOE, Ng~V1, DMg~V2, DMr~V3, Nr~V4, DMc~V5, Nc~V6, DMin.3~V7, Nin.3~V8, DMin.2~V9, Nin.2~V10, DMin.1~V11, Nin.1~V12, DMin~V13, Nin~V14, DMln.3~V15, Nln.3~V16, DMln.2~V17, Nln.2~V18, DMln.1~V19, Nln.1~V20, DMln~V21, Nln~V22, DMp~V23, Np~V24, DMsn.3~V25, Nsn.3~V26, DMsn.2~V27, Nsn.2~V28, DMsn.1~V29, Nsn.1~V30, DMsn~V31, Nsn~V32, d1~V33, d2~V34, d3~V35, d4~V36, d5~V37, d6~V38, d7~V39, Ngrain~V40, Nlamina~V41, Nstem~V42, Nchaff~V43, Nmobile~V44, Nculm~V45, Photoarea~V46, Totalphoto~V47, DMculm~V48, DMlamina~V49, DMchaff~V50, DMstem~V51, DMgrain~V52)
PBNEMA
}


resul <- function(){	
EFFECTSNgrain <- rsm(Ngrain ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = cod())
summary(EFFECTSNgrain)
}


#res <- function(){	
#sink("C:/DOETEST/OUTPUT5/EFFECTS.Ngrain.txt")
#print(resul())
#sink()
#}

res <- function (folder) {
folder=paste(folder,"/",sep="")
sink(paste(folder,"EFFECTS.Ngrain.txt",sep=""))
print(resul())
sink()
}



##tes <- function(){	
##EFFECTSNgrain <- rsm(V40 ~ FO(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, V21, V22, V23, V24, V25, V26, V27, V28, V29, V30, V31, V32), data = DOE)
##summary(EFFECTSNgrain)
##}

##res <- function(){	
##sink("C:/DOETEST/OUTPUT5/EFFECTS.Ngrain.txt")
##print(tes())
##sink()
##}


