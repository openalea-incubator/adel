library("rsm")

input.matrix <- function (folder) {

DOE <<- read.table(folder,header=FALSE)
}


rename <- function(){	
PBNEMA<-coded.data(DOE, Ng~V1, DMg~V2, DMr~V3, Nr~V4, DMc~V5, Nc~V6, DMin.3~V7, Nin.3~V8, DMin.2~V9, Nin.2~V10, DMin.1~V11, Nin.1~V12, DMin~V13, Nin~V14, DMln.3~V15, Nln.3~V16, DMln.2~V17, Nln.2~V18, DMln.1~V19, Nln.1~V20, DMln~V21, Nln~V22, DMp~V23, Np~V24, DMsn.3~V25, Nsn.3~V26, DMsn.2~V27, Nsn.2~V28, DMsn.1~V29, Nsn.1~V30, DMsn~V31, Nsn~V32, d1~V33, d2~V34, d3~V35, d4~V36, d5~V37, d6~V38, d7~V39, Ngrain~V40, Nlamina~V41, Nstem~V42, Nchaff~V43, Nmobile~V44, Nculm~V45, Photoarea~V46, Totalphoto~V47, DMculm~V48, DMlamina~V49, DMchaff~V50, DMstem~V51, DMgrain~V52)
PBNEMA
}


sensitivity.Ngrain <- function(){	
EFFECTS.Ngrain <- rsm(Ngrain ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = rename())
summary(EFFECTS.Ngrain)
}

sensitivity.Nlamina <- function(){	
EFFECTS.Nlamina<- rsm(Nlamina ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = rename())
summary(EFFECTS.Nlamina)
}

sensitivity.Nstem <- function(){	
EFFECTS.Nstem<- rsm(Nstem ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = rename())
summary(EFFECTS.Nstem)
}

sensitivity.Nchaff <- function(){	
EFFECTS.Nchaff<- rsm(Nchaff ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = rename())
summary(EFFECTS.Nchaff)
}

sensitivity.Nmobile <- function(){	
EFFECTS.Nmobile <- rsm(Nmobile ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = rename())
summary(EFFECTS.Nmobile)
}

sensitivity.Nculm <- function(){	
EFFECTS.Nculm <- rsm(Nculm ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = rename())
summary(EFFECTS.Nculm)
}

sensitivity.Photoarea <- function(){	
EFFECTS.Photoarea <- rsm(Photoarea ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = rename())
summary(EFFECTS.Photoarea)
}

sensitivity.Totalphoto <- function(){	
EFFECTS.Totalphoto <- rsm(Totalphoto ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = rename())
summary(EFFECTS.Totalphoto)
}

sensitivity.DMculm <- function(){	
EFFECTS.DMculm <- rsm(DMculm ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = rename())
summary(EFFECTS.DMculm)
}

sensitivity.DMlamina <- function(){	
EFFECTS.DMlamina <- rsm(DMlamina ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = rename())
summary(EFFECTS.DMlamina)
}

sensitivity.DMchaff <- function(){	
EFFECTS.DMchaff <- rsm(DMchaff ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = rename())
summary(EFFECTS.DMchaff)
}

sensitivity.DMstem <- function(){	
EFFECTS.DMstem <- rsm(DMstem ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = rename())
summary(EFFECTS.DMstem)
}

sensitivity.DMgrain <- function(){	
EFFECTS.DMgrain <- rsm(DMgrain ~ FO(Ng, DMg, DMr, Nr, DMc, Nc, DMin.3, Nin.3, DMin.2, Nin.2, DMin.1, Nin.1, DMin, Nin, DMln.3, Nln.3, DMln.2, Nln.2, DMln.1, Nln.1, DMln, Nln, DMp, Np, DMsn.3, Nsn.3, DMsn.2, Nsn.2, DMsn.1, Nsn.1, DMsn, Nsn), data = rename())
summary(EFFECTS.DMgrain)
}


effects.results <- function (folder1) {
folder1=paste(folder1,"/",sep="")

sink(paste(folder1,"EFFECTS.Ngrain.txt",sep=""))
print(sensitivity.Ngrain())
sink()

sink(paste(folder1,"EFFECTS.Nlamina.txt",sep=""))
print(sensitivity.Nlamina())
sink()

sink(paste(folder1,"EFFECTS.Nstem.txt",sep=""))
print(sensitivity.Nstem())
sink()

sink(paste(folder1,"EFFECTS.Nchaff.txt",sep=""))
print(sensitivity.Nchaff())
sink()

sink(paste(folder1,"EFFECTS.Nmobile.txt",sep=""))
print(sensitivity.Nmobile())
sink()

sink(paste(folder1,"EFFECTS.Nculm.txt",sep=""))
print(sensitivity.Nculm())
sink()

sink(paste(folder1,"EFFECTS.Photoarea.txt",sep=""))
print(sensitivity.Photoarea())
sink()

sink(paste(folder1,"EFFECTS.Totalphoto.txt",sep=""))
print(sensitivity.Totalphoto())
sink()

sink(paste(folder1,"EFFECTS.DMculm.txt",sep=""))
print(sensitivity.DMculm())
sink()

sink(paste(folder1,"EFFECTS.DMlamina.txt",sep=""))
print(sensitivity.DMlamina())
sink()

sink(paste(folder1,"EFFECTS.DMchaff.txt",sep=""))
print(sensitivity.DMchaff())
sink()

sink(paste(folder1,"EFFECTS.DMstem.txt",sep=""))
print(sensitivity.DMstem())
sink()

sink(paste(folder1,"EFFECTS.DMgrain.txt",sep=""))
print(sensitivity.DMgrain())
sink()

sink(paste(folder1,"EFFECTS.All.txt",sep=""))
print(sensitivity.Ngrain())
print(sensitivity.Nlamina())
print(sensitivity.Nstem())
print(sensitivity.Nchaff())
print(sensitivity.Nmobile())
print(sensitivity.Nculm())
print(sensitivity.Photoarea())
print(sensitivity.Totalphoto())
print(sensitivity.DMculm())
print(sensitivity.DMlamina())
print(sensitivity.DMchaff())
print(sensitivity.DMstem())
print(sensitivity.DMgrain())
sink()


}



