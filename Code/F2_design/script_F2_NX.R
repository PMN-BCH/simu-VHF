#--------------------#
# Simulations for F2 #
#--------------------#

source("/.../script_F2_fct.R",encoding="native.enc")
chemin.details <- "/.../details_F2_NX/"

#Simulation specificities
#--------------------------

#Sample size
spec <- "F2_NX"
N.F2 <- X 

#Simulation sumber
K <- 93639

#Importation de la table de resultat
#------------------------------------
res <- read.table("/home/pauline/M2_PMN/base/source/res_F2.csv",header=T,na.strings="NA",sep=";")
res$nc.F2 <- NA


#Automatic process
#--------------------

for (m in 1:nrow(res)){
  print(paste(spec,"> Scenario",m))
  # Identification du scenario et graine
  scenario <- res$yx[m]
  set.seed( res$graine[which(res$yx==scenario)] ) #  the same for all scenarios/designs
  
  #Parameters defined by the scenario
  pC <- res$pC.yx[m]
  pI <- res$pI.yx[m]
  
  #Execution 
  res<-F2.func()
}
#Output
#------
write.table(res,paste("/home/pauline/M2_PMN/base/res_",spec,".csv",sep=""),sep=";",row.names=FALSE)