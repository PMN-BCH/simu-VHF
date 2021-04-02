#--------------------#
# Simulations for S2 #
#--------------------#

source("/.../script_S2_fct.R",encoding="native.enc")
chemin.details <- "/.../details_S2_NX/"

#Simulation specifities
#-----------------------

#Sample size
spec <- "S2_NX"
N.S2 <- X

delta <- 20 #Interim analyses every delta=20 patients 

#Simulation number
K <- 93639


#Result importing
#-----------------
res <- read.table("/home/pauline/M2_PMN/base/source/res_S2.csv",header=T,na.strings="NA",sep=";")


#automatic process
#-------------------

for (m in 1:nrow(res)){
  print(paste(spec,"> Scenario",m))
  # Scenario identification
  scenario <- res$yx[m]
  set.seed( res$graine[which(res$yx==scenario)] )#  the same for all scenarios/designs
  
  #Parameters defined by the scenario
  pC <- res$pC.yx[m]
  pI <- res$pI.yx[m]
  
  #Execution 
  res<-S2.func()
}
#Output
#-------
write.table(res,paste("/home/pauline/M2_PMN/base/res_",spec,".csv",sep=""),sep=";",row.names=FALSE)