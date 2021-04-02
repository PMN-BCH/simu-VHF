#--------------------#
# Simulations for F1 #
#--------------------#

source("/.../script_F1_fct.R",encoding="native.enc")
chemin.details <- "/.../details_F1_NX/"

#Simulation specificities
#--------------------------

#Sample size
spec <- "F1_NX"
N.F1 <- X

#Simulation number
K <- 93639

#Result importing
#-----------------
res <- read.table("/home/pauline/M2_PMN/base/source/res_F1.csv",header=T,na.strings="NA",sep=";")


#Automatic processs
#----------------------

for (m in 1:nrow(res)){
  print(paste(spec,"> Scenario",m))
  # Scenario identification 
  scenario <- res$yx[m]
  set.seed( res$graine[which(res$yx==scenario)] ) #  the same for all scenarios/designs
  
  #Parameters defined by the scenario
  pC <- 0.5
  pI <- res$pI.yx[m]
  
  #Execution
  res<-F1.func()
}
#Output
#------
write.table(res,paste("/.../res_",spec,".csv",sep=""),sep=";",row.names=FALSE)
