#--------------------#
# Simulations for S1 #
#--------------------#

source("/.../script_S1_fct.R",encoding="native.enc")
chemin.details <- "/.../details_S1_NX/"


#Simulation specificities
#--------------------------

#Sample size
spec <- "S1_NX"
N.S1 <- X

delta <- 20 #Interim analyses every delta=20 patients 

#Simulation number
K <- 93639

#Result importing
#-----------------
res <- read.table("/.../res_S1.csv",header=T,na.strings="NA",sep=";")

#Automatic process
#------------------

for (m in 1:nrow(res)){
  print(paste(spec,"> Scenario",m))
  # Scenario identification 
  scenario <- res$yx[m]
  set.seed( res$graine[which(res$yx==scenario)] ) #  the same for all scenarios/designs
  
  #Parameters defined by the scenario
  pC <- 0.5
  pI <- res$pI.yx[m]
  
  #Execution 
  res<-S1.func()
}
#Output
#------
write.table(res,paste("/.../res_",spec,".csv",sep=""),sep=";",row.names=FALSE)