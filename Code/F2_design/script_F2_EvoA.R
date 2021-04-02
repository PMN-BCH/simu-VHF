#-----------------------------------------#
# Simulations for F2 - Changing with time #
#-----------------------------------------#

source("/.../script_F2_Evo_fct.R",encoding="native.enc")
chemin.details <- "/.../details_F2_EvoA/"

#Simulation specifities
#------------------------------
spec <- "F2_EvoA"

#Sample size
N.F2 <- 248  #NSN of fix design

#Simulation number
K <- 93639

#Evolution
augmentation <- 0.03
borne <- 0.1


#Result importing
#------------------
res <- read.table("/.../res_F2.csv",header=T,na.strings="NA",sep=";")
res$nc.F2 <- NA


#Automatic process
#-------------------

for (m in 1:nrow(res)){
  print(paste(spec,": Augmentation de",augmentation,", bornée à",borne,"> Scenario",m))
  # Scenario identification
  scenario <- res$yx[m]
  set.seed( res$graine[which(res$yx==scenario)] )#  the same for all scenarios/designs
  
  #Parameters defined by the scenario
  pC.init <- res$pC.yx[m]
  pI.init <- res$pI.yx[m]
  
  #Execution
  res<-F2.func()
}
#Output
#-------
write.table(res,paste("/.../res_",spec,".csv",sep=""),sep=";",row.names=F)