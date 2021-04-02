#-----------------------------------------#
# Simulations for F1 - Changing with time #
#-----------------------------------------#


source("/.../script_F1_Evo_fct.R",encoding="native.enc")
chemin.details <- "/.../details_F1_EvoA/"

#Simulation specificities
#--------------------------
spec <- "F1_EvoA"

#Sample size
N.F1 <- 60

#Simulation number
K <- 93639

#Evolution
augmentation <- 0.03
borne <- 0.1

#Result importing
#------------------
res <- read.table("/.../res_F1.csv",header=T,na.strings="NA",sep=";")


# Automatic process
#----------------------

for (m in 1:nrow(res)){
  print(paste(spec,": Increase of",augmentation,", limited to",borne,"> Scenario",m))
  # Scenario identification
  scenario <- res$yx[m]
  set.seed( res$graine[which(res$yx==scenario)] )# the same for all scenarios/designs
  
  #Parameters defined by the scenario
  pC <- 0.5
  pI.init <- res$pI.yx[m]
  
  #Execution
  res<-F1.func()
}
#Output
#------
write.table(res,paste("/.../res_",spec,".csv",sep=""),sep=";",row.names=F)
