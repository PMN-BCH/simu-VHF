#-----------------------------------------#
# Simulations for S2 - Changing with time #
#-----------------------------------------#

source("/.../script_S2_Evo_fct.R",encoding="native.enc")
chemin.details <- "/.../details_S2_EvoA/"

#Simulation specificities
#-----------------------------

#Sample size
spec <- "S2_EvoA"
N.S2 <- 378   #NSN of fix design

delta <- 20 #Interim analyses every 20 patients included

#Simulation number
K <- 93639


#Evolution
augmentation <- 0.03
borne <- 0.1


#Result importing
#-----------------
res <- read.table("/.../res_S2.csv",header=T,na.strings="NA",sep=";")

#Automatic process
#------------------

for (m in 1:nrow(res)){
  print(paste(spec,": Augmentation de",augmentation,", bornée à",borne,"> Scenario",m))
  # Scenario identification
  scenario <- res$yx[m]
  set.seed( res$graine[which(res$yx==scenario)] )
  
  #Parameters defined by the scenario
  pC.init <- res$pC.yx[m]
  pI.init <- res$pI.yx[m]
  
  #Execution 
  res<-S2.func()
}
#Output
#-------
write.table(res,paste("/.../res_",spec,".csv",sep=""),sep=";",row.names=F)