#-----------------------------------------#
# Simulations for S1 - Changing with time #
#-----------------------------------------#

source(".../script_S1_Evo_fct.R",encoding="native.enc")
chemin.details <- "/.../details_S1_EvoA/"


#Simulation specifities
#------------------------

#Sample size
spec <- "S1_EvoA"
N.S1 <- 91   #NSN of fix design

delta <- 20 #Interim analyses every 20 patients included

NF1 <- 60

#Simulation number
K <- 93639

#Evolution
augmentation <- 0.03
borne <- 0.1

#Result importing
#------------------
res <- read.table("/.../res_S1.csv",header=T,na.strings="NA",sep=";")


#Automatic process
#------------------

for (m in 1:nrow(res)){
  print(paste(spec,": Increase of",augmentation,", limited to",borne,"> Scenario",m))
  # Scenario identification
  scenario <- res$yx[m]
  set.seed( res$graine[which(res$yx==scenario)] )#  the same for all scenarios/designs
  
  #Parametes defined by the scenario
  pC <- 0.5
  pI.init <-res$pI.yx[m]
  
  #Execution
  res<-S1.func()
}
#Output
#-------
write.table(res,paste("/.../res_",spec,".csv",sep=""),sep=";",row.names=F)
