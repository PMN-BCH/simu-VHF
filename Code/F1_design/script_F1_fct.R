#-----------------#
# Function for F1 #
#-----------------#

F1.func<-function(){
  #Table of simulations
  F1.simu<-data.frame(array(NA,c(K,3)))
  colnames(F1.simu)<-c("nBras.F1","pI.F1","pval.F1")
  
  #Simulations
  for (k in 1:K){
    
    #Intervention gorup
    I <- factor( rbinom(N.F1,1,pI), levels=c("0","1") )
    
    #Data saving
    F1.simu$nBras.F1[k]<-length(I)
    F1.simu$pI.F1[k]<- prop.table(table(I))[2]
    
    E <- chisq.test(table(I),p=c(1-pC,pC))$expected
    #Statistical test
    if (E[1]<=5 | E[2]<=5) F1.simu$pval.F1[k] <- prop.test(table(I)[2],n=N.F1,p=pC,alternative="greater",conf.level=0.975)$p.value 
    if (E[1]>5 & E[2]>5)  F1.simu$pval.F1[k] <- prop.test(table(I)[2],n=N.F1,p=pC,alternative="greater",conf.level=0.975,correct=F)$p.value
    
  }
  write.table(F1.simu,paste(chemin.details,spec,"_",scenario,".csv",sep=""),sep=";",row.names=FALSE)
  
  res$N.F1[which(res$yx==scenario)] <- N.F1
  res$pI.F1[which(res$yx==scenario)] <- round(mean(F1.simu$pI.F1),3)
  res$p.F1[which(res$yx==scenario)] <- round((sum(F1.simu$pval.F1<0.025)/K),4)
  res$pMed.F1[which(res$yx==scenario)]  <- median(F1.simu$pval.F1)
  res$p5.F1[which(res$yx==scenario)]    <- quantile(F1.simu$pval.F1,0.05)
  res$p95.F1[which(res$yx==scenario)]   <- quantile(F1.simu$pval.F1,0.95)
  
  return(res)
}