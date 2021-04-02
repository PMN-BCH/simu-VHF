#-----------------#
# Function for F2 #
#-----------------#

F2.func<-function(){
  #Table of simulations
  
  F2.simu<-data.frame(array(NA,c(K,5)))
  colnames(F2.simu)<-c("nBras.F2","pC.F2","pI.F2","pval.F2","nc.F2")
  
  #Simulations
  
  for (k in 1:K){
    
    #Intervention group
    I <- factor( rbinom(N.F2/2,1,pI), levels=c("0","1") )
    #Control group (two-arm trials)
    C <- factor( rbinom(N.F2/2,1,pC), levels=c("0","1") )
    
    #Data saving
    F2.simu$nBras.F2[k]<-N.F2/2
    F2.simu$pC.F2[k]<-prop.table(table(C))[2]
    F2.simu$pI.F2[k]<-prop.table(table(I))[2]
    tab.C <- table(C)
    tab.I <- table(I)
    mt <- matrix(c(tab.C,tab.I),ncol=2)
    

    E <- chisq.test(mt)$expected
    #Statistical test
    if (E[1,1]==0 | E[1,2]==0 | E[2,1]==0 | E[2,2]==0) { F2.simu$pval.F2[k] <- NA ; F2.simu$nc.F2[k] <- 1 
    } else 
      if (E[1,1]<=5 | E[1,2]<=5 | E[2,1]<=5 | E[2,2]<=5) { F2.simu$pval.F2[k] <- prop.test(mt,alternative="greater",conf.level=0.975)$p.value }
    if (E[1,1]>5 & E[1,2]>5 & E[2,1]>5 & E[2,2]>5)  F2.simu$pval.F2[k] <- prop.test(mt,alternative="greater",conf.level=0.975,correct=F)$p.value
  }
  write.table(F2.simu,paste(chemin.details,spec,"_",scenario,".csv",sep=""),sep=";",row.names=FALSE)
  
  res$N.F2[which(res$yx==scenario)]     <- N.F2
  res$pC.F2[which(res$yx==scenario)]    <- round(mean(F2.simu$pC.F2),3)
  res$pI.F2[which(res$yx==scenario)]    <- round(mean(F2.simu$pI.F2),3)
  res$p.F2[which(res$yx==scenario)]     <- round((sum(F2.simu$pval.F2<0.025,na.rm=T)/K),4)
  res$pMed.F2[which(res$yx==scenario)]  <- median(F2.simu$pval.F2,na.rm=T)
  res$p5.F2[which(res$yx==scenario)]    <- quantile(F2.simu$pval.F2,0.05,na.rm=T)
  res$p95.F2[which(res$yx==scenario)]   <- quantile(F2.simu$pval.F2,0.95,na.rm=T)
  res$nc.F2[which(res$yx==scenario)]   <- sum(F2.simu$nc.F2,na.rm=T)/K
  
  return(res)
}
