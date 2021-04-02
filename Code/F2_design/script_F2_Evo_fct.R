#--------------------------------------#
# Function for F2 - Changing with time #
#--------------------------------------#

F2.func<-function(){

  #Table of simulations
  F2.simu<-data.frame(array(NA,c(K,10)))
  colnames(F2.simu)<-c("pI1.S1","pI2.S1","pI3.S1","pI4.S1","pI5.S1","pI6.S1","pI7.S1","nBras.F2","pI.F2","pval.F2")
  
  #Simulations
  
  for (k in 1:K){
    
    j<-1
    #Intervention & control groups
    pI <- pI.init 
	pC <- pC.init
    I <- rbinom(10,1,pI.init)
    C <- rbinom(10,1,pC.init)
    F2.simu[k,j] <- sum(I)/length(I) 
    sum(I)/length(I) ; sum(C)/
	
    while ( (length(I)+10)<=N.F2/2 ) {
      if((pI+augmentation)<1 & (pI<pI.init+borne) )
        {
        pI <- pI.init+augmentation*j
		pC <- pC.init+augmentation*j
        }
       I.new <- rbinom(10,1,pI)
       C.new <- rbinom(10,1,pC)
       j <- j+1
       I <- c(I, I.new ) 
       C <- c(C, C.new ) 
	   F2.simu[k,j] <- sum(I)/length(I)
	   sum(I)/length(I)
	   sum(C)/length(C)
    }
    I <- factor(I, levels=c("0","1"))
    C <- factor(C, levels=c("0","1"))
	
   
   #Data saving
    F2.simu$nBras.F2[k]<-N.F2/2
    F2.simu$pC.F2[k]<-prop.table(table(C))[2]
    F2.simu$pI.F2[k]<-prop.table(table(I))[2]
    tab.C <- table(C)
    tab.I <- table(I)
    mt <- matrix(c(tab.C,tab.I),ncol=2)
    

    E <- chisq.test(mt)$expected
    #statistical test
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
  res$nc.F2[which(res$yx==scenario)]    <- sum(F2.simu$nc.F2,na.rm=T)/K
  
  return(res)
}
