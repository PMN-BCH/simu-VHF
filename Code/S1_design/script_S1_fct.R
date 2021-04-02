#-----------------#
# Function for S1 #
#-----------------#

#Parameters for group-sequential design
NF1 <- 60 #NSN of fix design

alpha <- 0.025
beta <- 0.1
  
z1<- -qnorm(alpha,0,1)
z2<- -qnorm(beta,0,1)
  
theta.nc<-log( (0.7*(1-0.5)) / ((0.5*(1-0.7))) )
theta<-log( (0.7*(1-0.5)) / ((0.5*(1-0.7))) ) * ((2*z1)/(z1+z2))

a<-(2/theta)*log(1/(2*alpha))
c<-0.25*theta


S1.func<-function()
  {

  #Table of simulations
  S1.simu<-data.frame(array(NA,c(K,7)))
  colnames(S1.simu)<-c("NJ.S1","pIJ.S1","ZJ.S1","BJinf.S1","BJsup.S1","ccl.S1","sig_level.S1")
  
  for (k in 1:K){
    
    #Intervention group
    I <- rbinom(N.S1,1,pI)
    
    N<-0
    S<-0
    V<-0
    
    zsup<-0
    zinf<-0
    Z<-0
    indic<-1 #Indicator, if = 0 : trial end
    nc <- NULL #Inconclusive trial (yes/no)
    
    #Analysis j
    while (indic==1) 
      {
        #Sampling without replacement in I for delta patients included in j
        if ((length(I)-delta)<0)
          {
          index <- sample(1:length(I), length(I), replace=F)
        } else {
          index <- sample(1:length(I), delta, replace=F)
        }
        Ij <- I[index]
        I <- I[-index]
        deltaS <- sum(Ij)
      
        N<-N+length(index)
        S<-S+deltaS
        
        #Stopping rules
        deltaV<-delta*pC*(1-pC)
        V<-V+deltaV
        zsup<- a + c*V - 0.583*sqrt(deltaV)
        zinf<- -a + 3*c*V+0.583*sqrt(deltaV)
        
        #Test stat
        Z<-S-N*pC
        
        #Conclusion (ccl)
        ccl <- ifelse(Z>=zsup,1,0)
        
        #Indicator for While Loop
        indic <- ifelse((length(I)>0) & Z>=zinf & Z<zsup,1,0)
        
        #Inconclusive yes/no
        nc <- ifelse(ccl==0 & Z>zinf,1,0)
    }
    
    S1.simu$NJ.S1[k] <- N
    S1.simu$pIJ.S1[k] <- S/N
    S1.simu$ZJ.S1[k] <- Z
    S1.simu$BJinf.S1[k] <- zinf
    S1.simu$BJsup.S1[k] <- zsup
    S1.simu$ccl.S1[k] <- ccl
    S1.simu$nc.S1[k] <- nc
    S1.simu$sig_level.S1[k] <- ifelse(nc==1,(1-pnorm(abs(Z)/sqrt(V))),NA)
  }
  write.table(S1.simu,paste(chemin.details,spec,"_",scenario,".csv",sep=""),sep=";",row.names=F)
  
  res$NJMed.S1[which(res$yx==scenario)] <- median(S1.simu$NJ.S1) 
  res$NJ5.S1[which(res$yx==scenario)]   <- quantile(S1.simu$NJ.S1,0.05)
  res$NJ95.S1[which(res$yx==scenario)]  <- quantile(S1.simu$NJ.S1,0.95) 
  res$pct.supNF1[which(res$yx==scenario)]    <- round((sum(S1.simu$NJ.S1>NF1)/K),4)
  
  res$sig_levelMed.S1[which(res$yx==scenario)] <- median(S1.simu$sig_level.S1,na.rm=TRUE) 
  res$sig_level5.S1[which(res$yx==scenario)]   <- quantile(S1.simu$sig_level.S1,0.05,na.rm=TRUE)
  res$sig_level95.S1[which(res$yx==scenario)]  <- quantile(S1.simu$sig_level.S1,0.95,na.rm=TRUE) 
  res$sig_levelNA.S1[which(res$yx==scenario)]  <- round(sum(is.na(S1.simu$sig_level.S1))/K,5)
  res$sig_level_alpha.S1[which(res$yx==scenario)] <- round((sum(S1.simu$sig_level.S1<0.025,na.rm=TRUE)/K),4)
  
  res$pIJ.S1[which(res$yx==scenario)]   <- round(mean(S1.simu$pIJ.S1),3)
  res$p.S1[which(res$yx==scenario)]     <- round((sum(S1.simu$ccl.S1)/K),4)
  res$nc.S1[which(res$yx==scenario)]    <- round((sum(S1.simu$nc.S1)/K),4)
  res$futil.S1[which(res$yx==scenario)] <- round(nrow(subset(S1.simu,S1.simu$nc.S1==0 & S1.simu$ccl.S1==0))/K,4)
  
  res$NJMed.S1_Eff[which(res$yx==scenario)] <- median(S1.simu$NJ.S1[S1.simu$ccl.S1==1]) 
  res$NJ5.S1_Eff[which(res$yx==scenario)]   <- quantile(S1.simu$NJ.S1[S1.simu$ccl.S1==1],0.05)
  res$NJ95.S1_Eff[which(res$yx==scenario)]  <- quantile(S1.simu$NJ.S1[S1.simu$ccl.S1==1],0.95)
  res$NJMed.S1_NC[which(res$yx==scenario)]  <- median(S1.simu$NJ.S1[S1.simu$nc.S1==1]) 
  res$NJ5.S1_NC[which(res$yx==scenario)]    <- quantile(S1.simu$NJ.S1[S1.simu$nc.S1==1],0.05)
  res$NJ95.S1_NC[which(res$yx==scenario)]   <- quantile(S1.simu$NJ.S1[S1.simu$nc.S1==1],0.95)
  res$NJMed.S1_Futil[which(res$yx==scenario)]  <- median(S1.simu$NJ.S1[S1.simu$nc.S1==0 & S1.simu$ccl.S1==0]) 
  res$NJ5.S1_Futil[which(res$yx==scenario)]    <- quantile(S1.simu$NJ.S1[S1.simu$nc.S1==0 & S1.simu$ccl.S1==0],0.05)
  res$NJ95.S1_Futil[which(res$yx==scenario)]   <- quantile(S1.simu$NJ.S1[S1.simu$nc.S1==0 & S1.simu$ccl.S1==0],0.95)
  
  return(res)
}
