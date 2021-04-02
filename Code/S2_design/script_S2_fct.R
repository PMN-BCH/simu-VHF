#-----------------#
# Function for S2 #
#-----------------#

#Parameters for group-sequential design
NF2 <- 248  #NSN of fix design

alpha <- 0.025
beta <- 0.1

z1<- -qnorm(alpha,0,1)
z2<- -qnorm(beta,0,1)

theta.nc<-log( (0.7*(1-0.5)) / ((0.5*(1-0.7))) )
theta<-log( (0.7*(1-0.5)) / ((0.5*(1-0.7))) ) * ((2*z1)/(z1+z2))

a<-(2/theta)*log(1/(2*alpha))
c<-0.25*theta


S2.func<-function(){
  
  #Table of simulations
  S2.simu<-data.frame(array(NA,c(K,8)))
  colnames(S2.simu)<-c("NJ.S2","pCJ.S2","pIJ.S2","ZJ.S2","BJinf.S2","BJsup.S2","ccl.S2","sig_level.S2")
  
  for (k in 1:K){
    
    #Intervention group
    I <- rbinom(N.S2/2,1,pI)
    #Control group
    C <- rbinom(N.S2/2,1,pC)

    S<-0
    S.I<-0
    S.C <-0
    
    N<-0
    N.I<-0
    N.C<-0
    
    zsup<-0
    zinf<-0
    Z<-0
    V<-Vant<-0
    indic<-1 #Indicator = 0 : trial end
    nc <- NULL #Inconclusive trial (yes/no)
    
    #Analysis j
    while (indic==1) {
      
        #Sampling without replacement in I for delta patients included in j
        if ((length(I)-delta/2)<0){
          index.I <- sample(1:length(I), (length(I)), replace=F)
          index.C <- sample(1:length(C), (length(C)), replace=F)
        } else {
          index.I <- sample(1:length(I), delta/2, replace=F)
          index.C <- sample(1:length(C), delta/2, replace=F)
        }
        Ij <- I[index.I]
        I <- I[-index.I]
        deltaS.I <- sum(Ij)
        
        Cj <- C[index.C]
        C <- C[-index.C]
        deltaS.C <- sum(Cj)
        
        S.I<-S.I+deltaS.I
        S.C<-S.C+deltaS.C
        N.I<-N.I+length(index.I)
        N.C<-N.C+length(index.C)
        N<-N+length(index.I)+length(index.C)
        S<-S+deltaS.I+deltaS.C
        E<-N-S 
        
        #Test stat
        Z<- ((N.C*S.I-N.I*S.C)/N)
        V <- ( (N.I*N.C*S*E)/(N**3))
        deltaV <- V-Vant
        
        #Stopping rules
        zsup<- a + c*V - 0.583*sqrt(deltaV)
        zinf<- -a + 3*c*V + 0.583*sqrt(deltaV)
        
        #Conclusion (ccl)
        ccl <- ifelse(Z>=zsup,1,0)
        
        #Indicator for While Loop
        indic <- ifelse((length(I)>0) & Z>=zinf & Z<zsup,1,0)
        
        #Inconclusive (yes/no)
        nc <- ifelse(ccl==0 & Z>zinf,1,0)
        
        # V anterior
        Vant <- V
    }
    
    S2.simu$NJ.S2[k] <- N
    S2.simu$pCJ.S2[k] <- S.C/N.C
    S2.simu$pIJ.S2[k] <- S.I/N.I
    S2.simu$ZJ.S2[k] <- Z
    S2.simu$BJinf.S2[k] <- zinf
    S2.simu$BJsup.S2[k] <- zsup
    S2.simu$ccl.S2[k] <- ccl
    S2.simu$nc.S2[k] <- nc
    S2.simu$sig_level.S2[k] <- ifelse(nc==1,2*(1-pnorm(abs(Z)/sqrt(V))),NA)
  }
  write.table(S2.simu,paste(chemin.details,spec,"_",scenario,".csv",sep=""),sep=";",row.names=F)
 
  res$NJMed.S2[which(res$yx==scenario)] <- median(S2.simu$NJ.S2) 
  res$NJ5.S2[which(res$yx==scenario)]   <- quantile(S2.simu$NJ.S2,0.05)
  res$NJ95.S2[which(res$yx==scenario)]  <- quantile(S2.simu$NJ.S2,0.95)
  res$pct.supNF2[which(res$yx==scenario)]    <- round((sum(S2.simu$NJ.S2>NF2)/K),4)
  
  res$sig_levelMed.S2[which(res$yx==scenario)] <- median(S2.simu$sig_level.S2,na.rm=TRUE) 
  res$sig_level5.S2[which(res$yx==scenario)]   <- quantile(S2.simu$sig_level.S2,0.05,na.rm=TRUE)
  res$sig_level95.S2[which(res$yx==scenario)]  <- quantile(S2.simu$sig_level.S2,0.95,na.rm=TRUE) 
  res$sig_levelNA.S2[which(res$yx==scenario)]  <- round(sum(is.na(S2.simu$sig_level.S2))/K,5)
  res$sig_level_alpha.S2[which(res$yx==scenario)] <- round((sum(S2.simu$sig_level.S2<0.025,na.rm=TRUE)/K),4)
  
  res$pCJ.S2[which(res$yx==scenario)]   <- round(mean(S2.simu$pCJ.S2),3)
  res$pIJ.S2[which(res$yx==scenario)]   <- round(mean(S2.simu$pIJ.S2),3)
  res$p.S2[which(res$yx==scenario)]     <- round((sum(S2.simu$ccl.S2)/K),4)
  res$nc.S2[which(res$yx==scenario)]    <- round((sum(S2.simu$nc.S2)/K),4)
  res$futil.S2[which(res$yx==scenario)] <- round(nrow(subset(S2.simu,S2.simu$nc.S2==0 & S2.simu$ccl.S2==0))/K,4)
  
  res$NJMed.S2_Eff[which(res$yx==scenario)] <- median(S2.simu$NJ.S2[S2.simu$ccl.S2==1]) 
  res$NJ5.S2_Eff[which(res$yx==scenario)]   <- quantile(S2.simu$NJ.S2[S2.simu$ccl.S2==1],0.05)
  res$NJ95.S2_Eff[which(res$yx==scenario)]  <- quantile(S2.simu$NJ.S2[S2.simu$ccl.S2==1],0.95)
  res$NJMed.S2_NC[which(res$yx==scenario)]  <- median(S2.simu$NJ.S2[S2.simu$nc.S2==1]) 
  res$NJ5.S2_NC[which(res$yx==scenario)]    <- quantile(S2.simu$NJ.S2[S2.simu$nc.S2==1],0.05)
  res$NJ95.S2_NC[which(res$yx==scenario)]   <- quantile(S2.simu$NJ.S2[S2.simu$nc.S2==1],0.95)
  res$NJMed.S2_Futil[which(res$yx==scenario)]  <- median(S2.simu$NJ.S2[S2.simu$nc.S2==0 & S2.simu$ccl.S2==0]) 
  res$NJ5.S2_Futil[which(res$yx==scenario)]    <- quantile(S2.simu$NJ.S2[S2.simu$nc.S2==0 & S2.simu$ccl.S2==0],0.05)
  res$NJ95.S2_Futil[which(res$yx==scenario)]   <- quantile(S2.simu$NJ.S2[S2.simu$nc.S2==0 & S2.simu$ccl.S2==0],0.95)
  
  return(res)
}
