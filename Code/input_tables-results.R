#
# Simu-VHF
# --------
#

chemin.out <- "/.../"



# ------------------------------------
# Creation of databases for results 
# ------------------------------------


# References of scenarios with [yx]
#-----------------------------------
X <- c(0.75,0.65,0.5,0.45,0.40,0.35)
Y <- c(0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.4)
id <- data.frame(array(NA,c(length(Y),length(X))))
colnames(id)<-X ; rownames(id)<-Y

for (y in (1:length(Y))){
  for (x in (1:length(X))){
    id[y,x]<-paste(y,x)
  }
}

#M Scenarii
M<-length(X)*length(Y)


# Dataframe for results 
#-----------------------

#Scenarii
yx<-as.character(unique(unlist(id)))
pC.yx <- c(rep(X[1],length(Y)),rep(X[2],length(Y)),rep(X[3],length(Y)),
           rep(X[4],length(Y)),rep(X[5],length(Y)),rep(X[6],length(Y)))
DELTA.yx <- rep(Y,length(X))
pI.yx <- pC.yx+DELTA.yx
pI.yx[pI.yx<=0]<-0.01
pI.yx[pI.yx>1]<-1
res <- data.frame(yx,pC.yx,DELTA.yx,pI.yx)


# One random seed
res$graine <- round(runif(1,1234,453211234),0)


# Table for each design
res_F1 <- res
res_F2 <- res
res_S1 <- res
res_S2 <- res

setwd(chemin.out)

res_F1$N.F1<-rep(NA,M)    ; res_F1$pI.F1<-rep(NA,M) ; 
res_F1$p.F1<-rep(NA,M);
res_F1$pMed.F1<-rep(NA,M) ; res_F1$p5.F1<-rep(NA,M) ; res_F1$p95.F1<-rep(NA,M)
write.table(res_F1,"res_F1.csv",sep=";",row.names=F)

res_F2$N.F2<-rep(NA,M)    ; res_F2$pC.F2<-rep(NA,M) ;  res_F2$pI.F2<-rep(NA,M) ; 
res_F2$p.F2<-rep(NA,M)    ; 
res_F2$pMed.F2<-rep(NA,M); res_F2$p5.F2<-rep(NA,M) ; res_F2$p95.F2<-rep(NA,M)
write.table(res_F2,"res_F2.csv",sep=";",row.names=F)

res_S1$NJMed.S1<-rep(NA,M);       res_S1$NJ5.S1<-rep(NA,M);       res_S1$NJ95.S1<-rep(NA,M) ; 
res_S1$pct.supNF1<-rep(NA,M);     res_S1$pct.supNminS1<-rep(NA,M); 
res_S1$sig_levelMed.S1<-rep(NA,M) ; res_S1$sig_level5.S1<-rep(NA,M);   res_S1$sig_level95.S1<-rep(NA,M) ;   res_S1$sig_levelNA.S1<-rep(NA,M) ; 
res_S1$pIJ.S1<-rep(NA,M)  ;       
res_S1$p.S1<-rep(NA,M) ;          res_S1$nc.S1<-rep(NA,M) ;       res_S1$futil.S1<-rep(NA,M) ;
res_S1$NJMed.S1_Eff<-rep(NA,M);   res_S1$NJ5.S1_Eff<-rep(NA,M);   res_S1$NJ95.S1_Eff<-rep(NA,M) ; 
res_S1$NJMed.S1_NC<-rep(NA,M);    res_S1$NJ5.S1_NC<-rep(NA,M);    res_S1$NJ95.S1_NC<-rep(NA,M) ; 
res_S1$NJMed.S1_Futil<-rep(NA,M); res_S1$NJ5.S1_Futil<-rep(NA,M); res_S1$NJ95.S1_Futil<-rep(NA,M) ; 
write.table(res_S1,"res_S1.csv",sep=";",row.names=F)

res_S2$NJMed.S2<-rep(NA,M);      res_S2$NJ5.S2<-rep(NA,M);      res_S2$NJ95.S2<-rep(NA,M) ;
res_S1$pct.supNF2<-rep(NA,M);    res_S1$pct.supNminS2<-rep(NA,M); 
res_S1$sig_levelMed.S2<-rep(NA,M); res_S1$sig_level5.S2<-rep(NA,M);  res_S1$sig_level95.S2<-rep(NA,M) ;  res_S2$sig_levelNA.S2<-rep(NA,M) ; 
res_S2$pCJ.S2<-rep(NA,M);        res_S2$pIJ.S2<-rep(NA,M); 
res_S2$p.S2<-rep(NA,M) ;         res_S2$nc.S2<-rep(NA,M) ;      res_S2$futil.S2<-rep(NA,M) ;
res_S2$NJMed.S2_Eff<-rep(NA,M);  res_S2$NJ5.S2_Eff<-rep(NA,M);  res_S2$NJ95.S2_Eff<-rep(NA,M) ; 
res_S2$NJMed.S2_NC<-rep(NA,M);   res_S2$NJ5.S2_NC<-rep(NA,M);   res_S2$NJ95.S2_NC<-rep(NA,M) ; 
res_S2$NJMed.S2_Futil<-rep(NA,M);res_S2$NJ5.S2_Futil<-rep(NA,M);res_S2$NJ95.S2_Futil<-rep(NA,M) ; 
write.table(res_S2,"res_S2.csv",sep=";",row.names=F)
