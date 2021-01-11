
library(bnlearn)
bayesnet<-empty.graph(colnames(datm))

bayesnet<- set.arc(bayesnet,'PASI','psmj')
bayesnet<- set.arc(bayesnet,'BSA','psmj')

bayesnet<- set.arc(bayesnet,'SCC','syzb')
bayesnet<- set.arc(bayesnet,'BT','syzb')
bayesnet<- set.arc(bayesnet,'BJS','syzb')
bayesnet<- set.arc(bayesnet,'TNF','syzb')

bayesnet<- set.arc(bayesnet,'C3','BT')
bayesnet<- set.arc(bayesnet,'C4','BT')

bayesnet<- set.arc(bayesnet,'IL10','BJS')
bayesnet<- set.arc(bayesnet,'IL17','BJS')
bayesnet<- set.arc(bayesnet,'IL22','BJS')
bayesnet<- set.arc(bayesnet,'IL23','BJS')

bayesnet<- set.arc(bayesnet,'DLQI','shzl')
bayesnet<- set.arc(bayesnet,'SAS','shzl')
bayesnet<- set.arc(bayesnet,'SDS','shzl')

bayesnet<- set.arc(bayesnet,'XQ','bszz')
bayesnet<- set.arc(bayesnet,'CCS','bszz')
bayesnet<- set.arc(bayesnet,'PSQI','bszz')


bayesnet<- set.arc(bayesnet,'psmj','assessment')
bayesnet<- set.arc(bayesnet,'syzb','assessment')
bayesnet<- set.arc(bayesnet,'shzl','assessment')
bayesnet<- set.arc(bayesnet,'bszz','assessment')

library(Rgraphviz)
graphviz.plot(bayesnet)

######


scale<-read.table("./newdata2/score_weight.txt",T,row.names = 1,"\t")

library(easyAHP)
weight<-data.frame()
a<-easyAHP(as.matrix(rowMeans(scale)[c("PASI","BSA")]))
b<-data.frame(a$Makers$Weights)
b$CI<-a$Makers$CI
b$CR<-a$Makers$CR
weight<-rbind(weight,b)


a<-easyAHP(as.matrix(rowMeans(scale)[c("C3","C4")]))
b<-data.frame(a$Makers$Weights)
b$CI<-a$Makers$CI
b$CR<-a$Makers$CR
weight<-rbind(weight,b)


a<-easyAHP(as.matrix(rowMeans(scale)[c("SCC","BT","BJS","TNF")]))
b<-data.frame(a$Makers$Weights)
b$CI<-a$Makers$CI
b$CR<-a$Makers$CR
weight<-rbind(weight,b)

a<-easyAHP(as.matrix(rowMeans(scale)[c("IL10", "IL17","IL22","IL23")]))
b<-data.frame(a$Makers$Weights)
b$CI<-a$Makers$CI
b$CR<-a$Makers$CR
weight<-rbind(weight,b)


a<-easyAHP(as.matrix(rowMeans(scale)[c("DLQI","SAS","SDS")]))
b<-data.frame(a$Makers$Weights)
b$CI<-a$Makers$CI
b$CR<-a$Makers$CR
weight<-rbind(weight,b)

a<-easyAHP(as.matrix(rowMeans(scale)[c( "XQ","CCS","PSQI")]))
b<-data.frame(a$Makers$Weights)
b$CI<-a$Makers$CI
b$CR<-a$Makers$CR
weight<-rbind(weight,b)


a<-easyAHP(as.matrix(rowMeans(scale)[c("psmj","syzb","shzl","bszz")]))
b<-data.frame(a$Makers$Weights)
b$CI<-a$Makers$CI
b$CR<-a$Makers$CR
weight<-rbind(weight,b)
colnames(weight)[1]<-"weight"

#####

data<-read.csv("./newdata/data.csv",T)

data<-data[,-c(1:5)]

dat<-data.frame(apply(data,2,function(x){(x-min(x))/(max(x)-min(x))*10}))

dat$psmj<-dat$PASI*weight["PASI","weight"]+dat$BSA*weight["BSA","weight"]

dat$BT<-dat$C3*weight["C3","weight"]+dat$C4*weight["C4","weight"]

dat$BJS<-dat$IL22*weight["IL22","weight"]+dat$IL10*weight["IL10","weight"]+
  dat$IL17*weight["IL17","weight"]+dat$IL23*weight["IL23","weight"]

dat$syzb<-dat$SCC*weight["SCC","weight"]+dat$BT*weight["BT","weight"]+
  dat$BJS*weight["BJS","weight"]+dat$TNF*weight["TNF","weight"]

dat$shzl<-dat$DLQI*weight["DLQI","weight"]+
  dat$SAS*weight["SAS","weight"]+dat$SDS*weight["SDS","weight"]

dat$bszz<-dat$XQ*weight["XQ","weight"]+
  dat$CCS*weight["CCS","weight"]+
  dat$PSQI*weight["PSQI","weight"]

dat$all<-dat$psmj*weight["psmj","weight"]+ 
  dat$shzl*weight["shzl","weight"]+
  dat$syzb*weight["syzb","weight"]+
  dat$bszz*weight["bszz","weight"]


datm<- discretize(dat,method='quantile',breaks = 2)

fitted <- bn.fit(bayesnet, datm,method='mle')


w.new<-data.frame(weight=weight[,-c(2:3)])
row.names(w.new)<-row.names(weight)
e<--sum((w.new*log(w.new)))

del<-c()
delta<-1


all.weight<-data.frame(weight=weight$weight) 
row.names(all.weight)<-row.names(weight)
i<-1
all.weight$i<-i
all.weight$variable<-row.names(all.weight)

  while (abs(delta)>0.00000001) {
  E1<--sum((w.new*log(w.new)))
  
  #PASI+BAS
  prob.pas<-data.frame(fitted$PASI$prob)
  names(prob.pas)[1]<-"PASI"
  prob.pas<-merge(data.frame(PASI=datm$PASI),prob.pas,by = "PASI")
  
  w1<-w.new[row.names(w.new)=="PASI",1]
  
  td<-sum(prob.pas$Freq/w1)
  w1<-w1+td*0.02
  
  prob.BAS<-data.frame(fitted$BSA$prob)
  names(prob.BAS)[1]<-"BAS"
  prob.BAS<-merge(data.frame(BAS=datm$BSA),prob.BAS,by = "BAS")
  
  w2<-psmj.V[2]
  td<-sum(prob.BAS$Freq/w2)
  w2<-w2+td*0.02
  
  w.new[row.names(w.new)=="PASI",1]<-w1/sum(w1,w2)
  w.new[row.names(w.new)=="BSA",1]<-w2/sum(w1,w2)
  
  ##C3+C4
  prob.C3<-data.frame(fitted$C3$prob)
  names(prob.C3)[1]<-"C3"
  prob.C3<-merge(data.frame(C3=datm$C3),prob.C3,by = "C3")
  
  w1<-w.new[row.names(w.new)=="C3",1]
  
  td<-sum(prob.C3$Freq/w1)
  w1<-w1+td*0.02
  
  
  prob.C4<-data.frame(fitted$C4$prob)
  names(prob.C4)[1]<-"C4"
  prob.C4<-merge(data.frame(C4=datm$C4),prob.C4,by = "C4")
  
  w2<-w.new[row.names(w.new)=="C4",1]
  
  td<-sum(prob.C4$Freq/w2)
  w2<-w2+td*0.02
  
  
  w.new[row.names(w.new)=="C3",1]<-w1/sum(w1,w2)
  w.new[row.names(w.new)=="C4",1]<-w2/sum(w1,w2)
  
  ##IL2+IL17+IL22+IL10
  
  prob.IL10<-data.frame(fitted$IL10$prob)
  names(prob.IL10)[1]<-"IL10"
  prob.IL10<-merge(data.frame(IL10=datm$IL10),prob.IL10,by = "IL10")
  
  w1<-w.new[row.names(w.new)=="IL10",1]
  
  td<-sum(prob.IL10$Freq/w1)
  w1<-w1+td*0.02
  
  prob.IL17<-data.frame(fitted$IL17$prob)
  names(prob.IL17)[1]<-"IL17"
  prob.IL17<-merge(data.frame(IL17=datm$IL17),prob.IL17,by = "IL17")
  
  w2<-w.new[row.names(w.new)=="IL17",1]
  
  td<-sum(prob.IL17$Freq/w2)
  w2<-w2+td*0.02
  
  
  
  prob.IL22<-data.frame(fitted$IL22$prob)
  names(prob.IL22)[1]<-"IL22"
  prob.IL22<-merge(data.frame(IL22=datm$IL22),prob.IL22,by = "IL22")
  
  w3<-w.new[row.names(w.new)=="IL22",1]
  
  td<-sum(prob.IL22$Freq/w3)
  w3<-w3+td*0.02
  
  
  prob.IL23<-data.frame(fitted$IL23$prob)
  names(prob.IL23)[1]<-"IL23"
  prob.IL23<-merge(data.frame(IL23=datm$IL23),prob.IL23,by = "IL23")
  
  w4<-w.new[row.names(w.new)=="IL23",1]
  
  td<-sum(prob.IL23$Freq/w4)
  w4<-w4+td*0.02
  
  w.new[row.names(w.new)=="IL10",1]<-w1/sum(w1,w2,w3,w4)
  w.new[row.names(w.new)=="IL17",1]<-w2/sum(w1,w2,w3,w4)
  w.new[row.names(w.new)=="IL22",1]<-w3/sum(w1,w2,w3,w4)
  w.new[row.names(w.new)=="IL23",1]<-w4/sum(w1,w2,w3,w4)
  
  
  ##SCC+BT+BJS+TNF
  
  prob.SCC<-data.frame(fitted$SCC$prob)
  names(prob.SCC)[1]<-"SCC"
  prob.SCC<-merge(data.frame(SCC=datm$SCC),prob.SCC,by = "SCC")
  
  w1<-w.new[row.names(w.new)=="SCC",1]
  
  td<-sum(prob.SCC$Freq/w1)
  w1<-w1+td*0.02
  
  prob.TNF<-data.frame(fitted$TNF$prob)
  names(prob.TNF)[1]<-"TNF"
  prob.TNF<-merge(data.frame(TNF=datm$TNF),prob.TNF,by = "TNF")
  
  w2<-w.new[row.names(w.new)=="TNF",1]
  
  td<-sum(prob.TNF$Freq/w2)
  w2<-w2+td*0.02
  
  
  prob.BT<-data.frame(fitted$BT$prob)
  prob.BT<-merge(data.frame(BT=datm$BT,C3=datm$C3,C4=datm$C4),prob.BT,by=c("BT","C3","C4"))
  w3<-w.new[row.names(w.new)=="BT",1]
  td<-sum(prob.BT$Freq/w3)
  w3<-w3+td*0.02
  
  
  prob.BJS<-data.frame(fitted$BJS$prob)
  prob.BJS<-merge(data.frame(BJS=datm$BJS,IL10=datm$IL10,IL17=datm$IL17,IL22=datm$IL22,IL23=datm$IL23),
                  prob.BJS,by=c("BJS","IL10","IL17","IL22","IL23"))
  w4<-w.new[row.names(w.new)=="BJS",1]
  td<-sum(prob.BJS$Freq/w4)
  w4<-w4+td*0.02
  
  w.new[row.names(w.new)=="SCC",1]<-w1/sum(w1,w2,w3,w4)
  w.new[row.names(w.new)=="TNF",1]<-w2/sum(w1,w2,w3,w4)
  w.new[row.names(w.new)=="BT",1]<-w3/sum(w1,w2,w3,w4)
  w.new[row.names(w.new)=="BJS",1]<-w4/sum(w1,w2,w3,w4)
  
  
  ####psmj+syzb+shzl+bszz
  prob.psmj<-data.frame(fitted$psmj$prob)
  prob.psmj<-merge(data.frame(psmj=datm$psmj,PASI=datm$PASI,BSA=datm$BSA),
                   prob.psmj,by=c("psmj","PASI","BSA"))
  w1<-w.new[row.names(w.new)=="psmj",1]
  td<-sum(prob.psmj$Freq/w1)
  w1<-w1+td*0.02
  
  
  prob.syzb<-data.frame(fitted$syzb$prob)
  prob.syzb<-merge(data.frame(syzb=datm$syzb,SCC=datm$SCC,BT=datm$BT,BJS=datm$BJS,TNF=datm$TNF),
                   prob.syzb,by=c("syzb","SCC","BT","BJS","TNF"))
  w2<-w.new[row.names(w.new)=="syzb",1]
  td<-sum(prob.syzb$Freq/w2)
  w2<-w2+td*0.02
  
  prob.shzl<-data.frame(fitted$shzl$prob)
  prob.shzl<-merge(data.frame(shzl=datm$shzl,DLQI=datm$DLQI,SAS=datm$SAS,SDS=datm$SDS),
                   prob.shzl,by=c("shzl","DLQI","SAS","SDS"))
  w3<-w.new[row.names(w.new)=="shzl",1]
  td<-sum(prob.shzl$Freq/w3)
  w3<-w3+td*0.02
  
  prob.bszz<-data.frame(fitted$bszz$prob)
  prob.bszz<-merge(data.frame(bszz=datm$bszz,XQ=datm$XQ,CCS=datm$CCS,PSQI=datm$PSQI),
                   prob.bszz,by=c("bszz","XQ","CCS","PSQI"))
  w4<-w.new[row.names(w.new)=="bszz",1]
  td<-sum(prob.bszz$Freq/w4)
  w4<-w4+td*0.02
  
  
  w.new[row.names(w.new)=="psmj",1]<-w1/sum(w1,w2,w3,w4)
  w.new[row.names(w.new)=="syzb",1]<-w2/sum(w1,w2,w3,w4)
  w.new[row.names(w.new)=="shzl",1]<-w3/sum(w1,w2,w3,w4)
  w.new[row.names(w.new)=="bszz",1]<-w4/sum(w1,w2,w3,w4)
  ###DLQI+SAS+SDS
  
  
  prob.DLQI<-data.frame(fitted$DLQI$prob)
  names(prob.DLQI)[1]<-"DLQI"
  prob.DLQI<-merge(data.frame(DLQI=datm$DLQI),prob.DLQI,by = "DLQI")
  
  w1<-w.new[row.names(w.new)=="DLQI",1]
  
  td<-sum(prob.DLQI$Freq/w1)
  w1<-w1+td*0.02
  
  
  prob.SAS<-data.frame(fitted$SAS$prob)
  names(prob.SAS)[1]<-"SAS"
  prob.SAS<-merge(data.frame(SAS=datm$SAS),prob.SAS,by = "SAS")
  
  w2<-w.new[row.names(w.new)=="SAS",1]
  
  td<-sum(prob.SAS$Freq/w2)
  w2<-w2+td*0.02
  
  
  prob.SDS<-data.frame(fitted$SDS$prob)
  names(prob.SDS)[1]<-"SDS"
  prob.SDS<-merge(data.frame(SDS=datm$SDS),prob.SDS,by = "SDS")
  
  w3<-w.new[row.names(w.new)=="SDS",1]
  
  td<-sum(prob.SDS$Freq/w3)
  w3<-w3+td*0.02
  
  w.new[row.names(w.new)=="DLQI",1]<-w1/sum(w1,w2,w3)
  w.new[row.names(w.new)=="SAS",1]<-w2/sum(w1,w2,w3)
  w.new[row.names(w.new)=="SDS",1]<-w3/sum(w1,w2,w3)
  
  
  E2<--sum((w.new*log(w.new)))
  
  delta<- E2-E1
  
  rm(list = c("w1","w2","w3","w4"))
  
  i<-i+1
  all.weight<-rbind(all.weight,data.frame(weight=w.new[,1],i=i,variable=row.names(w.new)))
  e<-c(e,E2)
  del<-c(del,delta)
  
  print(c(i,E2,abs(delta)))
  
  
}


plot(e,type = "l")

colnames(all.weight)[2]<-"i"

library(ggplot2)

ggplot(aes(y=e,x=i),data = data.frame(e=e,i=1:length(e)))+
  geom_point(col="gray50")+geom_line(col="blue")+
  ylab(label = "Entropy")+xlab(label="Learn number")+theme_bw()

ggplot(aes(x=i,y=weight,group=variable),data=all.weight)+
  geom_line(aes(col=variable),lwd=0.01)+
  geom_point(aes(col=variable),size=0.1)+
  facet_wrap(~variable,scales = "free_y")+
  theme(legend.position = "none")+theme_bw()+theme(legend.position = "none")


ggplot()+theme_bw()+
  geom_col(aes(y=variable,x=weight,fill=method),position = "dodge",width = 0.6,
           data= rbind(data.frame(weight=weight[,1],method="AHP",variable=row.names(weight)),
                       data.frame(weight=w.new[,1],method="self learn",variable=row.names(w.new))))


##
self<- function(x){bn.fit(bayesnet, x,method='mle')
  
  
  fitted <- bn.fit(bayesnet,x,method='mle')
  
  ######
  
  
  
  w.new<-data.frame(weight=weight[,-c(2:3)])
  row.names(w.new)<-row.names(weight)
  e<--sum((w.new*log(w.new)))
  
  del<-c()
  delta<-1
  
  
  all.weight<-data.frame(weight=weight$weight) 
  row.names(all.weight)<-row.names(weight)
  i<-1
  all.weight$i<-i
  all.weight$variable<-row.names(all.weight)
  
  
  
  while (abs(delta)>0.00000001) {
    
    E1<--sum((w.new*log(w.new)))
    
    #PASI+BAS
    prob.pas<-data.frame(fitted$PASI$prob)
    
    
    names(prob.pas)[1]<-"PASI"
    prob.pas<-merge(data.frame(PASI=x$PASI),prob.pas,by = "PASI")
    
    w1<-w.new[row.names(w.new)=="PASI",1]
    
    td<-sum(prob.pas$Freq/w1)
    w1<-w1+td*0.02
    
    prob.BAS<-data.frame(fitted$BSA$prob)
    names(prob.BAS)[1]<-"BAS"
    prob.BAS<-merge(data.frame(BAS=x$BSA),prob.BAS,by = "BAS")
    
    w2<-psmj.V[2]
    td<-sum(prob.BAS$Freq/w2)
    w2<-w2+td*0.02
    
    w.new[row.names(w.new)=="PASI",1]<-w1/sum(w1,w2)
    w.new[row.names(w.new)=="BSA",1]<-w2/sum(w1,w2)
    
    ##C3+C4
    prob.C3<-data.frame(fitted$C3$prob)
    names(prob.C3)[1]<-"C3"
    prob.C3<-merge(data.frame(C3=x$C3),prob.C3,by = "C3")
    
    w1<-w.new[row.names(w.new)=="C3",1]
    
    td<-sum(prob.C3$Freq/w1)
    w1<-w1+td*0.02
    
    
    prob.C4<-data.frame(fitted$C4$prob)
    names(prob.C4)[1]<-"C4"
    prob.C4<-merge(data.frame(C4=x$C4),prob.C4,by = "C4")
    
    w2<-w.new[row.names(w.new)=="C4",1]
    
    td<-sum(prob.C4$Freq/w2)
    w2<-w2+td*0.02
    
    
    w.new[row.names(w.new)=="C3",1]<-w1/sum(w1,w2)
    w.new[row.names(w.new)=="C4",1]<-w2/sum(w1,w2)
    
    ##IL2+IL17+IL22+IL10
    
    prob.IL10<-data.frame(fitted$IL10$prob)
    names(prob.IL10)[1]<-"IL10"
    prob.IL10<-merge(data.frame(IL10=x$IL10),prob.IL10,by = "IL10")
    
    w1<-w.new[row.names(w.new)=="IL10",1]
    
    td<-sum(prob.IL10$Freq/w1)
    w1<-w1+td*0.02
    
    prob.IL17<-data.frame(fitted$IL17$prob)
    names(prob.IL17)[1]<-"IL17"
    prob.IL17<-merge(data.frame(IL17=x$IL17),prob.IL17,by = "IL17")
    
    w2<-w.new[row.names(w.new)=="IL17",1]
    
    td<-sum(prob.IL17$Freq/w2)
    w2<-w2+td*0.02
    
    
    
    prob.IL22<-data.frame(fitted$IL22$prob)
    names(prob.IL22)[1]<-"IL22"
    prob.IL22<-merge(data.frame(IL22=x$IL22),prob.IL22,by = "IL22")
    
    w3<-w.new[row.names(w.new)=="IL22",1]
    
    td<-sum(prob.IL22$Freq/w3)
    w3<-w3+td*0.02
    
    
    prob.IL23<-data.frame(fitted$IL23$prob)
    names(prob.IL23)[1]<-"IL23"
    prob.IL23<-merge(data.frame(IL23=x$IL23),prob.IL23,by = "IL23")
    
    w4<-w.new[row.names(w.new)=="IL10",1]
    
    td<-sum(prob.IL23$Freq/w4)
    #print(td)
    w4<-w4+td*0.02
    
    w.new[row.names(w.new)=="IL10",1]<-w1/sum(w1,w2,w3,w4)
    w.new[row.names(w.new)=="IL17",1]<-w2/sum(w1,w2,w3,w4)
    w.new[row.names(w.new)=="IL22",1]<-w3/sum(w1,w2,w3,w4)
    w.new[row.names(w.new)=="IL23",1]<-w4/sum(w1,w2,w3,w4)
    
    
    ##SCC+BT+BJS+TNF
    
    prob.SCC<-data.frame(fitted$SCC$prob)
    names(prob.SCC)[1]<-"SCC"
    prob.SCC<-merge(data.frame(SCC=x$SCC),prob.SCC,by = "SCC")
    
    w1<-w.new[row.names(w.new)=="SCC",1]
    
    td<-sum(prob.SCC$Freq/w1)
    #print(td)
    w1<-w1+td*0.02
    
    prob.TNF<-data.frame(fitted$TNF$prob)
    names(prob.TNF)[1]<-"TNF"
    prob.TNF<-merge(data.frame(TNF=x$TNF),prob.TNF,by = "TNF")
    
    w2<-w.new[row.names(w.new)=="TNF",1]
    
    td<-sum(prob.TNF$Freq/w2)
    #print(td)
    w2<-w2+td*0.02
    
    
    prob.BT<-data.frame(fitted$BT$prob)
    prob.BT<-merge(data.frame(BT=x$BT,C3=x$C3,C4=x$C4),prob.BT,by=c("BT","C3","C4"))
    w3<-w.new[row.names(w.new)=="BT",1]
    td<-sum(prob.BT$Freq/w3)
    w3<-w3+td*0.02
    
    
    prob.BJS<-data.frame(fitted$BJS$prob)
    prob.BJS<-merge(data.frame(BJS=x$BJS,IL10=x$IL10,IL17=x$IL17,IL22=x$IL22,IL23=x$IL23),
                    prob.BJS,by=c("BJS","IL10","IL17","IL22","IL23"))
    w4<-w.new[row.names(w.new)=="BJS",1]
    td<-sum(prob.BJS$Freq/w4)
    w4<-w4+td*0.02
    
    w.new[row.names(w.new)=="SCC",1]<-w1/sum(w1,w2,w3,w4)
    w.new[row.names(w.new)=="TNF",1]<-w2/sum(w1,w2,w3,w4)
    w.new[row.names(w.new)=="BT",1]<-w3/sum(w1,w2,w3,w4)
    w.new[row.names(w.new)=="BJS",1]<-w4/sum(w1,w2,w3,w4)
    
    
    ####psmj+syzb+shzl+bszz
    prob.psmj<-data.frame(fitted$psmj$prob)
    prob.psmj<-merge(data.frame(psmj=x$psmj,PASI=x$PASI,BSA=x$BSA),
                     prob.psmj,by=c("psmj","PASI","BSA"))
    w1<-w.new[row.names(w.new)=="psmj",1]
    td<-sum(prob.psmj$Freq/w1)
    w1<-w1+td*0.02
    
    
    prob.syzb<-data.frame(fitted$syzb$prob)
    prob.syzb<-merge(data.frame(syzb=x$syzb,SCC=x$SCC,BT=x$BT,BJS=x$BJS,TNF=x$TNF),
                     prob.syzb,by=c("syzb","SCC","BT","BJS","TNF"))
    w2<-w.new[row.names(w.new)=="syzb",1]
    td<-sum(prob.syzb$Freq/w2)
    w2<-w2+td*0.02
    
    prob.shzl<-data.frame(fitted$shzl$prob)
    prob.shzl<-merge(data.frame(shzl=x$shzl,DLQI=x$DLQI,SAS=x$SAS,SDS=x$SDS),
                     prob.shzl,by=c("shzl","DLQI","SAS","SDS"))
    w3<-w.new[row.names(w.new)=="shzl",1]
    td<-sum(prob.shzl$Freq/w3)
    w3<-w3+td*0.02
    
    prob.bszz<-data.frame(fitted$bszz$prob)
    prob.bszz<-merge(data.frame(bszz=x$bszz,XQ=x$XQ,CCS=x$CCS,PSQI=x$PSQI),
                     prob.bszz,by=c("bszz","XQ","CCS","PSQI"))
    w4<-w.new[row.names(w.new)=="bszz",1]
    td<-sum(prob.bszz$Freq/w4)
    w4<-w4+td*0.02
    
    
    w.new[row.names(w.new)=="psmj",1]<-w1/sum(w1,w2,w3,w4)
    w.new[row.names(w.new)=="syzb",1]<-w2/sum(w1,w2,w3,w4)
    w.new[row.names(w.new)=="shzl",1]<-w3/sum(w1,w2,w3,w4)
    w.new[row.names(w.new)=="bszz",1]<-w4/sum(w1,w2,w3,w4)
    ###DLQI+SAS+SDS
    
    
    prob.DLQI<-data.frame(fitted$DLQI$prob)
    names(prob.DLQI)[1]<-"DLQI"
    prob.DLQI<-merge(data.frame(DLQI=x$DLQI),prob.DLQI,by = "DLQI")
    
    w1<-w.new[row.names(w.new)=="DLQI",1]
    
    td<-sum(prob.DLQI$Freq/w1)
    w1<-w1+td*0.02
    
    
    prob.SAS<-data.frame(fitted$SAS$prob)
    names(prob.SAS)[1]<-"SAS"
    prob.SAS<-merge(data.frame(SAS=x$SAS),prob.SAS,by = "SAS")
    
    w2<-w.new[row.names(w.new)=="SAS",1]
    
    td<-sum(prob.SAS$Freq/w2)
    w2<-w2+td*0.02
    
    
    prob.SDS<-data.frame(fitted$SDS$prob)
    names(prob.SDS)[1]<-"SDS"
    prob.SDS<-merge(data.frame(SDS=x$SDS),prob.SDS,by = "SDS")
    
    w3<-w.new[row.names(w.new)=="SDS",1]
    
    td<-sum(prob.SDS$Freq/w3)
    w3<-w3+td*0.02
    
    w.new[row.names(w.new)=="DLQI",1]<-w1/sum(w1,w2,w3)
    w.new[row.names(w.new)=="SAS",1]<-w2/sum(w1,w2,w3)
    w.new[row.names(w.new)=="SDS",1]<-w3/sum(w1,w2,w3)
    
    
    E2<--sum((w.new*log(w.new)))
    
    delta<- E2-E1
    
    rm(list = c("w1","w2","w3","w4"))
    
    i<-i+1
    all.weight<-rbind(all.weight,data.frame(weight=w.new[,1],i=i,variable=row.names(w.new)))
    e<-c(e,E2)
    del<-c(del,delta)
    
    
    
  } 
  
  return(list(entropy=e,weight=all.weight))
  
}  

weight.iter<-data.frame()
entropy.iter<-data.frame()

set.seed(1)
for(l in seq(5,100,5)){
  entropy<-data.frame()
  wei<-data.frame()
  for(i in 1:20){print(c(l,i))
    x<-datm[sample(1:dim(datm)[1],l,F),]
    dim(x)
    res<-self(x)
    plot(res$entropy,main=paste("sample_size=",l,i))
    entro<-data.frame(i=1:length(res$entropy),entropy=res$entropy)
    entro$iter_sample<-i
    entro$size<-l
    entropy<-rbind(entropy,entro)
    
    w<-res$weight
    w$iter_sample<-i
    w$size<-l
    wei<-rbind(wei,w)
    
  }
  
  weight.iter<-rbind(weight.iter,wei)
  entropy.iter<-rbind(entropy.iter,entropy)
  
  
}

library(psych)
library(reshape2)


learn_num<-aggregate(entropy.iter$i,by=list(size=entropy.iter$size,iter_sample=entropy.iter$iter_sample),max)
colnames(learn_num)[3]<-"learn_numbers"


for(i in 1:dim(learn_num)[1]){
  learn_num$'entropy'[i]<-entropy.iter[entropy.iter$size==learn_num$size[i] &
                                         entropy.iter$iter_sample == learn_num$iter_sample[i] &
                                         entropy.iter$i==learn_num$`learn_numbers`[i],'entropy']
}


m<-aggregate(learn_num[,3:4],by = list(size=learn_num$size),mean)
s<-aggregate(learn_num[,3:4],by = list(size=learn_num$size),sd)
colnames(s)[2:3]<-paste0(colnames(s)[2:3],".sd")

et.sta<-cbind(m,s)
et.sta<-et.sta[,-1]

ggplot(aes(x=size,y=entropy),data=et.sta)+geom_line(col="gray50")+
  geom_errorbar(aes(ymin=entropy-entropy.sd,ymax=entropy+entropy.sd),width=2,col="gray30")+
  geom_point()+
  theme_b


ggplot(aes(x=size,y=learn_numbers),data=et.sta)+
  geom_errorbar(aes(ymin=learn_numbers-learn_numbers.sd,ymax=learn_numbers+learn_numbers.sd),width=2,col="gray30")+
  geom_point()+
  theme_bw()+geom_smooth(method = "lm",se=F)+
  geom_text(aes(x=25,y=150),label=paste( "r = ",format(cor(et.sta[,c("size","learn_numbers")])[1,2],digits = 3)))


ggplot(aes(x=learn_numbers,y=entropy),data=et.sta)+
  geom_errorbar(aes(xmin=learn_numbers-learn_numbers.sd,xmax=learn_numbers+learn_numbers.sd),col="gray30")+
  geom_errorbar(aes(ymin=entropy-entropy.sd,ymax=entropy+entropy.sd),width=2,col="gray30")+
  geom_point(aes(col=size))+scale_color_gradient(low="blue",high = "red")+
  theme_bw()




weight.sta<-data.frame()
for ( i in unique(weight.iter$size)){
  for(l in unique(weight.iter$iter_sample)){
    a<- learn_num[learn_num$size==i & learn_num$iter_sample == l ,"learn_numbers"]
    wt<-weight.iter[weight.iter$size==i & weight.iter$iter_sample==l &  weight.iter$i==a,]
    weight.sta<-rbind(weight.sta,wt)
  }
}


w.m<-aggregate(weight.sta[,'weight'],by=list(variable=weight.sta$variable,size=weight.sta$size),mean)
colnames(w.m)[3]<-"mean"
w.s<-aggregate(weight.sta[,'weight'],by=list(variable=weight.sta$variable,size=weight.sta$size),sd)

wt.ms<-cbind(w.m,w.s[,-c(1:2)])
colnames(wt.ms)[4]<-"sd"


ggplot(aes(x=size,y=mean,col=variable),data=wt.ms)+
  geom_line(aes(group=variable))+
  geom_point(aes(group=variable))+
  geom_errorbar(aes(group=variable,ymin=mean-sd,ymax=mean+sd),width=1.5)+
  theme_bw()+ylab("Weight")

ggplot(aes(x=size,y=mean,col=variable),data=wt.ms)+
  geom_line(aes(group=variable))+
  geom_point(aes(group=variable))+
  geom_errorbar(aes(group=variable,ymin=mean-sd,ymax=mean+sd),width=1.5)+
  theme_bw()+ylab("Weight")+facet_wrap(~variable,scales = "free_y")+theme(legend.position = "none")



###############ESPA

weight<-all.weight[all.weight$i==142,]
row.names(weight)<-weight$variable


weight<-weight[colnames(data),"weight"]

dat<-data
row.names(dat)<-paste0("l",1:100)


ESPA<-function(dat){
  u<-apply(dat,2,max)
  v<-apply(dat,2,min)
  
  aa<-data.frame()
  bb<-data.frame()
  cc<-data.frame()
  for(i in 1:dim(dat)[1]){  
    a<-(u- dat[i,])*(u+v-dat[i,])/(u*(u-v))
    b<-2*(u- dat[i,])*(dat[i,]-v)/(u*(u-v))
    c<- dat[i,]*(dat[i,]-v)/(u*(u-v))
    
    aa<-rbind(aa,a)
    bb<-rbind(bb,b)
    cc<-rbind(cc,c)
  }
  
  
  AA<-data.frame()
  BB<-data.frame()
  CC<-data.frame()
  for(i in 1:dim(dat)[1]){ 
    a<-dat[i,]*(dat[i,]-v)/(u*(u-v))
    b<-2*(u- dat[i,])*(dat[i,]-v)/(u*(u-v))
    c<- (u- dat[i,])*(u+v-dat[i,])/(u*(u-v))
    
    AA<-rbind(AA,a)
    BB<-rbind(BB,b)
    CC<-rbind(CC,c)
  }
  
  aa$IL10<-AA$IL10  
  bb$IL10<-BB$IL10
  cc$IL10<-CC$IL10
  
  
  CDC<-aa-cc
  
  DD_CDC<-abs(aa-matrix(rep(1,100*16),nrow=100))
  
  DD_a<-abs(aa-1)
  DD_b<-abs(bb-0)
  DD_c<-abs(cc-0)
  
  SD_CD<-1-(DD_CDC+DD_a+DD_b+DD_c)/4
  
  ED_CD<-c()
  for(i in 1:dim(SD_CD)[1]){
    
    s<-sqrt(sum((weight*(1-SD_CD[i,]))^2))
    ED_CD<-c(ED_CD,s)
  }
  
  names(ED_CD)<-row.names(SD_CD)
  
  
  return(ED_CD)
}


ESPA(dat = dat )




library(foreach)
library(doParallel)

cl<- makeCluster(50)   
registerDoParallel(cl)

doespa<-function(dat,i,l){
  da<-dat
  da[l,]<-da[l,]*i
  ed<-ESPA(dat = da )
  ed<-data.frame(ed)
  ed$index<-row.names(ed)
  ed$change<-row.names(dat)[l]
  ed$raio<-i
  return(ed)
}

ED<-data.frame()
for(i in ratio){
  print(i)
  d<-foreach(l=1:dim(dat)[1],.combine=rbind) %dopar%
    doespa(dat,i,l)
  ED<-rbind(ED,d)
}


ED$index<-factor(ED$index,levels = paste0("l",1:100))
ED$change<-factor(ED$change,levels =paste0("l",1:100) )
ggplot(aes(x=raio,y=ed),data=ED[ED$index %in% paste0("l",1:20) & ED$change%in%paste0("l",1:20), ])+
  geom_line(aes(group=index,col=index))+
  geom_point(aes(group=index,col=index),size=0.8)+
  facet_wrap(~change,nrow = 5)+
  theme_bw()+theme(legend.title = element_text(colour = "white") )+
  ylab("ED")+xlab("Ratio of input values")















