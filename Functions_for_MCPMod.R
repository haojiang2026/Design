
get.pTE.true=function(ndose,pT0.true,pT1.true,pE0.true,pE1.true,rho0,rho1){
  
  pTE0.true=matrix(0,nrow=4,ncol=ndose)
  
  pTE1.true=matrix(0,nrow=4,ncol=ndose)
  
  for(j in 1:ndose){
    
    mean1=c(qnorm(pT0.true[j]),qnorm(pE0.true[j]))
    
    mean2=c(qnorm(pT1.true[j]),qnorm(pE1.true[j]))
    
    sigma1=matrix(c(1,rho0,rho0,1),2,2)
    
    sigma2=matrix(c(1,rho1,rho1,1),2,2)
    
    pTE0.true[1,j]=pmvnorm(lower=c(0,0),upper=c(Inf,Inf),mean=mean1,sigma=sigma1)[1] #P(I=0,T=1,E=1)
    
    pTE0.true[2,j]=pmvnorm(lower=c(0,-Inf),upper=c(Inf,0),mean=mean1,sigma=sigma1)[1] #P(I=0,T=1,E=0)
    
    pTE0.true[3,j]=pmvnorm(lower=c(-Inf,0),upper=c(0,Inf),mean=mean1,sigma=sigma1)[1] #P(I=0,T=0,E=1)
    
    pTE0.true[4,j]=pmvnorm(lower=c(-Inf,-Inf),upper=c(0,0),mean=mean1,sigma=sigma1)[1] #P(I=0,T=0,E=0)
    
    pTE1.true[1,j]=pmvnorm(lower=c(0,0),upper=c(Inf,Inf),mean=mean2,sigma=sigma2)[1] #P(I=1,T=1,E=1)
    
    pTE1.true[2,j]=pmvnorm(lower=c(0,-Inf),upper=c(Inf,0),mean=mean2,sigma=sigma2)[1] #P(I=1,T=1,E=0)
    
    pTE1.true[3,j]=pmvnorm(lower=c(-Inf,0),upper=c(0,Inf),mean=mean2,sigma=sigma2)[1] #P(I=1,T=0,E=1)
    
    pTE1.true[4,j]=pmvnorm(lower=c(-Inf,-Inf),upper=c(0,0),mean=mean2,sigma=sigma2)[1] #P(I=1,T=0,E=0)
  }
  
  return(list(pTE0.true=pTE0.true,pTE1.true=pTE1.true))
  
}


get.data=function(stage,ndose,n,pI.ture,pTE.true,ob_data){
  
  pTE0.true=pTE.true$pTE0.true
  
  pTE1.true=pTE.true$pTE1.true
  
  data.current=data.frame(ID=integer(0),stage=integer(0),dose=integer(0),yI=integer(0),yT=integer(0),yE=integer(0))
  
  for (j in 1:ndose){
    
    nj=n[j]
    
    id.start=nrow(ob_data)+nrow(data.current)+1
    
    id.end=nrow(ob_data)+nrow(data.current)+nj
    
    if(nj==0){next}
    
    data=data.frame(ID=integer(nj),stage=integer(nj),dose=integer(nj),yI=integer(nj),yT=integer(nj),yE=integer(nj))
    
    data$ID=id.start:id.end
    
    data$stage=rep(stage,nj)
    
    data$dose=rep(j,nj)
    
    data$yI=rbinom(nj,1,pI.true[j])
    
    for(i in 1:nj){
      
      if(data[i,'yI']==0){ind=which(rmultinom(1,1,prob=pTE0.true[,j])==1)}
      
      if(data[i,'yI']==1){ind=which(rmultinom(1,1,prob=pTE1.true[,j])==1)}
      
      if(ind==1){data[i,c('yT','yE')]=c(1,1)}
      
      if(ind==2){data[i,c('yT','yE')]=c(1,0)}
      
      if(ind==3){data[i,c('yT','yE')]=c(0,1)}
      
      if(ind==4){data[i,c('yT','yE')]=c(0,0)}
      
    }
    
    data.current=rbind(data.current,data)
    
  }
  
  ob_data=rbind(ob_data,data.current)
  
  return(ob_data)
  
}


monitor.T=function(ob_data,ndose){
  
  try({
    
    pT.data=table(ob_data$dose,factor(ob_data$yT,levels=c(0,1)))
    
    pT.data=prop.table(pT.data,margin=1)[,2]
    
    dose=0:(ndose-1)
    
    dat=data.frame(dose=dose,resp=pT.data)
    
    suppressMessages({
      
      fit.emax=fitMod(dose,resp,data=dat,model='emax')
      
      fit.sigEmax=fitMod(dose,resp,data=dat,model='sigEmax')
      
      fit.linlog=fitMod(dose,resp,data=dat,model='linlog')
      
      fit.exponential=fitMod(dose,resp,data=dat,model='exponential')
      
      fit.logistic=fitMod(dose,resp,data=dat,model='logistic')
      
      fit.linear=fitMod(dose,resp,data=dat,model='linear')
      
    })
    
    models=try({
      
      Mods(emax=as.numeric(fit.emax$coefs),
                  
                  sigEmax=as.numeric(fit.sigEmax$coefs),
                  
                  linlog=as.numeric(fit.linlog$coefs),
                  
                  exponential=as.numeric(fit.exponential$coefs),
                  
                  logistic=as.numeric(fit.logistic$coefs),
                  
                  linear=as.numeric(fit.linear$coefs),
                  
                  doses=dose,fullMod=TRUE)
      
    },silent=TRUE)
    
    if(inherits(models,'try-error')){
      
      model.selection='emax'
      
      suppressMessages({model.best=fitMod(dose,resp,data=dat,model=model.selection)})
      
      pT.model=predict(model.best,newdata=dat,predType="full-model")
      
    }else{
      
      model.logistic=glm(yT~factor(dose)-1,data=ob_data,family=binomial)
      
      S=vcov(model.logistic)
      
      suppressWarnings({model.test=MCTtest(dose,resp,dat,S,models=models,type='general',df=Inf)})
      
      ttest=model.test[["tStat"]]
      
      model.selection=names(ttest)[which.max(ttest)]
      
      suppressMessages({model.best=fitMod(dose,resp,data=dat,model=model.selection)})
      
      pT.model=predict(model.best,newdata=dat,predType="full-model")
      
    }
    
    return(list(model.selection=model.selection,model.parameter=model.best$coefs,pT=pT.model))
    
  },silent=TRUE)
}
  

monitor.E=function(ob_data,ndose){
  
  try({
    
    pE.data=table(ob_data$dose,factor(ob_data$yE,levels=c(0,1)))
    
    pE.data=prop.table(pE.data,margin=1)[,2]
    
    dose=0:(ndose-1)
    
    dat=data.frame(dose=dose,resp=pE.data)
    
    suppressMessages({
      
      fit.emax=fitMod(dose,resp,data=dat,model='emax')
      
      fit.sigEmax=fitMod(dose,resp,data=dat,model='sigEmax')
      
      fit.linlog=fitMod(dose,resp,data=dat,model='linlog')
      
      fit.exponential=fitMod(dose,resp,data=dat,model='exponential')
      
      fit.logistic=fitMod(dose,resp,data=dat,model='logistic')
      
      fit.linear=fitMod(dose,resp,data=dat,model='linear')
      
    })
    
    models=try({
      
      Mods(emax=as.numeric(fit.emax$coefs),
           
           sigEmax=as.numeric(fit.sigEmax$coefs),
           
           linlog=as.numeric(fit.linlog$coefs),
           
           exponential=as.numeric(fit.exponential$coefs),
           
           logistic=as.numeric(fit.logistic$coefs),
           
           linear=as.numeric(fit.linear$coefs),
           
           doses=dose,fullMod=TRUE)
      
    },silent=TRUE)
    
    if(inherits(models,'try-error')){
      
      model.selection='emax'
      
      suppressMessages({model.best=fitMod(dose,resp,data=dat,model=model.selection)})
      
      pE.model=predict(model.best,newdata=dat,predType="full-model")
      
    }else{
      
      model.logistic=glm(yE~factor(dose)+0,data=ob_data,family=binomial)
      
      S=vcov(model.logistic)
      
      suppressWarnings({model.test=MCTtest(dose,resp,dat,S,models=models,type='general',df=Inf)})
      
      ttest=model.test[["tStat"]]
      
      model.selection=names(ttest)[which.max(ttest)]
      
      suppressMessages({model.best=fitMod(dose,resp,data=dat,model=model.selection)})
      
      pE.model=predict(model.best,newdata=dat,predType="full-model")
      
    }
    
    return(list(model.selection=model.selection,model.parameter=model.best$coefs,pE=pE.model))
    
  },silent=TRUE)
}


monitor.U=function(adm.set,ndose,toxicity,efficacy,utable){
  
  colnames(utable)=c('yI','yT','yE','omega')
  
  utable.MCPMod=aggregate(omega~yT+yE,data=utable,sum)
  
  utable.MCPMod=arrange(utable.MCPMod,yT,yE)
  
  omega=utable.MCPMod$omega
  
  pT=toxicity$pT
  
  pE=efficacy$pE
  
  utility=omega[1]*(1-pT)*(1-pE)+
          omega[2]*(1-pT)*pE+
          omega[3]*pT*(1-pE)+
          omega[4]*pT*pE
  
  return(list(utility=utility))
  
}


get.adm.set=function(ndose,phi.pT,phi.pE,cf.pT,cf.pE,toxicity,efficacy){
  
  pT=toxicity$pT
  
  pE=efficacy$pE
  
  adm.toxicity=rep(0,ndose)
  
  adm.efficacy=rep(0,ndose)
  
  for(j in 1:ndose){
    
    if(pT[j]<phi.pT*(2-cf.pT)){adm.toxicity[j]=1}
    
    if(pE[j]>phi.pE*cf.pE){adm.efficacy[j]=1}
  }
  
  if(any(adm.toxicity==0)){adm.toxicity[min(which(adm.toxicity==0)):ndose]=0}
  
  if(any(adm.efficacy==0)){adm.efficacy[1:max(which(adm.efficacy==0))]=0}
  
  adm.set=ifelse(adm.toxicity==1 & adm.efficacy==1,1,0)
  
  return(adm.set)
  
}


get.true=function(ndose,pI.true,pTE.true,utable){
  
  colnames(utable)=c('yI','yT','yE','omega')
  
  utable=arrange(utable,yI,yT,yE)
  
  omega=utable$omega
  
  pT0.true=pTE.true$pTE0.true[1,]+pTE.true$pTE0.true[2,]
  
  pE0.true=pTE.true$pTE0.true[1,]+pTE.true$pTE0.true[3,]
  
  pT1.true=pTE.true$pTE1.true[1,]+pTE.true$pTE1.true[2,]
  
  pE1.true=pTE.true$pTE1.true[1,]+pTE.true$pTE1.true[3,]
  
  pT.true=pI.true*pT1.true+(1-pI.true)*pT0.true
  
  pE.true=pI.true*pE1.true+(1-pI.true)*pE0.true
  
  U.true=omega[1]*(1-pI.true)*(1-pT0.true)*(1-pE0.true)+
         omega[2]*(1-pI.true)*(1-pT0.true)*pE0.true+
         omega[3]*(1-pI.true)*pT0.true*(1-pE0.true)+
         omega[4]*(1-pI.true)*pT0.true*pE0.true+
         omega[5]*pI.true*(1-pT1.true)*(1-pE1.true)+
         omega[6]*pI.true*(1-pT1.true)*pE1.true+
         omega[7]*pI.true*pT1.true*(1-pE1.true)+
         omega[8]*pI.true*pT1.true*pE1.true
  
  return(list(pI.true=pI.true,pT.true=pT.true,pE.true=pE.true,U.true=U.true))
  
}


get.oc=function(pI.true,pT0.true,pT1.true,pE0.true,pE1.true,rho0,rho1,
                       phi.pT,phi.pE,cf.pT,cf.pE,ndose,ntrial,ndraw,nstage,nsample,
                       utable,seed_number){
  
  library(DoseFinding)
  library(mvtnorm)
  library(dplyr)
  
  set.seed(seed_number)
  
  mud.vec=rep(0,ndose+1)
  
  patient.vec=matrix(0,nrow=ntrial,ncol=ndose)
  
  pTE.true=get.pTE.true(ndose,pT0.true,pT1.true,pE0.true,pE1.true,rho0,rho1)
  
  #simulation
  for(m in 1:ntrial){
    
    print(m)
    
    ob_data=data.frame(ID=integer(0),stage=integer(0),dose=integer(0),yI=integer(0),yT=integer(0),yE=integer(0))
    
    stage=1
    
    n=rep(nsample[1]/ndose,ndose)
    
    adm.set=rep(1,ndose)
    
    while(stage<=nstage){
      
      ob_data=get.data(stage,ndose,n,pI.ture,pTE.true,ob_data)
      
      toxicity=monitor.T(ob_data,ndose)
      
      efficacy=monitor.E(ob_data,ndose)
      
      if(inherits(toxicity,'try-error') | inherits(efficacy,'try-error')){
        
        mud.vec[ndose+1]=mud.vec[ndose+1]+1
        
        break
      }
      
      adm.set=get.adm.set(ndose,phi.pT,phi.pE,cf.pT,cf.pE,toxicity,efficacy)
      
      if(sum(adm.set)==0){
        
        mud.vec[ndose+1]=mud.vec[ndose+1]+1
        
        break
        
      }
      
      utility=monitor.U(adm.set,ndose,toxicity,efficacy,utable)
      
      if(stage==nstage){
        
        utility=as.numeric(utility$utility)
        
        mud=which.max(utility[which(adm.set==1)])
        
        mud.vec[mud]=mud.vec[mud]+1
        
        break
        
      }
      
      stage=stage+1
      
      utility=as.numeric(utility$utility)
      
      n=rep(0,ndose)
      
      n[which.max(utility[which(adm.set==1)])]=nsample[stage]
      
    }
    
    patient.vec[m,]=as.vector(table(factor(ob_data$dose,levels=1:ndose)))
    
  }
  
  #result
  
  mud.per=mud.vec/ntrial
  
  pat.mean=apply(patient.vec,2,mean)
  
  pat.per=pat.mean/sum(pat.mean)
  
  pat=sum(pat.mean)
  
  result.true=get.true(ndose,pI.true,pTE.true,utable)
  
  return(list(
    
    pI.true=result.true$pI.true,
    
    pT.true=result.true$pT.true,
    
    pE.true=result.true$pE.true,
    
    U.true=result.true$U.true,
    
    mud.per=mud.per,
    
    pat.per=pat.per,
    
    pat.mean=pat.mean,
    
    pat=pat
    
  ))
  
}



