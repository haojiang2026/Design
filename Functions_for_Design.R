
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


monitor.I=function(ob_data,ndose,ndraw){
  
  alpha.pI=rep(0,ndose)
  
  beta.pI=rep(0,ndose)
  
  pI.draw=matrix(0,nrow=ndraw,ncol=ndose)
  
  pI.draw.pava=matrix(0,nrow=ndraw,ncol=ndose)
  
  for(j in 1:ndose){
    
    vector.I=ob_data[ob_data$dose==j,'yI']
    
    alpha.pI[j]=1+sum(vector.I==1)
    
    beta.pI[j]=1+sum(vector.I==0)
    
    pI.draw[,j]=rbeta(n=ndraw,shape1=alpha.pI[j],shape2=beta.pI[j])
    
  }
  
  var.pI=(alpha.pI*beta.pI)/(((alpha.pI+beta.pI)^2)*(alpha.pI+beta.pI+1))
  
  for(i in 1:ndraw){pI.draw.pava[i,]=pava(y=pI.draw[i,],w=1/var.pI)}
  
  return(list(pI.draw.pava=pI.draw.pava))
  
}


monitor.T=function(ob_data,ndose,ndraw,immunity){
  
  pI.draw.pava=as.matrix(immunity$pI.draw.pava)
  
  alpha.pT=matrix(0,nrow=2,ncol=ndose)
  
  beta.pT=matrix(0,nrow=2,ncol=ndose)
  
  pT0.draw=matrix(0,nrow=ndraw,ncol=ndose)
  
  pT1.draw=matrix(0,nrow=ndraw,ncol=ndose)
  
  pT0.draw.pava=matrix(0,nrow=ndraw,ncol=ndose)
  
  pT1.draw.pava=matrix(0,nrow=ndraw,ncol=ndose)
  
  pT.draw.pava=matrix(0,nrow=ndraw,ncol=ndose)
  
  for(i in 0:1){
    
    for(j in 1:ndose){
      
      vector.T=ob_data[ob_data$dose==j & ob_data$yI==i,'yT']
      
      alpha.pT[i+1,j]=1+sum(vector.T==1)
      
      beta.pT[i+1,j]=1+sum(vector.T==0)
      
      if(i==0){pT0.draw[,j]=rbeta(n=ndraw,shape1=alpha.pT[i+1,j],shape2=beta.pT[i+1,j])}
        
      if(i==1){pT1.draw[,j]=rbeta(n=ndraw,shape1=alpha.pT[i+1,j],shape2=beta.pT[i+1,j])}
        
    }
  }
  
  var.pT=(alpha.pT*beta.pT)/(((alpha.pT+beta.pT)^2)*(alpha.pT+beta.pT+1))
  
  for(i in 1:ndraw){
    
    d=rbind(pT0.draw[i,],pT1.draw[i,])
  
    d=biviso(y=d,w=1/var.pT)

    pT0.draw.pava[i,]=d[1,]
    
    pT1.draw.pava[i,]=d[2,]
    
  }
  
  pT.draw.pava=pI.draw.pava*pT1.draw.pava+(1-pI.draw.pava)*pT0.draw.pava
  
  return(list(pT0.draw.pava=pT0.draw.pava,pT1.draw.pava=pT1.draw.pava,pT.draw.pava=pT.draw.pava))
  
}


monitor.E=function(ob_data,ndose,ndraw,immunity){
  
  pI.draw.pava=as.matrix(immunity$pI.draw.pava)
  
  alpha.pE=matrix(0,nrow=2,ncol=ndose)
  
  beta.pE=matrix(0,nrow=2,ncol=ndose)
  
  pE0.draw=matrix(0,nrow=ndraw,ncol=ndose)
  
  pE1.draw=matrix(0,nrow=ndraw,ncol=ndose)
  
  pE0.draw.pava=matrix(0,nrow=ndraw,ncol=ndose)
  
  pE1.draw.pava=matrix(0,nrow=ndraw,ncol=ndose)
  
  pE.draw.pava=matrix(0,nrow=ndraw,ncol=ndose)
  
  for(i in 0:1){
    
    for(j in 1:ndose){
      
      vector.E=ob_data[ob_data$dose==j & ob_data$yI==i,'yE']
      
      alpha.pE[i+1,j]=1+sum(vector.E==1)
      
      beta.pE[i+1,j]=1+sum(vector.E==0)
      
      if(i==0){pE0.draw[,j]=rbeta(n=ndraw,shape1=alpha.pE[i+1,j],shape2=beta.pE[i+1,j])}
      
      if(i==1){pE1.draw[,j]=rbeta(n=ndraw,shape1=alpha.pE[i+1,j],shape2=beta.pE[i+1,j])}

    }
  }
  
  var.pE=(alpha.pE*beta.pE)/(((alpha.pE+beta.pE)^2)*(alpha.pE+beta.pE+1))
  
  for(i in 1:ndraw){
    
    d=rbind(pE0.draw[i,],pE1.draw[i,])
    
    d=biviso(y=d,w=1/var.pE)
    
    pE0.draw.pava[i,]=d[1,]
    
    pE1.draw.pava[i,]=d[2,]
    
  }
  
  pE.draw.pava=pI.draw.pava*pE1.draw.pava+(1-pI.draw.pava)*pE0.draw.pava
  
  return(list(pE0.draw.pava=pE0.draw.pava,pE1.draw.pava=pE1.draw.pava,pE.draw.pava=pE.draw.pava))
}


monitor.U=function(immunity,toxicity,efficacy,utable){
  
  colnames(utable)=c('yI','yT','yE','omega')
  
  utable=arrange(utable,yI,yT,yE)
  
  omega=utable$omega
  
  pI=as.matrix(immunity$pI.draw.pava)
  
  pT0=as.matrix(toxicity$pT0.draw.pava)
  
  pT1=as.matrix(toxicity$pT1.draw.pava)
  
  pE0=as.matrix(efficacy$pE0.draw.pava)
  
  pE1=as.matrix(efficacy$pE1.draw.pava)
  
  utility.draw=omega[1]*(1-pI)*(1-pT0)*(1-pE0)+
               omega[2]*(1-pI)*(1-pT0)*pE0+
               omega[3]*(1-pI)*pT0*(1-pE0)+
               omega[4]*(1-pI)*pT0*pE0+
               omega[5]*pI*(1-pT1)*(1-pE1)+
               omega[6]*pI*(1-pT1)*pE1+
               omega[7]*pI*pT1*(1-pE1)+
               omega[8]*pI*pT1*pE1  
  
  return(list(utility.draw=utility.draw))
  
}


get.adm.set=function(ndose,phi.pT,phi.pE,cf.pT,cf.pE,toxicity,efficacy){
  
  pT.draw.pava=as.matrix(toxicity$pT.draw.pava)
  
  pE.draw.pava=as.matrix(efficacy$pE.draw.pava)
  
  adm.toxicity=rep(0,ndose)
  
  adm.efficacy=rep(0,ndose)
  
  for(j in 1:ndose){
    
    if(mean(pT.draw.pava[,j]>phi.pT)<cf.pT){adm.toxicity[j]=1}
    
    if(mean(pE.draw.pava[,j]<phi.pE)<cf.pE){adm.efficacy[j]=1}
  }
  
  adm.set=ifelse(adm.toxicity==1 & adm.efficacy==1,1,0)
  
  return(adm.set)
  
}


get.allocation=function(adm.set,stage,nsample,ndraw,utility){
  
  utility.draw=as.matrix(utility$utility.draw)
  
  utility.draw[,which(adm.set==0)]=0
  
  umean=colMeans(utility.draw)
  
  umax.ids=max.col(utility.draw,ties.method='first')
  
  umax=tabulate(umax.ids,nbins=ncol(utility.draw))
  
  pmud=umax/ndraw
  
  n=as.numeric(rmultinom(1,nsample[stage],pmud))
  
  return(list(n=n,pmud=pmud,umean=umean))
  
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
  
  library(Iso)
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
      
      immunity=monitor.I(ob_data,ndose,ndraw)
      
      toxicity=monitor.T(ob_data,ndose,ndraw,immunity)
      
      efficacy=monitor.E(ob_data,ndose,ndraw,immunity)
      
      adm.set=get.adm.set(ndose,phi.pT,phi.pE,cf.pT,cf.pE,toxicity,efficacy)
      
      if(sum(adm.set)==0){mud.vec[ndose+1]=mud.vec[ndose+1]+1;break}
      
      utility=monitor.U(immunity,toxicity,efficacy,utable)
      
      if(stage==nstage){
        
        allocation=get.allocation(adm.set,stage,nsample,ndraw,utility)
        
        mud=which.max(allocation$umean)
        
        mud.vec[mud]=mud.vec[mud]+1
        
        break
        
      }
      
      stage=stage+1
      
      allocation=get.allocation(adm.set,stage,nsample,ndraw,utility)
      
      n=allocation$n
      
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


















