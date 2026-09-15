
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

