#ManuaLoaChangingsig
#____________________________________________________________________________________________________
#THIS FUNCTION WILL RUN EITHER EM OR GD with GLOBAL AR(order) AND GLOBAL REGRESSION Xg
cptswglobreg.Changingsig=function(data,Xc,Xg, order=1, maxit,pen0=NULL,phi=0,
                                  intercept=TRUE,GD=TRUE,L=NULL,tol=0.0001){
  run_cptregAR=function(y, Xc, phi, pen0, minseglen) {
    
    # cpt.regAR expects response + regressors in one matrix
    regdata=cbind(y, Xc)
    
    fit=cpt.regAR(
      data=regdata,
      phi=phi,
      penalty="Manual",
      pen.value=pen0,
      method="PELT",
      class=TRUE,
      param.estimates=TRUE,
      minseglen=minseglen
    )
    cpts <- cpts(fit)
    # Calculate fitted changing-regression component
    muhatc=numeric(length(y))
    tch=c(0, cpts, length(y))
    for (i in seq_len(length(tch)-1)){
      seg=(tch[i]+1):tch[i+1]
      fitseg=lm(y[seg]~0+Xc[seg, , drop=FALSE])
      muhatc[seg]=fitted(fitseg)
    }
    list(
      cpts=cpts,
      muhatc=as.vector(muhatc)
    )
  }
  
  
  y=data
  n=length(y)
  lam=eigen(t(Xg)%*%Xg)$values
  if(is.null(L)) L=2/(min(lam)+max(lam))
  ng=ncol(Xg)
  nc=ncol(Xc)
  numchangingpars=qr(Xc)$rank 
  if(is.null(pen0)) pen0=log(n)*(numchangingpars+2) # changing SIC penalty
  Xcorig=Xc
  
  
  #Xcperp=Xc-Xg%*%solve(t(Xg)%*%Xg)%*%t(Xg)%*%Xc
  
  
  #FULL FIT NO CHANGEPOINTS
  X=cbind(Xg,Xc)
  fitnc=lm(y~0+X)
  bglob=as.vector(fitnc$coef[1:ng])
  ytilde=y-Xg%*%bglob
  fit1=lm(ytilde~0+Xc)
  enc=fit1$resid
  nregpars=qr(X)$rank
  if(is.null(phi))
  {
    
    phi=arima(enc,order=c(order,0,0),include.mean=FALSE)$coef[1:order]
  }
  #if(phi==0) phi=rep(0,order)
  signc=arima(enc,order=c(order,0,0),include.mean=FALSE,fixed=phi)$sigma2
  objnc=n*log(2*pi)+n+log(n)*(nregpars+1+order*(max(abs(phi)>0)))+ n*log(signc)
  print(objnc)
  
  
  sig=signc
  
  globalpars=c(phi,bglob)
  loglik=objnc
  nchanges=0
  
  #Initial phi
  out=run_cptregAR(
    y=ytilde,
    Xc=Xc,
    phi=phi,
    pen0=pen0,
    minseglen=5+ncol(Xc)+ncol(Xg)
  )
  cpts=out$cpts
  nchanges=c(nchanges,length(cpts))
  iterations=1
  
  grad=matrix(NA,ncol=ng)
  error=Inf
  while(iterations<=maxit & error>tol){
    
    if(length(cpts)==0){
      Xcpts=Xc
    }
    else{  
      #GIVEN CURRENT SET OF CHANGES FIND LS ESTIMATES OF ALL MEAN PARAMETERS CHANGING & GLOBAL
      #UPDATE Ystar and bglob based on current changes
      tch=c(0,cpts,n)
      nc=ncol(Xc)
      ncpts=length(cpts)+1
      Xcpts=matrix(NA,nrow=n,ncol=ncpts*nc)
      seg1=1:cpts[1]
      nch=length(tch)-1
      
      
      Xcpts[,1:(nc)]=rbind(Xc[seg1, ,drop=FALSE],matrix(nrow=n-cpts[1],ncol=(nc),0))
      #CREATE SEGMENT X
      for(i in 2:(length(tch)-1)){
        seg=(tch[i]+1):tch[i+1]
        Xtmp=Xc[seg, ,drop=FALSE]
        Xcpts[,(nc*(i-1)+1):(nc*i)]=rbind(matrix(nrow=tch[i],ncol=nc,0),Xtmp,matrix(nrow=n-tch[i+1],ncol=nc,0))
        
      } 
    }
    
    
    #NOW UPDATE GLOBAL PARAMETERS.  EITHER NR(EM) OR GD
    if(!GD){ 
      #FULL FIT
      X=cbind(Xg,Xcpts)
      lmfullfit=lm(y~0+X)
      ng=ncol(Xg)
      #UPDATE bglob regression
      bglob=lmfullfit$coef[(1):(ng)]
      ytilde=as.vector(y-Xg%*%bglob)
      #UPDATE PHI
      lsres=lmfullfit$resid
      fit=arima(lsres,order=c(order,0,0),include.mean=FALSE)
      res=fit$resid
      phi=fit$coef[1:order]
      
      
      
      seg=1:n
      if(length(cpts)>0) seg=1:cpts[1]
      nseg=length(seg)
      loss=1:(length(cpts)+1)    
      ytildeseg=ytilde[seg]
      Xseg=Xc[seg,]
      lmseg=lm(ytildeseg~0+Xseg)
      sigseg=arima(lmseg$resid,order=c(order,0,0),include.mean=FALSE,fixed=phi)$sigma2 
      loss=nseg*log(sigseg) 
      if(length(cpts)>0){
        for(i in 2:(length(tch)-1)){
          seg=(tch[i]+1):tch[i+1]
          nseg=length(seg)
          ytildeseg=ytilde[seg]
          Xseg=Xc[seg,]
          lmseg=lm(ytildeseg~0+Xseg)
          sigseg=arima(lmseg$resid,order=c(order,0,0),include.mean=FALSE,fixed=phi)$sigma2 
          loss=loss+nseg*log(sigseg)        
        } }
      loss=loss+n*log(2*pi)+n
      
      
      X=Xcpts
      nc=ncol(X)
      muhatc=lm(ytilde~0+X)$fitted
      ystar=y-muhatc    		
    }
    else{
      #GRADIENT STEP FOR bglob 
      
      grad=rbind(grad,as.vector(t(Xg)%*%(y-out$muhatc-Xg%*%bglob)))
      bglob=bglob+L*grad[iterations+1,]	
      muglob=Xg%*%bglob
      ytilde=y-muglob
      #UPDATE MUHATC based on new bglob and then update phi
      muhatc=lm(ytilde~0+Xcpts)$fitted
      yhat=muglob+muhatc
      e=y-yhat
      fit=arima(e,order=c(order,0,0),include.mean=FALSE)
      res=e
      phi=fit$coef[1:order]
      seg=1:n
      if(length(cpts)>0) seg=1:cpts[1]
      nseg=length(seg)
      loss=1:(length(cpts)+1)    
      ytildeseg=ytilde[seg]
      Xseg=Xc[seg,]
      lmseg=lm(ytildeseg~0+Xseg)
      sigseg=arima(lmseg$resid,order=c(order,0,0),include.mean=FALSE,fixed=phi)$sigma2 
      loss=nseg*log(sigseg) 
      if(length(cpts)>0){
        for(i in 2:(length(tch)-1)){
          seg=(tch[i]+1):tch[i+1]
          nseg=length(seg)
          ytildeseg=ytilde[seg]
          Xseg=Xc[seg,]
          lmseg=lm(ytildeseg~0+Xseg)
          sigseg=arima(lmseg$resid,order=c(order,0,0),include.mean=FALSE,fixed=phi)$sigma2 
          loss=loss+nseg*log(sigseg)        
        } }
      loss=loss+n*log(2*pi)+n
      
      
    } 		
    globalpars=rbind(globalpars,c(phi,bglob))
    
    
    #UPDATE LIKELIHOOD BASED ON CHANGEPOINTS AND GLOBAL PARS	
    penobj=loss+log(n)*(nregpars+1+order*(max(abs(phi)>0)))+(pen0)*(length(cpts))
    loglik=c(loglik,penobj)
    print(penobj)
    
    #REESTIMATE CHANGES BASED ON UPDATED GLOBAL PARS    	
    out=run_cptregAR(
      y=ytilde,
      Xc=Xc,
      phi=phi,
      pen0=pen0,
      minseglen=5+ncol(Xc)+ncol(Xg)
    )
    cpts=out$cpts
    print(cpts)
    ncpts=length(cpts)	
    nchanges=c(nchanges,length(cpts))
    
    error=abs(loglik[iterations]-loglik[iterations+1])/loglik[iterations]
    if(iterations==maxit){warning("Maximum iterations reached, consider increasing iterations")}
    
    
    iterations=iterations+1
  }
  
  return(
    list(
      changes=cpts,
      loglik=loglik,
      nchanges=nchanges,
      globpars=globalpars,
      grad=grad,
      muhatc=muhatc
    )
  )
  
}