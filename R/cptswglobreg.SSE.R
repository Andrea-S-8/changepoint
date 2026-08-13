# ManuaLoaScaledSSE
# ____________________________________________________________________________________________________
# THIS FUNCTION WILL RUN EITHER EM OR GD with GLOBAL AR(order) AND GLOBAL REGRESSION Xg
cptswglobreg.SSE=function(data,Xc,Xg,order=1,maxit,pen0=NULL,GD=TRUE,
                          L=NULL,tol=0.0001,bginit=NULL,verbose=FALSE){
  run_cptregAR=function(y,Xc,phi,pen0,minseglen){
    # cpt.regAR expects response + regressors in one matrix
    regdata=cbind(y,Xc)
    
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
    
    cpts=cpts(fit)
    
    # Calculate fitted changing-regression component
    muhatc <- numeric(length(y))
    tch <- c(0, cpts, length(y))
    
    for (i in seq_len(length(tch) - 1)) {
      seg <- (tch[i] + 1):tch[i + 1]
      fitseg <- lm(
        y[seg] ~ 0 + Xc[seg, , drop = FALSE]
      )
      muhatc[seg] <- fitted(fitseg)
    }
    
    muhatc <- as.vector(muhatc)
    
    list(
      cpts=cpts,
      muhatc=muhatc
    )
  }
  

  y=data
  n=length(y)
  
  lam=eigen(t(Xg)%*%Xg)$values
  
  if(is.null(L))
    L=2/(min(lam)+max(lam))
  
  ng=ncol(Xg)
  nc=ncol(Xc)
  
  numchangingpars=qr(Xc)$rank
  
  if(is.null(pen0))
    pen0=log(n)*(numchangingpars+1)

  #FULL FIT NO CHANGEPOINTS
  X=cbind(Xg,Xc)
  fitnc=lm(y~0+X)
  bglob=as.vector(fitnc$coef[1:ng])
  ytilde=as.vector(y-Xg%*%bglob)
  enc=lm(ytilde~0+Xc)$resid
  nregpars=qr(X)$rank
  
  phi=arima(enc,order=c(order,0,0),include.mean=FALSE)$coef[1:order]
  
  signc=arima(enc,order=c(order,0,0),include.mean=FALSE,fixed=phi)$sigma2
  sig=signc  
  objnc=n+log(n)*(nregpars+1+order*(max(abs(phi)>0)))
  print(objnc)
  globalpars=c(phi,bglob)
  penloss=objnc
  nchanges=0
  
  #Initial Phi 
  ytilde_scaled=ytilde/sqrt(signc)
  out=run_cptregAR(
    y=ytilde_scaled,
    Xc=Xc,
    phi=phi,
    pen0=pen0,
    minseglen=5+ncol(Xc)+ncol(Xg)
  )
  cpts=out$cpts
  if(verbose){print(cpts)}
  nchanges=c(nchanges,length(cpts))
  iterations=1
  
  grad=matrix(NA,ncol=ng)
  error=Inf
  
  while(
    iterations<=maxit &
    error>tol
  ){
  
    if(length(cpts)==0){
      Xcpts=Xc
      segport=1
      tch=c(0,n)
    }
    else{
      #GIVEN CURRENT SET OF CHANGES FIND CREATE MATRIX to find LS ESTIMATES OF ALL MEAN PARAMETERS CHANGING & GLOBAL
      #UPDATE Ystar and bglob based on current changes
      tch=c(0,cpts,n)
      nc=ncol(Xc)
      ncpts=length(cpts)+1
      Xcpts=matrix(
        nrow=n,
        ncol=ncpts*nc
      )
      seg1=1:cpts[1]
      nch=length(tch)-1
      segport=seg1/n
      
      
      Xcpts[,1:nc]=rbind(Xc[seg1,],matrix(nrow=n-cpts[1],ncol=nc,0))
      # Create remaining segment-specific regressors
      for(i in 2:(length(tch)-1)){
        seg=(tch[i]+1):tch[i+1]
        segport=c(
          segport,
          seg/n
        )
        Xtmp=Xc[seg,]
        Xcpts[,(nc*(i-1)+1):(nc*i)]=rbind(matrix(nrow=tch[i],ncol=nc,0),Xtmp,
          matrix(nrow=n-tch[i+1],ncol=nc,0))
      }
    }
    
    
    #NOW UPDATE GLOBAL PARAMETERS.  EITHER NR(EM) OR GD
    if(!GD){
      # FULL FIT
      X=cbind(Xg,Xcpts)
      lmfullfit=lm(y~0+X)
      ng=ncol(Xg)
      # Update global regression parameters
      bglob=lmfullfit$coef[1:ng]
      ytilde=as.vector(y-Xg%*%bglob)
  
      # Update Phi
      lsres=lm(ytilde~0+Xcpts)$resid
      fit=arima(lsres,order=c(order,0,0),include.mean=FALSE)
      phi=fit$coef[1:order]
      sig=fit$sigma2
    }
    else{
      #GRADIENT STEP FOR bglob 
      
      grad=rbind(grad,as.vector(t(Xg)%*%(y-
              out$muhatc*sqrt(signc)-Xg%*%bglob)))
      bglob=as.vector(
        bglob+
          L*grad[iterations+1,]
      )
      muglob=Xg%*%bglob
      ytilde=as.vector(y-muglob)
      
      # ------------------------------------------------------------
      # SCALE THE RESPONSE
      # Match the original manual ScaledSSE implementation.
      # ------------------------------------------------------------
      
      ytilde_scaled=
        ytilde/sqrt(signc)
      
      # UPDATE MUHATC
      muhatc=lm(ytilde_scaled~0+Xcpts)$fitted
      # Update phi using fitted values on the ORIGINAL scale
      yhat=muglob+muhatc*sqrt(signc)
      e=y-yhat
      fit=arima(e,order=c(order,0,0),include.mean=FALSE)
      phi=fit$coef[1:order]
      sig=fit$sigma2
    }
    
    globalpars=rbind(globalpars,c(phi,bglob))
    
    #UPDATE LIKELIHOOD BASED ON CHANGEPOINTS AND GLOBAL PARS
    
    ytilde_scaled=ytilde/sqrt(signc)
    loss=-pen0
    for(i in 1:(length(tch)-1)){
      seg=(tch[i]+1):tch[i+1]
      nseg=length(seg)
      ytmp=ytilde_scaled[seg]
      X=Xc[seg,]
      ehat=lm(ytmp~0+X)$resid
      sighat=arima(
        ehat,
        order=c(order,0,0),
        include.mean=FALSE,
        fixed=phi
      )$sigma2
      loss=loss+nseg*sighat+pen0
    }
    
    penobj=loss+log(n)*(nregpars+1+
          order*(max(abs(phi)>0)))
    penloss=c(penloss,penobj)
    print(penobj)
    
    #REESTIMATE CHANGES BASED ON UPDATED GLOBAL PARS 
    out=run_cptregAR(
      y=ytilde_scaled,
      Xc=Xc,
      phi=phi,
      pen0=pen0,
      minseglen=5+ncol(Xc)+ncol(Xg)
    )
    
    cpts=out$cpts
    if(verbose){print(cpts)}
    ncpts=length(cpts)
    nchanges=c(nchanges,length(cpts))

    error=abs(penloss[iterations]-penloss[iterations+1]
    )/penloss[iterations]
    if(iterations==maxit){
      warning("Maximum iterations reached, consider increasing iterations")
    }
    iterations=iterations+1
  }
  
  # Recreate Xcpts using final changepoints
  if(length(cpts)==0){
    Xcpts=Xc
  }else{
    tch=c(0,cpts,n)
    nc=ncol(Xc)
    ncpts=length(cpts)+1
    Xcpts=matrix(
      nrow=n,
      ncol=ncpts*nc
    )
    seg1=1:cpts[1]
    Xcpts[,1:nc]=rbind(
      Xc[seg1,],
      matrix(
        nrow=n-cpts[1],
        ncol=nc,
        0
      )
    )
    for(i in 2:(length(tch)-1)){
      seg=(tch[i]+1):tch[i+1]
      Xtmp=Xc[seg,]
      Xcpts[
        ,
        (nc*(i-1)+1):(nc*i)
      ]=rbind(
        matrix(
          nrow=tch[i],
          ncol=nc,
          0
        ),
        Xtmp,
        matrix(
          nrow=n-tch[i+1],
          ncol=nc,
          0
        )
      )
    }
  }
  
  # FINAL FIT
  X=cbind(Xg,Xcpts)
  lsfit=lm(y~X+0)
  arfit=arima(
    lsfit$resid,
    order=c(order,0,0),
    include.mean=FALSE
  )
  finallike=
    n*log(arfit$sigma2)+
    n*log(2*pi)+
    n+
    pen0*length(cpts)+
    log(n)*
    (
      nregpars+
        1+
        order*(max(abs(phi)>0))
    )
  
  return(
    list(
      changes=cpts,
      penloss=penloss,
      nchanges=nchanges,
      globpars=globalpars,
      grad=grad,
      muhatc=out$muhatc,
      finallike=finallike
    )
  )
}