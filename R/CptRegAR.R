####################


check_data <- function(data, minseglen=3){
  if(!is.array(data) || !is.numeric(data))
    stop("Argument 'data' must be a numerical matrix.")
  if(length(dim(data))!=2)
    stop("Argument 'data' must be a numerical matrix.")
  if(!is.numeric(minseglen) || length(minseglen)>1)
    stop("Argument 'minseglen' is invalid.")
  
  n <- nrow(data)
  p <- ncol(data)-1
  if(p==0) stop("Dimension of data is 1, no regressors found.")
  if(n<p) stop("More regressors than observations.")
  ##Check to see if there is only 1 intercpt regressor
  intercept <- apply(as.matrix(data[,2]),2,function(a){all(a[1]==a)})
  if(sum(intercept)<1){
    warning("Missing intercept regressor. Append 1's to right of data.")
    data <- cbind(data,1)
    p <- p+1
  }else if(sum(intercept)>1){
    i <- which(intercept)[-1]+1
    warning("Multiple intercepts found. Keeping only first instance.")
    data <- data[,-i]
  }
  return(data)
}



ChangepointRegressionAR <- function(data, phi=0,penalty="MBIC", penalty.value=0,
                                    method="AMOC", dist="Normal", minseglen=3, cpts.only=TRUE, shape=0, MBIC=0, tol=1e-07){
  ##Assume all input arguments are entered correctly
  if(!is.logical(cpts.only) && length(cpts.only)>1)
    stop("Argument 'cpts.only' is invalid.")
  
  # create an extended dataset including the lagged data and lagged regressors - in the columns to the right!
  # datalag <- cbind(
  #   data,
  #   rbind(rep(0,ncol(data)), data[-nrow(data),])
  # )
  
  if(method=="AMOC" && dist=="Normal"){
    out <- CptRegAR_AMOC_Normal(data, phi,penalty, penalty.value, minseglen, shape, MBIC, tol)
  }else if(method=="PELT" && dist=="Normal"){
    out <- CptRegAR_PELT_Normal(data, phi,penalty.value, minseglen, shape,MBIC, tol)
  }else{
    stop("Changepoint in regression method not recognised.")
  }
  if(cpts.only){
    return(sort(out$cpts))
  }else{
    return(out)
  }
}


CptRegAR_AMOC_Normal <- function(data, phi=0,penalty="MBIC", penalty.value=0, minseglen=3,
                                 shape=0, MBIC=0, tol=1e-07){
  n <- as.integer(nrow(data))
  p <- as.integer(ncol(data) -1)
  if(p<1 || n<p) stop("Invalid data dimensions.")
  #tol <- 1e-07 #Rank tolerance (see lm.fit)
  #shape <- -1 #-1=RSS,0=-2logLik,>0=-2logLik with this fixed variance
  if(!is.numeric(shape) || length(shape)!=1)
    stop("Argument 'shape' is invalid.")
  
  
  ##Clean-up on exit
  answer=list()
  answer[[5]]=1
  on.exit(.C("Free_CptRegAR_Normal_AMOC",answer[[6]]))
  
  
  answer <- .C("CptRegAR_Normal_AMOC", data=as.double(data), phi=as.double(phi),
               n=as.integer(n),m=as.integer(ncol(data)), pen=as.double(penalty.value), err=0L,
               shape=as.double(shape), minseglen=as.integer(minseglen), tol=as.double(tol),
               tau=0L, nulllike=vector("double",1), taulike=vector("double",1), 
               tmplike=vector("double",n), MBIC=as.integer(MBIC))
  
  
  #Check if error has occured
  if(answer$err!=0) stop("C code error:",answer$err,call.=F)
  tmp <- c(answer$tau,answer$nulllike,answer$taulike)
  #Following line in cpt.mean (RSS) but not in old cpt.reg (-2Loglik)
  #if(penalty=="MBIC") tmp[3] = tmp[3] + log(tmp[1]) + log(n-tmp[1]+1)
  out <- changepoint::decision(tau = tmp[1], null = tmp[2], alt = tmp[3], 
                               penalty = penalty, n=n, diffparam=p, pen.value = penalty.value)
  names(out) <- c("cpts","pen.value")
  return(out)  #return list of cpts & pen.value
}


CptRegAR_PELT_Normal <- function(data, phi=0,penalty.value=0, minseglen=3, shape=0, 
                                 MBIC=0, tol=1e-07){
  n <- as.integer(nrow(data))
  p <- as.integer(ncol(data) -1)
  #tol <- 1e-07 #Rank tolerance (see lm.fit)
  #shape <- -1 #-1=RSS,0=-2logLik,>0=-2logLik with this fixed variance
  if(!is.numeric(shape) || length(shape)!=1)
    stop("Argument 'shape' is invalid.")
  
  #Check if error has occurred
  answer=list()
  answer[[6]]=1
  on.exit(.C("Free_CptRegAR_PELT_new",answer[[6]]))
  
  y=data[,1]
  design=data[,-1,drop=FALSE]
  p=as.integer(ncol(design))
  answer <- .C("CptRegAR_PELT", y=y,design=design, phi=as.double(phi),n=as.integer(n),
               p=as.integer(ncol(design)), pen=as.double(penalty.value), cptsout=integer(n),
               error=0L, shape=as.double(shape), minseglen=as.integer(minseglen), 
               tol=as.double(tol), lastchangelike=numeric(n+1), 
               lastchangecpts=integer(n+1), numchangecpts=integer(n+1), 
               MBIC=as.integer(MBIC),ncheck=integer(n+1))
  
  
  if(answer$error!=0){
    stop("C code error:",answer$error,call.=F)
  }
  return(list(lastchangecpts=answer$lastchangecpts,
              cpts=sort(answer$cpt[answer$cpt>0]), 
              lastchangelike=answer$lastchangelike, 
              ncpts=answer$numchangecpts,
              ncheck=answer$ncheck))
}






########################################################################