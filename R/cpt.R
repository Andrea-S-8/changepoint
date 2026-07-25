cpt.mean=function(data,penalty="MBIC",pen.value=0,method="PELT",Q=5,test.stat="Normal",class=TRUE,param.estimates=TRUE,minseglen=1){
  checkData(data)
  if(method=="SegNeigh" & minseglen>1){stop("minseglen not yet implemented for SegNeigh method, use PELT instead.")}
  if(minseglen<1){minseglen=1;warning('Minimum segment length for a change in mean is 1, automatically changed to be 1.')}
  if(!((test.stat=="Normal")||(test.stat=="CUSUM"))){ stop("Invalid test statistic, must be Normal or CUSUM") }

  if(penalty == "CROPS"){
    # browser()
    if(is.numeric(pen.value)){
      if(length(pen.value) == 2){
        if(pen.value[2] < pen.value[1]){
          pen.value = rev(pen.value)
        }
        #run range of penalties
        return(CROPS(data=data, method=method, pen.value=pen.value, test.stat=test.stat, class=class, param.est=param.estimates, minseglen=minseglen, func="mean"))
      }else{
        stop('The length of pen.value must be 2')
      }
    }else{
      stop('For CROPS, pen.value must be supplied as a numeric vector and must be of length 2')
    }
  }

  if(test.stat=="Normal"){
    if(method=="AMOC"){
      return(single.mean.norm(data,penalty,pen.value,class,param.estimates,minseglen))
    }
    else if(method=="PELT" || method=="BinSeg"){

      return(multiple.mean.norm(data,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else if(method=="SegNeigh"){
      warning("SegNeigh is computationally slow, use PELT instead")
      return(multiple.mean.norm(data,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else{
      stop("Invalid Method, must be AMOC, PELT, SegNeigh or BinSeg")
    }
  }
  else if(test.stat=="CUSUM"){
    if(method=="AMOC"){
      tmp=single.mean.cusum(data,penalty,pen.value,class,param.estimates,minseglen)
      warning('Traditional penalty values are not appropriate for the CUSUM test statistic')
      return(tmp)
    }
    else if(method=="SegNeigh" || method=="BinSeg"){
      tmp=multiple.mean.cusum(data,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen)
      warning('Traditional penalty values are not appropriate for the CUSUM test statistic')
      return(tmp)
    }
    else{
      stop("Invalid Method, must be AMOC, SegNeigh or BinSeg")
    }
  }
}

#cpt.reg=function(data,penalty="MBIC",pen.value=0,method="AMOC",Q=5,test.stat="Normal",class=TRUE,param.estimates=TRUE){
#	if(test.stat !="Normal"){ stop("Invalid test statistic, must be Normal") }
#	if(method=="AMOC"){
#		return(single.reg.norm(data,penalty,pen.value,class,param.estimates))
#	}
#	else if(method=="PELT" || method=="BinSeg"){
#		return(multiple.reg.norm(data,mul.method=method,penalty,pen.value,Q,class,param.estimates))
#	}
#	else if(method=="SegNeigh"){
#		warning("SegNeigh is computationally slow, use PELT instead")
#		return(multiple.reg.norm(data,mul.method=method,penalty,pen.value,Q,class,param.estimates))
#	}
#	else{
#		stop("Invalid Method, must be AMOC, PELT, SegNeigh or BinSeg")
#	}
#}

cpt.var=function(data,penalty="MBIC",pen.value=0,know.mean=FALSE, mu=NA,method="PELT",Q=5,test.stat="Normal",class=TRUE,param.estimates=TRUE,minseglen=2){
  checkData(data)
  if(method=="SegNeigh" & minseglen>2){stop("minseglen not yet implemented for SegNeigh method, use PELT instead.")}
  if(minseglen<2){minseglen=2;warning('Minimum segment length for a change in variance is 2, automatically changed to be 2.')}
  if(penalty == "CROPS"){
    # browser()
    if(is.numeric(pen.value)){
      if(length(pen.value) == 2){
        if(pen.value[2] < pen.value[1]){
          pen.value = rev(pen.value)
        }
        #run range of penalties
        return(CROPS(data=data, method=method,pen.value=pen.value, test.stat=test.stat, class=class, param.est=param.estimates, minseglen=minseglen, func="var"))
      }else{
        stop('The length of pen.value must be 2')
      }
    }else{
      stop('For CROPS, pen.value must be supplied as a numeric vector and must be of length 2')
    }
  }

  if(test.stat =="Normal"){

    if(method=="AMOC"){
      return(single.var.norm(data,penalty,pen.value,know.mean,mu,class,param.estimates,minseglen))
    }
    else if(method=="PELT" || method=="BinSeg"){

      return(multiple.var.norm(data,mul.method=method,penalty,pen.value,Q,know.mean,mu,class,param.estimates,minseglen))
    }
    else if(method=="SegNeigh"){
      warning("SegNeigh is computationally slow, use PELT instead")
      return(multiple.var.norm(data,mul.method=method,penalty,pen.value,Q,know.mean,mu,class,param.estimates,minseglen))
    }
    else{
      stop("Invalid Method, must be AMOC, PELT, SegNeigh or BinSeg")
    }
  }
  else if(test.stat=="CSS"){
    warning('Traditional penalty values are not appropriate for the CSS test statistic')
    if(method=="AMOC"){
      return(single.var.css(data,penalty,pen.value,class,param.estimates,minseglen))
    }
    else if(method=="PELT" || method=="SegNeigh" || method=="BinSeg"){
      return(multiple.var.css(data,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else{
      stop("Invalid Method, must be AMOC, SegNeigh or BinSeg")
    }
  }
  else{
    stop("Invalid test statistic, must be Normal or CSS")
  }
}

cpt.meanvar=function(data,penalty="MBIC",pen.value=0,method="PELT",Q=5,test.stat="Normal",class=TRUE,param.estimates=TRUE,shape=1,minseglen=2){
  checkData(data)
  if(method=="SegNeigh" & minseglen>2){stop("minseglen not yet implemented for SegNeigh method, use PELT instead.")}
  if(minseglen<2){
    if(!(minseglen==1 & (test.stat=="Poisson"|test.stat=="Exponential"))){
      minseglen=2;warning('Minimum segment length for a change in mean and variance is 2, automatically changed to be 2.')}
  }
  if(penalty == "CROPS"){
    if(is.numeric(pen.value)){
      if(length(pen.value) == 2){
        if(pen.value[2] < pen.value[1]){
          pen.value = rev(pen.value)
        }
        #run range of penalties
        return(CROPS(data=data, method=method,pen.value=pen.value, test.stat=test.stat, class=class, param.est=param.estimates, minseglen=minseglen, shape=shape, func="meanvar"))
      }else{
        stop('The length of pen.value must be 2')
      }
    }else{
      stop('For CROPS, pen.value must be supplied as a numeric vector and must be of length 2')
    }
  }
  if(test.stat=="Normal"){

    if(method=="AMOC"){
      return(single.meanvar.norm(data,penalty,pen.value,class,param.estimates,minseglen))
    }
    else if(method=="PELT" || method=="BinSeg"){

      return(multiple.meanvar.norm(data,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else if(method=="SegNeigh"){
      warning("SegNeigh is computationally slow, use PELT instead")
      return(multiple.meanvar.norm(data,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else{
      stop("Invalid Method, must be AMOC, PELT, SegNeigh or BinSeg")
    }
  }
  else if(test.stat=="Gamma"){
    if(method=="AMOC"){
      return(single.meanvar.gamma(data,shape,penalty,pen.value,class,param.estimates,minseglen))
    }
    else if(method=="PELT" || method=="BinSeg"){
      return(multiple.meanvar.gamma(data,shape,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else if(method=="SegNeigh"){
      warning("SegNeigh is computationally slow, use PELT instead")
      return(multiple.meanvar.gamma(data,shape,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else{
      stop("Invalid Method, must be AMOC, PELT, SegNeigh or BinSeg")
    }
  }
  else if(test.stat=="Exponential"){
    if(method=="AMOC"){
      return(single.meanvar.exp(data,penalty,pen.value,class,param.estimates,minseglen))
    }
    else if(method=="PELT" || method=="BinSeg"){

      return(multiple.meanvar.exp(data,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else if(method=="SegNeigh"){
      warning("SegNeigh is computationally slow, use PELT instead")
      return(multiple.meanvar.exp(data,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else{
      stop("Invalid Method, must be AMOC, PELT, SegNeigh or BinSeg")
    }
  }
  else if(test.stat=="Poisson"){
    if(method=="AMOC"){
      return(single.meanvar.poisson(data,penalty,pen.value,class,param.estimates,minseglen))
    }
    else if(method=="PELT" || method=="BinSeg"){

      return(multiple.meanvar.poisson(data,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else if(method=="SegNeigh"){
      warning("SegNeigh is computationally slow, use PELT instead")
      return(multiple.meanvar.poisson(data,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else{
      stop("Invalid Method, must be AMOC, PELT, SegNeigh or BinSeg")
    }
  }
  else{
    stop("Invalid test statistic, must be Normal, Gamma, Exponential or Poisson")
  }
}

checkData = function(data){
  if(!is.numeric(data)){
    stop("Only numeric data allowed")
  }
  if(anyNA(data)){stop("Missing value: NA is not allowed in the data as changepoint methods are only sensible for regularly spaced data.")}

}


cpt.regAR <- function(data, phi=0, penalty="MBIC", pen.value=0, method="AMOC", dist="Normal",
                      class=TRUE, param.estimates=TRUE, minseglen=3, shape = 0, tol = 1e-07){
  MBIC=0
  ##Check arguments are valid
  if(!is.array(data) || !is.numeric(data))  ##Further checks applied later
    stop("Argument 'data' must be a numerical matrix/array.")
  if(!is.character(penalty) || length(penalty)>1)
    stop("Argument 'penalty' is invalid.")
  #value of 'penalty' & 'pen.value' checked within changepoint::penalty_decision
  if(!is.character(method) || length(method)>1)
    stop("Argument 'method' is invalid.") 
  if(method!="AMOC" && method != "PELT") ##RESTRICTION IN USE
    stop("Invalid method, must be AMOC or PELT.")
  if(!is.character(dist) || length(dist)>1)
    stop("Argument 'dist' is invalid.")
  if(dist != "Normal"){  ##RESTRICTION IN USE
    warning(paste0("dist = ",dist," is not supported. Converted to dist='Normal'"))
    dist <- "Normal"
  }
  if(!is.logical(class) || length(class)>1)
    stop("Argument 'class' is invalid.")
  if(!is.logical(param.estimates) || length(param.estimates)>1)
    stop("Argument 'param.estimates' is invalid.")
  if(!is.numeric(minseglen) || length(minseglen)>1)
    stop("Argument 'minseglen' is invalid.")
  if(minseglen <= 0 || minseglen%%1 != 0) ##Further checks applied later
    stop("Argument 'minseglen' must be positive integer.")
  if(!is.numeric(tol) || length(tol)!=1)
    stop("Argument 'tol' is invalid.")
  if(tol<0) stop("Argument 'tol' must be positive.")
  ##Argument shape is assessed by the command where it is to be used.
  
  
  #Single data set? convert to multiple data set array.
  if(length(dim(data)) == 2) data <- array(data,dim=c(1,dim(data)))
  
  
  ans <- vector("list",dim(data)[1])
  for(i in 1:dim(data)[1]){  ##To each data set.
    #Check format for ith data set
    datai <- check_data(data[i, , ])
    
    
    #Evaluate penalty value
    pen.value <- changepoint::penalty_decision(penalty=penalty,
                                               pen.value=pen.value, n=nrow(datai), diffparam=ncol(datai),
                                               asymcheck="cpt.reg", method=method)
    if(penalty=="MBIC"){MBIC=1}
    
    
    #Check value of minseglen
    if(minseglen < (ncol(datai)-1)){
      warning(paste("minseglen is too small, set to:",ncol(datai)))
      minsegleni <- ncol(datai)
    }else if(nrow(datai) < (2*minseglen)){
      stop("Minimum segment length is too large to include a change in this data.")
    }else{
      minsegleni <- minseglen
    }
    
    
    CPTS <- ChangepointRegressionAR(data=datai, phi=phi, penalty=penalty,
                                    penalty.value=pen.value, method=method, dist=dist, minseglen=minseglen,
                                    shape = shape, tol=tol, cpts.only=class, MBIC=MBIC)
    
    
    if(class){
      #Convert to cpt.reg object
      ansi <- new("cpt.reg")
      data.set(ansi) <- datai
      cpttype(ansi) <- "regressionAR"
      method(ansi) <- method
      distribution(ansi) <- dist
      pen.type(ansi) <- penalty
      pen.value(ansi) <- pen.value
      cpts(ansi) <- CPTS
      if(method=="PELT") ncpts.max(ansi) <- Inf
      if(param.estimates) ansi = param(ansi)
      param.est(ansi)=list(beta=param.est(ansi)$beta,sig2=param.est(ansi)$sig2,phi=phi) # add phi to the parameter list
      ans[[i]] <- ansi
    }else{
      #store what is returned from calculation
      ans[[i]] <- CPTS
    }
  }
  
  
  ##Return result, only first case if there is only a singe data set
  if(dim(data)[1]==1){ return(ans[[1]]) }else{ return(ans) }  
}
