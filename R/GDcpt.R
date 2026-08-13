GDcpt <- function(data,Xc,Xg,method=c("ScaledSSE", "Changingsig"),...){
  method <- match.arg(method)
  if(method=="ScaledSSE"){
    message("Calling ScaledSSE")
    cptswglobreg.SSE(
      data = data,
      Xc = Xc,
      Xg = Xg,
      ...
    )
  }
  else {
    message("Calling Changingsig")
    cptswglobreg.Changingsig(
      data = data,
      Xc = Xc,
      Xg = Xg,
      ...
    )
  }
}
