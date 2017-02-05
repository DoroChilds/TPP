robustNlsPredict <- function(model, newdata){
  ## A wrapper function for predict.nls() that checks whether model actually 
  ## converged before starting the prediction.
  if (class(model)!="try-error" && !is.null(model)){
    if (is.null(newdata)){
      prediction <- predict(model)
    } else {
      prediction <- predict(model, newdata=newdata)
    }
  } else {
    prediction <- rep(NA_real_, times=length(newdata[[1]]))
  }
  return(prediction)
}