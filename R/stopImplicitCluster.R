stopImplicitCluster <- function(){
  .options <- doParallel:::.options 
  if(exists(".revoDoParCluster", where=.options) && 
       !is.null(.options[['.revoDoParCluster']]))
  {
    stopCluster(.options[['.revoDoParCluster']])
    remove('.revoDoParCluster', envir=.options)
  }
}
