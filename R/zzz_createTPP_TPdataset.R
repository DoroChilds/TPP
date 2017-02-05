# to do

# createTPP_TPdataset <- function(data.list, configTable, fcStr, idVar){
#   # problem: function not necessary any more? (nobody calls it)
#   # extract different concentrations from configTable
#   concentrations <- extractConc(configTable)
#   temperatures <- configTable$Temperature
#   # pre-define new list
#   conc.list <- list()
#   # find out which proteins are consistantly found
#   protein.df.list <- lapply(names(data.list), function(nam) return(as.character(data.list[[nam]][,idVar])))
#   consist.prot.ids <- Reduce(intersect, protein.df.list)
#   
#   # throw out rows which are not consistantly found
#   consist.data.list <- lapply(names(data.list), function(df.name){
#     consist.prot <- which(data.list[[df.name]][,idVar] %in% consist.prot.indizes)
#     return(data.list[[df.name]][consist.prot,])
#   })
#   
#   # write all dfs next to each other rearrange data.frame list
#   big.df <- do.call(cbind, consist.data.list)
#   # loop over all concentrations and 
#   restructured.list <- lapply(unique(concentrations), function(conc){
#     # get colmnames of concentrations by ids
#     inds <- which(concentrations %in% conc)
#     colns <- sapply(inds, function(i){
#       return(names(concentrations[i]))
#     })
#     column.ids <- c(1, sort(as.vector(sapply(colns, function(n){
#       return(which(colnames(big.df) %in% paste(fcStr,n, sep="")))
#     }))))
#     conc.df <- big.df[,column.ids]
#     colnames(conc.df) <- c("clustername", temperatures)
#     return(conc.df)
#   })
#   names(restructured.list) <- paste("conc", gsub("\\.", "_", unique(concentrations)), sep="")
#   return(restructured.list)
# }
