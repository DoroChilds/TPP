% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tppccrNormalizeToReference.R
\name{tppccrNormalizeToReference}
\alias{tppccrNormalizeToReference}
\title{Normalize fold changes of TPP-CCR experiment to a reference column}
\usage{
tppccrNormalizeToReference(data, refCol = NULL)
}
\arguments{
\item{data}{expressionSet object containing the data to be normalized}

\item{refCol}{column number to use as a reference. Will contain only 1s after 
the normalization.}
}
\value{
List of expressionSet objects storing the normalized fold changes, 
as well as row and column metadata. In each expressionSet \code{S}, the fold 
changes can be accessed by \code{Biobase::exprs(S)}. Protein expNames can be accessed 
by \code{featureNames(S)}. Isobaric labels and the corresponding 
concentrations are returned by \code{S$label} and \code{S$concentration}.
}
\description{
Normalize fold changes of TPP-CCR experiment to a reference 
column (usually that with the lowest concentration) to ensure that the 
transformation by \link{tppccrTransform} yields values between 0 and 1.
}
\examples{
data(hdacCCR_smallExample)
tppccrData <- tppccrImport(configTable=hdacCCR_config, data = hdacCCR_data)
tppccrNorm <- tppccrNormalize(data=tppccrData)
# Normalize to lowest concentration (in the first column):
tppccrNormToRef <- tppccrNormalizeToReference(data=tppccrNorm, refCol=1)
# Obtain results per replicate:
refTransf_replicate1 <- tppccrNormToRef$Panobinostat_1
head(Biobase::exprs(refTransf_replicate1))
# Perform transformation:
tppccrTransformed <- tppccrTransform(data=tppccrNormToRef)
# Obtain transformed measurements per replicate:
transf_replicate1 <- tppccrTransformed$Panobinostat_1
transf_replicate2 <- tppccrTransformed$Panobinostat_2
# Inspect transformed data in replicate 1:
effects_replicate1 <- Biobase::featureData(transf_replicate1)$compound_effect
newData_repl1 <- data.frame(Biobase::exprs(transf_replicate1), 
                              Type=effects_replicate1)[!is.na(effects_replicate1),]
                              
}
