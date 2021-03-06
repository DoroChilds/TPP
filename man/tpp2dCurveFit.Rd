% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tpp2dCurveFit.R
\name{tpp2dCurveFit}
\alias{tpp2dCurveFit}
\title{Run TPP-CCR analysis for 2D-TPP experiment}
\usage{
tpp2dCurveFit(
  configFile = NULL,
  data,
  nCores = 1,
  naStrs = NULL,
  fcStr = NULL,
  idVar = NULL,
  nonZeroCols = NULL,
  r2Cutoff = 0.8,
  fcCutoff = 1.5,
  slopeBounds = c(1, 50),
  fcTolerance = 0.1
)
}
\arguments{
\item{configFile}{DEPCRECATED}

\item{data}{data frame that contains the data of the 2D-TPP 
experiment for each temperature.}

\item{nCores}{numeric value stating how many cores are to be used for computation}

\item{naStrs}{DEPCRECATED}

\item{fcStr}{DEPCRECATED}

\item{idVar}{DEPCRECATED}

\item{nonZeroCols}{DEPCRECATED}

\item{r2Cutoff}{Quality criterion on dose response curve fit.}

\item{fcCutoff}{Cutoff for highest compound concentration fold change.}

\item{slopeBounds}{Bounds on the slope parameter for dose response curve 
fitting.}

\item{fcTolerance}{tolerance for the fcCutoff parameter. See details.}
}
\value{
A data frames in which the fit results are stored row-wise for each
  protein.
}
\description{
Performs analysis of a TPP-CCR experiment by invoking the routine 
  for TPP-CCR curve fitting for each temperature of the sample.
}
\examples{
# Preparation:
data(panobinostat_2DTPP_smallExample)

# Import data:
datIn <- tpp2dImport(configTable = panobinostat_2DTPP_config,
                      data = panobinostat_2DTPP_data,
                      idVar = "representative",
                      addCol = "clustername",
                      intensityStr = "sumionarea_protein_",
                      nonZeroCols = "qusm")

# Compute fold changes:
datFC <- tpp2dComputeFoldChanges(data = datIn)

# Perform median normalization:
datNorm <- tpp2dNormalize(data = datFC)

# View updated attributes. Now contain field 'fcStrNorm' indicating prefix
# of the fold change columns after normalization.
attr(datNorm, "importSettings")["fcStrNorm"]

# Perform dose response curve fitting and pEC50 calculation:
datFit <- tpp2dCurveFit(data = datNorm)
  
}
