% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tpptrPlotSplines.R
\name{tpptrPlotSplines}
\alias{tpptrPlotSplines}
\title{Plot spline fits per protein}
\usage{
tpptrPlotSplines(
  data,
  factorsH1 = NULL,
  factorsH0 = NULL,
  fittedModels,
  testResults,
  resultPath = NULL,
  individual = TRUE,
  overview = FALSE,
  returnPlots = FALSE,
  control = list(nCores = "max", maxRank = 500, highlightBelow = 0.05),
  maxRank = NULL,
  highlightBelow = NULL,
  plotIndividual = NULL,
  plotAlphabetical = NULL
)
}
\arguments{
\item{data}{long table of proteins measurements that were used for spline fitting.}

\item{factorsH1}{DEPRECATED}

\item{factorsH0}{DEPRECATED}

\item{fittedModels}{long table of fitted models. 
Output of \code{\link{tpptrFitSplines}}.}

\item{testResults}{long table of p-values per protein. 
Output of \code{\link{tpptrFTest}}.}

\item{resultPath}{an optional character vector with the name of the path where the plots should be saved.}

\item{individual}{logical. Export each plot to individual files?}

\item{overview}{logical. Generate summary pdfs?}

\item{returnPlots}{logical. Should the ggplot objects be returned as well?}

\item{control}{a list of general settings.}

\item{maxRank}{DEPRECATED}

\item{highlightBelow}{DEPRECATED}

\item{plotIndividual}{DEPRECATED}

\item{plotAlphabetical}{DEPRECATED
Contains the following fields:
\itemize{ 
\item{\code{nCores}: number of CPUs for parallel production of plots per 
protein if \code{individual = TRUE} (default: "max")}
\item{\code{maxRank}: how many of the top hits should 
be plotted if \code{overview = TRUE} (default: 500)} 
 \item{\code{highlightBelow}: maximum adjusted p-value 
for which a protein is highlighted by a different background color if 
\code{overview = TRUE} (default: 0.05)}
}}
}
\value{
None
}
\description{
Plot spline fits per protein
}
\details{
Plots of the natural spline fits will be stored in a subfolder with 
name \code{Spline_Fits} at the location specified by \code{resultPath}.

Exporting each plot to individual files (individual = TRUE) can 
cost runtime and the resulting files can be tedious to browse. 
If you just want to browse the results, use \code{overview = TRUE} 
instead.

If \code{overview = TRUE}, two summary PDFs are created that enable quick 
browsing through all results.  They contain the plots in alphacetical order 
(\code{splineFit_alphabetical.pdf}), or ranked by p-values 
(\code{splineFit_top_xx.pdf}, where xx is the maximum rank defined by 
\code{overviewSettings$maxRank}).
}
\examples{
data(hdacTR_smallExample)
tpptrData <- tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
tidyData <- tpptrTidyUpESets(tpptrData)
splineFits <- tpptrFitSplines(data = tidyData, nCores = 1, splineDF = 4:5,
                              factorsH1 = "condition", returnModels = TRUE)
testResults <- tpptrFTest(fittedModels = splineFits, doPlot = FALSE)
tpptrPlotSplines(data = tidyData, fittedModels = splineFits,
                 individual = FALSE,
                 testResults = testResults, resultPath = getwd())

}
\seealso{
\code{\link{ns}, \link{AICc}, 
\link{tpptrFitSplines}, \link{tpptrFTest}}
}
