\name{curatedBreastDataExprSetList}
\alias{curatedBreastDataExprSetList}
\docType{data}
\title{
curatedBreastDataExprSetList
}
\description{
A list of ExpressionSet objects, one for each curated study, containing study-specific gene expression and phenotype data. FeatureNames are gene symbols. Data is already quantile normalized according to standard protocols for 1 and 2-channel arrays, depending on the platform used for this study.
}
\usage{data("curatedBreastDataExprSetList")}

\format{
A list, with each index containing an ExpressionSet object from a specific study, and potentially a specific batch. 
}
\value{
No return value as this is not a function but rather a data object.
}
\details{
Batches from studies are treated as individual datasets, as the signal can differ between batches. Thus, an expression object named using a GSE study number followed by an underscore means this ExpressionSet contains samples either from a distinct platform (and the study used >1 platforms), or from a distinct batch or tissue site.  An "all" tag means that there were no batches for this study. Raw data files downloaded from GEO oftentimes have clear batch/site information appended to sample names; this was often the source of batch identification and how the package developer chose to create the batch name string.
}

\references{
Planey, Butte. Database integration of 4923 publicly-available samples of breast cancer molecular and clinical data. AMIA Joint Summits Translational Science Proceedings. (2003) PMC3814460
}
\examples{
data(curatedBreastDataExprSetList)
#what are all the names of the studies?
names(curatedBreastDataExprSetList)
#what is the dimension of the gene 
#expression matrix for study GSE17705 from the JBI 
#(as opposed to MDACC) site?
dim(exprs(curatedBreastDataExprSetList$study_17705_GPL96_JBI_Tissue_BC_Tamoxifen))
}
\keyword{datasets}
