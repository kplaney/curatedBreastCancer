\name{processExpressionSet}
\alias{processExpressionSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Post-process a normalized assayData in an ExpressionSet object
}
\description{
This function Post-processes a normalized assayData in an ExpressionSet object. Is is assumed the assay data is already baseline normalized (for example, for microarray data, this could mean quantile normalized and then logged.)
}
\usage{
processExpressionSet(exprSet, outputFileDirectory = "./",numTopVarGenes,
minVarPercentile, maxVarPercentile = 1, minVar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprSet}{
expressionSet S4 object with expression (assay) data, featureData and phenoData.
}
  \item{outputFileDirectory}{
Output file directory for messages that print status of post-processing the ExpressionSet.  
}
  \item{minVarPercentile}{
Minimum variance percentile. Must be provided in conjunction with maxVarPercentile to use percentiles to threshold genes.
}
  \item{maxVarPercentile}{
Maximum variance percentile. Defaul is 1, i.e. 1\%. Must be provided in conjunction with minVarPercentile to use percentiles to threshold genes.
}

  \item{minVar}{
If maxVar is provided, as opposed to minVarPercentile and maxVarPercentile, genes are removed that are below a certain variance magnitude. This is helpful before running certain algorithms, such as the popular Combat batch normalization technique, that can throw errors if genes with extremely low variances are in the data matrix. May be used in conjunction with maxVar or in isolation.
}
  \item{numTopVarGenes}{
A numeric value indicating the number of genes (features) to select; the function will only take this number of genes that have the highest variance across all genes.
}
}
\details{
This function performs several post-processing tasks: filtering out genes and samples with high NA rates, imputing missing values, collapsing duplicated features/genes to make a unique feature list, removing any samples for which there is already a sample with the sample patient ID (duplicated samples), and filtering genes by variance.  This function is a wrapper for the functions: filterAndImputeSamples(), collapseDupProbes(), removeDuplicatedPatients(), and filterGenesByVariance(). It is is run after initial dataset normalization, such as quantile normalization on microarray datasets.
}
\value{
A post-processed S4 expressionSet.  Tests are run to confirm the final S4 object is a valid ExpressionObject before it is returned.
}

\author{
Katie Planey <katie.planey@gmail.com>
}

\examples{
#load up our datasets
data(curatedBreastDataExprSetList);

#just perform on one dataset as an example, GSE9893. 
#This dataset does have NA values, so
#you'll see the impute.knn progress printed to the screen.
#also take only genes that fall in 
#the variance percentiles between .75 and 1 
#(i.e. top 75th percentile genes by variance.)

post_procExprSet <- processExpressionSet(exprSet=
curatedBreastDataExprSetList[[5]], 
outputFileDirectory = "./",
minVarPercentile=.75, maxVarPercentile = 1)

}