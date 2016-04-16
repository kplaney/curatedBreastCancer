# curatedBreastCancer
#NOTE: as of 2016, I have an updated script for the gene expression processing. Please see the breastProcessAndGeneFeatures_script.R in the Concide repo for accessing this new updated processing function. The same funcitons are indeed on the curatedBreastData repo, but not in the Bioconductor package, as I am in the process of updating this package submission. I plan to remove the processing functions from the curatedBreastData package and just point users to the Coincide package, as Bioconductor encourages developers to keep code out of database packages exactly to avoid this issue of trying to update code quicker.

The only difference is a very small bug that crashes on one dataset (so it won't produce incorrect output, it just stops) as it was a minor ,drop=FALSE indexing issue.

The GSE.... R files are mainly very early records of how I processed these files, but are not intended to be directly re-run.  However, I do have notes on the different datasets if you need more specifics than is provided in the corresponding publications, https://www.ncbi.nlm.nih.gov/pubmed/24303324 and  https://www.ncbi.nlm.nih.gov/pubmed/26961683 (the latter has some more details in the supplementary methods section).

If anyone is ever interested in adding more datasets, we can perhaps collaborate to introduce a second database package - please email me at catherine.planey@gmail.com. I am not currently in academia, but will do my best to respond promptly.
