# curatedBreastCancer

I released a slightly updated processing code in Bioconductor 3.4 in the fall of 2016. The main difference is a very small bug that crashes on one dataset (so it won't produce incorrect output, it just stops) as it was a minor ,drop=FALSE indexing issue. I am in the process of fixing the DFS_months_or_MIN months variables for GSE16446, as it is in days, not months (you can just divide the days by 28 to get months, but this  is only needed for this one specific study.)

The GSE.... R files are mainly very early records of how I processed these files, but are not intended to be directly re-run.  However, I do have notes on the different datasets if you need more specifics than is provided in the corresponding publications, https://www.ncbi.nlm.nih.gov/pubmed/24303324 and  https://www.ncbi.nlm.nih.gov/pubmed/26961683 (the latter has some more details in the supplementary methods section).

If anyone is ever interested in adding more datasets, we can perhaps collaborate to introduce a second database package - please email me at catherine.planey@gmail.com. I am not currently in academia, but will do my best to respond promptly.
