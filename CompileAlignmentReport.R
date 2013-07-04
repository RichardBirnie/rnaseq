library(knitr)
library(markdown)

#grab command line arguments
args <- commandArgs(trailingOnly=TRUE)

#calls the Rmarkdown AlignmentReport.Rmd and compiles it into an html document with 
#figures summarising the QC metrics for the current batch of samples
fp <- file.path(args[3])
knit(file.path(args[2]), output=fp)
htmlFile <- sub('.md', '.html', fp)
markdownToHTML(file=fp, output=htmlFile)
file.remove(fp)