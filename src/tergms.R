# Racialization in Testosterone Research: ERGMs
# Brandon L. Kramer 

# loading packages 
rm(list = ls())
knitr::opts_knit$set(root.dir = "C:/Users/bkram/CloudStation/Biomedical MultipliciTs/Racial Differences in T/racialization-in-testosterone-research/data")
for (pkg in c("tidyverse", "statnet", "htmlwidgets", "latticeExtra", "parallel")) {library(pkg, character.only = TRUE)}

# loading data 
setwd("C:/Users/bkram/CloudStation/Biomedical MultipliciTs/Racial Differences in T/racialization-in-testosterone-research/data")
nodelist <- read_csv("racialization-nodelist.csv")
edgelist <- read_csv("racialization-edgelist.csv")

# converting to network 
el=as.matrix(edgelist)            
el[,1]=as.character(el[,1])  
el[,2]=as.character(el[,2])
tnet=network(edgelist,matrix.type="edgelist",directed=TRUE); rm(el) 
plot(tnet)

detectCores() # 12 cores




References 

#http://statnet.org/tergm_tutorial.html
#http://repository.essex.ac.uk/25129/1/v83i06.pdf
#https://eprints.gla.ac.uk/121525/1/
#https://arxiv.org/pdf/1506.06696.pdf
#https://arxiv.org/pdf/1512.06141.pdf
#https://static1.squarespace.com/static/5951237c725e258922c1a8b5/t/5ab918bcf950b759a809d0fb/1522079932793/BWC_egoTERGM_Manuscript_FINAL03132018.pdf
#https://eprints.gla.ac.uk/145136/7/145136.pdf (uses xergm for GOF)
# https://www.cambridge.org/core/journals/political-analysis/article/predicting-network-events-to-assess-goodness-of-fit-of-relational-event-models/DB71F72569196DD4F54C115353734A19
# https://sites.cs.ucsb.edu/~xyan/papers/REM%20Theory.pdf (REM vs TERGMs)
# https://cran.r-project.org/web/packages/rem/rem.pdf
# https://journals.sagepub.com/doi/abs/10.1177/0018726719895322
# https://www.youtube.com/watch?v=XwcoWKj-URo 
# https://ethz.ch/content/dam/ethz/special-interest/gess/social-networks-dam/documents/goldfish_EUSN2017/Butts_2008.pdf




