
# Racialization in Testosterone Research: ERGMs
# Brandon L. Kramer 

# loading packages 
rm(list = ls())
knitr::opts_knit$set(root.dir = "C:/Users/bkram/CloudStation/git/racialization-in-testosterone-research")
for (pkg in c("tidyverse", "igraph", "googlesheets", "networkD3", 
              "ergm", "sna", "statnet", "btergm", "parallel")) {library(pkg, character.only = TRUE)}

# loading data 
setwd("C:/Users/bkram/CloudStation/git/data")
nodelist <- read_csv("racialization-nodelist.csv")
edgelist <- read_csv("racialization-edgelist.csv")

# converting to network 
el=as.matrix(edgelist)            
el[,1]=as.character(el[,1])  
el[,2]=as.character(el[,2])
tnet=network(edgelist,matrix.type="edgelist",directed=TRUE); rm(el) 
plot(tnet)

detectCores() # 12 cores

# network descriptives 
network.size(tnet) # 309
network.edgecount(tnet) #595
network.dyadcount(tnet) #95172

# applying nodal attributes to network 

tnet%v%"evidence" <- nodelist$evidence 
tnet%v%"year" <- nodelist$year
tnet%v%"groupcode" <- nodelist$groupcode
tnet%v%"rttcode" <- nodelist$rttcode
tnet%v%"gstimescited" <- nodelist$gstimescited

# BASELINE AND NODAL ATTRIBUTE TESTING 

# baseline (AIC: 6706    BIC: 6715)
t <- ergm(tnet ~ edges, 
          control=control.ergm(parallel=8, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t)

# evidence (AIC: 6166    BIC: 6185) 
t_ev <- ergm(tnet ~ edges + nodematch("evidence"), 
          control=control.ergm(parallel=8, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_ev)

# year (AIC: 6643    BIC: 6662)
t_yr <- ergm(tnet ~ edges + nodecov("year"),
             control=control.ergm(parallel=8, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_yr)

# population grouping (AIC: 5899    BIC: 5947)  
t_pop <- ergm(tnet ~ edges + nodematch("groupcode"), 
             control=control.ergm(parallel=8, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_pop)
# nodefactor works better here 

# racial testosterone tested (AIC: 6651    BIC: 6670)    
t_rtt <- ergm(tnet ~ edges + nodematch("rttcode"), 
              control=control.ergm(parallel=8, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_rtt)

# times cited  (AIC: 6707    BIC: 6726) # ns effect   
t_rtt <- ergm(tnet ~ edges + nodecov("gstimescited"), 
              control=control.ergm(parallel=8, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_rtt)




# PRELIMINARY STRUCTURAL FACTORS 

















# dont test slippage (must finish coding)