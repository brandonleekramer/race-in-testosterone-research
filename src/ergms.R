
# Racialization in Testosterone Research: ERGMs
# Brandon L. Kramer 

# loading packages 
rm(list = ls())
#knitr::opts_knit$set(root.dir = "C:/Users/bkram/CloudStation/git/racialization-in-testosterone-research")
setwd("~/racialization-data")

for (pkg in c("tidyverse", "igraph", "networkD3", 
              "ergm", "sna", "statnet", "parallel")) {library(pkg, character.only = TRUE)}

# loading data 
#setwd("C:/Users/bkram/CloudStation/Biomedical MultipliciTs/Racial Differences in T/racialization-in-testosterone-research")
nodelist <- read_csv("racialization-nodelist.csv")
edgelist <- read_csv("racialization-edgelist.csv")

# converting to network 
el=as.matrix(edgelist)            
el[,1]=as.character(el[,1])  
el[,2]=as.character(el[,2])
tnet=network(edgelist,matrix.type="edgelist",directed=TRUE); rm(el) 
plot(tnet)

# setting the number of cores for parallelizing 
#(num_cores <- detectCores() - 1) 
num_cores <- 9

# network descriptives 
network.size(tnet) # 309
network.edgecount(tnet) #595
network.dyadcount(tnet) #95172

# applying nodal attributes to the network 

# article attributes 
tnet%v%"year" <- nodelist$year
tnet%v%"evidence" <- nodelist$evidence 
tnet%v%"outcome" <- nodelist$outcomecode
tnet%v%"domain" <- nodelist$domaincode
tnet%v%"journal" <- nodelist$journalcode
tnet%v%"gstimescited" <- nodelist$gstimescited
tnet%v%"region" <- nodelist$region
tnet%v%"drd" <- nodelist$datarecyclingbinary
tnet%v%"drm" <- nodelist$datarecyclingmulti
# population difference tests 
tnet%v%"group" <- nodelist$groupcode 
tnet%v%"testedmen" <- nodelist$testedmen 
tnet%v%"rttcode" <- nodelist$rttcode
tnet%v%"wbcomp" <- nodelist$wbcompcode
tnet%v%"bacomp" <- nodelist$bacompcode
tnet%v%"wacomp" <- nodelist$wacompcode
tnet%v%"hisp" <- nodelist$hisp
tnet%v%"wvsnw" <- nodelist$wvsnw
# control tests 
tnet%v%"timediff" <- nodelist$timediff
tnet%v%"agemeasured" <- nodelist$agemeasured	
tnet%v%"agecontrol" <- nodelist$agecontrol	
tnet%v%"agediff" <- nodelist$agediff	
tnet%v%"bmimeasured" <- nodelist$bmimeasured	
tnet%v%"bmicontrolled" <- nodelist$bmicontrolled	
tnet%v%"bmidiff" <- nodelist$bmidiff	
tnet%v%"illnesscontrolled" <- nodelist$illnesscontrolled	
tnet%v%"sescontrolled" <- nodelist$sescontrolled	
tnet%v%"medscontrolled" <- nodelist$medscontrolled

nodelist %>% 
  count(region)


################################################### BASELINE AND STUDY-LEVEL ATTRIBUTES 

# baseline (AIC: 6706    BIC: 6715)
t <- ergm(tnet ~ edges, 
          control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t)
gof_t <- gof(t)


# year 
t_yr <- ergm(tnet ~ edges + nodecov("year"),
             control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_yr)

# evidence 
t_ev <- ergm(tnet ~ edges + nodefactor("evidence"), 
          control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_ev)

# outcome 
t_oc <- ergm(tnet ~ edges + nodefactor("outcome"), 
             control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_oc)

# domain
t_dom <- ergm(tnet ~ edges + nodefactor("domain"), 
                  control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom)

# journal
t_jnl <- ergm(tnet ~ edges + nodefactor("journal"), 
              control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_jnl)

# times cited     
t_tc <- ergm(tnet ~ edges + nodecov("gstimescited"), 
                control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_tc)

# region
t_rgn <- ergm(tnet ~ edges + nodefactor("region"), 
              control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_rgn)

# data recycling (dummy)
t_drd <- ergm(tnet ~ edges + nodefactor("drd"), 
              control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_drd)

# data recycling (multi) 
t_drm <- ergm(tnet ~ edges + nodefactor("drm"), 
              control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_drm)
plot(gof(t_drm))
# domain, region, and data recycling had biggest initial impact 

# domain + drd
t_dom_drd <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drd"), 
              control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drd)

# domain + drm
t_dom_drm <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm"), 
                  control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm)

# domain,region
t_dom_rgn <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("region"), 
                      control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_rgn)

# region,drm
t_rgn_drm <- ergm(tnet ~ edges + nodefactor("region") + nodefactor("drm"), 
                  control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_rgn_drm)

# region,drm
t_rgn_drd <- ergm(tnet ~ edges + nodefactor("region") + nodefactor("drd"), 
                  control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_rgn_drd)

# domain,drd,region = best combination so far 
t_dom_drm_rgn <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodefactor("region"), 
                      control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn)

# double-checking to see if nsig covariates matter again 

# domain,drd,region,year
t_dom_drm_rgn_yr <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + 
                        nodefactor("region") + nodecov("year"), 
                      control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_yr)

# domain,drd,region,outcome
t_dom_drm_rgn_oc <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + 
                           nodefactor("region") + nodefactor("outcome"), 
                         control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_oc)

# domain,drd,region,journal
t_dom_drm_rgn_jnl <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + 
                           nodefactor("region") + nodefactor("journal"), 
                         control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_jnl)

####################################################################################################################
# domain,drd,region = best combination so far 
t_dom_drm <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm"), 
                      control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm)
gof_dom_drm <- gof(t_dom_drm) 
plot(gof_dom_drm); gof_dom_drm


####################################################################################################################

################################################### SPECIFIC FORMS OF RACIAL DIFFERENCE TESTING 

# tested_men  
t_dom_drm_rgn_men <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodefactor("region") + 
                           nodefactor("testedmen"),
                         control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_men)

# population grouping 
t_dom_drm_rgn_gp <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodefactor("region") + 
                      nodefactor("group"),
                      control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_gp)

# racial testosterone tested     
t_dom_drm_rgn_gp_rtt <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodefactor("region") + 
                           nodefactor("group") + nodematch("rttcode"),
                         control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_gp_rtt)

# white-black differences 
t_dom_drm_rgn_gp_wb <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodefactor("region") + 
                               nodefactor("group") + nodematch("wbcomp"),
                             control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_gp_wb)

# white-asian differences
t_dom_drm_rgn_gp_wa <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodefactor("region") + 
                              nodefactor("group") + nodematch("wacomp"),
                            control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_gp_wa)

# black-asian differences
t_dom_drm_rgn_gp_ba <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodefactor("region") + 
                              nodefactor("group") + nodematch("bacomp"),
                            control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_gp_ba)

# black-asian differences
t_dom_drm_rgn_gp_hsp <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodefactor("region") + 
                              nodefactor("group") + nodematch("hisp"),
                            control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_gp_hsp)

# western-nonwestern comps 
t_dom_drm_rgn_gp_wvsnw <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodefactor("region") + 
                               nodefactor("group") + nodematch("wvsnw"),
                             control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_gp_wvsnw)


################################################### PRELIMINARY STRUCTURAL FACTORS 

# baseline (AIC: 6706    BIC: 6715)
t <- ergm(tnet ~ edges, 
          control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t)

# controlling for isolates
t_iso <- ergm(tnet ~ edges + isolates, 
          control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_iso)

# controlling for isolates
summary(tnet~edges + triangle)
t_tri <- ergm(tnet ~ edges + triangle, 
              control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_tri)

# looking at degree distributions 
summary(tnet ~ edges + idegree(1:10))

# degree = 1
t_id1 <- ergm(tnet ~ edges + idegree(1), 
              control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_id1)

# degree = 1:2
t_id1_2 <- ergm(tnet ~ edges + idegree(1:2), 
              control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_id1_2)

# degree = 1 + iso
t_id1_iso <- ergm(tnet ~ edges + idegree(1) + isolates, 
              control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_id1_iso)

# adding in isolates + degree term  
t_dom_drm_rgn_iso_dg1 <- ergm(tnet ~ edges + nodematch("domain") + nodefactor("drm") + nodefactor("region") +
                                idegree(1) + isolates, 
                              control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_iso_dg1)
gof_dom_drm_rgn_iso_dg1 <- gof(t_dom_drm_rgn_iso_dg1) 
plot(gof_dom_drm_rgn_iso_dg1); gof_dom_drm_rgn_iso_dg1

# comparison with previous best model = slightly better but GOF still bad 
gof_dom_drm_rgn <- gof(t_dom_drm_rgn) 
plot(gof_dom_drm_rgn); gof_dom_drm_rgn

# gwesp(0.25) = model degeneracy
t_gw25 <- ergm(tnet ~ edges + gwesp(0.25),
                  control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_gw25)

# gwesp(0.50) = model degeneracy
t_gw50 <- ergm(tnet ~ edges + gwesp(0.50),
               control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_gw50)

# gwesp(0.10) = model degeneracy
t_gw10 <- ergm(tnet ~ edges + gwesp(0.10),
               control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_gw10)

# gwesp(0.01) = model degeneracy
t_gw01 <- ergm(tnet ~ edges + gwesp(0.01),
               control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_gw01)



# must finish coding slippage variable ######################


references 

https://statnet.org/trac/raw-attachment/wiki/Sunbelt2016/ergm_tutorial.html

http://statnet.csde.washington.edu/workshops/EUSN/current/ergm/ergm_tutorial.html




