
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
edgelist <- edgelist %>% select(-WorkSpace)

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
nodelist$indegree <- degree(tnet, cmode = "indegree")
nodelist$outdegree <- degree(tnet, cmode = "outdegree")
tnet%v%"medscontrolled" <- nodelist$indegree
tnet%v%"medscontrolled" <- nodelist$outdegree

################################################### BASELINE AND STUDY-LEVEL ATTRIBUTES 

# baseline (AIC: 6706    BIC: 6715)
t <- ergm(tnet ~ edges, 
          control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t)

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
t_dom_drm_rgn <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodematch("region"), 
                      control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK", eval.loglik = TRUE))
summary(t_dom_drm_rgn)
gof_dom_drm_rgn <- gof(t_dom_drm_rgn) 
plot(gof_dom_drm_rgn); gof_dom_drm_rgn

####################################################################################################################

t_dom_drd_rgn <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("drd") + nodematch("region"), 
                      control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drd_rgn)
gof_dom_drd_rgn <- gof(t_dom_drd_rgn) 
plot(gof_dom_drd_rgn); gof_dom_drd_rgn

################################################### SPECIFIC FORMS OF RACIAL DIFFERENCE TESTING 

# tested_men  
t_dom_drm_rgn_men <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodematch("region") + 
                            nodematch("testedmen"),
                          control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_men)

# racial testosterone tested = sig     
t_dom_drm_rgn_men_rtt <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodematch("region") + 
                               nodematch("testedmen") + nodematch("rttcode"),
                             control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_men_rtt)

# western-nonwestern comps = sig 
t_dom_drm_rgn_men_wvnw <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodematch("region") + 
                                 nodematch("testedmen") + nodematch("wvsnw"),
                               control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_men_wvnw)

# western-nonwestern comps = sig 
t_dom_drm_rgn_men_rtt_wvnw <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodematch("region") + 
                                 nodematch("testedmen") + nodematch("rttcode") + nodematch("wvsnw"),
                               control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dom_drm_rgn_men_rtt_wvnw)

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

# degree = 1 + iso
t_id1_iso <- ergm(tnet ~ edges + idegree(1) + isolates, 
              control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_id1_iso)

# gwesp(0.60) = model degeneracy
t_gw60 <- ergm(tnet ~ edges + gwesp(0.60),
               control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_gw60)

# gwesp(0.60) = model degeneracy
t_gw60_iso <- ergm(tnet ~ edges + gwesp(0.60) + isolates,
                   control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_gw60_iso)

summary(tnet ~ edges + esp(1:20))
summary(tnet ~ edges + dsp(1:20))
summary(tnet ~ edges + nsp(1:20))

# esp(0.10) = didnt converge
t_esp10 <- ergm(tnet ~ edges + esp(0.10),
               control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_esp10)

# esp(0.25) = didnt converge
t_esp25 <- ergm(tnet ~ edges + esp(0.25),
                control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_esp25)

# esp(0.5) = didnt converge
t_esp50 <- ergm(tnet ~ edges + esp(0.50),
                control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_esp50)

# esp(0.75) = didnt converge
t_esp75 <- ergm(tnet ~ edges + esp(0.75),
                control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_esp75)


# dsp(0.10) = didnt converge
t_dsp10 <- ergm(tnet ~ edges + dsp(0.10),
                control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dsp10)

# dsp(0.25) = didnt converge
t_dsp25 <- ergm(tnet ~ edges + dsp(0.25),
                control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dsp25)

# dsp(0.5) = didnt converge
t_dsp50 <- ergm(tnet ~ edges + dsp(0.50),
                control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dsp50)

# dsp(0.75) = didnt converge
t_dsp75 <- ergm(tnet ~ edges + dsp(0.75),
                control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK"))
summary(t_dsp75)

# read over the statnet proposal so i'm going to try that suggestion

# gwesp0 + dm/drm/rgn = converged 
t_dm_drm_rgn_gwe0 <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodematch("region") + 
                  gwesp(0, fixed=TRUE),
                control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK", seed=123))
summary(t_dm_drm_rgn_gwe0)

# gwesp10 + dm/drm/rgn =  
t_dm_drm_rgn_gwe10 <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodematch("region") + 
                            gwesp(0.10, fixed=TRUE),
                          control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK", seed=123))
summary(t_dm_drm_rgn_gwe10)

# gwesp20 + dm/drm/rgn =  
t_dm_drm_rgn_gwe20 <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodematch("region") + 
                             gwesp(0.20, fixed=TRUE),
                           control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK", seed=123))
summary(t_dm_drm_rgn_gwe20)

# gwesp30 + dm/drm/rgn =  
t_dm_drm_rgn_gwe30 <- ergm(tnet ~ edges + nodefactor("domain") + nodefactor("drm") + nodematch("region") + 
                             gwesp(0.30, fixed=TRUE),
                           control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, parallel.type="PSOCK", seed=123))
summary(t_dm_drm_rgn_gwe30)

c(t_dm_drm_rgn_gwe0$mle.lik, t_dm_drm_rgn_gwe10$mle.lik, t_dm_drm_rgn_gwe20$mle.lik, t_dm_drm_rgn_gwe30$mle.lik)
gof_t_dm_drm_rgn_gwe20 <- gof(t_dm_drm_rgn_gwe20)
plot(gof_t_dm_drm_rgn_gwe20); gof_t_dm_drm_rgn_gwe20

# going with the gwesp(0.2) model but the GOF still is not great but the ideg + odeg do not fit well yet 

# let's bring everything together that we know has an impact and makes sense theoretically 
# domain and region def matters
# data recycling is going to be a 1/0 variable (using nodematch)
# no rtt code but will include testing men and western v non-western samples 

t_dom_rgn_drd_men_wvnw_gwe20 <- ergm(tnet ~ edges + 
                                       nodefactor("domain") + 
                                       nodematch("drd") + 
                                       nodematch("region") + 
                                       nodematch("testedmen") + 
                                       nodematch("wvsnw") +
                                       gwesp(0.20, fixed=TRUE),
                                     control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                          parallel.type="PSOCK", seed=123))
summary(t_dom_rgn_drd_men_wvnw_gwe20)

# converged and interpretable 
t_dom_rgn_men_wvnw_gwe20 <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                       nodematch("testedmen") + nodematch("wvsnw") +
                                       gwesp(0.20, fixed=TRUE),
                                     control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                          parallel.type="PSOCK", seed=123, eval.loglik = TRUE))
summary(t_dom_rgn_men_wvnw_gwe20)
gof_dom_rgn_men_wvnw_gwe20 <- gof(t_dom_rgn_men_wvnw_gwe20)
plot(gof_dom_rgn_men_wvnw_gwe20); gof_dom_rgn_men_wvnw_gwe20
bgof_dom_rgn_men_wvnw_gwe20 <- btergm::gof(t_dom_rgn_men_wvnw_gwe20)
plot(bgof_dom_rgn_men_wvnw_gwe20); bgof_dom_rgn_men_wvnw_gwe20

summary(tnet ~ edges + idegree(1:20))
summary(tnet ~ edges + odegree(1:20))

t_dom_rgn_men_wvnw_gwe20_gwi2_iso <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                      nodematch("testedmen") + nodematch("wvsnw") +
                                      gwidegree(decay=2, fixed=TRUE) + isolates + gwesp(0.20, fixed=TRUE),
                                      control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                           parallel.type="PSOCK", seed=123))
summary(t_dom_rgn_men_wvnw_gwe20_gwi2_iso)

t_dom_rgn_men_wvnw_gwe20_gwi2 <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                            nodematch("testedmen") + nodematch("wvsnw") +
                                            gwidegree(decay=2, fixed=TRUE) + gwesp(0.20, fixed=TRUE),
                                          control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                               parallel.type="PSOCK", seed=123))
summary(t_dom_rgn_men_wvnw_gwe20_gwi2)

t_dom_rgn_men_wvnw_gwe20_in1 <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                        nodematch("testedmen") + nodematch("wvsnw") +
                                        idegree(1) + gwesp(0.20, fixed=TRUE),
                                      control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                           parallel.type="PSOCK", seed=123))
summary(t_dom_rgn_men_wvnw_gwe20_in1)

t_dom_rgn_men_wvnw_gwe20_in1_iso <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                       nodematch("testedmen") + nodematch("wvsnw") +
                                       idegree(1) + isolates + gwesp(0.20, fixed=TRUE),
                                     control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                          parallel.type="PSOCK", seed=123))
summary(t_dom_rgn_men_wvnw_gwe20_in1_iso)

t_dom_rgn_men_wvnw_gwe20_iso <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                           nodematch("testedmen") + nodematch("wvsnw") +
                                           isolates + gwesp(0.20, fixed=TRUE),
                                         control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                              parallel.type="PSOCK", seed=123))
summary(t_dom_rgn_men_wvnw_gwe20_iso)

t_dom_rgn_wvnw_gwe20_in1men <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                       #nodematch("testedmen") + 
                                       nodematch("wvsnw") +
                                       idegree(1, by="testedmen",homophily=TRUE) + gwesp(0.20, fixed=TRUE),
                                       control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                            parallel.type="PSOCK", seed=123))
summary(t_dom_rgn_wvnw_gwe20_in1men)

t_dom_rgn_wvnw_gwe20_in1wvsnw <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                      nodematch("testedmen") + 
                                      #nodematch("wvsnw") +
                                      idegree(1, by="wvsnw",homophily=TRUE) + gwesp(0.20, fixed=TRUE),
                                    control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                         parallel.type="PSOCK", seed=123, eval.loglik = TRUE))
summary(t_dom_rgn_wvnw_gwe20_in1wvsnw)
anova(t_dom_drm_rgn, t_dom_rgn_wvnw_gwe20_in1wvsnw)


t_dom_rgn_wvnw_gwe20_i12o14 <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                        nodematch("testedmen") + nodematch("wvsnw") +
                                        idegree(1:2) + odegree(1:4) +gwesp(0.20, fixed=TRUE, cutoff=5),
                                      control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                           parallel.type="PSOCK", seed=123))
summary(t_dom_rgn_wvnw_gwe20_i12o14)

t_dom_rgn_wvnw_gwe20_i12o12 <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                      nodematch("testedmen") + nodematch("wvsnw") +
                                      idegree(1:2) + odegree(1:2) +gwesp(0.20, fixed=TRUE, cutoff=5),
                                    control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                         parallel.type="PSOCK", seed=123))
summary(t_dom_rgn_wvnw_gwe20_i12o12)
gof_dom_rgn_wvnw_gwe20_i12o12 <- btergm::gof(t_dom_rgn_wvnw_gwe20_i12o12)
plot(gof_dom_rgn_wvnw_gwe20_i12o12); gof_dom_rgn_wvnw_gwe20_i12o12

# lowered the esp term and it worked! #################################################################
t_dom_rgn_wvnw_gwe10_i12o12 <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                      nodematch("testedmen") + nodematch("wvsnw") +
                                      idegree(1:2) + odegree(1:2) +gwesp(0.10, fixed=TRUE, cutoff=5),
                                    control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                         parallel.type="PSOCK", seed=123))
summary(t_dom_rgn_wvnw_gwe10_i12o12)
gof_t_dom_rgn_wvnw_gwe10_i12o12 <- btergm::gof(t_dom_rgn_wvnw_gwe10_i12o12)
plot(gof_t_dom_rgn_wvnw_gwe10_i12o12); gof_t_dom_rgn_wvnw_gwe10_i12o12
#######################################################################################################

# trying it again - didnt work so 0.1 it is 
t_dom_rgn_wvnw_gwe05_i12o12 <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                      nodematch("testedmen") + nodematch("wvsnw") +
                                      idegree(1:2) + odegree(1:2) +gwesp(0.05, fixed=TRUE, cutoff=5),
                                    control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                         parallel.type="PSOCK", seed=123))
summary(t_dom_rgn_wvnw_gwe05_i12o12)
gof_dom_rgn_wvnw_gwe05_i12o12 <- btergm::gof(t_dom_rgn_wvnw_gwe05_i12o12)
plot(gof_dom_rgn_wvnw_gwe05_i12o12); gof_dom_rgn_wvnw_gwe05_i12o12

# going to adjust the cutoff now 
t_dom_rgn_wvnw_gwe102_i12o12 <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                      nodematch("testedmen") + nodematch("wvsnw") +
                                      idegree(1:2) + odegree(1:2) +gwesp(0.10, fixed=TRUE, cutoff=2),
                                    control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                         parallel.type="PSOCK", seed=123))
summary(t_dom_rgn_wvnw_gwe102_i12o12)
gof_dom_rgn_wvnw_gwe102_i12o12 <- btergm::gof(t_dom_rgn_wvnw_gwe102_i12o12)
plot(gof_dom_rgn_wvnw_gwe102_i12o12); gof_dom_rgn_wvnw_gwe102_i12o12
# no affect 

t_dom_rgn_wvnw_gwe10_gwi10_gwo10 <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                      nodematch("testedmen") + nodematch("wvsnw") +
                                      gwidegree(0.10, fixed=TRUE, cutoff=4) + 
                                      gwodegree(0.10, fixed=TRUE, cutoff=10) + 
                                      gwesp(0.10, fixed=TRUE, cutoff=5),
                                    control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                         parallel.type="PSOCK", seed=123))
summary(t_dom_rgn_wvnw_gwe10_gwi10_gwo10)
gof_dom_rgn_wvnw_gwe10_gwi10_gwo10 <- btergm::gof(t_dom_rgn_wvnw_gwe10_gwi10_gwo10)
plot(gof_dom_rgn_wvnw_gwe10_gwi10_gwo10); gof_dom_rgn_wvnw_gwe10_gwi10_gwo10
#this didnt work


# best available model (ROC 0.69) #####################################################################
t_dom_rgn_wvnw_gwe10_i12o12 <- ergm(tnet ~ edges + nodefactor("domain") + nodematch("region") + 
                                      nodematch("testedmen") + nodematch("wvsnw") +
                                      idegree(1:2) + odegree(1:2) +gwesp(0.10, fixed=TRUE, cutoff=5),
                                    control=control.ergm(parallel=num_cores, MCMLE.maxit = 25, 
                                                         parallel.type="PSOCK", seed=123))
summary(t_dom_rgn_wvnw_gwe10_i12o12)
gof_t_dom_rgn_wvnw_gwe10_i12o12 <- btergm::gof(t_dom_rgn_wvnw_gwe10_i12o12)
plot(gof_t_dom_rgn_wvnw_gwe10_i12o12); gof_t_dom_rgn_wvnw_gwe10_i12o12
anova(t_dom_drm_rgn, t_dom_rgn_wvnw_gwe20_in1wvsnw, t_dom_rgn_wvnw_gwe20_i12o14)
#######################################################################################################










summary(tnet ~ edges + ostar(1:10))
summary(tnet ~ edges + istar(1:10))
summary(tnet ~ edges + idegree(1:20))
summary(tnet ~ edges + odegree(1:20))




# must finish coding slippage variable ######################


references 

#https://statnet.org/trac/raw-attachment/wiki/Sunbelt2016/ergm_tutorial.html

#http://statnet.csde.washington.edu/workshops/EUSN/current/ergm/ergm_tutorial.html




