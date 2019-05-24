# This is R replication code for a paper under review named "Racialization in Testosterone Research." The data include 149 papers that have examined population differences in testosterone from 1966-2017. I plan to make this data available at a later date, largely because I plan to write to some additional papers using these data. 

#In the first section, the code outlines descriptive network statistics using statnet, including degree and betweenness centrality. In the second section, I use igraph and visNetwork to visualize the data based on two cuts. The first looks at age/gender groupings in my data and the second looks at outcomes, which is the main focus of the paper I currently have under review. 

#clearing old data
rm(list = ls())
gc()

########################################################## network stats ##################

library(sna)
library(statnet)

#creating a statnet object 
dat <- gs_url("https://docs.google.com/spreadsheets/d/1Ofq18UN39QH-8OPxT-gzF0_J_dVrckpG0th2X5ktbIg/edit#gid=922546532")

#setting path (if applicable)
#setwd("C:/Users/soren/Google Drive/Biomedical MultipliciTs/Racial Differences in T/Network")

dat <- dat %>% 
  gs_read(ws = "Edgelist")
#dat=read.csv(file.choose(),header=TRUE)
el=as.matrix(dat)            
el[,1]=as.character(el[,1])  
el[,2]=as.character(el[,2])
testosterone=network(dat,matrix.type="edgelist",directed=TRUE) 
plot(testosterone)

#network size, number of edges, number of dyads 
summary(testosterone)
network.size(testosterone)
network.edgecount(testosterone)
network.dyadcount(testosterone)
as.sociomatrix(testosterone)

#output 
netstats <- data.frame(
  ego=network.vertex.names(testosterone),
  indegree=degree(testosterone, cmode="indegree"),             
  outdegree=degree(testosterone, cmode="outdegree"), 
  betweenness=betweenness(testosterone, cmode="directed"))
netstats
setwd("C:/Users/soren/Google Drive/Biomedical MultipliciTs/Racial Differences in T/Network")
write.csv((netstats),file = 'RacializationxTestosterone_NetworkStats.csv')
# I then manually added in the applicable network statistics into my original data file.

detach("package:statnet")
detach("package:sna")
detach("package:igraph")

#clearing old data
rm(list = ls())
gc()

########################################################## visualizations ##################

library('RColorBrewer')
library('extrafont')
library('igraph')
library('plyr') 
library('visNetwork')

#detach("package:igraph")
#detach("package:visNetwork")

# creating nodes & edges from google sheet
network_data <- gs_url("https://docs.google.com/spreadsheets/d/1Ofq18UN39QH-8OPxT-gzF0_J_dVrckpG0th2X5ktbIg/edit#gid=922546532")

nodes <- network_data %>% 
  gs_read(ws = "Nodelist")
edges <- network_data %>% 
  gs_read(ws = "Edgelist")

#creating nodes & edges from local csv files (if applicable)
#setwd("C:/Users/soren/Google Drive/Biomedical MultipliciTs/Racial Differences in T/Network")
#setwd("C:/Users/bkram/CloudStation/Biomedical MultipliciTs/Racial Differences in T/Network")
#nodes <- read.csv(file.choose(),header=TRUE, as.is=TRUE) 
#edges <- read.csv(file.choose(),header=TRUE, as.is=TRUE)

#remove all the self-loops (required for igraph, not visNetwork)
edges <- edges[which(edges$loop=='0'),]

#checking the data out 
#note nodelist needs "id", edgelist needs "from", "to" in visNetwork
head(nodes)
head(edges)
nrow(nodes); length(unique(nodes$id))
nrow(edges); nrow(unique(edges[,c("from", "to")]))

#making the network with igraph and mapping with visIgraph
#net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T) 
#net <- simplify(g, remove.multiple = F, remove.loops = T)
#visIgraph(net)

#mapping data with visNetwork
visNetwork(nodes, edges)

#creating color palettes for the network graphs  
palette1 <- colorRampPalette(brewer.pal(11,"Spectral"))
palette2 <- colorRampPalette(brewer.pal(11,"PiYG"))
palette3 <- colorRampPalette(brewer.pal(11,"BrBG"))
palette4 <- colorRampPalette(brewer.pal(4,"Greys"))
palette5 <- colorRampPalette(brewer.pal(11,"RdYlGn"))
display.brewer.pal(11, "Spectral")
display.brewer.pal(11, "PiYG")
display.brewer.pal(11, "BrBG")
display.brewer.pal(4, "Greys")
display.brewer.pal(11, "RdYlGn")

#applying color palettes to age/gender graph 
nodes$color = nodes$group
nodes$color = gsub("AdultMale",palette1(11)[10],nodes$color) #blue 
nodes$color = gsub("AdultFemale",palette2(11)[9],nodes$color) #green
nodes$color = gsub("Children", palette1(11)[2],nodes$color) #red
nodes$color = gsub("MixedSample",palette3(11)[3],nodes$color) #orange
nodes$color = gsub("NonEvidence",palette4(4)[1],nodes$color) #light grey 

#nodes based on indegree
nodes$size <- nodes$indegree*7
nodes$label <- nodes$id

#graph with random seed, force atlas, no overlap 
visNetwork(nodes, edges) %>% 
  visOptions(selectedBy = "group") %>%
  visLayout(randomSeed = 12345) %>% 
  visPhysics(solver = "forceAtlas2Based",
            forceAtlas2Based = list(avoidOverlap = 0.01))

#output graph to html file 
netviz1 <- visNetwork(nodes, edges) %>% 
  visOptions(selectedBy = "group") %>%
  visLayout(randomSeed = 12345) %>% 
  visPhysics(solver = "forceAtlas2Based",
             forceAtlas2Based = list(avoidOverlap = 0.0001))
setwd("C:/Users/bkram/CloudStation/Biomedical MultipliciTs/Racial Differences in T/Network")
visSave(netviz1, file = "by-groups.html")

#applying color palettes to outcomes graph 
nodes$color=nodes$outcome
nodes$color=gsub("Mixed",palette1(11)[3],nodes$color) #orange
nodes$color=gsub("Mixed",palette5(11)[5],nodes$color) #orange 
nodes$color=gsub("Difference",palette2(11)[10],nodes$color) #green
nodes$color=gsub("Null",palette1(11)[1],nodes$color) #magenta
nodes$color=gsub("NoComp",palette4(4)[1],nodes$color) #light grey

#nodes based on indegree
nodes$size <- nodes$indegree*7
nodes$label <- nodes$id

visNetwork(nodes, edges) %>% 
  visOptions(selectedBy = "outcome", nodesIdSelection = TRUE) %>%
  visLayout(randomSeed = 12345) %>% 
  visPhysics(solver = "forceAtlas2Based",
             forceAtlas2Based = list(avoidOverlap = 0.01))

#output graph to html file
netviz2 <- visNetwork(nodes, edges) %>% 
  visOptions(selectedBy = "outcome") %>%
  visLayout(randomSeed = 12345) %>% 
  visPhysics(solver = "forceAtlas2Based",
             forceAtlas2Based = list(avoidOverlap = 0.0001))
visSave(netviz2, file = "by-outcomes.html")





























