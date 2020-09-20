rm(list = ls())

# Libraries
#detach("package:statnet", unload=TRUE)
#detach("package:sna", unload=TRUE)
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)

setwd("~/racialization-data")
nodelist <- read_csv("racialization-nodelist.csv")
edgelist <- read_csv("racialization-edgelist.csv")
# Create a graph object with the igraph library
mygraph <- graph_from_data_frame( edgelist, vertices=nodelist )

mean_distance(mygraph, directed = TRUE)

# now, we need to construct a nodelist 
for_cents <- data.frame(id = c(1:(igraph::vcount(mygraph))), label = igraph::V(mygraph)$name)
for_cents$deg_cent <- degree(mygraph)
for_cents$indeg_cent <- sna::degree(get.adjacency(mygraph,sparse=F), cmode=="indegree")
for_cents$outdeg_cent <- sna::degree(get.adjacency(mygraph,sparse=F), cmode = "outdegree")
for_cents$btw_cent <- betweenness(mygraph)
for_cents$wtd_deg_cent <- strength(mygraph)
for_cents$eigen_cent <- eigen_centrality(mygraph)$vector
for_cents$page_rank <- page_rank(mygraph)$vector
for_cents$auth_score <- authority.score(mygraph)$vector
for_cents$hub_score <- hub.score(mygraph)$vector
for_cents$info_cent <- sna::infocent(get.adjacency(mygraph,sparse=F))
for_cents$subgraph <- igraph::subgraph_centrality(mygraph)
for_cents$k_core <- coreness(mygraph)
for_cents <- for_cents %>% select(-id)

nodelist <- nodelist %>% 
  full_join(for_cents, by = "label") %>% 
  rename(name = id, pop_tested = group)

rm(mygraph)

test <- nodelist %>% filter(evidence != 0)
cor(test$indeg_cent, test$outdeg_cent)
cor(test$btw_cent, test$eigen_cent)

# create a correlation table of select centrality measures 
cor_table <- for_cents %>% select(indeg_cent, outdeg_cent, btw_cent, eigen_cent, auth_score, hub_score)
cor(cor_table)

# creating first level of hierarchy
unique(nodelist$outcome)

first_level <- data.frame(from="Origin", to=c("Null", "Mixed", "NoComp", "Difference"))

# creating second level of hierarchy
second_level <- nodelist
second_level <- second_level %>%
  #mutate(domain_viz = recode(outcome, `0` = "NonEvidence", `1` = "Biosocial", `2` = "Difference", 
  #                  `3` = "Prostate Cancer", `4` = "PCOS", `5` = "Menopause", `6` = "Other")) %>%
  rename(from = outcome, to = name) %>%
  arrange(from, -deg_cent) %>%
  select(from, to)

# bind levels into a hierarchy, remove levels
hierarchy <- rbind(first_level, second_level)
rm(first_level, second_level)

# create a vertices data.frame. One line per object of our hierarchy, giving features of nodes.
vertices  <-  data.frame(
  name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to)))
)
vertices$group  <-  hierarchy$from[ match( vertices$name, hierarchy$to ) ]
#nodelist <- nodelist %>% rename(name = label)
vertices <- vertices %>% full_join(nodelist, by = "name")
vertices$na_label <- NA

# Create a graph object with the igraph library 
mygraph <- graph_from_data_frame( hierarchy, vertices=vertices )

# create color palette 
palette <- colorRampPalette(brewer.pal(11,"Spectral"))
display.brewer.pal(11, "Spectral")

# graph as a dendrogram
ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
  geom_edge_diagonal() +
  theme_void()

# get the edgelist and cut it down
setwd("~/racialization-data")
connect <- read_csv("racialization-edgelist.csv")
connect <- connect %>%
  rename(from = source, to = target) %>%
  filter(from != to) 
  #top_n(5000, weight) %>%
  #top_frac(0.75, weight) %>%
  #mutate(weight = weight) %>%
  #mutate(weight = (weight-min(weight))/(max(weight)-min(weight))) %>%
  #rename(value = weight)

#set degree
deg <- 360/length(unique(nodelist$name))
#first_level (origin, 4 outcomes)
first_level <- c(0,0,0,0,0)
#find 90 on the right side 
# that was Kim2013B in the 64th row (56th row after removing top 8 randoms)
subset1 = c(1:56)
for (i in 1:56) {subset1[i] = (56*deg) - (i * deg)}
#second
subset2 = c(1:77)
for (i in 1:77) {subset2[i] = 0 - (i * deg)}
#third
subset3 = c(1:77)
for (i in 1:77) {subset3[i] = 90 - (i * deg)}
#fourth
subset4 = c(1:77)
for (i in 1:77) {subset4[i] = 0 - (i * deg)}
#fifth
subset5 = c(1:21)
for (i in 1:21) {subset5[i] = 90 - (i * deg)}
## new experiment
vertices$new_angle <- c(first_level,subset1,subset2,89,subset3,subset4,subset5)

## create hjust
hjust1 = c(1:133)
for (i in 1:133) {hjust1[i] = 0}
hjust2 = c(1:154)
for (i in 1:154) {hjust2[i] = 1}
hjust3 = c(1:21)
for (i in 1:21) {hjust3[i] = 0}
vertices$hjust <- c(first_level,hjust1,89,hjust2,hjust3)

# customize label information 
vertices <- vertices %>%
  mutate(heb_label = ifelse(test = str_detect(string = group, pattern = "NoComp"), yes = NA, no = label)) 

heb_labels <- vertices %>%
  # indeg
  mutate(indeg_colors = ifelse(vertices$indeg_cent > mean(na.omit(vertices$indeg_cent)), "black", "#A9A9A9")) %>%
  mutate(indeg_sizes = ifelse(vertices$indeg_cent > mean(na.omit(vertices$indeg_cent)), 3, 2)) %>%
  # outdeg 
  mutate(outdeg_colors = ifelse(vertices$outdeg_cent > mean(na.omit(vertices$outdeg_cent)), "black", "#A9A9A9")) %>%
  mutate(outdeg_sizes = ifelse(vertices$outdeg_cent > mean(na.omit(vertices$outdeg_cent)), 3, 2)) %>%
  # btw 
  mutate(btw_colors = ifelse(vertices$btw_cent > mean(na.omit(vertices$btw_cent)), "black", "#A9A9A9")) %>%
  mutate(btw_sizes = ifelse(vertices$btw_cent > mean(na.omit(vertices$btw_cent)), 3, 2)) %>%
  # eigen
  mutate(eigen_colors = ifelse(vertices$eigen_cent > mean(na.omit(vertices$eigen_cent)), "black", "#A9A9A9")) %>%
  mutate(eigen_sizes = ifelse(vertices$eigen_cent > mean(na.omit(vertices$eigen_cent)), 3, 2)) %>%
  # hub
  mutate(hub_colors = ifelse(vertices$hub_score > mean(na.omit(vertices$hub_score)), "black", "#A9A9A9")) %>%
  mutate(hub_sizes = ifelse(vertices$hub_score > mean(na.omit(vertices$hub_score)), 3, 2)) %>%
  #auth 
  mutate(auth_colors = ifelse(vertices$auth_score > mean(na.omit(vertices$auth_score)), "black", "#A9A9A9")) %>%
  mutate(auth_sizes = ifelse(vertices$auth_score > mean(na.omit(vertices$auth_score)), 3, 2)) %>%
  #times cited
  mutate(cited_colors = ifelse(vertices$gstimescited > mean(na.omit(vertices$gstimescited)), "black", "#A9A9A9")) %>%
  mutate(cited_sizes = ifelse(vertices$gstimescited > mean(na.omit(vertices$gstimescited)), 3, 2)) %>%
  drop_na(indeg_colors, indeg_sizes) %>% 
  select(indeg_colors, indeg_sizes, outdeg_colors, outdeg_sizes, btw_colors, btw_sizes, 
         eigen_colors, eigen_sizes, hub_colors, hub_sizes, auth_colors, auth_sizes, cited_colors, cited_sizes)

# connect the datasets
from <- match( connect$from, vertices$name)
to <- match( connect$to, vertices$name)

connect$test <- recode(connect$valence, `0` = "blue", `1` = "red", `2` = "orange", `3` = "green")

# domain graphs ###############
edgelist %>% 
  filter(source != target) %>% 
  count()

259/544

# indegree 
ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
  geom_conn_bundle(data = get_con(from = from, to = to), 
                   tension = 0.6, aes(colour=..index..), show.legend = FALSE) +
  scale_edge_colour_gradient(low = "#000000" , high = "#D3D3D3") +
  #scale_edge_linetype_manual(values = connect$test, guide = FALSE) +
  #scale_edge_colour_identity() +
  geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf, label=vertices$heb_label, 
                     hjust=vertices$hjust, angle = vertices$new_angle), 
                 size=heb_labels$indeg_sizes, alpha=1, colour=heb_labels$indeg_colors) +
  theme_void() +
  theme(legend.position=c(0.115,0.185), legend.title=element_blank(), 
        legend.text = element_text(size = 9), plot.margin=unit(c(0,0,0,0),"cm")) +
  guides(size = FALSE, alpha = FALSE, color = guide_legend(override.aes = list(size=4))) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3)) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group, alpha = 0.8, size=vertices$indeg_cent)) + 
  scale_size_continuous( range = c(3,10) ) +
  scale_colour_manual(values= c(  palette(11)[10], "#FFDB58","#D3D3D3", palette(11)[2]))

# outdegree 
ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
  geom_conn_bundle(data = get_con(from = from, to = to), 
                   tension = 0.6, aes(colour=..index..), show.legend = FALSE) +
  scale_edge_colour_gradient(low = "#000000" , high = "#D3D3D3") +
  geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf, label=vertices$heb_label, 
                     hjust=vertices$hjust, angle = vertices$new_angle), 
                 size=heb_labels$outdeg_sizes, alpha=1, colour=heb_labels$outdeg_colors) +
  theme_void() +
  theme(legend.position=c(0.115,0.185), legend.title=element_blank(), 
        legend.text = element_text(size = 9), plot.margin=unit(c(0,0,0,0),"cm")) +
  guides(size = FALSE, alpha = FALSE, color = guide_legend(override.aes = list(size=4))) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3)) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group, alpha = 0.8, size=vertices$outdeg_cent)) + 
  scale_size_continuous( range = c(3,10) ) +
  scale_colour_manual(values= c(  palette(11)[10], "#FFDB58","#D3D3D3", palette(11)[2]))

# btw 
ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
  geom_conn_bundle(data = get_con(from = from, to = to), 
                   tension = 0.6, aes(colour=..index..), show.legend = FALSE) +
  scale_edge_colour_gradient(low = "#000000" , high = "#D3D3D3") +
  geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf, label=vertices$heb_label, 
                     hjust=vertices$hjust, angle = vertices$new_angle), 
                 size=heb_labels$btw_sizes, alpha=1, colour=heb_labels$btw_colors) +
  theme_void() +
  theme(legend.position=c(0.115,0.185), legend.title=element_blank(), 
        legend.text = element_text(size = 9), plot.margin=unit(c(0,0,0,0),"cm")) +
  guides(size = FALSE, alpha = FALSE, color = guide_legend(override.aes = list(size=4))) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3)) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group, alpha = 0.8, size=vertices$btw_cent)) + 
  scale_size_continuous( range = c(3,10) ) +
  scale_colour_manual(values= c(  palette(11)[10], "#FFDB58","#D3D3D3", palette(11)[2]))

# eigen 
ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
  geom_conn_bundle(data = get_con(from = from, to = to), 
                   tension = 0.6, aes(colour=..index..), show.legend = FALSE) +
  scale_edge_colour_gradient(low = "#000000" , high = "#D3D3D3") +
  geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf, label=vertices$heb_label, 
                     hjust=vertices$hjust, angle = vertices$new_angle), 
                 size=heb_labels$eigen_sizes, alpha=1, colour=heb_labels$eigen_colors) +
  theme_void() +
  theme(legend.position=c(0.115,0.185), legend.title=element_blank(), 
        legend.text = element_text(size = 9), plot.margin=unit(c(0,0,0,0),"cm")) +
  guides(size = FALSE, alpha = FALSE, color = guide_legend(override.aes = list(size=4))) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3)) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group, alpha = 0.8, size=vertices$eigen_cent)) + 
  scale_size_continuous( range = c(3,10) ) +
  scale_colour_manual(values= c(  palette(11)[10], "#FFDB58","#D3D3D3", palette(11)[2]))

# times cited 
ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
  geom_conn_bundle(data = get_con(from = from, to = to), 
                   tension = 0.6, aes(colour=..index..), show.legend = FALSE) +
  scale_edge_colour_gradient(low = "#000000" , high = "#D3D3D3") +
  geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf, label=vertices$heb_label, 
                     hjust=vertices$hjust, angle = vertices$new_angle), 
                 size=heb_labels$cited_sizes, alpha=1, colour=heb_labels$cited_colors) +
  theme_void() +
  theme(legend.position=c(0.115,0.185), legend.title=element_blank(), 
        legend.text = element_text(size = 9), plot.margin=unit(c(0,0,0,0),"cm")) +
  guides(size = FALSE, alpha = FALSE, color = guide_legend(override.aes = list(size=4))) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3)) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group, alpha = 0.8, size=vertices$gstimescited)) + 
  scale_size_continuous( range = c(3,10) ) +
  scale_colour_manual(values= c(  palette(11)[10], "#FFDB58","#D3D3D3", palette(11)[2]))


nodelist %>% 
  group_by(outcome) %>% 
  summarise(indeg = mean(indeg_cent),
            cites = mean(gstimescited))











# References
# https://www.r-graph-gallery.com/311-add-labels-to-hierarchical-edge-bundling.html


 



