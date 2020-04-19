rm(list = ls())

# Libraries
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)

setwd("~/racialization-data")
nodelist <- read_csv("racialization-nodelist.csv")
edgelist <- read_csv("racialization-edgelist.csv")
# Create a graph object with the igraph library
mygraph <- graph_from_data_frame( edgelist, vertices=nodelist )

# now, we need to construct a nodelist 
for_cents <- data.frame(id = c(1:(igraph::vcount(mygraph))), label = igraph::V(mygraph)$name)
for_cents$deg_cent <- degree(mygraph)
for_cents$btw_cent <- betweenness(mygraph)
for_cents$wtd_deg_cent <- strength(mygraph)
for_cents$eigen_cent <- eigen_centrality(mygraph)$vector
for_cents$page_rank <- page_rank(mygraph)$vector
for_cents$auth_score <- authority.score(mygraph)$vector
for_cents$hub_score <- hub.score(mygraph)$vector
for_cents$k_core <- coreness(mygraph)
components <- components(mygraph)
for_cents$component <- components$membership
for_cents <- for_cents %>% select(-id)

nodelist <- nodelist %>% 
  full_join(for_cents, by = "label") %>% 
  rename(name = id, pop_tested = group)

rm(mygraph)

# creating first level of hierarchy
unique(nodelist$domaincode)

first_level <- data.frame(from="Origin",
                          to=c("NonEvidence", "Biosocial", "Difference",
                               "Prostate Cancer", "PCOS", "Menopause",  "Other"))

# creating second level of hierarchy
second_level <- nodelist
second_level <- second_level %>%
  mutate(domain_viz = recode(domaincode, `0` = "NonEvidence", `1` = "Biosocial", `2` = "Difference", 
                    `3` = "Prostate Cancer", `4` = "PCOS", `5` = "Menopause", `6` = "Other")) %>%
  rename(from = domain_viz, to = name) %>%
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
#first_level (origin, 7 domains)
first_level <- c(0,0,0,0,0,0,0,0)
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
  mutate(heb_label = ifelse(test = str_detect(string = group, pattern = "NonEvidence"), yes = NA, no = label)) 
heb_labels <- vertices %>%
  mutate(label_colors = ifelse(vertices$btw_cent > mean(na.omit(vertices$btw_cent)), "black", "#A9A9A9")) %>%
  mutate(label_sizes = ifelse(vertices$btw_cent > mean(na.omit(vertices$btw_cent)), 3, 2)) %>%
  drop_na(label_colors, label_sizes) %>% select(label_colors, label_sizes)

mean(na.omit(vertices$btw_cent))

# connect the datasets
from <- match( connect$from, vertices$name)
to <- match( connect$to, vertices$name)

# graph connections
ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
  geom_conn_bundle(data = get_con(from = from, to = to), tension = 0.6, aes(colour=..index..)) +
  scale_edge_colour_gradient(low = "#000000" , high = "#D3D3D3") +
  geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf, label=vertices$heb_label, 
                     hjust=vertices$hjust, angle = vertices$new_angle
                     ), size=heb_labels$label_sizes, alpha=1, colour=heb_labels$label_colors) +
  theme_void() +
  theme(legend.position="none", plot.margin=unit(c(0,0,0,0),"cm")) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3)) +
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group, size=vertices$btw_rank, 
                      alpha = 0.8)) +
  scale_size_continuous( range = c(3,10) ) +
  scale_colour_manual(values= c(palette(11)[9], #Menopause 
                                palette(11)[5], #Difference 
                                palette(11)[3], 
                                "#D3D3D3",     #NonEvidence
                                palette(11)[1], 
                                palette(11)[11], #PCOS 
                                palette(11)[10]  #ProstateCancer
                                ))














# References
# https://www.r-graph-gallery.com/311-add-labels-to-hierarchical-edge-bundling.html


 



