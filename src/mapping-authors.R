rm(list = ls())

library(tidyverse)
library(maps)
library(RColorBrewer)
library(geosphere)

setwd("~/racialization-data")
nodelist <- read_csv("racialization-nodelist.csv")
edgelist <- read_csv("racialization-edgelist.csv")

testlist1 <- edgelist %>% 
  filter(source != target) %>% 
  select(source) %>% 
  inner_join(nodelist %>% 
             rename(source = id) %>% 
             select(source, authorlon, authorlat), 
             by = "source") %>% 
  #distinct(source, .keep_all = TRUE) %>% 
  select(-source) 

testlist2 <- edgelist %>%
  filter(source != target) %>% 
  select(target) %>% 
  inner_join(nodelist %>% 
               rename(target = id) %>% 
               select(target, authorlon, authorlat), 
             by = "target") %>% 
  #distinct(target, .keep_all = TRUE) %>% 
  select(-target) 


testlist3 <- cbind(testlist1, testlist2) 
testlist3 <- testlist3 %>% rename


  
  
  ###3
data <- rbind(
    Buenos_aires=c(-58,-34),
    Paris=c(2,49),
    Melbourne=c(145,-38),
    Saint.Petersburg=c(30.32, 59.93),
    Abidjan=c(-4.03, 5.33),
    Montreal=c(-73.57, 45.52),
    Nairobi=c(36.82, -1.29),
    Salvador=c(-38.5, -12.97)
  )  %>% as.data.frame()
  colnames(data)=c("long","lat")
  
all_pairs <- cbind(t(combn(data$long, 2)), t(combn(data$lat, 2))) %>% as.data.frame()
colnames(all_pairs) <- c("long1","long2","lat1","lat2")
  
  
  
  
  
  
  
  

# A function to plot connections
plot_my_connection=function( dep_lon, dep_lat, arr_lon, arr_lat, ...){
  inter <- gcIntermediate(c(dep_lon, dep_lat), c(arr_lon, arr_lat), n=50, addStartEnd=TRUE, breakAtDateLine=F)             
  inter=data.frame(inter)
  diff_of_lon=abs(dep_lon) + abs(arr_lon)
  if(diff_of_lon > 180){
    lines(subset(inter, lon>=0), ...)
    lines(subset(inter, lon<0), ...)
  }else{
    lines(inter, ...)
  }
}

# https://www.r-graph-gallery.com/how-to-draw-connecting-routes-on-map-with-r-and-great-circles.html


# Background map
map('world',col="#f2f2f2", fill=TRUE, bg="white", lwd=0.05,mar=rep(0,4),border=0, ylim=c(-80,80) )

# Circles for cities
points(x=nodelist$authorlon, y=nodelist$authorlat, col="slateblue", cex=3, pch=20)

# Connections
plot_my_connection(testlist1[1], testlist1[2], testlist2[1], testlist2[2], col="slateblue", lwd=2)


plot_my_connection(Buenos_aires[1], Buenos_aires[2], Melbourne[1], Melbourne[2], col="slateblue", lwd=2)
plot_my_connection(Buenos_aires[1], Buenos_aires[2], Paris[1], Paris[2], col="slateblue", lwd=2)


testlist1[1]


map('world',
    col="#f2f2f2", fill=TRUE, bg="white", lwd=0.05,
    mar=rep(0,4),border=0, ylim=c(-80,80) 
)
points(x=nodelist$authorlon, y=nodelist$authorlat, col="blue", cex=1, pch=20)

# Compute the connection between Buenos Aires and Paris
inter <- gcIntermediate(Paris,  Buenos_aires, n=50, addStartEnd=TRUE, breakAtDateLine=F)

# Show this connection
lines(inter, col="slateblue", lwd=2)

Buenos_aires <- c(-58,-34)
Paris <- c(2,49)
Melbourne <- c(145,-38)
data <- rbind(Buenos_aires, Paris, Melbourne) %>% 
  as.data.frame()
colnames(data) <- c("long","lat")

# Compute the connection between Buenos Aires and Paris
inter <- gcIntermediate(Paris,  Buenos_aires, n=50, addStartEnd=TRUE, breakAtDateLine=F)

# Show this connection
lines(inter, col="slateblue", lwd=2)



# graph the communities on the world map to see if they need consolidating
world_map <- map_data("world") 

world_map %>% 
  ggplot() +
  theme_minimal() +
  geom_polygon(aes(x = long, y = lat, fill = "", group = group), color = "grey") +
  geom_point(data = nodelist, aes(x = authorlon, y = authorlat)) 



  coord_fixed(1.2) +
  guides(fill=FALSE) + 
  scale_fill_manual(values = c("#E57200", # dark red    
                               "#628ed8", # light blue
                               "#990000", # orange 
                               "#232D4B", # dark blue
                               "#eaaa31"  # yellow
  )) +
  theme(plot.title = element_text(size = 13, hjust = 0.5), 
        plot.caption = element_text(size = 10, vjust = 8, hjust = 0.85),
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank()) +
  labs(title= " Modularity Groupings for Country-Country \n Collaboration Networks (GitHub, 2008-2019)",
       caption = "Note: Groupings Calculated with Loops in Network")