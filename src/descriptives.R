# Racialization in Testosterone Research 

# This project examines how testosterone has been used to racialize populations in scientific research. For a synopsis of the methods, please see Kramer's "Racialization in Testosterone Research." In the first section of this document, I outline the basic descriptives of population comparisons in testosterone research. 

# Loading data from Google Sheets 

library(googlesheets)
suppressMessages(library(dplyr))

racialization_raw <- gs_url("https://docs.google.com/spreadsheets/d/1Ofq18UN39QH-8OPxT-gzF0_J_dVrckpG0th2X5ktbIg/edit#gid=922546532")

racialization <- racialization_raw %>% 
  gs_read(ws = "Nodelist")
racialization

# Subsetting data based on evidence only 
# (Removing cited studies that do not contain population comparisons)

evidence <- subset(racialization, evidence == 1)
evidence

# Calculating "racial testosterone theory" frequencies 

wbcomp <- evidence$wbcomp
bacomp <- evidence$bacomp
wacomp <- evidence$wacomp

wb <- table(wbcomp)
wb
149 - 45
104 / 149
prop.table(wb)

wa <- table(wacomp)
wa
149 - 103
46 / 149
prop.table(wa)

ba <- table(bacomp)
ba
149 - 130
19 / 149 
prop.table(ba)

# Breaking down racial testosterone theory across men, women and children 

table(racialization$group_stats)

male_only <- subset(racialization, group_stats == "AdultMale" | 
                      group_stats == "MaleFemale" | group_stats == "MaleFemaleChildren")
male_wbonly <- subset(male_only, wbcomp != "NoComp")
male_waonly <- subset(male_only, wacomp != "NoComp")
male_baonly <- subset(male_only, bacomp != "NoComp")
female_only <- subset(racialization, group_stats == "AdultFemale" | 
                        group_stats == "MaleFemale" | group_stats == "MaleFemaleChildren")
female_wbonly <- subset(female_only, wbcomp != "NoComp")
female_waonly <- subset(female_only, wacomp != "NoComp")
female_baonly <- subset(female_only, bacomp != "NoComp")
children_only <- subset(racialization, group_stats == "Children" | group_stats == "MaleFemaleChildren")
children_wbonly <- subset(children_only, wbcomp != "NoComp")
children_waonly <- subset(children_only, wacomp != "NoComp")
children_baonly <- subset(children_only, bacomp != "NoComp")

wbmale <- table(male_wbonly$wbcomp)
wbmale; sum(wbmale); prop.table(wbmale); 
wamale <- table(male_waonly$wacomp)
wamale; (wamale); prop.table(wamale)
bamale <- table(male_baonly$bacomp)
bamale; (bamale); prop.table(bamale)

wbfemale <- table(female_wbonly$wbcomp)
wbfemale; sum(wbfemale); prop.table(wbfemale)
wafemale <- table(female_waonly$wacomp)
wafemale; sum(wafemale); prop.table(wafemale)
bafemale <- table(female_baonly$bacomp)
bafemale; sum(bafemale); prop.table(bafemale)

wbchild <- table(children_wbonly$wbcomp)
wbchild; sum(wbchild); prop.table(wbchild)
wachild <- table(children_waonly$wacomp)
wachild; sum(wachild); prop.table(wachild)
bachild <- table(children_baonly$bacomp)
bachild; sum(bachild); prop.table(bachild)

ev_table <- table(evidence$outcome)
prop.table(ev_table)

test <- table(evidence$datarecycling); 56/149
table(evidence$timediff); 149-93; 56/149
table(evidence$agediff); 49/149
table(evidence$bmidiff); 75/149
table(evidence$illnesscontrolled); 62/149
table(evidence$sescontrolled); 11/149
table(evidence$medscontrolled); 54/149

# Graphing publications that racialize populations in testosterone research over time 

library(googlesheets)
library(ggplot2)

timeseries_raw <- gs_url("https://docs.google.com/spreadsheets/d/1fkHvv-Hq0s_qzJDzyh24Tjp06ulsUcoX7REeG44xSgE/edit#gid=1154073267")

time_data <- timeseries_raw %>% 
  gs_read(ws = "Table 3 (Over Time)")

time_data

date <- time_data$YearPublished
Publications <- time_data$TotalPublications

time_series <- ggplot(data = time_data, aes(x = date, y = Publications)) +
  geom_line(size = 1.25) + 
  theme_light() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_continuous(name="Publications Per Year") +
  theme(axis.title.x=element_blank())

time_series

# Creating a six figure panel of studies that conduct population comparisons in testosterone research across four time periods (1966-1989, 1966-1999, 1966-2009, 1966-2017) using SVG visualizations created in Gephi 0.9.2. The last two panels add in the breakdown of the full network with men, women and children and the final panel includes the graph of the publications in the network over time. 

library("magick")
library("multipanelfigure")

#setwd("C:/Users/soren/Google Drive/Biomedical MultipliciTs/Racial Differences in T/Network/Over Time")
setwd("C:/Users/bkram/CloudStation/Biomedical MultipliciTs/Racial Differences in T/Network/Over Time")

svg1989 <- image_read_svg('1966-1989.svg')
svg1999 <- image_read_svg('1966-1999.svg')
svg2009 <- image_read_svg('1966-2009.svg')
svg2017 <- image_read_svg('1966-2017_NoLabels.svg')
svggender <- image_read_svg('1966-2017_Gender.svg')

print(svg1989)
print(svg1999)
print(svg2009)
print(svg2017)
print(svggender)

svg1989R <- image_scale(svg1989, "350")
svg1999R <- image_scale(svg1999, "350")
svg2009R <- image_scale(svg2009, "350")
svg2017R <- image_scale(svg2017, "350")
svggenderR <- image_scale(svggender, "350")

print(svg1989R)
print(svg1999R)
print(svg2009R)
print(svg2017R)
print(svggenderR)

image_write(svg1989R, path = "svg1989R.svg", format = "svg")
image_write(svg1999R, path = "svg1999R.svg", format = "svg")
image_write(svg2009R, path = "svg2009R.svg", format = "svg")
image_write(svg2017R, path = "svg2017R.svg", format = "svg")
image_write(svggenderR, path = "svggenderR.svg", format = "svg")

#creating figure 

figure <- multi_panel_figure(width = 300, height = 150, 
                             row_spacing = 4.2, column_spacing = 3.4, 
                             columns = 3, rows = 2,
                             panel_label_type = "upper-alpha")
figure %<>% fill_panel(paste0("svg1989R.svg"), rows = 1, columns = 1)
figure %<>% fill_panel(paste0("svg1999R.svg"), rows = 1, columns = 2)
figure %<>% fill_panel(paste0("svg2009R.svg"), rows = 1, columns = 3)
figure %<>% fill_panel(paste0("svg2017R.svg"), rows = 2, columns = 1)
figure %<>% fill_panel(paste0("svggenderR.svg"), rows = 2, columns = 2)
figure %<>% fill_panel(time_series, rows = 2, columns = 3)
figure

figure %>% save_multi_panel_figure(filename = "BioSocieties.svg")
figure %>% save_multi_panel_figure(filename = "BioSocieties.png")
