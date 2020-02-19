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
