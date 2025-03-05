setwd("/home/mrna/R/x86_64-pc-linux-gnu-library/4.4/")

#install.packages("tidyr")
#install.packages("devtools", dependencies = "TRUE")
#devtools::install_github("yulab-smu/ggacmap")
#install.packages("tidyverse")
#Load the Racmacs package
library(Racmacs)
library(ggplot2)
library(ggacmap)
#Set an option for the number of computer cores to run in parallel when optimizing maps
#The default when running on CRAN is to use 2
options(Racoptimizer.num_cores = 1)
#However you can also set the number of cores to the maximum number like this
#options(RacOptimizer.num_cores = parallel::detectCores())

#Read in the titer table
path_to_titre_file <- system.file("extdata/lab _titer_data.csv", package = "Racmacs")
titre_table <- read.titerTable(path_to_titre_file)
print(titre_table[1:8,1:9])

#Create the acmap object, specifying the titer table
map <- acmap(
  titer_table = titre_table
)

#Perform some optimization runs on the map object to try and determine a best map
map <- optimizeMap(
  map = map,
  number_of_dimensions = 2,
  number_of_optimizations = 1000,
  minimum_column_basis = "none"
)
#> a
#> Took 51.6 secs
#> Warning in optimizeMap(map = map, number_of_dimensions = 2,
#> number_of_opytimizations = 500, : There is some variation (9.57 AU for one point) in the top runs, 
#> this may be an indication that more optimization runs could help achieve a better optimum. If this 
#> still fails to help see ?unstableMaps for further causes.
#> )

plot(map)
view(map)

#setwd("/home/mrna/R/x86_64-pc-linux-gnu-library/4.4/")
#f <- system.file("extdata/Day1/files/lab_titer_data.ace", package = "ggacmap")
#lab_map <- read.acmap(f)

#p <- ggacmap(lab_map) 
#geom_point(aes(color = I(color)))
#p + annotate("text", x = 5, y = 5)
#print(p)

plot(
  map,
  optimization_number = 1,
  xlim = NULL,
  ylim = NULL,
  plot_ags = TRUE,
  plot_sr = TRUE,
  plot_labels = FALSE,
  plot_blobs = TRUE,
  point_opacity = "automatic",
  show_procrustes = TRUE,
  show_error_lines = FALSE,
  plot_stress = FALSE,
  indicate_outliers = "arrowheads",
  grid.col = "grey90",
  grid.margin.col = "grey50",
  outlier.arrow.col = grid.col,
  fill.alpha = 0.8,
  outline.alpha = 0.8,
  procrustes.lwd = 2,
  procrustes.col = "black",
  procrustes.arr.type = "triangle",
  procrustes.arr.length = 0.2,
  procrustes.arr.width = 0.15,
  label.offset = 0,
  padding = 1,
  cex = 1,
  margins = rep(0.5, 4)
)







