# NOTE -- COMPATIBILITY
# This app was built using the following software versions:
#
# R version 3.3.1 (2016-06-21)
# igraph_1.0.1 
# bio3d_2.2-4
# shinyjs_0.7  
# shinyBS_0.61   
# shiny_0.13.2

# Libraries

library(bio3d)
library(igraph)
library(shinyBS)
library(shinyjs)
library(ggplot2)

# Source dependencies

source("cluster_function_gui.R")
source("file_choose2.R")
source("create_graph_from_membership_vector.R")
source("convert_to_kinetic_membership_vector.R")
