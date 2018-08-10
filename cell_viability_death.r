library(tidyverse)

## Load CSV file ##
cell_death_data <- read_csv("cell_viability_cell_death.csv")
# Get rid of blank columns and well number column
cell_death_data <- dplyr::select(cell_death_data, -contains("X"), -contains("Well"))
# Change some variables to factor type
cell_death_data$Repeat <- as.factor(cell_death_data$Repeat)
cell_death_data$Compound <- as.factor(cell_death_data$Compound)
cell_death_data$CellType <- as.factor(cell_death_data$CellType)
