library(tidyverse)
library(ggplot2)

## Load CSV file ##
cell_death_data <- read_csv("cell_viability_cell_death.csv")
# Get rid of blank columns and well number column
cell_death_data <- dplyr::select(cell_death_data, -contains("X"), -contains("Well"))
# Change some variables to factor type
cell_death_data$Repeat <- as.factor(cell_death_data$Repeat)
cell_death_data$Compound <- as.factor(cell_death_data$Compound)
cell_death_data$CellType <- as.factor(cell_death_data$CellType)

cell_death_data <- filter(cell_death_data, Compound != "n/a")
glimpse(cell_death_data)

# Box plots, see variation
ggplot(data = cell_death_data, aes(x = Compound, y = FluorescencePrestoBlue, color = Compound)) +
  geom_boxplot() + 
  facet_grid(~CellType)
ggplot(data = cell_death_data, aes(x = Compound, y = FluorescenceCyQUANT, color = Compound)) +
  geom_boxplot() + 
  facet_grid(~CellType)

## Presto blue
## Group by variables, then take mean values
mean_presto_data <- group_by(cell_death_data, Compound, TimePoint, CellType) %>%
  summarise(FluorescencePrestoBlue = mean(FluorescencePrestoBlue, na.rm = TRUE))

ggplot(data = mean_presto_data, aes(x = TimePoint, y = FluorescencePrestoBlue, color = Compound)) + 
  geom_line() +
  facet_grid(~CellType)
ungroup(cell_death_data)


## CyQuant
## Group by variables, then take mean values
mean_cyquant_data <- group_by(cell_death_data, Compound, TimePoint, CellType) %>%
  summarise(FluorescenceCyQUANT = mean(FluorescenceCyQUANT, na.rm = TRUE))

ggplot(data = mean_cyquant_data, aes(x = TimePoint, y = FluorescenceCyQUANT, color = Compound)) + 
  geom_line() +
  geom_point() +
  facet_grid(~CellType)
ungroup(cell_death_data)
