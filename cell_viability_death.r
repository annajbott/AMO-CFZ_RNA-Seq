library(tidyverse)
library(ggplot2)
source("fun.R")

## Load CSV file ##
cell_death_data <- read_csv("cell_viability_cell_death.csv")
# Get rid of blank columns and well number column
cell_death_data <- dplyr::select(cell_death_data, -contains("X"), -contains("Well"))
# Change some variables to factor type
cell_death_data$Repeat <- as.factor(cell_death_data$Repeat)
cell_death_data$Compound <- as.factor(cell_death_data$Compound)
cell_death_data$CellType <- as.factor(cell_death_data$CellType)

medium_data <- filter(cell_death_data, str_detect(SampleName, "medium"))

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


## Just  ##
cell_death_data_cfz <- filter(cell_death_data, CellType != "WT")


max_cells <- filter(cell_death_data_cfz, TimePoint == 6)
max_cells$TimePoint <- 0
# zero time equal to average dmso at 6 hrs
max_cells$FluorescencePrestoBlue <- 190743
cell_death_data_cfz <- rbind(max_cells,cell_death_data_cfz)
cell_death_data_cfz <- filter(cell_death_data_cfz, Compound != "DMSO")


cell_death_data_cfz_summary <- summarySE(cell_death_data_cfz, measurevar="FluorescencePrestoBlue", groupvars=c("Compound", "TimePoint"))

pd <- position_dodge(0) # move them .05 to the left and right
plot1 <- ggplot(cell_death_data_cfz_summary, aes(x=TimePoint, y=FluorescencePrestoBlue, colour=Compound))  +
  geom_errorbar(aes(ymin=FluorescencePrestoBlue-se, ymax=FluorescencePrestoBlue+se), width=3, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) + xlim(0,50)

plot1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +ylab("Fluorescence") +xlab("Time (hours)")
