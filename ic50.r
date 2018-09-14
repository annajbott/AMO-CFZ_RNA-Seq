library("ic50")
library(tidyverse)
library("nplr")

###############
## IC50 data ##
###############

IC50_data <- read_csv("temp/IC50s.csv")

IC50_data <- filter(IC50_data, Compound != "DMSO")
IC50_data <- rename(IC50_data, `Cell Type` = "CellType", `Concentration (uM)` = "Concentration", `Fluorescence 1` = "F1", `Fluorescence 2` = "F2", `Fluorescence 3` = "F3", `Average Fluorescence (minus background)` = "Average")

IC50_data$CellType <- as.factor(IC50_data$CellType)
IC50_data$Compound <- as.factor(IC50_data$Compound)

IC50_WT <- filter(IC50_data, CellType == "WT")
IC50_CFZ <- filter(IC50_data, CellType == "CFZ")

### CFZ ###
IC50_tidy_CFZ <- gather(IC50_CFZ, key = Replicate, value = Fluorescence, F1, F2, F3)

mean_dmso <- mean(c(130877,135524,132845,130621,128336,122588,125774,122523,131636,130491))
max_dmso <- max(c(130877,135524,132845,130621,128336,122588,125774,122523,131636,130491))
mean_medium <- mean(33412.1,35434.7,42388.8)

## NCP 26
IC50_tidy_CFZ_26 <- filter(IC50_tidy_CFZ, Compound == "26")
# Use Time zero for minimal value to consider, just using average of medium with no cells. Ctrl is treated cells with control so DMSO
IC50_prop_CFZ_26 <- convertToProp(IC50_tidy_CFZ_26$Fluorescence, T0 = mean_medium, Ctrl = max_dmso)
np1 <- nplr(x=IC50_tidy_CFZ_26$Concentration, y=IC50_prop_CFZ_26)
plot(np1, pcol="grey40", lcol="skyblue1", showEstim=0.5, showGOF = FALSE, showInfl=FALSE, showCI = FALSE, main="NCP 26 treated AMO-CFZ cells", cex.main=1.5,xlab=expression(Log[10](Concentration) (uM)), ylab="Proportion")

## NCP 22
IC50_tidy_CFZ_22 <- filter(IC50_tidy_CFZ, Compound == "22")
# Use Time zero for minimal value to consider, just using average of medium with no cells. Ctrl is treated cells with control so DMSO
IC50_prop_CFZ_22 <- convertToProp(IC50_tidy_CFZ_22$Fluorescence, T0 = mean_medium, Ctrl = max_dmso)
np1 <- nplr(x=IC50_tidy_CFZ_22$Concentration, y=IC50_prop_CFZ_22)
plot(np1, pcol="grey40", lcol="skyblue1", showEstim=0.5, showGOF = FALSE, showInfl=FALSE, showCI = FALSE, main="NCP 22 treated AMO-CFZ cells", cex.main=1.5,xlab=expression(Log[10](Concentration) (uM)), ylab="Proportion")

## NCP 13
IC50_tidy_CFZ_13 <- filter(IC50_tidy_CFZ, Compound == "13")
# Use Time zero for minimal value to consider, just using average of medium with no cells. Ctrl is treated cells with control so DMSO
IC50_prop_CFZ_13 <- convertToProp(IC50_tidy_CFZ_13$Fluorescence, T0 = mean_medium, Ctrl = max_dmso)
np1 <- nplr(x=IC50_tidy_CFZ_13$Concentration, y=IC50_prop_CFZ_13)
plot(np1, pcol="grey40", lcol="skyblue1", showEstim=0.5, showGOF = FALSE, showInfl=FALSE, showCI = FALSE, main="NCP 13 treated AMO-CFZ cells", cex.main=1.5,xlab=expression(Log[10](Concentration) (uM)), ylab="Proportion")


## NCP 18
IC50_tidy_CFZ_18 <- filter(IC50_tidy_CFZ, Compound == "18")
# Use Time zero for minimal value to consider, just using average of medium with no cells. Ctrl is treated cells with control so DMSO
IC50_prop_CFZ_18 <- convertToProp(IC50_tidy_CFZ_18$Fluorescence, T0 = mean_medium, Ctrl = max_dmso)
np1 <- nplr(x=IC50_tidy_CFZ_18$Concentration, y=IC50_prop_CFZ_18)
plot(np1, pcol="grey40", lcol="skyblue1", showEstim=0.5, showCI = FALSE, showGOF = FALSE, showInfl=FALSE,  main="NCP 18 treated AMO-CFZ cells", cex.main=1.5,xlab=expression(Log[10](Concentration) (uM)), ylab="Proportion")

####
### WT ###
IC50_tidy_WT <- gather(IC50_WT, key = Replicate, value = Fluorescence, F1, F2, F3)

mean_dmso <- mean(c(96842,116872,117621,113390,121062,111754,129983,111803,118298,128105))
max_dmso <- max(c(96842,116872,117621,113390,121062,111754,129983,111803,118298,128105))
mean_medium <- mean(21258.7,21707.4,19920.7)

## NCP 26
IC50_tidy_WT_26 <- filter(IC50_tidy_WT, Compound == "26")
# Use Time zero for minimal value to consider, just using average of medium with no cells. Ctrl is treated cells with control so DMSO
IC50_prop_WT_26 <- convertToProp(IC50_tidy_WT_26$Fluorescence, T0 = mean_medium, Ctrl = max_dmso)
np1 <- nplr(x=IC50_tidy_WT_26$Concentration, y=IC50_prop_WT_26)
plot(np1, pcol="grey40", lcol="skyblue1", showEstim=0.5, showGOF = FALSE, showInfl=FALSE, showCI = FALSE, main="NCP 26 treated AMO-WT cells", cex.main=1.5,xlab=expression(Log[10](Concentration) (uM)), ylab="Proportion")

## NCP 22
IC50_tidy_WT_22 <- filter(IC50_tidy_WT, Compound == "22")
# Use Time zero for minimal value to consider, just using average of medium with no cells. Ctrl is treated cells with control so DMSO
IC50_prop_WT_22 <- convertToProp(IC50_tidy_WT_22$Fluorescence, T0 = mean_medium, Ctrl = max_dmso)
np1 <- nplr(x=IC50_tidy_WT_22$Concentration, y=IC50_prop_WT_22)
plot(np1, pcol="grey40", lcol="skyblue1", showEstim=0.5, showGOF = FALSE, showInfl=FALSE, showCI = FALSE, main="NCP 22 treated AMO-WT cells", cex.main=1.5,xlab=expression(Log[10](Concentration) (uM)), ylab="Proportion")

## MAZ 13
IC50_tidy_WT_13 <- filter(IC50_tidy_WT, Compound == "13")
# Use Time zero for minimal value to consider, just using average of medium with no cells. Ctrl is treated cells with control so DMSO
IC50_prop_WT_13 <- convertToProp(IC50_tidy_WT_13$Fluorescence, T0 = mean_medium, Ctrl = max_dmso)
np1 <- nplr(x=IC50_tidy_WT_13$Concentration, y=IC50_prop_WT_13)
plot(np1, pcol="grey40", lcol="skyblue1", showEstim=0.5, showGOF = FALSE, showInfl=FALSE, showCI = FALSE, main="MAZ 13 treated AMO-WT cells", cex.main=1.5,xlab=expression(Log[10](Concentration) (uM)), ylab="Proportion")

## MAZ 18
IC50_tidy_WT_18 <- filter(IC50_tidy_WT, Compound == "18")
# Use Time zero for minimal value to consider, just using average of medium with no cells. Ctrl is treated cells with control so DMSO
IC50_prop_WT_18 <- convertToProp(IC50_tidy_WT_18$Fluorescence, T0 = mean_medium, Ctrl = max_dmso)
np1 <- nplr(x=IC50_tidy_WT_18$Concentration, y=IC50_prop_WT_18)
plot(np1, pcol="grey40", lcol="skyblue1", showEstim=0.5, showGOF = FALSE, showInfl=FALSE, showCI = FALSE, main="MAZ 18 treated AMO-WT cells", cex.main=1.5,xlab=expression(Log[10](Concentration) (uM)), ylab="Proportion")

