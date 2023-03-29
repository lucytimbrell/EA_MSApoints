# Project: Eastern African Middle Stone Age point variability
# Process: Metric analysis
# Author: L. J. Timbrell
# Date: 27/03/2023

setwd("/Volumes/Middle Stone Age points project/East outline analysis")

database <- read.csv("final_easternafrican_msa_points.csv") # this file is created in the GMM script

library(tidyverse)
library(ggpubr)

# Summary statistics for raw materials

length_sum <- database %>%
  group_by(Raw.material) %>%
  get_summary_stats(Length..mm., type = "mean_sd")

width_sum <- database%>%
  group_by(Raw.material) %>%
  get_summary_stats(Width..mm., type = "mean_sd")

thickness_sum <- database %>%
  group_by(Raw.material) %>%
  get_summary_stats(Thickness..mm., type = "mean_sd")

# Boxplots

length_bxp <- ggboxplot(database,
                        x = "Raw.material",
                        y = "Length..mm.", 
                        xlab = "",
                        ylab = "Length (mm)",
                        bxp.errorbar = TRUE)
length_bxp

width_bxp <- ggboxplot(database,
                       x = "Raw.material",
                       y = "Width..mm.", 
                       xlab = "",
                       ylab = "Width (mm)",
                       bxp.errorbar = TRUE)
width_bxp

thickness_bxp <- ggboxplot(database,
                           x = "Raw.material",
                           y = "Thickness..mm.", 
                           xlab = "Raw material",
                           ylab = "Thickness (mm)",
                           bxp.errorbar = TRUE)

thickness_bxp

length_bxp/width_bxp/thickness_bxp +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 20))) 
# ANOVA and Tukey HSD
res <- aov(Length..mm.~factor(Raw.material), database)
summary(res)
TukeyHSD(res)

res1 <- aov(Width..mm. ~ factor(Raw.material), database)
summary(res1)
TukeyHSD(res1)

res2 <- aov(Thickness..mm. ~ factor(Raw.material), database)
summary(res2)
TukeyHSD(res2)

# Summary statistics for site

length_sum <- database %>%
  group_by(Site) %>%
  get_summary_stats(Length..mm., type = "mean_sd")


width_sum <- database%>%
  group_by(Site) %>%
  get_summary_stats(Width..mm., type = "mean_sd")


thickness_sum <- database %>%
  group_by(Site) %>%
  get_summary_stats(Thickness..mm., type = "mean_sd")

# Boxplots
cols = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888", "#F8766D")

length_bxp <- ggboxplot(database,
                        x = "Site",
                        y = "Length..mm.", 
                        xlab = "",
                        ylab = "Length (mm)",
                        fill  = "Site",
                        palette = cols,
                        bxp.errorbar = TRUE)+ theme(legend.position = "none")
length_bxp

width_bxp <- ggboxplot(database,
                       x = "Site",
                       y = "Width..mm.", 
                       xlab = "",
                       ylab = "Width (mm)",
                       fill  = "Site",
                       palette = cols,
                       bxp.errorbar = TRUE)+ theme(legend.position = "none")
width_bxp

thickness_bxp <- ggboxplot(database,
                           x = "Site",
                           y = "Thickness..mm.", 
                           ylab = "Thickness (mm)",
                           fill  = "Site",
                           palette = cols,
                           bxp.errorbar = TRUE)+ theme(legend.position = "bottom")

thickness_bxp

length_bxp/width_bxp/thickness_bxp +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 20))) 


# ANOVA and Tukey HSD
res <- aov(Length..mm.~factor(Site), database)
summary(res)
TukeyHSD(res)

res1 <- aov(Width..mm. ~ factor(Site), database)
summary(res1)
TukeyHSD(res1)

res2 <- aov(Thickness..mm. ~ factor(Site), database)
summary(res2)
TukeyHSD(res2)

# Summary statistics for assemblage

length_sum <- database %>%
  group_by(Assemblage_2) %>%
  get_summary_stats(Length..mm., type = "mean_sd")

width_sum <- database%>%
  group_by(Assemblage_2) %>%
  get_summary_stats(Width..mm., type = "mean_sd")

thickness_sum <- database %>%
  group_by(Assemblage_2) %>%
  get_summary_stats(Thickness..mm., type = "mean_sd")

# Boxplots
database <- as.tibble(database)
database$Assemblage_2[database$Assemblage_2==""] <- NA

database$Assemblage <- database$Assemblage_2

a <- ggboxplot(subset(database,!is.na(Assemblage_2)), "Site", "Length..mm.", fill = 'Assemblage',  ylab = "Length (mm)",
               bxp.errorbar = TRUE)+ theme(legend.position = "none")
b <- ggboxplot(subset(database, !is.na(Assemblage_2)), "Site", "Width..mm.", fill = 'Assemblage',ylab = "Width (mm)",
               bxp.errorbar = TRUE)+ theme(legend.position = "none")
c <- ggboxplot(subset(database, !is.na(Assemblage_2)), "Site", "Thickness..mm.", fill = 'Assemblage',  ylab = "Thickness (mm)",
               bxp.errorbar = TRUE)+ theme(legend.position = "bottom")

a/b/c +
  plot_annotation(theme = theme(plot.title = element_text( size = 20))) 

# ANOVA and Tukey HSD
res <- aov(Length..mm.~factor(Assemblage_2), database)
summary(res)
thd1 <- TukeyHSD(res)

p1 <- thd1$`factor(Assemblage_2)`[,4]

res1 <- aov(Width..mm. ~ factor(Assemblage_2), database)
summary(res1)
thd2 <- TukeyHSD(res1)

p2 <- thd2$`factor(Assemblage_2)`[,4]

res2 <- aov(Thickness..mm. ~ factor(Assemblage_2), database)
summary(res2)
thd3 <- TukeyHSD(res2)

p3 <- thd3$`factor(Assemblage_2)`[,4]

#Create table of results
ptable <- cbind(p1, p2, p3)
colnames(ptable) <- c("Length","Width", "Thickness")

# Write files
write.csv(ptable, "metric_tukey_results.csv")
write.csv(length_sum, "length_summary_assemblages.csv")
write.csv(width_sum, "width_summary_assemblages.csv")
write.csv(thickness_sum, "thickness_summary_assemblages.csv")


