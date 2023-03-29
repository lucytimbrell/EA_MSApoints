# Project: Eastern African Middle Stone Age point variability
# Process: Outline geometric morphometric analysis
# Author: L. J. Timbrell
# Date: 27/03/2023

# packages
library(Momocs)
library(rio)
library(patchwork)
library(ggplot2)
library(tidyverse)


setwd("/YOUR/WORKING/DIRECTORY")

# Import data
tpslines <- import_tps("East_data.tps")
database <- read.csv("East_spreadsheet.csv")

shape <- Out(tpslines$coo, fac = database)
names(shape) <- database$Edited.photo.name

panel(shape, main = "Outline data") # Visualization of points in their original orientation

shapenorm <- shape %>% 
  coo_centre() %>% 
  coo_scale() %>% 
  coo_slidedirection("right") %>% 
  coo_close() 

panel(shapenorm, main = "Normalised outline data") # Visualization of points in their original orientation
stack(shapenorm)

saveRDS(shapenorm, file = "Normalised_outlinesEMSA.rds") # save normalised and transformed landmarks

#### Reload and subset the data  ####

outlines <- import("Normalised_outlinesEMSA.rds")
database <- import("East_MSA_database.rds")  

outlines$fac$Site <- as.factor(outlines$fac$Site)
panel(outlines, main = "Outline data", fac ="Site") 

#### Filter out tip damage ####
complete_outlines <- Momocs:: filter(outlines, Preservation=="Complete")
database <- Momocs:: filter(database, Preservation=="Complete")


#### Filter out Acheulean assemblages from Kaputhurin ####
complete_outlines1 <- Momocs:: filter(complete_outlines, !Assemblage.layer.no.=="GnJh 16")
database1 <- Momocs:: filter(database, !Assemblage.layer.no.=="GnJh 16")

complete_outlines1 <- Momocs:: filter(complete_outlines1,  !Assemblage.layer.no.=="GnJh 17")
database1 <- Momocs:: filter(database1,  !Assemblage.layer.no.=="GnJh 17")

complete_outlines1 <- Momocs:: filter(complete_outlines1,  !Assemblage.layer.no.=="GnJh 51")
database1 <- Momocs:: filter(database1,  !Assemblage.layer.no.=="GnJh 51")

complete_outlines1 <- Momocs:: filter(complete_outlines1,  !Assemblage.layer.no.=="GnJh 81")
database1 <- Momocs:: filter(database1,  !Assemblage.layer.no.=="GnJh 81")

complete_outlines1 <- Momocs:: filter(complete_outlines1,  !Assemblage.layer.no.=="GnJh40")
database1 <- Momocs:: filter(database1,  !Assemblage.layer.no.=="GnJh40")

complete_outlines1 <- Momocs:: filter(complete_outlines1,  !Assemblage_2=="GsJi 13")
database1 <- Momocs:: filter(database1,  !Assemblage_2=="GsJi 13")

table(complete_outlines1$fac$Assemblage.layer.no.)
table(complete_outlines1$fac$Assemblage_2)

## Calculate mean number of points representing outlines
meanpoints <- matrix(0, nrow = length(complete_outlines1), ncol = 1)
for(i in 1:length(complete_outlines1)){
  artefact <- complete_outlines[i]
  nlandmarks <- length(unlist(artefact))/2
  meanpoints[i,] <- nlandmarks
}

mean(meanpoints)


#  Data tidying

complete_outlines1$fac$Raw.material <- as.factor(complete_outlines1$fac$Raw.material)
complete_outlines1$fac$Raw.material[complete_outlines1$fac$Raw.material=="Basalt "] <- "Basalt"
complete_outlines1$fac$Raw.material[complete_outlines1$fac$Raw.material=="Basalt?"] <- "Basalt"
complete_outlines1$fac$Raw.material[complete_outlines1$fac$Raw.material=="Chert "] <- "Chert"
complete_outlines1$fac$Raw.material[complete_outlines1$fac$Raw.material=="Chert  "] <- "Chert"
complete_outlines1$fac$Raw.material[complete_outlines1$fac$Raw.material=="Obsidian "] <- "Obsidian"
complete_outlines1$fac$Raw.material[complete_outlines1$fac$Raw.material=="Rhyolite?"] <- "Rhyolite"


complete_outlines1$fac$Assemblage_2 <- as.factor(complete_outlines1$fac$Assemblage_2)
complete_outlines1$fac$Assemblage_2[complete_outlines1$fac$Assemblage_2=="Complex II, Iic-a"] <- "Complex II, IIc-a"
complete_outlines1$fac$Assemblage_2[complete_outlines1$fac$Assemblage_2=="Complex II, Iid-f"] <- "Complex II, IId-f"
complete_outlines1$fac$Assemblage_2[complete_outlines1$fac$Assemblage_2=="06N-14W iiic"] <- "06N-14W IIIc"
complete_outlines1$fac$Assemblage_2[complete_outlines1$fac$Assemblage_2=="08N-07W iiib"] <- "08N-07W IIIb"
complete_outlines1$fac$Assemblage_2[complete_outlines1$fac$Assemblage_2=="GnJh-15 "] <- "GnJh-15"

complete_outlines1$fac$Raw.material <- droplevels(complete_outlines1$fac$Raw.material)
table(complete_outlines1$fac$Raw.material) 

complete_outlines1$fac$Assemblage_2 <- droplevels(complete_outlines1$fac$Assemblage_2)
table(complete_outlines1$fac$Assemblage_2)

database1 <- complete_outlines1$fac

write.csv(complete_outlines1$fac, "final_easternafrican_msa_points.csv")

#### EFA
calibrate_harmonicpower_efourier(complete_outlines1, nb.h = 20, plot = FALSE) # 11 harmonics

efashape <- efourier(complete_outlines1, nb.h = 10, norm = FALSE)

####  PCA  ####

pcashape <- PCA(efashape) 


# Plot

p1 <- scree_plot(pcashape, nax =1:7) # PC1-7 gives over 95% of cum variance in the data, PC1 = 63% of variance
p1 <- p1  + theme_minimal()
p1

plot.new()
gg <- PCcontrib(pcashape, nax = 1:6, plot = FALSE)
gg$gg + 
  geom_polygon(fill="gray", col="black") 


## Build new database with PCs and centroid size
tidy.table <- cbind(as.tibble(database1), as.tibble(pcashape$x[,1:8]))
tidy.table$Site <- as.factor(tidy.table$Site)

centroidsize <- as_tibble(coo_centsize(complete_outlines1))
centroidsize <- rename(centroidsize, CS = "value")

tidy.table <- cbind(tidy.table, centroidsize)
write.csv(tidy.table, "East_analysis_shapedata.csv")


pcascores <- as_tibble(pcashape$x)
databasedata <- cbind(complete_outlines1$fac, centroidsize, pcascores) # new database with PCs and centroid size

## PC1
library(ggpubr)
p1 <- ggscatter(databasedata, x = "PC1", y = "CS", 
                add = "reg.line", conf.int = TRUE,
                 cor.method = "pearson",
                xlab = "PC1", ylab = "Centroid size")

cor.test(databasedata$PC1, databasedata$CS)
cor(databasedata$PC1, databasedata$CS)

summary(lm(CS~PC1, data = databasedata))

## PC2
p2 <- ggscatter(databasedata, x = "PC2", y = "CS", 
                add = "reg.line", conf.int = TRUE,
                cor.method = "pearson",
                xlab = "PC2", ylab = "Centroid size")

cor.test(databasedata$PC2, databasedata$CS)
cor(databasedata$PC2, databasedata$CS)

summary(lm(CS~PC2, data = databasedata))

## PC3
p3 <- ggscatter(databasedata, x = "PC3", y = "CS", 
                add = "reg.line", conf.int = TRUE, 
                 cor.method = "pearson",
                xlab = "PC3", ylab = "Centroid size")

cor.test(databasedata$PC3, databasedata$CS)
cor(databasedata$PC3, databasedata$CS)

summary(lm(CS~PC3, data = databasedata))

(p1| p2| p3)+
  plot_annotation(title = "Allometry", tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 24))) 



## PC4
p4 <- ggscatter(databasedata, x = "PC4", y = "CS", 
                add = "reg.line", conf.int = TRUE,
                cor.method = "pearson",
                xlab = "PC4", ylab = "Centroid size")

cor.test(databasedata$PC4, databasedata$CS)
cor(databasedata$PC4, databasedata$CS)

summary(lm(CS~PC4, data = databasedata))

## PC2
p5 <- ggscatter(databasedata, x = "PC5", y = "CS", 
                add = "reg.line", conf.int = TRUE,
                cor.method = "pearson",
                xlab = "PC5", ylab = "Centroid size")

cor.test(databasedata$PC5, databasedata$CS)
cor(databasedata$PC5, databasedata$CS)

summary(lm(CS~PC5, data = databasedata))

## PC3
p6 <- ggscatter(databasedata, x = "PC6", y = "CS", 
                add = "reg.line", conf.int = TRUE, label = "Edited.photo.name",
                cor.method = "pearson",
                xlab = "PC6", ylab = "Centroid size")

cor.test(databasedata$PC6, databasedata$CS)
cor(databasedata$PC6, databasedata$CS)

summary(lm(CS~PC6, data = databasedata))

(p1| p2| p3)/(p4| p5| p6)+
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 24))) 



## Raw materials

a <- ggplot(tidy.table, aes(PC1, PC2, colour = Raw.material)) + 
  geom_point(size = 2) +
  stat_ellipse() +
  theme_minimal()+ 
  theme(legend.position = "none")

b <- ggplot(tidy.table, aes(PC1, PC3, colour = Raw.material)) + 
  geom_point(size = 2) +
  stat_ellipse() +
  theme_minimal()+ 
  theme(legend.position = "none")

c <- ggplot(tidy.table, aes(PC2, PC3, colour = Raw.material)) + 
  geom_point(size = 2) +
  stat_ellipse() +
  theme_minimal()+ 
  theme(legend.position = "none")


d <- tidy.table %>% 
  dplyr::select(PC1, PC2, PC3, Raw.material) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Raw.material = as.factor(Raw.material)) %>%
  ggplot(aes(name, value, fill = Raw.material)) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

(a | b | c) / d +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 30))) 

a <- ggplot(tidy.table, aes(PC4, PC5, colour = Raw.material)) + 
  geom_point(size = 2) +
  stat_ellipse() +
  theme_minimal()+ 
  theme(legend.position = "none")

b <- ggplot(tidy.table, aes(PC4, PC6, colour = Raw.material)) + 
  geom_point(size = 2) +
  stat_ellipse() +
  theme_minimal()+ 
  theme(legend.position = "none")

c <- ggplot(tidy.table, aes(PC5, PC6, colour = Raw.material)) + 
  geom_point(size = 2) +
  stat_ellipse() +
  theme_minimal()+ 
  theme(legend.position = "none")


d <- tidy.table %>% 
  dplyr::select(PC4, PC5, PC6, Raw.material) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Raw.material = as.factor(Raw.material)) %>%
  ggplot(aes(name, value, fill = Raw.material)) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

(a | b | c) / d +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 30))) 

## LDA
table(tidy.table$Raw.material)

dashape95 <- LDA(pcashape, ~Raw.material, prior = c(46,61, 1, 84, 1, 3)/nrow(tidy.table), retain = 0.90, cv = TRUE) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

## ANOVA and MANOVA -
pcashape %>% MANOVA(~Raw.material, retain = 0.95)

res <- aov(PC1~Raw.material, tidy.table)
res1 <- aov(PC2~Raw.material, tidy.table)
res2 <- aov(PC3~Raw.material, tidy.table)
res3 <- aov(PC4~Raw.material, tidy.table)
res4 <- aov(PC5~Raw.material, tidy.table)
res5 <- aov(PC6~Raw.material, tidy.table)
res6 <- aov(PC7~Raw.material, tidy.table)
res7 <- aov(CS~Raw.material, databasedata)


TukeyHSD(res)
TukeyHSD(res1)
TukeyHSD(res2)
TukeyHSD(res3)
TukeyHSD(res4)
TukeyHSD(res5)
TukeyHSD(res6)
TukeyHSD(res7)

n = 6
cols = gg_color_hue(n)

p <- CLUST(pcashape, ~Raw.material, dist_method = "euclidean", retain = 0.95, hclust_method = "complete", type = "horizontal", palette = pal_manual(cols))
p + theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.title.x = element_blank())

## Site
tidy.table$Site <- as.factor(tidy.table$Site)

a <- ggplot(tidy.table, aes(PC1, PC2, colour = Site)) + 
  geom_point(size = 2) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888", "#F8766D"))+
  theme_minimal()+ 
  theme(legend.position = "none")

b <- ggplot(tidy.table, aes(PC1, PC3, colour = Site)) + 
  geom_point(size = 2) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888", "#F8766D"))+
  theme_minimal()+ 
  theme(legend.position = "none")

c <- ggplot(tidy.table, aes(PC2, PC3, colour = Site)) + 
  geom_point(size = 2) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888", "#F8766D"))+
  theme_minimal()+ 
  theme(legend.position = "none")


d <- tidy.table %>% 
  dplyr::select(PC1, PC2, PC3, Site) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Site = as.factor(Site)) %>%
  mutate(Site = fct_relevel(Site, "Goda Buticha ", "Porc Epic", "Kapthurin Formation", "Prolonged Drift", "Prospect Farm", "Omo Kibish")) %>%
  ggplot(aes(name, value, fill = Site)) +
  scale_fill_manual(guide = guide_legend(nrow = 2, byrow = TRUE), limits = c("Goda Buticha ", "Porc Epic", "Kapthurin Formation", "Prolonged Drift", "Prospect Farm", "Omo Kibish"),values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888", "#F8766D"))+
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

  (a | b | c) / d +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.title = element_text(size = 30))) 

a <- ggplot(tidy.table, aes(PC4, PC5, colour = Site)) + 
  geom_point(size = 2) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888", "#F8766D"))+
  theme_minimal()+ 
  theme(legend.position = "none")

b <- ggplot(tidy.table, aes(PC4, PC6, colour = Site)) + 
  geom_point(size = 2) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888", "#F8766D"))+
  theme_minimal()+ 
  theme(legend.position = "none")

c <- ggplot(tidy.table, aes(PC5, PC6, colour = Site)) + 
  geom_point(size = 2) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888", "#F8766D"))+
  theme_minimal()+ 
  theme(legend.position = "none")


d <- tidy.table %>% 
  dplyr::select(PC4, PC5, PC6, Site) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Site = as.factor(Site)) %>%
  mutate(Site = fct_relevel(Site, "Goda Buticha ", "Porc Epic", "Kapthurin Formation", "Prolonged Drift", "Prospect Farm", "Omo Kibish")) %>%
  ggplot(aes(name, value, fill = Site)) +
  scale_fill_manual(guide = guide_legend(nrow = 2, byrow = TRUE), limits = c("Goda Buticha ", "Porc Epic", "Kapthurin Formation", "Prolonged Drift", "Prospect Farm", "Omo Kibish"),values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888", "#F8766D"))+
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

(a | b | c) / d +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 30))) 
## LDA
table(tidy.table$Site)

dashape95 <- LDA(pcashape, ~Site, prior = c(32, 22, 12, 73, 8, 56)/nrow(tidy.table), retain = 0.90, cv = TRUE) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

## ANOVA and MANOVA -
pcashape %>% MANOVA(~Site, retain = 0.95)

res <- aov(PC1~Site, tidy.table)
res1 <- aov(PC2~Site, tidy.table)
res2 <- aov(PC3~Site, tidy.table)
res3 <- aov(PC4~Site, tidy.table)
res4 <- aov(PC5~Site, tidy.table)
res5 <- aov(PC6~Site, tidy.table)
res6 <- aov(PC7~Site, tidy.table)


TukeyHSD(res)
TukeyHSD(res1)
TukeyHSD(res2)
TukeyHSD(res3)
TukeyHSD(res4)
TukeyHSD(res5)


n = 6
cols = gg_color_hue(n)

p <- CLUST(pcashape, ~Site, dist_method = "euclidean", retain = 0.90, hclust_method = "complete", type = "horizontal", palette = pal_manual(cols))
p + theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())


## Assemblage
tidy.table$Assemblage_2 <- as.factor(tidy.table$Assemblage_2)
tidy.table$Assemblage_2[tidy.table$Assemblage_2==""] <- NA

a <- ggboxplot(subset(tidy.table, !is.na(Assemblage_2)), "Site", "PC1", fill = 'Assemblage_2',xlab = "", 
               bxp.errorbar = TRUE)+ theme(legend.position = "none")
b <- ggboxplot(subset(tidy.table, !is.na(Assemblage_2)), "Site", "PC2", fill = 'Assemblage_2',xlab = "", 
               bxp.errorbar = TRUE)+ theme(legend.position = "none")
c <- ggboxplot(subset(tidy.table, !is.na(Assemblage_2)), "Site", "PC3", fill = 'Assemblage_2',xlab = "", 
               bxp.errorbar = TRUE)+ theme(legend.position = "none")

d <- ggboxplot(subset(tidy.table, !is.na(Assemblage_2)), "Site", "PC4", fill = 'Assemblage_2', xlab = "", 
               bxp.errorbar = TRUE)+ theme(legend.position = "none")
e <- ggboxplot(subset(tidy.table, !is.na(Assemblage_2)), "Site", "PC5", fill = 'Assemblage_2',xlab = "", 
               bxp.errorbar = TRUE)+ theme(legend.position = "none")
f <- ggboxplot(subset(tidy.table, !is.na(Assemblage_2)), "Site", "PC6", fill = 'Assemblage_2',xlab = "", 
               bxp.errorbar = TRUE)+ theme(legend.position = "bottom")

a/b/c/d/e/f +
  plot_annotation(theme = theme(plot.title = element_text(size = 20))) 


pcashape %>% MANOVA(~Assemblage_2, retain = 0.95)

res <- aov(PC1~Assemblage_2, tidy.table)
res1 <- aov(PC2~Assemblage_2, tidy.table)
res2 <- aov(PC3~Assemblage_2, tidy.table)
res3 <- aov(PC4~Assemblage_2, tidy.table)
res4 <- aov(PC5~Assemblage_2, tidy.table)
res5 <- aov(PC6~Assemblage_2, tidy.table)
res6 <- aov(CS~Assemblage_2, tidy.table)


t1 <- TukeyHSD(res)
t2 <- TukeyHSD(res1)
t3 <- TukeyHSD(res2)
t4 <- TukeyHSD(res3)
t5<- TukeyHSD(res4)
t6 <- TukeyHSD(res5)
t7 <- TukeyHSD(res6)


p1 <- t1$Assemblage_2[,4]
p2 <- t2$Assemblage_2[,4]
p3 <- t3$Assemblage_2[,4]
p4 <- t4$Assemblage_2[,4]
p5 <- t5$Assemblage_2[,4]
p6 <- t6$Assemblage_2[,4]
p7 <- t7$Assemblage_2[,4]

ptable <- cbind(p1, p2, p3, p4, p5, p6, p7)
colnames(ptable) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "CS")

write.csv(ptable, "pcs_tukey_results.csv")
library(tidyverse)
library(rstatix)


mean_CS <- tidy.table %>%
  group_by(Assemblage_2) %>%
  get_summary_stats(CS, type = "mean_sd")

meanshape <- MSHAPES(efashape, ~ Assemblage_2)

p <- plot_MSHAPES(meanshape, size = 1,  palette = pal_manual(c("blue", "red")))
p + theme(plot.margin = unit(c(1,1,1,1), "cm"))

meanshape <- MSHAPES(pcashape, ~Assemblage_2)

test <- MSHAPES(efashape, ~Assemblage_2)
ms <- test$shp
par(mfrow=c(3,5))
coo_plot(ms$`06N-14W 60-200`, col = "light blue", main = "06N-14W 60-200cm")
coo_plot(ms$`08N-07W IIIa`, col = 'light blue', main = "08N-07W IIIa")
coo_plot(ms$`08N-07W  IIIb`, col = "light blue", main = "08N-07W  IIIb")
coo_plot(ms$`08N-07W  IIIc`, col = "light blue", main = "08N-07W  IIIc")
coo_plot(ms$`08N-07W IIId`, col = "light blue", main = "08N-07W IIId")
coo_plot(ms$'GsJi 7',  col = "pink", main = "GsJi 7")
coo_plot(ms$'GsJi 8', col = "pink", main = "GsJi 8")
coo_plot(ms$AHS, col = "dark green", main = "AHS")
coo_plot(ms$` BNS `, col = "dark green", main = "BNS")
coo_plot(ms$`Complex II, IIc-a`, col = "red", main = "Complex II, IIc-a")
coo_plot(ms$`Complex II, IId-f`, col = "red", main = "Complex II, IId-f")
coo_plot(ms$`GnJh 78`, col = "gold", main = "GnJh 78")
coo_plot(ms$`GnJh-15`, col = "gold", main = "GnJh 15")
coo_plot(ms$`GrJi 11`, col = "darkviolet", main = "GrJi 11")

write.csv(meanshape, "final_meanshape_pcs.csv")

table(tidy.table$Assemblage_2, tidy.table$Raw.material)

KMEDOIDS(pcashape, k = 2, metric = "euclidean")

n = 14
cols = gg_color_hue(n)CLUST(pcashape, ~Assemblage_2,lwd = 0.75,cex =0.5,  dist_method = "euclidean", retain = 0.90, hclust_method = "complete", type = "horizontal", palette = pal_manual(cols))

 pcs <- tidy.table[,29:34]
 row.names(pcs) <- tidy.table$Edited.photo.name
 
 d <- dist(pcs, method = "euclidean") 
 fit <- hclust(d, method="complete") 
 plot(fit) # display dendogram
 groups <- cutree(fit, k=6) # cut tree into 5 clusters
 # draw dendogram with red borders around the 5 clusters 
 rect.hclust(fit, k=6, border="red")

 groups <- cutree(fit, k=14) # cut tree into 5 clusters
 tidy.table$cluster2 <- groups
 table(tidy.table$cluster2, tidy.table$Assemblage_2)
 
 groups <- cutree(fit, k=14) # cut tree into 5 clusters
 tidy.table$cluster <- groups
 table(tidy.table$cluster, tidy.table$Assemblage_2)
 
 library(ggpubr)
 p1 <- ggscatter(final_dataset, x = "PC2", y = "Mean_precipiation", 
                 add = "reg.line", conf.int = TRUE,
                 cor.method = "pearson",
                 xlab = "PC1", ylab = "Precipitation")
 
 cor.test(databasedata$PC1, databasedata$CS)
 cor(databasedata$PC1, databasedata$CS)
 
 summary(lm(Mid.age~PC1, data = tidy.table))