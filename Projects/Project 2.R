# NOTES
# A
# B

# 8 Jun 2026
# Code: Paola Barajas: paola.barajas@ufz.de, paola.barajas@idiv.de

R.Version()   #  R version 4.5.2 (2025-10-31 ucrt)

# 1. Load packages -------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggExtra)

library(hypervolume)
library(BAT)
library(fundiversity)

# 2. Set working directory and load data ---------------------------------------

setwd("G:/My Drive/_Teaching/2026_Class_functional ecology/MSc/Projects")
data_Tenerife <- read.csv("Barajasetal2023_data_22052023.csv")

# 3. Visualizing trait data ---------------------------------------

# Tenerife data

data_Tenerife %>% group_by(Biogeo_status) %>% summarise(n = n())
data_Tenerife %>% group_by(Ende_status) %>% summarise(n = n())

data_Tenerife %>% group_by(Ende_status) %>% summarise(n = n())


# # PCA 
PCA             <- prcomp(data_Tenerife[, c("Leaf_area","LMA","Leaf_N","Leaf_th","LDMC","Seed_mass","Stem_density","Height")] )
summary(PCA)

PCAvalues       <- data.frame(Species = data_Tenerife$Species1, Biogeo_status = data_Tenerife$Biogeo_status, 
                              Ende_status = data_Tenerife$Ende_status, PCA$x)
PCAvalues$PC1   <- PCAvalues$PC1*-1 
#PCAvalues$PC2 <- PCAvalues$PC2*-1 # for visualization purposes

PCAloadings     <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)  
PCAloadings$PC1 <- PCAloadings$PC1*-1 
#PCAloadings$PC2 <- PCAloadings$PC2*-1

# # 
# 1. Plots Figure 1  
Tenerife <- PCAvalues %>%
  ggplot(aes(PC1, PC2)) +                               # Plot Tenerife PCA
  stat_density_2d(
    geom = "polygon",
    contour = TRUE,
    aes(fill = after_stat(level)),
    colour = "gray",
    bins = 34  ) +
  scale_fill_distiller(palette = "BuGn", direction = 1) +
  geom_jitter(alpha = 0.6, size = 2, colour = "turquoise4") +  # Display species as points
  geom_text(
    data = PCAloadings,
    aes(x = PC1 * 4.7, y = PC2 * 4.7, label = Variables),
    size = 2.3  ) +
  geom_segment(
    data = PCAloadings,
    linewidth = 0.2,                                     
    aes(x = 0, xend = PC1 * 4.2, y = 0, yend = PC2 * 4.2),
    arrow = arrow(length = unit(0.1, "cm")),
    colour = "black"  ) +
  xlab("PC1") + ylab("PC2 (25%)") +
  xlim(-5, 5) + ylim(-5, 4) +
  theme_minimal()

Tenerife
  
# Plot each group separately
NNE <- dplyr::filter(PCAvalues, Biogeo_status == "NNE")
NNE<- NNE  %>% ggplot(aes(PC1,  PC2))+ 
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 10) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5,  size = 2  , colour = "black") +    
  xlim(-5 , 5) + ylim(-5, 4) +
  xlab("n = 54") + ylab("")   + theme_minimal()
NNE

TE <- dplyr::filter(PCAvalues, Biogeo_status == "TE")
TE <-  TE  %>% ggplot(aes(PC1,  PC2))+  
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 15) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5,  size =2   , colour = "dodgerblue3") +     
  xlim(-5 , 5) + ylim(-5, 4) +
  xlab("PC1 (30%) n = 85") + ylab("PC2 (25%)")   + theme_minimal()
TE

CLADO <- dplyr::filter(PCAvalues, Ende_status == "Cla")
CLADO <- CLADO %>% ggplot(aes(PC1,  PC2))+ 
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 20) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5, size = 2  , colour = "mediumorchid2") +  
  xlim(-5 , 5) + ylim(-5, 4) +
  xlab("PC1 (30%) n = 195") + ylab("")   +theme_minimal()
CLADO

Ten<- ggExtra:: ggMarginal(TENERIFE, type = "density", fill="transparent", size = 15)
ns <- ggExtra:: ggMarginal(NNE, type = "density", fill="transparent", size = 15) 
sie<- ggExtra:: ggMarginal(TE, type = "density", fill="transparent", size = 15)
cla<- ggExtra:: ggMarginal(CLADO, type = "density", fill="transparent", size = 15)

ggpubr::ggarrange(Tenerife, ns,sie, cla, ncol =2, nrow = 2) 

# 4. Estimating functional diversity ------------------------------------------

# We calculate the three following indices: kernel.alpha, kernel.dispersion and kernel.evennes. 
# Indices based on hypervolume

?estimate_bandwidth

all_bw <- estimate_bandwidth(PCAvalues[, c("PC1", "PC2", "PC3")], # Estimate bandwidth for the island
                             method = "cross-validation")

ae  <- PCAvalues[which(PCAvalues$Biogeo_status == "CE"),  c("PC1", "PC2", "PC3")] # Archipelago endemic species
mac <- PCAvalues[which(PCAvalues$Biogeo_status == "MAC"), c("PC1", "PC2", "PC3")] # Macaronesia endemic species
ns  <- PCAvalues[which(PCAvalues$Biogeo_status == "NNE"), c("PC1", "PC2", "PC3")] # Non-endemic Native species
te <- PCAvalues[which(PCAvalues$Biogeo_status == "TE"),  c("PC1", "PC2", "PC3")] # Single island endemic species - Tenerife endemics

all_sp <- PCAvalues[, c("PC1", "PC2", "PC3")] # All Tenerife Island species

groups  <- list(AE = ae, MAC = mac, NS = ns, TE = te, Island = all_sp)

hv_list <- list()

hv_list$hv_ae <-  hypervolume::hypervolume_gaussian(
  groups$AE, kde.bandwidth = all_bw,
  quantile.requested = 0.95, quantile.requested.type = "probability",
  verbose = FALSE)

hv_list$hv_mac <-  hypervolume_gaussian(
  groups$MAC, kde.bandwidth = all_bw,
  quantile.requested = 0.95, quantile.requested.type = "probability",
  verbose = FALSE)

hv_list$hv_ns <-  hypervolume_gaussian(
  groups$NS, kde.bandwidth = all_bw,
  quantile.requested = 0.95, quantile.requested.type = "probability",
  verbose = FALSE)

hv_list$hv_te <-  hypervolume_gaussian(
  groups$TE, kde.bandwidth = all_bw,
  quantile.requested = 0.95, quantile.requested.type = "probability",
  verbose = FALSE)

hv_list$hv_island <-  hypervolume_gaussian(
  groups$Island, kde.bandwidth = all_bw,
  quantile.requested = 0.95, quantile.requested.type = "probability",
  verbose = FALSE)

# Check here for an example of a robust hypervolume estimation:
# https://github.com/paolabarajas/Assembly-of-FD-on-the-oceanic-island-Tenerife/blob/main/Code/Figure_3_Richness_Div_Eve.R

# Hypervolume is an S4 object (with slots like @volume). We convert it into df for plotting
hv_metrics_data <- data.frame(
  cat = c("AE", "MAC", "NS", "TE", "Island"),  
  hv = c(
    hv_list$hv_ae@Volume,      
    hv_list$hv_mac@Volume,
    hv_list$hv_ns@Volume,
    hv_list$hv_te@Volume,
    hv_list$hv_island@Volume) )

coco = c("grey","gold2", "mediumorchid2","green4","black","dodgerblue3")

ggplot(hv_metrics_data) +
  geom_point(aes( factor(cat, level = c ("AE", "MAC", "NS", "TE", "Island")  ), 
                  y = hv, color = cat), size = 5) +
  scale_colour_manual(values = coco) +
  theme_minimal()




