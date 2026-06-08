
# NOTES
# the script prepares biodiversity data into a species-by-plot abundance matrix. 
# Calculates taxonomic and functional diversity for each assemblage with iNEXT.3D. 
# Functional diversity estimation uses a distance matrix.

# 8 Jun 2026
# Code: M. Paola Barajas Barbosa: paola.barajas@ufz.de, paola.barajas@idiv.de

# 1. Load packages -------------------------------------------------------------

library(dplyr)
library(tidyr)
library(tibble)
library(cluster)
library(iNEXT.3D)
library(ggplot2)

# 2. Set working directory and load data ---------------------------------------

setwd("G:/My Drive/_Teaching/2026_Class_functional ecology/MSc/Projects")

load("datt_20230419.RData")
load("data_Trait_imputed_20230418.RData")

# 3. Prepare abundance matrix for iNEXT.3D -------------------------------------

# Minimum number of individuals required per plot
min_ind <- 15

# Summarise abundance for each species within each plot
tmp <- datt %>%
  group_by(Scientific_name, PlotID) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") # Result: one row per species-plot combination

# Convert data from long to wide format:
Abun_p <- tmp %>%
  pivot_wider(           # rows = species, columns = plots, values = abundances
    names_from = PlotID,
    values_from = Abundance,
    values_fill = 0
  ) %>%
  column_to_rownames("Scientific_name")

# iNEXT.3D requires a numeric matrix
Abun_p <- as.matrix(Abun_p)
storage.mode(Abun_p) <- "numeric"

# Remove plots with too few individuals
Abun_p <- Abun_p[, colSums(Abun_p) > min_ind, drop = FALSE]

# Remove plots with fewer than 5 observed species
Abun_p <- Abun_p[, colSums(Abun_p > 0) >= 5, drop = FALSE] # # This is important because iNEXT.3D needs sufficient species information

# 4. Estimate taxonomic diversity (TD) -----------------------------------------

# function help page 
?estimate3D

# TD = taxonomic diversity
# q = Hill numbers of order 0, 1, and 2
# base = "coverage" means standardization by sample coverage
# level = target coverage
# nboot = number of bootstrap replications for uncertainty estimates
est_TD <- estimate3D(
  data = Abun_p,
  diversity = "TD",
  q = c(0, 1, 2),
  datatype = "abundance",
  base = "coverage",
  level = 0.8,
  nboot = 10
)

# Extract the taxonomic diversity estimates in wide format
est_out_TD <- est_TD %>%
  pivot_wider(
    id_cols = Assemblage,
    names_from = Order.q,
    values_from = qTD    )

est_out_TD <- data.frame(
  PlotID = as.factor(est_out_TD$Assemblage) ,
  TD_0 = est_out_TD$`0` ,
  TD_1 = est_out_TD$`1`,
  TD_2 = est_out_TD$`2`  )


# 5. Prepare trait data for functional diversity (FD) --------------------------

# Keep only species that are present in the abundance matrix
species_keep <- data.frame(Scientific_name = rownames(Abun_p))

traits <- species_keep %>%
  left_join(impute_out, by = "Scientific_name")

# Set species names as row names
rownames(traits) <- traits$Scientific_name
traits$Scientific_name <- NULL

# Compute euclidean distance among species based on traits
distM <- cluster::daisy(x = traits, metric = "euclidean") %>%
  as.matrix()

# 6. Estimate functional diversity (FD) ----------------------------------------

est_FD <- estimate3D(
  data = Abun_p,
  diversity = "FD",
  q = c(0, 1, 2),
  datatype = "abundance",
  base = "coverage",
  level = 0.8,
  nboot = 10,
  FDdistM = distM  )

est_out_FD <- est_FD %>%
  pivot_wider(
    id_cols = Assemblage,
    names_from = Order.q,
    values_from = qFD    )

est_out_FD <- data.frame(
  PlotID = as.factor(est_out_FD$Assemblage) ,
  FD_0 = est_out_FD$`0` ,
  FD_1 = est_out_FD$`1`,
  FD_2 = est_out_FD$`2`  )

est_out_FD_TD <- left_join(est_out_FD, est_out_TD, by = "PlotID")

# 7. Relating diversity estimates to plot attributes ---------------------------

# Summary of plot attributes to relate them to plot level diversity 

names(datt)

plot_attributes <- datt %>%
  group_by(geo_entity2, PlotID) %>%
  summarise(Elev_m  = mean(Elev_m), 
            MAT     = mean(MAT),                 
            MAP     = mean(MAP),                   
            PET     = mean(PET),                  
            AridInd = mean(AridInd),               
            SubstrateAge_range =mean(SubstrateAge_range), 
            .groups = "drop") # Result: one row per species-plot combination

plot_attributes$PlotID <- as.factor(plot_attributes$PlotID)

plot_div_attrib <- left_join( est_out_FD_TD, plot_attributes, by = "PlotID")

coco = c(
  "#00AFBB", 
  "gray60", 
  "#E7B800")

# Example: diversity vs MAP

ggplot(plot_div_attrib, aes(x = MAP, y = FD_0)) +
  geom_point(aes(color = geo_entity2), size = 1.5, alpha = 0.7) +
  geom_smooth(
    aes(color = geo_entity2, fill = geo_entity2),
    method = glm,
    method.args = list(family = Gamma(link = "log")),
    se = TRUE,
    linewidth = 0.8
  ) +
  scale_color_manual(values = coco) +
  scale_fill_manual(values = coco) +
  labs(
    x = "Mean annual precipitation (MAP)",
    y = "Functional diversity (q = 0)"
  ) +
  theme_minimal(base_size = 8)

# Example: diversity vs Substrate Age 

ggplot(plot_div_attrib, aes(x = SubstrateAge_range, y = FD_0)) +
  geom_point(aes(color = geo_entity2), size = 1.5, alpha = 0.7) +
  geom_smooth(
    aes(color = geo_entity2, fill = geo_entity2),
    method = glm,
    method.args = list(family = Gamma(link = "log")),
    se = TRUE,
    linewidth = 0.8
  ) +
  scale_color_manual(values = coco) +
  scale_fill_manual(values = coco) +
  labs(
    x = "Mean annual precipitation (MAP)",
    y = "Functional diversity (q = 0)"
  ) +
  theme_minimal(base_size = 8)
