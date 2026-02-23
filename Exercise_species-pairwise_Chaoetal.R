# _____________________________________________________________
# Functional diversity Quantification I
# Species-pairwise functional distances sensu Anne Chao et al
# Code: by M. Paola Barajas: paola.barajas@ufz.de -  paolabarajas@gmail.com 
# R version 4.5.2 (2025-10-31 ucrt)
# Last code update: 23 Feb 2026
# _____________________________________________________________

# Data visualization
library(ggplot2)
library(GGally)

# Data handling
library(tidyverse)

# Biodiversity estimation
library(cluster)
library(iNEXT.3D) # https://github.com/KaiHsiangHu/iNEXT.3D

setwd("C:/Users/barajasb/Documents/Class_2026_functional ecology")

# Load data ----
traits <- read.csv("Trait_imputed_144.csv") # Species by trait matrix
Abun_p <- read.csv("Abun_p.csv")            # Species by site matrix - columns are plots

Abun_p <- tibble::column_to_rownames(Abun_p, "Scientific_name")
Abun_p <- Abun_p[, colSums(Abun_p > 0) >= 5 &   # minimum amount of species
                   colSums(Abun_p) >= 20,       # minimum amount of individuals
                          drop = FALSE ]

# 1. Biodiversity Estimation ----

# 1.1 - Taxonomic diversity ----
estimate_TD <- iNEXT.3D::estimate3D(Abun_p,                  # species by site matrix
                                    diversity = 'TD',        # biodiversity facet
                                    q = c(0, 1, 2),          # hill number
                                    datatype = "abundance", 
                                    base = "coverage", 
                                    level = 0.93)

# 1.2 - Functional diversity ----
mask                 <- as.data.frame(rownames(Abun_p))
mask$Scientific_name <- mask$`rownames(Abun_p)`

traits_1 <- mask %>%
  dplyr::left_join(traits, by = "Scientific_name") %>%
  tibble::column_to_rownames("Scientific_name") %>%
  dplyr::select(-`rownames(Abun_p)`)

# - First: Exploring trait data   
GGally::ggpairs(traits_1, columns = 1:4) #  ?GGally

distM <- cluster::daisy(x = traits_1,   # calculate distance matrix
                        metric = "euclidean") %>% as.matrix()

# - Second: Estimating diversity 
estimate_FD <- iNEXT.3D::estimate3D(Abun_p, 
                                    diversity = 'FD', 
                                    q = c(0, 1, 2), 
                                    datatype = 'abundance', 
                                    base = 'coverage',
                                    level = 0.93, 
                                    nboot = 10, 
                                    FDdistM = distM, 
                                    FDtype = "tau_values", 
                                    FDtau = NULL)            # Uses Tau dmean

estimate_FD_min <- iNEXT.3D::estimate3D(Abun_p, 
                                        diversity = 'FD', 
                                        q = c(0, 1, 2), 
                                        datatype = 'abundance', 
                                        base = 'coverage',
                                        level = 0.93, 
                                        nboot = 10, 
                                        FDdistM = distM, 
                                        FDtype = "tau_values", 
                                        FDtau = min(distM) )  # Uses Tau dmin

# Let's see the results
View(estimate_FD) 
str(estimate_FD) 
?estimate3D

plot(estimate_FD$Order.q, estimate_FD$qFD)
plot(estimate_FD_min$Order.q, estimate_FD_min$qFD)
plot(estimate_TD$Order.q, estimate_TD$qTD)

# Plot functional diversity estimation as a function of level of threshold distinctiveness 
estimate_FD_tmp <- rbind(estimate_FD, estimate_FD_min)

ggplot(estimate_FD_tmp, aes( x= Tau , y= qFD ) ) + 
  geom_point()  + 
  xlab("Tau (functional threshold of distinctiveness)") + 
  ylab("Functional Diversity") +
  theme_minimal()

# Plot functional diversity estimation for different hill numbers
ggplot(estimate_FD_tmp, aes( x= Order.q , y= qFD ) ) + 
  geom_point()  + 
  xlab("Hill number") + 
  ylab("Functional Diversity") +
  theme_minimal()

#_____________________________________________________________________________
# Task: Let's talk about what is the ecological meaning of the graph.

# What are the different Hill number (q) telling us?
# Why is q=0 values higher than q=2
# Why functional diversity decreases with increasing Tau value?
