rm(list=ls())
# First execute 'ODDS_ratio_calculation.m' to generate the odds-ratio
setwd("U:/SMART-MC/Simulation and Real Data Analysis")
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(patchwork)
coord_ratio <- 0.3
dd_quartiles_raw <- read.table("Output_dd_quartiles_raw.csv", header = FALSE, sep = ",") 
OR_q25 <- read.table("Output_OR_q25.csv", header = FALSE, sep = ",") 
OR_q50 <- read.table("Output_OR_q50.csv", header = FALSE, sep = ",") 
OR_q75 <- read.table("Output_OR_q75.csv", header = FALSE, sep = ",") 

heatmap.rownames <- c("IB to Fin", "IB to DF", "Fin to BcD", "IB to Nat", "DF to BcD", "Nat to BcD", "Nat to Fin")
heatmap.colnames <- c("30,F,W", "30,F,B", "30,F,O", "30,M,W", "30,M,B", "30,M,O",
              "60,F,W", "60,F,B", "60,F,O", "60,M,W", "60,M,B", "60,M,O")

dd_quartiles_raw

OR_q25.heat <- OR_q25[, 3:14]
OR_q50.heat <- OR_q50[, 3:14]
OR_q75.heat <- OR_q75[, 3:14]

########### Q1 #################################################################

matrix_data <- matrix(unlist(OR_q25.heat), ncol = length(OR_q25.heat), byrow = FALSE)
rownames(matrix_data) <- heatmap.rownames; colnames(matrix_data) <- heatmap.colnames
melted_data <- melt(matrix_data, varnames = c("Var1", "Var2"), value.name = "value")

# Plot with transposed axes
plot.Q1 <- ggplot(melted_data, aes(Var2, Var1, fill = log10(value))) +  # Swap Var1 and Var2
  geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), na.value = "white", 
                       name = "log10(OR)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Patients", 
    y = "Treatment transitions",
    title = "Disease duration = 9 months (25-th percentile)"
  ) +
  coord_fixed(ratio = coord_ratio) 
#plot.Q1

########### Q2 #################################################################

matrix_data <- matrix(unlist(OR_q50.heat), ncol = length(OR_q50.heat), byrow = FALSE)
rownames(matrix_data) <- heatmap.rownames; colnames(matrix_data) <- heatmap.colnames
melted_data <- melt(matrix_data, varnames = c("Var1", "Var2"), value.name = "value")

# Plot with transposed axes
plot.Q2 <- ggplot(melted_data, aes(Var2, Var1, fill = log10(value))) +  # Swap Var1 and Var2
  geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), na.value = "white", 
                       name = "log10(OR)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Patients", 
    y = "Treatment transitions",
    title = "Disease duration = 14 months (median)"
  ) +
  coord_fixed(ratio = coord_ratio) 
#plot.Q2


########### Q3 #################################################################

matrix_data <- matrix(unlist(OR_q75.heat), ncol = length(OR_q75.heat), byrow = FALSE)
rownames(matrix_data) <- heatmap.rownames; colnames(matrix_data) <- heatmap.colnames
melted_data <- melt(matrix_data, varnames = c("Var1", "Var2"), value.name = "value")

# Plot with transposed axes
plot.Q3 <- ggplot(melted_data, aes(Var2, Var1, fill = log10(value))) +  # Swap Var1 and Var2
  geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), na.value = "white", 
                       name = "log10(OR)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Patients", 
    y = "Treatment transitions",
    title = "Disease duration = 20 months (75-th percentile)"
  ) +
  coord_fixed(ratio = coord_ratio) 
#plot.Q3
#ggsave("Odds_ratio_Q3.png", plot = plot.Q3, width = 12, height = 4, dpi = 400)

combined_plot <- plot.Q1 / plot.Q2 / plot.Q3
ggsave("Plots/Odds_ratio_ALL_quartile.png", plot = combined_plot, width = 12, height = 12, dpi = 400)

