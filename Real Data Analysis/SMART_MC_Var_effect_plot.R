rm(list=ls())

# Load libraries
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(cowplot) 
library(reshape2)
library(stringr)
library(gridGraphics)

# Define helper function to generate individual plots
plot_transition <- function(filename, xvar, xgrid, xlab, title) {
  df <- read_csv(filename, show_col_types = FALSE)
  df_trans <- tail(df, 9)
  df_long <- df_trans %>%
    pivot_longer(cols = starts_with("grid"), names_to = "Grid", values_to = "Probability") %>%
    mutate(Grid = as.integer(gsub("grid", "", Grid)),
           X = xgrid[Grid])
  
  palette <- brewer.pal(9, "Set1")
  
  p <- ggplot(df_long, aes(x = X, y = Probability, color = Transition)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = palette) +
    labs(title = title, x = xlab, y = "Trans. Prob.", color = "Transition") +
    theme_minimal() +
    theme(
      text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  return(p)
}

# Grids
min_age <- 12; max_age <- 67; Numgrids <- 100
ages <- seq(min_age, max_age, length.out = Numgrids)
min_dd <- 0; max_dd <- 63
DDs <- seq(min_dd, max_dd, length.out = Numgrids)

# Generate 8 plots
p1 <- plot_transition("Effect_age_MW.csv", "Age", ages, "Age", "Trans. Prob. by Age (Male, White)")
p2 <- plot_transition("Effect_age_MB.csv", "Age", ages, "Age", "Trans. Prob. by Age (Male, Black)")
p3 <- plot_transition("Effect_age_FW.csv", "Age", ages, "Age", "Trans. Prob. by Age (Female, White)")
p4 <- plot_transition("Effect_age_FB.csv", "Age", ages, "Age", "Trans. Prob. by Age (Female, Black)")
p5 <- plot_transition("Effect_DD_MW.csv", "DD", DDs, "Disease Duration (months)", "Trans. Prob. by Dis. Dur. (Male, White)")
p6 <- plot_transition("Effect_DD_MB.csv", "DD", DDs, "Disease Duration (months)", "Trans. Prob. by Dis. Dur. (Male, Black)")
p7 <- plot_transition("Effect_DD_FW.csv", "DD", DDs, "Disease Duration (months)", "Trans. Prob. by Dis. Dur. (Female, White)")
p8 <- plot_transition("Effect_DD_FB.csv", "DD", DDs, "Disease Duration (months)", "Trans. Prob. by Dis. Dur. (Female, Black)")

# First collect legend from one of the plots
legend_plot <- (p1 | p2) + 
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "top", 
    legend.direction = "horizontal",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16, face = "bold")
  )
legend <- get_legend(legend_plot)

# Now remove legends from individual plots
p_list <- list(p1, p2, p3, p4, p5, p6, p7, p8)
p_list <- lapply(p_list, function(p) p + theme(legend.position = "none"))

# Build rows without legends
top_row <- p_list[[1]] | p_list[[2]]
second_row <- p_list[[3]] | p_list[[4]]
third_row <- p_list[[5]] | p_list[[6]]
fourth_row <- p_list[[7]] | p_list[[8]]

# Combine everything: legend on top
final_plot <- plot_grid(
  legend,
  top_row,
  second_row,
  third_row,
  fourth_row,
  ncol = 1,
  rel_heights = c(0.4, 1, 1, 1, 1) # adjust legend height here
)

# Show plot
print(final_plot)

# Save as full page high-res jpg
ggsave("Plots/SMART_MC_Transitions.jpg", final_plot, width = 12, height = 15, dpi = 600)

################################################################################
### Initial probability Heatmap
################################################################################

ages <- seq(min_age, max_age, length.out = Numgrids)

# Target ages
target_ages <- c(30, 60)

# Subpopulations
populations <- expand.grid(
  Gender = c("M", "F"),
  Race = c("W", "B")
) %>%
  mutate(ID = paste(Gender, Race, sep = "_"))

# Map filenames
filename_map <- list(
  "M_W" = "Effect_age_MW.csv",
  "M_B" = "Effect_age_MB.csv",
  "F_W" = "Effect_age_FW.csv",
  "F_B" = "Effect_age_FB.csv"
)

# Expand for both ages
pop_expanded <- populations %>%
  slice(rep(1:n(), each = length(target_ages))) %>%
  mutate(Age = rep(target_ages, times = nrow(populations)))

# Initialize
heatmap_data <- data.frame()

for (i in 1:nrow(pop_expanded)) {
  
  gender <- pop_expanded$Gender[i]
  race <- pop_expanded$Race[i]
  age_val <- pop_expanded$Age[i]
  
  id <- paste(gender, race, sep = "_")
  
  file <- filename_map[[id]]
  
  df <- read_csv(file, show_col_types = FALSE)
  df_init <- df[1:7,]
  
  df_long <- df_init %>%
    pivot_longer(cols = starts_with("grid"), names_to = "Grid", values_to = "Probability") %>%
    mutate(Grid = as.integer(gsub("grid", "", Grid)),
           AgeGrid = ages[Grid]) 
  
  interp_probs <- df_long %>%
    group_by(Transition) %>%
    summarise(Prob = approx(AgeGrid, Probability, xout = age_val)$y) %>%
    ungroup()
  
  subpop_label <- paste0(age_val, " (age)\n", gender, ", ", race)
  
  interp_probs <- interp_probs %>%
    mutate(Population = subpop_label)
  
  heatmap_data <- bind_rows(heatmap_data, interp_probs)
}

# Treatment order
treatment_order <- c("BcD", "GA", "IB", "DF", "Nat", "S1P", "AL")
heatmap_data <- heatmap_data %>%
  mutate(Transition = factor(Transition, levels = treatment_order))

heatmap_data <- heatmap_data %>%
  mutate(Transition = factor(Transition, levels = rev(treatment_order)))


# Reshape for heatmap
heatmap_wide <- heatmap_data %>%
  pivot_wider(names_from = Population, values_from = Prob)

heatmap_long <- heatmap_wide %>%
  pivot_longer(cols = -Transition, names_to = "Subpopulation", values_to = "Probability")

# Fix order of subpopulations
heatmap_long <- heatmap_long %>%
  mutate(
    Age = as.numeric(str_extract(Subpopulation, "^[0-9]+")),
    Gender = str_extract(Subpopulation, "\n([MF]),") %>% str_remove_all("\n|,"),
    Race = str_extract(Subpopulation, "[WB]$"),
    Subpopulation = factor(Subpopulation, levels = heatmap_long %>%
                             arrange(Age, Gender, Race) %>%
                             pull(Subpopulation) %>% unique())
  )

# Plot
init_prob_plot <- ggplot(heatmap_long, aes(x = Subpopulation, y = Transition, fill = Probability)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    title = "SMART-MC Estimated Initial Treatment Probabilities by Age, Sex, and Race",
    x = NULL, #"Subpopulations",
    y = "Initial Treatment"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold", vjust = -0.5),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.margin = unit(c(2, 1, 1, 1), "lines")
  ) +
  scale_x_discrete(position = "top")

print(init_prob_plot)

ggsave("Plots/SMART_MC_Initial_probs.jpg", init_prob_plot, width = 12, height = 4, dpi = 600)

################################################################################
### Initial probability Heatmap (Long format, subpopulations on y-axis)
################################################################################

ages <- seq(min_age, max_age, length.out = Numgrids)

# Target ages
target_ages <- c(30, 60)

# Subpopulations
populations <- expand.grid(
  Gender = c("M", "F"),
  Race = c("W", "B")
) %>%
  mutate(ID = paste(Gender, Race, sep = "_"))

# Map filenames
filename_map <- list(
  "M_W" = "Effect_age_MW.csv",
  "M_B" = "Effect_age_MB.csv",
  "F_W" = "Effect_age_FW.csv",
  "F_B" = "Effect_age_FB.csv"
)

# Expand for both ages
pop_expanded <- populations %>%
  slice(rep(1:n(), each = length(target_ages))) %>%
  mutate(Age = rep(target_ages, times = nrow(populations)))

# Initialize
heatmap_data <- data.frame()

for (i in 1:nrow(pop_expanded)) {
  
  gender <- pop_expanded$Gender[i]
  race <- pop_expanded$Race[i]
  age_val <- pop_expanded$Age[i]
  
  id <- paste(gender, race, sep = "_")
  file <- filename_map[[id]]
  
  df <- read_csv(file, show_col_types = FALSE)
  df_init <- df[1:7,]
  
  df_long <- df_init %>%
    pivot_longer(cols = starts_with("grid"), names_to = "Grid", values_to = "Probability") %>%
    mutate(Grid = as.integer(gsub("grid", "", Grid)),
           AgeGrid = ages[Grid]) 
  
  interp_probs <- df_long %>%
    group_by(Transition) %>%
    summarise(Prob = approx(AgeGrid, Probability, xout = age_val)$y) %>%
    ungroup()
  
  subpop_label <- paste0(age_val, " (age)\n", gender, ", ", race)
  
  interp_probs <- interp_probs %>%
    mutate(Population = subpop_label)
  
  heatmap_data <- bind_rows(heatmap_data, interp_probs)
}

# Treatment order
treatment_order <- c("BcD", "GA", "IB", "DF", "Nat", "S1P", "AL")
heatmap_data <- heatmap_data %>%
  mutate(Transition = factor(Transition, levels = treatment_order))

# Reshape for heatmap
heatmap_wide <- heatmap_data %>%
  pivot_wider(names_from = Population, values_from = Prob)

heatmap_long <- heatmap_wide %>%
  pivot_longer(cols = -Transition, names_to = "Subpopulation", values_to = "Probability")

# Set custom subpopulation order: top to bottom in plot
custom_order <- c(
  "30 (age)\nF, W",
  "30 (age)\nF, B",
  "30 (age)\nM, W",
  "30 (age)\nM, B",
  "60 (age)\nF, W",
  "60 (age)\nF, B",
  "60 (age)\nM, W",
  "60 (age)\nM, B"
)

heatmap_long <- heatmap_long %>%
  mutate(Subpopulation = factor(Subpopulation, levels = rev(custom_order)))  # reverse for top-to-bottom order

# Plot (x-axis moved to top)
init_prob_plot <- ggplot(heatmap_long, aes(x = Transition, y = Subpopulation, fill = Probability)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    x = "Initial Treatment",
    y = "Subpopulations"
  ) +
  theme_minimal() +
  theme(
    axis.text.x.top = element_text(angle = 0, vjust = 0.5, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x.top = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  ) +
  scale_x_discrete(position = "top")

print(init_prob_plot)

# Save with adjusted height
ggsave("Plots/SMART_MC_Initial_probs_long.jpg", init_prob_plot, width = 6, height = 8, dpi = 600)