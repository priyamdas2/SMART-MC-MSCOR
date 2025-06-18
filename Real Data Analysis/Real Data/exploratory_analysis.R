rm(list = ls())
setwd("U:/SMART-MC/Real Data Analysis/Real Data")
library(tidyr)
library(dplyr)
library(readr)
library(ggalluvial)
library(ggplot2)
library(purrr)
################################################################################
### 1 Patient summary 
################################################################################
data <- read.csv("MS covariates.csv", header = FALSE)
colnames(data) <- c("patient_id", "age", "disease_duration_months", "female", "white", "black")

sample_size <- nrow(data)

mean_age <- mean(data$age)
sd_age <- sd(data$age)

mean_duration <- mean(data$disease_duration_months)
sd_duration <- sd(data$disease_duration_months)

# Female count and percent
n_female <- sum(data$female)
perc_female <- 100 * n_female / sample_size

# White count and percent
n_white <- sum(data$white)
perc_white <- 100 * n_white / sample_size

# Black count and percent
n_black <- sum(data$black)
perc_black <- 100 * n_black / sample_size

# Print formatted output
cat("Sample size:", sample_size, "\n")
cat(sprintf("Mean age (SD): %.1f (%.1f)", mean_age, sd_age), "\n")
cat(sprintf("Mean disease duration (months) (SD): %.1f (%.1f)", mean_duration, sd_duration), "\n")
cat(sprintf("Female n (%%): %d (%.1f%%)", n_female, perc_female), "\n")
cat(sprintf("White n (%%): %d (%.1f%%)", n_white, perc_white), "\n")
cat(sprintf("Black n (%%): %d (%.1f%%)", n_black, perc_black), "\n")


################################################################################
### 2 Treatment transition summary
################################################################################


# Read the data
data <- read_csv("MS_treatment_sequences_collapsed.csv", col_names = FALSE)
colnames(data) <- c("patient_id", "treatment")

# Build full sequences per patient
sequences <- data %>%
  group_by(patient_id) %>%
  summarise(seq = list(treatment))

# Determine maximum sequence length
max_steps <- max(sapply(sequences$seq, length))

# Expand sequences to wide format (padding shorter sequences with NA)
sequence_df <- sequences %>%
  mutate(seq = lapply(seq, function(x) c(x, rep(NA, max_steps - length(x))))) %>%
  unnest_wider(seq, names_sep = "_")

# Rename columns as Step_1, Step_2, ...
colnames(sequence_df)[2:(max_steps+1)] <- paste0("Step_", 1:max_steps)

# Map treatment numbers to names
treatment_map <- c(
  "1" = "B-cell depletion",
  "2" = "Glatiramer acetate",
  "3" = "Interferon-beta",
  "4" = "Dimethyl fumarate",
  "5" = "Natalizumab",
  "6" = "S1P modulators",
  "7" = "Aggressive / Legacy therapies"
)

sequence_df <- sequence_df %>%
  mutate(across(starts_with("Step_"), ~treatment_map[as.character(.)]))

# Convert wide to long format for plotting
sequence_long <- sequence_df %>%
  pivot_longer(cols = starts_with("Step_"), names_to = "step", values_to = "treatment") %>%
  drop_na(treatment)

# Extract numeric step order
sequence_long <- sequence_long %>%
  mutate(step_num = as.numeric(sub("Step_", "", step)))

# Create proper treatment factor levels for ordered legend
treatment_levels <- c(
  "B-cell depletion",
  "Glatiramer acetate",
  "Interferon-beta",
  "Dimethyl fumarate",
  "Natalizumab",
  "S1P modulators",
  "Aggressive / Legacy therapies"
)

sequence_long$treatment <- factor(sequence_long$treatment, levels = treatment_levels)

# Aggregate counts
step_counts <- sequence_long %>%
  group_by(step_num, treatment) %>%
  summarise(count = n(), .groups = "drop")

# Plot stacked bar chart as approximate alluvial
ggplot(step_counts, aes(x = step_num, y = count, fill = treatment)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_continuous(breaks = 1:max_steps, labels = paste("Visit", 1:max_steps)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(7, "Set2")) +
  labs(x = "Visit Number", y = "Number of Patients", fill = "Treatment") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = c(0.74, 0.76),  # <<<<< THIS is your top-right legend positioning
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_blank()
  )

ggsave("Treatment_Sequences_Alluvial.jpeg", width = 7, height = 6, units = "in", dpi = 300)

################################################################################
### 3 Transition Matrix Plot
################################################################################

data <- read_csv("MS_treatment_sequences_collapsed.csv", col_names = FALSE)
colnames(data) <- c("patient_id", "treatment")

# Build transitions
transitions <- data %>%
  group_by(patient_id) %>%
  arrange(patient_id) %>%
  summarise(seq = list(treatment)) %>%
  pull(seq) %>%
  map_df(~{
    if(length(.) > 1) {
      data.frame(from = .[-length(.)], to = .[-1])
    } else {
      data.frame(from = integer(0), to = integer(0))
    }
  })

# Map treatment codes to full names first
treatment_full <- c(
  "1" = "B-cell depletion",
  "2" = "Glatiramer acetate",
  "3" = "Interferon-beta",
  "4" = "Dimethyl fumarate",
  "5" = "Natalizumab",
  "6" = "S1P modulators",
  "7" = "Aggressive / Legacy therapies"
)

transitions_named <- transitions %>%
  mutate(from = treatment_full[as.character(from)],
         to = treatment_full[as.character(to)])

# Now apply short name mapping
short_labels <- c(
  "B-cell depletion" = "BcD",
  "Glatiramer acetate" = "GA",
  "Interferon-beta" = "IB",
  "Dimethyl fumarate" = "DF",
  "Natalizumab" = "Nat",
  "S1P modulators" = "S1P",
  "Aggressive / Legacy therapies" = "AL"
)

transitions_named <- transitions_named %>%
  mutate(from = short_labels[from],
         to = short_labels[to])

# Create factor levels to preserve ordering as before
treatment_levels_short <- c("BcD", "GA", "IB", "DF", "Nat", "S1P", "AL")

# Compute empirical transition counts
transition_matrix <- transitions_named %>%
  count(from, to, name = "count") %>%
  mutate(from = factor(from, levels = treatment_levels_short),
         to = factor(to, levels = treatment_levels_short))

# Plot heatmap
ggplot(transition_matrix, aes(x = to, y = from, fill = count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = count), size = 3.5) +
  scale_fill_gradient(
    low = "white", high = "steelblue", trans = "sqrt",
    breaks = c(0, 50, 200, 500, 1000, 2000),
    limits = c(0, max(transition_matrix$count))
  ) +
  labs(x = "To (treatments)", y = "From (treatments)", fill = "Count") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12)
  )

ggsave("Treatment_trans_mat.jpeg", width = 7.5, height = 6, units = "in", dpi = 300)


################################################################################
### 4 Transition summary
################################################################################


# Read the data
data <- read_csv("MS_treatment_sequences_collapsed.csv", col_names = FALSE)
colnames(data) <- c("patient_id", "treatment")

# Map treatment codes to short names
treatment_map <- c(
  "1" = "BcD",
  "2" = "GA",
  "3" = "IB",
  "4" = "DF",
  "5" = "Nat",
  "6" = "S1P",
  "7" = "AL"
)

data <- data %>%
  mutate(treatment = treatment_map[as.character(treatment)])

# Create full sequences per patient
sequences <- data %>%
  group_by(patient_id) %>%
  summarise(seq = list(treatment))

# Calculate number of transitions per patient
sequences <- sequences %>%
  mutate(num_transitions = sapply(seq, function(x) max(length(x)-1, 0)))

# Total sample size
N <- nrow(sequences)

# Proportion who switched at least once
n_switched <- sum(sequences$num_transitions >= 1)
prop_switched <- round(100 * n_switched / N, 1)

# Summary of number of transitions
mean_transitions <- round(mean(sequences$num_transitions), 2)
median_transitions <- round(median(sequences$num_transitions), 0)
max_transitions <- max(sequences$num_transitions)

# Create full transition pairs
transitions <- data %>%
  group_by(patient_id) %>%
  arrange(patient_id) %>%
  summarise(seq = list(treatment)) %>%
  pull(seq) %>%
  map_df(~{
    if(length(.) > 1) {
      data.frame(from = .[-length(.)], to = .[-1])
    } else {
      data.frame(from = character(0), to = character(0))
    }
  })

# Count most frequent transitions
top_transitions <- transitions %>%
  count(from, to, name = "count") %>%
  arrange(desc(count)) %>%
  slice(1:20)

# Print everything
cat("Total patients:", N, "\n")
cat("Patients who switched at least once:", n_switched, "(", prop_switched, "%)\n")
cat("Mean transitions:", mean_transitions, "\n")
cat("Median transitions:", median_transitions, "\n")
cat("Max transitions:", max_transitions, "\n")
cat("\nTop transitions:\n")
print(top_transitions)