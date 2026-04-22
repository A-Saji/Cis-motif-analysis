#motifplots
setwd("~/angeo_C4/cismotifs/positivecontrocis/meme_positivecontrol")
library(tidyverse)

df <- read_tsv("combined_motif_table.tsv")

summary_df <- df %>%
  group_by(Gene) %>%
  summarise(
    Total_motifs = n(),
    Absent_in_C3 = sum(C3_promoters_with_motif == 0)
  )
colnames(summary_df)[2]=c("Enriched")
colnames(summary_df)[3]=c("Exclusive")


plot_df <- summary_df %>%
  pivot_longer(
    cols = c(Enriched, Exclusive),
    names_to = "Type",
    values_to = "Count"
  )

plot_df
testmotifs=ggplot(plot_df, aes(x = Gene, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("violet", "yellow")) +
  scale_y_continuous(
    breaks = seq(min(plot_df$Count),
                 max(plot_df$Count), 
                 by = 1),
    expand = c(0, 0)
  ) +
  labs(
    title = "Motif counts per gene",
    x = "Gene",
    y = "Number of motifs",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

setwd("~/angeo_C4/cismotifs")

# Analysis
enriched_analysis_c4 <- read_csv("all_motif_with_orthogroup.csv")
exclusive_analysis_c4 <- read.csv("Uniquemotif.csv")
enriched_analysis_c3 <- read_tsv("Moitfs_inC3_orthogroups.tsv")
exclusive_analysis_c3 <- read_tsv("unique_c3_pos_c4_neg_analysis.csv")


# Control
enriched_control_c4 <- read_tsv("combined_motif_table_control.csv")
exclusive_control_c4 <- read_tsv("Unique_motifs_control.csv")
enriched_control_c3 <- read_tsv("Motifs_in_c3_control.tsv")
exclusive_control_c3 <- read_tsv("unique_c3_pos_c4_neg.csv")

# ----------------------------
# BUILD SUMMARY TABLE
# ----------------------------

plot_df <- tibble(
  Dataset = c("Test","Test","Test","Test",
              "Control","Control","Control","Control"),
  
  Type = c("C4","C4","C3","C3",
           "C4","C4","C3","C3"),
  
  Category = c("Enriched","Exclusive","Enriched","Exclusive",
               "Enriched","Exclusive","Enriched","Exclusive"),
  
  Count = c(
    nrow(enriched_analysis_c4),
    nrow(exclusive_analysis_c4),
    nrow(enriched_analysis_c3),
    nrow(exclusive_analysis_c3),
    
    nrow(enriched_control_c4),
    nrow(exclusive_control_c4),
    nrow(enriched_control_c3),
    nrow(exclusive_control_c3)
  )
)

# Combine Dataset + Type for color control
plot_df <- plot_df %>%
  mutate(Group = paste(Dataset, Type, sep = "_"))

plot_df <- plot_df %>%
  group_by(Dataset, Type) %>%
  mutate(
    total = sum(Count),
    perc = ifelse(Category == "Exclusive", Count / total * 100, NA)
  ) %>%
  ungroup()

# ----------------------------
# COLOR SCHEME (your request)
# ----------------------------

my_colors <- c(
  "Test_C4" = "#1f78b4",  # dark blue
  "Test_C3" = "#a6cee3",  # light blue
  "Control_C4"  = "#e31a1c",  # dark red
  "Control_C3"  = "#fb9a99"   # light red
)


# ----------------------------
# PLOT
# ----------------------------
library(ggpattern)

c <- ggplot(plot_df, aes(x = Dataset, y = Count, fill = Group)) +
  
  geom_bar_pattern(
    aes(
      pattern = Category,
      group = interaction(Dataset, Type)
    ),
    stat = "identity",
    position = position_dodge(width = 0.7),
    width = 0.6,
    
    # pattern settings
    pattern_fill = "black",
    pattern_density = 0.3,
    pattern_spacing = 0.02
  ) +
  
  scale_fill_manual(values = my_colors) +
  
  scale_pattern_manual(values = c(
    "Enriched" = "none",
    "Exclusive" = "stripe"
  )) +
  
  geom_text(
    data = subset(plot_df, Category == "Exclusive"),
    aes(
      label = paste0(round(perc, 1), "%"),
      group = interaction(Dataset, Type)
    ),
    position = position_dodge(width = 0.7),
    vjust = -0.5,
    size = 5,
    fontface = "bold"
  ) +
  
  labs(
    title = "Motif comparison: Test vs Control",
    x = "Dataset",
    y = "Number of motifs",
    fill = "Dataset + Type",
    pattern = "Motif Category"   # 👈 legend for pattern
  ) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )

testmotifs+c
allgenes=read.csv("genenames.csv")
allcontrol=read.csv("genenamescontrol.csv")

analysis_og <- enriched_analysis_c4 %>%
  filter(Orthogroup != "NA") %>%
  count(Orthogroup)

control_og <- enriched_control_c4 %>%
  filter(Gene != "NA") %>%
  count(Gene)

colnames(control_og)[1] <- "Orthogroup"

# Analysis
analysis_full <- allgenes %>%
  rename(Orthogroup = 1) %>%   # adjust if needed
  left_join(analysis_og, by = "Orthogroup") %>%
  mutate(n = ifelse(is.na(n), 0, n))

# Control
control_full <- allcontrol %>%
  rename(Orthogroup = 1) %>%
  left_join(control_og, by = "Orthogroup") %>%
  mutate(n = ifelse(is.na(n), 0, n))
analysis_dist <- analysis_full %>%
  count(n)

control_dist <- control_full %>%
  count(n)
analysis_dist <- analysis_dist %>%
  mutate(freq = nn / sum(nn))

control_dist <- control_dist %>%
  mutate(freq = nn / sum(nn))

combined_dist <- bind_rows( analysis_dist %>% mutate(Dataset = "Test"), control_dist %>% mutate(Dataset = "Control") )
all_n <- 0:max(combined_dist$n)

combined_dist <- combined_dist %>%
  tidyr::complete(n, Dataset, fill = list(freq = 0))
combined=ggplot(combined_dist, aes(x = n, y = freq, fill = Dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  
  scale_fill_manual(values = c(
    "Test" = "steelblue",
    "Control"  = "red"
  )) +
  
  scale_x_continuous(
    breaks = seq(min(combined_dist$n),
                 max(combined_dist$n),
                 by = 1)
  ) +
  
  labs(
    title = "Distribution of motifs per orthogroup",
    x = "Motifs per orthogroup",
    y = "Frequency",
    fill = "Dataset"
  ) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )


set.seed(69)
sample_once <- function(df, size = 51) {
  df_sample <- df %>% sample_n(size)
  
  dist <- df_sample %>%
    count(n)
  
  list(
    raw = df_sample,
    dist = dist
  )
}

samples_list <- replicate(10, sample_once(analysis_full), simplify = FALSE)
dist_list <- lapply(samples_list, function(x) x$dist)
all_n <- 0:max(analysis_full$n)
dist_list_fixed <- lapply(dist_list, function(df) {
  df %>%
    right_join(tibble(n = all_n), by = "n") %>%
    mutate(nn = ifelse(is.na(nn), 0, nn)) %>%
    mutate(freq = nn / sum(nn)) %>%
    select(n, freq)
})
combined_samples <- bind_rows(dist_list_fixed, .id = "replicate")

mean_freq <- combined_samples %>%
  group_by(n) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")
control_freq <- control_full %>%
  count(n) %>%
  right_join(tibble(n = all_n), by = "n") %>%
  mutate(nn = ifelse(is.na(nn), 0, nn)) %>%
  mutate(freq = nn / sum(nn)) %>%
  select(n, freq)

freq_stats <- combined_samples %>%
  group_by(n) %>%
  summarise(
    mean = mean(freq),
    sd = sd(freq),
    .groups = "drop"
  )


analysis_df <- freq_stats %>%
  select(n, mean, sd) %>%
  mutate(Group = "Test")

control_df <- control_freq %>%
  rename(mean = freq) %>%
  mutate(
    sd = NA,
    Group = "Control"
  )

plot_df <- bind_rows(analysis_df, control_df)
library(ggplot2)

colors <- c("Test" = "steelblue", "Control" = "red")

meanfreq=ggplot(plot_df, aes(x = n, y = mean, fill = Group)) +
  
  # --- Bars ---
  geom_col(position = position_dodge(width = 0.8), alpha = 0.7) +
  
  # --- Error bars (only for Analysis) ---
  geom_errorbar(
    data = subset(plot_df, Group == "Test"),
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.3,
    position = position_dodge(width = 0.8),
    color = "purple3"
  ) +
  
  # --- X axis: show ALL integers ---
  scale_x_continuous(
    breaks = seq(min(plot_df$n), max(plot_df$n), by = 1)
  ) +
  
  scale_fill_manual(values = colors) +
  
  labs(
    title = "Motif Distribution: Test vs Control",
    x = "Motifs per orthogroup",
    y = "Frequency",
    fill = "Dataset"
  ) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )


#evalue distribution
enriched_analysis_c4$Group <- "Test"
enriched_control_c4$Group  <- "Control"
analysis_df <- data.frame(
  Orthogroup = enriched_analysis_c4$Orthogroup,
  Motif      = enriched_analysis_c4$Consensus_Sequence,
  E_value    = enriched_analysis_c4$E_value,
  Group      = "Test"
)

control_df <- data.frame(
  Orthogroup = enriched_control_c4$Gene,
  Motif      = enriched_control_c4$Motif_Name,
  E_value    = enriched_control_c4$E_value,
  Group      = "Control"
)

df <- rbind(analysis_df, control_df)
df$logE <- -log10(df$E_value + 1e-300)
library(ggplot2)
colors <- c("Test" = "steelblue", "Control" = "red")

g=ggplot(df, aes(x = logE, color = Group, fill = Group)) +
  geom_density(alpha = 0.3, size = 1) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(
    title = "Distribution of motif E-values (C4)",
    x = "-log10(E-value)",
    y = "Density",
    color = "Group",
    fill  = "Group"
  )
h=ggplot(df, aes(x = Group, y = logE, fill = Group)) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(
    title = "E-value comparison",
    y = "-log10(E-value)",
    fill = "Group"
  )

g+h

combined+meanfreq
library(readr)

# --- Load data ---
analysis_tomtom <- read_tsv("combined_tomtom_0.01.tsv")
control_tomtom  <- read_tsv("combined_tomtom_control.tsv")

analysis_motifs <- read_csv("all_motif_with_orthogroup.csv")
control_motifs  <- read_tsv("combined_motif_table_control.csv")

# Fix column name
colnames(analysis_motifs)[6] <- "Gene"
colnames(analysis_tomtom)[2]="Motif_ID"
analysis_tomtom_clean <- analysis_tomtom %>%
  group_by(Gene, Motif_ID) %>%
  slice_min(p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

control_tomtom_clean <- control_tomtom %>%
  group_by(Gene, Motif_ID) %>%
  slice_min(p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

analysis_annotated <- analysis_tomtom_clean %>%
  filter(q_value < 0.05) %>%
  distinct(Gene, Motif_ID)

control_annotated <- control_tomtom_clean %>%
  filter(q_value < 0.05) %>%
  distinct(Gene, Motif_ID)

analysis_gene <- analysis_motifs %>%
  distinct(Gene, Motif_ID) %>%
  group_by(Gene) %>%
  summarise(Total = n()) %>%
  left_join(
    analysis_annotated %>%
      group_by(Gene) %>%
      summarise(Annotated = n()),
    by = "Gene"
  ) %>%
  mutate(
    Annotated = replace_na(Annotated, 0),
    Percent = (Annotated / Total) * 100,
    Group = "Test"
  )

control_gene <- control_motifs %>%
  distinct(Gene, Motif_ID) %>%
  mutate(Annotated = paste(Gene, Motif_ID) %in%
           paste(control_annotated$Gene, control_annotated$Motif_ID)) %>%
  group_by(Gene) %>%
  summarise(
    Total = n(),
    Annotated = sum(Annotated),
    Percent = (Annotated / Total) * 100
  ) %>%
  mutate(Group = "Control")

df <- bind_rows(analysis_gene, control_gene)
library(ggplot2)

colors <- c("Test" = "steelblue", "Control" = "red")

analysis_total <- analysis_motifs %>%
  distinct(Gene, Motif_ID) %>%
  nrow()

analysis_annotated_n <- analysis_annotated %>%
  nrow()

control_total <- control_motifs %>%
  distinct(Gene, Motif_ID) %>%
  nrow()

control_annotated_n <- control_annotated %>%
  nrow()
plot_df <- data.frame(
  Dataset = c("Test", "Control"),
  Total = c(analysis_total, control_total),
  Annotated = c(analysis_annotated_n, control_annotated_n)
) %>%
  mutate(
    Non_Annotated = Total - Annotated,
    Percent = Annotated / Total * 100
  )
plot_long <- plot_df %>%
  select(Dataset, Annotated, Non_Annotated) %>%
  pivot_longer(
    cols = c(Annotated, Non_Annotated),
    names_to = "Category",
    values_to = "Count"
  )

library(ggplot2)
plot_long <- plot_long %>%
  mutate(Fill = paste(Dataset, Category, sep = "_"))
my_colors <- c(
  "Test_Annotated" = "#1f4e79",        # dark blue
  "Test_Non_Annotated" = "#4f81bd",    # light blue
  "Control_Annotated" = "#8b0000",     # dark red
  "Control_Non_Annotated" = "#ff4d4d"  # light red
)
g=ggplot(plot_long,
       aes(x = Dataset, y = Count, fill = Fill)) +
  
  geom_bar(stat = "identity", width = 0.6) +
  
  # percentage labels
  geom_text(
    data = plot_df,
    aes(
      x = Dataset,
      y = Annotated,
      label = paste0(round(Percent, 1), "%")
    ),
    inherit.aes = FALSE,
    vjust = -7,
    size = 5,
    fontface = "bold"
  ) +
  
  scale_fill_manual(values = my_colors) +
  
  labs(
    title = "Annotated Motif Fraction (C4)",
    x = "",
    y = "Number of motifs",
    fill = "Dataset"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )


colors <- c("Test" = "steelblue", "Control" = "red")
h=ggplot(df, aes(x = Percent, fill = Group)) +
  geom_histogram(
    position = "dodge",
    bins = 20,
    alpha = 0.7
  ) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(
    title = "% Annotated motifs per gene",
    x = "% Annotated motifs",
    y = "Number of genes",
    fill = "Group"
  )
g+h

library(dplyr)
colnames(analysis_motifs)[6]=c("Gene")
analysis_best_hits <- analysis_tomtom %>%
  group_by(Gene, Motif_ID) %>%
  slice_min(order_by = E_value, n = 1, with_ties = FALSE) %>%
  ungroup()


control_best_hits <- control_tomtom  %>%
  group_by(Gene, Motif_ID) %>%
  slice_min(order_by = E_value, n = 1, with_ties = FALSE) %>%
  ungroup()

jaspar_meta <- read_tsv("jaspar_plants_metadata.csv")

jaspar_meta <- jaspar_meta %>%
  mutate(matrix_base = base_id)

analysis_mapped <- analysis_best_hits %>%
  left_join(jaspar_meta, by = c("Target_ID" = "matrix_id"))

control_mapped <- control_best_hits %>%
  left_join(jaspar_meta, by = c("Target_ID" = "matrix_id"))

analysis_tf_counts <- analysis_mapped %>%
  filter(!is.na(name)) %>%   # remove unmapped
  count(name, sort = TRUE)

control_tf_counts <- control_mapped %>%
  filter(!is.na(name)) %>%
  count(name, sort = TRUE)
library(ggplot2)

x=ggplot(analysis_tf_counts, aes(x = reorder(name, -n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  
  geom_text(
    aes(y = n, label = name),  # explicitly anchor to bar height
    angle = 90,
    vjust = 0.5,
    hjust= -0.3
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.15),),
                     breaks = seq(0, max(analysis_tf_counts$n), by = 1)) +  # adds space above
  
  coord_cartesian(clip = "off") +  # prevents cutting labels
  
  labs(
    title = "TF Binding Site Counts (Test)",
    x = "Transcription Factor",
    y = "Count"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)  # extra space
  )

y=ggplot(control_tf_counts, aes(x = reorder(name, -n), y = n)) +
  geom_bar(stat = "identity", fill = "red") +
  geom_text(
    aes(y = n, label = name),  # explicitly anchor to bar height
    angle = 90,
    vjust = 0.5,
    hjust= -0.3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25),),
                     breaks = seq(0, max(analysis_tf_counts$n), by = 1)) +
  labs(
    title = "TF Binding Site Counts (Control)",
    x = "Transcription Factor",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)  # extra space
  )
x/y


library(dplyr)
library(tidyr)

analysis_mapped$Motif_ID=as.character(analysis_mapped$Motif_ID)
library(dplyr)
library(tidyr)

combined_tf <- bind_rows(
  analysis_mapped %>% mutate(Dataset = "Test"),
  control_mapped %>% mutate(Dataset = "Control")
) 
family_counts <- combined_tf %>%
  filter(!is.na(family)) %>%
  count(family, Dataset) %>%
  pivot_wider(names_from = Dataset, values_from = n, values_fill = 0)

write.csv(family_counts, "family_counts.csv", row.names = FALSE)

# convert to long format
df_long <- family_counts %>%
  pivot_longer(
    cols = c(Test, Control),
    names_to = "Dataset",
    values_to = "Count"
  )

# ensure order
df_long$Dataset <- factor(df_long$Dataset, levels = c("Test", "Control"))

# reorder families by total abundance (important!)
df_long <- df_long %>%
  group_by(family) %>%
  mutate(total = sum(Count)) %>%
  ungroup() %>%
  mutate(family = reorder(family, -total))

# plot
ggplot(df_long,
       aes(x = family, y = Count, fill = Dataset)) +
  
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.7)
  ) +
  
  # optional labels
  geom_text(
    aes(label = Count),
    position = position_dodge(width = 0.7),
    vjust = -0.3,
    size = 3.5
  ) +
  
  scale_fill_manual(values = c(
    "Test" = "steelblue",
    "Control" = "red"
  )) +
  
  labs(
    title = "TF Family Distribution (Test vs Control)",
    x = "TF Family",
    y = "Count",
    fill = "Dataset"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )
df_norm <- family_counts %>%
  mutate(
    Test_prop = Test / sum(Test),
    Control_prop = Control / sum(Control)
  )

df_plot <- df_norm %>%
  select(family, Test_prop, Control_prop) %>%
  pivot_longer(
    cols = c(Test_prop, Control_prop),
    names_to = "Condition",
    values_to = "Proportion"
  ) %>%
  mutate(
    Condition = recode(Condition,
                       "Test_prop" = "Test",
                       "Control_prop" = "Control")
  )

library(ggplot2)
ggplot(df_plot, aes(x = Proportion, y = family, fill = Condition)) +
  
  geom_bar(stat = "identity", position = "dodge") +
  
  scale_fill_manual(values = c(
    "Test" = "steelblue",
    "Control" = "red"
  )) +
  
  labs(
    title = "TF Family Distribution (Normalized)",
    x = "Proportion of motifs",
    y = "TF Family",
    fill = "Condition"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )

df_plot <- family_counts %>%
  mutate(
    Test_prop = Test / sum(Test),
    Control_prop = Control / sum(Control)
  ) %>%
  mutate(
    Test_plot = Test_prop,
    Control_plot = -Control_prop   # 👈 flip to negative
  ) %>%
  select(family, Test_plot, Control_plot) %>%
  pivot_longer(
    cols = c(Test_plot, Control_plot),
    names_to = "Condition",
    values_to = "Proportion"
  ) %>%
  mutate(
    Condition = recode(Condition,
                       "Test_plot" = "Test",
                       "Control_plot" = "Control")
  )

ggplot(df_plot, aes(x = Proportion, y = family, fill = Condition)) +
  
  geom_bar(stat = "identity") +
  
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +  # 👈 center line
  
  scale_fill_manual(values = c(
    "Test" = "steelblue",
    "Control" = "red"
  )) +
  
  labs(
    title = "TF Family Distribution (Normalized, Diverging)",
    x = "Proportion (Test → | ← Control)",
    y = "TF Family",
    fill = "Condition"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )

results <- family_counts %>%
  rowwise() %>%
  mutate(
    total_test = sum(family_counts$Test),
    total_control = sum(family_counts$Control),
    
    a = Test,
    c = Control,
    b = total_test - a,
    d = total_control - c,
    
    fisher = list(fisher.test(matrix(c(a, b, c, d), nrow = 2))),
    
    p_value = fisher$p.value,
    odds_ratio = fisher$estimate
  ) %>%
  ungroup()

results <- results %>%
  mutate(
    Enriched = case_when(
      odds_ratio > 1 ~ "Test",
      odds_ratio < 1 ~ "Control",
      TRUE ~ "Neutral"
    ),
    
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

results <- results %>%
  mutate(
    prop_test = Test / total_test,
    prop_control = Control / total_control,
    log2FC = log2((prop_test + 1e-3) / (prop_control + 1e-6))
  )
write_csv(results,"results_fischer's.csv")
results <- results %>%
  arrange(log2FC) %>%
  mutate(family = factor(family, levels = family))

ggplot(results,
       aes(x = log2FC, y = family, fill = Enriched)) +
  
  geom_bar(stat = "identity") +
  
  # zero reference line
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # significance labels
  geom_text(
    aes(label = sig),
    hjust = ifelse(results$log2FC > 0, -0.2, 1.2),
    size = 4
  ) +
  
  scale_fill_manual(values = c(
    "Test" = "steelblue",
    "Control" = "tomato",
    "Neutral" = "grey70"
  )) +
  
  labs(
    title = "TF Family Enrichment (Test vs Control)",
    x = "log2 Fold Change (proportion-based)",
    y = "TF Family",
    fill = "Enriched in"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )
