setwd("~/angeo_C4/cismotifs")
library(readr)
library(dplyr)

analysis_tomtom <- read_tsv("tomtom_main.tsv")
control_tomtom  <- read_tsv("tomtom_contorl.tsv")

analysis_best_hits <- analysis_tomtom %>%
  group_by(Query_ID) %>%
  slice_min(order_by = `E-value`, n = 1, with_ties = FALSE) %>%
  ungroup()

control_best_hits <- control_tomtom %>%
  group_by(Query_ID) %>%
  slice_min(order_by = `E-value`, n = 1, with_ties = FALSE) %>%
  ungroup()

colnames(analysis_best_hits) <- make.names(colnames(analysis_best_hits))
colnames(control_best_hits)  <- make.names(colnames(control_best_hits))

jaspar_meta <- read_tsv("jaspar_plants_metadata.csv")

analysis_mapped <- analysis_best_hits %>%
  left_join(jaspar_meta, by = c("Target_ID" = "matrix_id"))
write.csv(analysis_mapped,"Test_mapped_not_syntenic1.csv")

control_mapped <- control_best_hits %>%
  left_join(jaspar_meta, by = c("Target_ID" = "matrix_id"))
write.csv(control_mapped,"control_mapped_not_syntenic1.csv")


analysis_tf_counts <- analysis_mapped %>%
  filter(!is.na(name)) %>%
  count(name, sort = TRUE)

control_tf_counts <- control_mapped %>%
  filter(!is.na(name)) %>%
  count(name, sort = TRUE)

analysis_tf_counts <- analysis_tf_counts %>%
  mutate(Dataset = "Test")

control_tf_counts <- control_tf_counts %>%
  mutate(Dataset = "Control")

combined_tf <- bind_rows(analysis_tf_counts, control_tf_counts)
library(ggplot2)
colors <- c("Test" = "steelblue", "Control" = "red")
a= ggplot(combined_tf,
       aes(x = reorder(name, -n), y = n, fill = Dataset)) +
  
  geom_bar(stat = "identity", position = "dodge") +
  
  scale_fill_manual(values = colors) +
  
  labs(
    title = "Top TF motifs (Analysis vs Control)",
    x = "Transcription Factor",
    y = "Count"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



analysis_family <- analysis_mapped %>%
  filter(!is.na(family)) %>%
  count(family, sort = TRUE) %>%
  mutate(Dataset = "Test")

control_family <- control_mapped %>%
  filter(!is.na(family)) %>%
  count(family, sort = TRUE) %>%
  mutate(Dataset = "Control")
colors <- c("Test" = "steelblue", "Control" = "red")
combined_family <- bind_rows(analysis_family, control_family)
ggplot(combined_family,
       aes(x = reorder(family, -n), y = n, fill = Dataset)) +
  
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = colors) +
  labs(
    title = "TF Family Distribution (Analysis vs Control)",
    x = "TF Family",
    y = "Count"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



# not nyntenic dataset coming from Motifstenyreal.r
not_syntenic_test

tf_df <- analysis_mapped %>%
  mutate(
    Gene = str_extract(Query_ID, "^[^_]+"),
    Motif_ID = str_extract(Query_ID, "MEME-[0-9]+"),
    TF = name
  )

tf_df <- tf_df %>%
  select(Gene, Motif_ID, TF)

network_df <- not_syntenic_test %>%
  inner_join(tf_df, by = c("Gene", "Motif_ID"))

edges <- network_df %>%
  distinct(Gene, TF)
library(igraph)

nodes <- data.frame(
  name = unique(c(edges$Gene, edges$TF)),
  type = ifelse(unique(c(edges$Gene, edges$TF)) %in% edges$Gene, "Gene", "TF")
)

g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

plot(
  g,
  vertex.color = ifelse(V(g)$type == "Gene", "skyblue", "orange"),
  vertex.size = 15,
  vertex.label.cex = 0.7
)

library(ggraph)
library(tidygraph)

g_tbl <- as_tbl_graph(g)

ggraph(g_tbl, layout = "fr") +
  geom_edge_link(alpha = 0.5) +
  geom_node_point(aes(color = type), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c("Gene" = "steelblue", "TF" = "red")) +
  theme_void()


analysis_tomtom_c4shifted <- read_tsv("tomtom_C4shifted.tsv")
control_tomtom_c4shifted  <- read_tsv("tomtom_c4shift.tsv")

analysis_best_hits_c4shifted <- analysis_tomtom_c4shifted %>%
  group_by(Query_ID) %>%
  slice_min(order_by = `E-value`, n = 1, with_ties = FALSE) %>%
  ungroup()

control_best_hits_c4shifted <- control_tomtom_c4shifted %>%
  group_by(Query_ID) %>%
  slice_min(order_by = `E-value`, n = 1, with_ties = FALSE) %>%
  ungroup()

colnames(analysis_best_hits_c4shifted) <- make.names(colnames(analysis_best_hits_c4shifted))
colnames(control_best_hits_c4shifted)  <- make.names(colnames(control_best_hits_c4shifted))

jaspar_meta <- read_tsv("jaspar_plants_metadata.csv")

analysis_mapped_c4shifted <- analysis_best_hits_c4shifted %>%
  left_join(jaspar_meta, by = c("Target_ID" = "matrix_id"))
write.csv(analysis_mapped_c4shifted,"Test_mapped__c4shifted.csv")

control_mapped_c4shifted <- control_best_hits_c4shifted %>%
  left_join(jaspar_meta, by = c("Target_ID" = "matrix_id"))


analysis_tf_counts_c4shifted <- analysis_mapped_c4shifted %>%
  filter(!is.na(name)) %>%
  count(name, sort = TRUE)

control_tf_counts_c4shifted <- control_mapped_c4shifted %>%
  filter(!is.na(name)) %>%
  count(name, sort = TRUE)

analysis_tf_counts_c4shifted <- analysis_tf_counts_c4shifted %>%
  mutate(Dataset = "Test")

control_tf_counts_c4shifted <- control_tf_counts_c4shifted %>%
  mutate(Dataset = "Control")

combined_tf_c4shifted <- bind_rows(analysis_tf_counts_c4shifted, control_tf_counts_c4shifted)
library(ggplot2)
colors <- c("Test" = "steelblue", "Control" = "red")
b=ggplot(combined_tf_c4shifted,
       aes(x = reorder(name, -n), y = n, fill = Dataset)) +
  
  geom_bar(stat = "identity", position = "dodge") +
  
  scale_fill_manual(values = colors) +
  
  labs(
    title = "Top TF motifs (Analysis vs Control)",
    x = "Transcription Factor",
    y = "Count"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
a/b


analysis_family_c4shifted <- analysis_mapped_c4shifted %>%
  filter(!is.na(family)) %>%
  count(family, sort = TRUE) %>%
  mutate(Dataset = "Test")

control_family_c4shifted <- control_mapped_c4shifted %>%
  filter(!is.na(family)) %>%
  count(family, sort = TRUE) %>%
  mutate(Dataset = "Control")
colors <- c("Test" = "steelblue", "Control" = "red")
combined_family_c4shifted <- bind_rows(analysis_family_c4shifted, control_family_c4shifted)
b=ggplot(combined_family_c4shifted,
       aes(x = reorder(family, -n), y = n, fill = Dataset)) +
  
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = colors) +
  labs(
    title = "TF Family Distribution (Analysis vs Control)",
    x = "TF Family",
    y = "Count"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# not nyntenic dataset coming from Motifstenyreal.r
not_syntenic_test

tf_df <- analysis_mapped %>%
  mutate(
    Gene = str_extract(Query_ID, "^[^_]+"),
    Motif_ID = str_extract(Query_ID, "MEME-[0-9]+"),
    TF = name
  )

tf_df <- tf_df %>%
  select(Gene, Motif_ID, TF)


tf_df_c4shifted <- analysis_mapped_c4shifted %>%
  mutate(
    Gene = str_extract(Query_ID, "^[^_]+"),
    Motif_ID = str_extract(Query_ID, "MEME-[0-9]+"),
    TF = name
  )

tf_df_c4shifted <- tf_df_c4shifted %>%
  select(Gene, Motif_ID, TF)


network_df_c4shifted <- C4_shifted_test %>%
  inner_join(tf_df_c4shifted, by = c("Gene", "Motif_ID"))

common_cols <- intersect(names(network_df), names(network_df_c4shifted))
network_df_all <- rbind(network_df[common_cols], network_df_c4shifted[common_cols])
write_csv(network_df_all,"networkraw.csv")
network_df_all=read.csv("networkraw.csv")
gene_class <- network_df_all %>%
  group_by(Gene) %>%
  summarise(
    Gene_Class = case_when(
      any(Position_Class == "c4_shifted") ~ "c4_shifted",
      any(Position_Class == "syntenic") ~ "syntenic",
      TRUE ~ "not_syntenic"
    ),
    .groups = "drop"
  )

# ---------------------------
# Directed edges (TF → Gene)
# ---------------------------
edges <- network_df_all %>%
  distinct(Gene, TF) %>%
  transmute(from = TF, to = Gene)

# ---------------------------
# Nodes (DO NOT overwrite later)
# ---------------------------
nodes <- data.frame(
  name = unique(c(edges$from, edges$to))
) %>%
  mutate(
    type = ifelse(name %in% edges$to, "Gene", "TF")
  ) %>%
  left_join(gene_class, by = c("name" = "Gene")) %>%
  mutate(
    color = case_when(
      type == "TF" ~ "purple",
      Gene_Class == "c4_shifted" ~ "#2196F3",   # blue
      Gene_Class == "syntenic" ~ "#4CAF50",     # green
      Gene_Class == "not_syntenic" ~ "#9E9E9E", # grey
      TRUE ~ "black"
    )
  )

nodes <- nodes %>%
  mutate(
    Node_Group = case_when(
      type == "TF" ~ "TF",
      Gene_Class == "c4_shifted" ~ "c4_shifted",
      Gene_Class == "syntenic" ~ "syntenic",
      Gene_Class == "not_syntenic" ~ "not_syntenic",
      TRUE ~ "other"
    )
  )
# ---------------------------
# Build graph (DIRECTED)
# ---------------------------
g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)

# ---------------------------
# igraph plot (quick view)
# ---------------------------
plot(
  g,
  vertex.color = V(g)$color,
  vertex.size = 18,
  vertex.label.cex = 0.7,
  edge.arrow.size = 0.5,
  edge.color = "grey60"
)

# ---------------------------
# ggraph plot (publication-ready)
# ---------------------------
g_tbl <- as_tbl_graph(g)

ggraph(g_tbl, layout = "fr") +
  geom_edge_link(
    arrow = arrow(length = unit(3, "mm")),
    end_cap = circle(3, "mm"),
    alpha = 1,
    color = "black"
  ) +
  geom_node_point(aes(color = Node_Group), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  
  scale_color_manual(
    name = "Node Type",
    values = c(
      "TF" = "yellow",
      "c4_shifted" = "#2196F3",
      "syntenic" = "#4CAF50",
      "not_syntenic" = "orange"
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right"
  )

