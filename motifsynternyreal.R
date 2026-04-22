library(dplyr)
library(readr)
library(ggplot2)

process_dataset <- function(base_path, dataset_name,
                            min_c4_species = 4,
                            min_c3_species = 1,
                            position_tolerance = 1,
                            support_threshold = 0.5) {
  
  setwd(base_path)
  
  meme_df <- read_tsv("meme_positions.csv")
  fimo_df <- read.delim("fimo_c3_positions.tsv", stringsAsFactors = FALSE)
  
  all_motifs <- bind_rows(
    meme_df %>% select(Gene, Sequence, Motif_ID, Start, Strand),
    fimo_df %>% select(Gene, Sequence, Motif_ID, Start, Strand)
  )
  
  # --- annotate groups + species ---
  all_motifs <- all_motifs %>%
    mutate(
      Group = case_when(
        grepl("^(Sevir|Urofu|Pavag|ELECO|GRMZM|Pahal)", Sequence) ~ "C4_PACMAD",
        grepl("^(OEL|Chala)", Sequence) ~ "C3_PACMAD",
        grepl("^(Bradi|LOC)", Sequence) ~ "C3_BEP",
        TRUE ~ "Unknown"
      ),
      Species = case_when(
        grepl("^GRMZM", Sequence) ~ "Maize",
        grepl("^Sevir", Sequence) ~ "Setaria_viridis",
        grepl("^Urofu", Sequence) ~ "Urochloa_fusca",
        grepl("^Pavag", Sequence) ~ "Paspalum_vaginatum",
        grepl("^ELECO", Sequence) ~ "Eleusine_coracana",
        grepl("^Pahal", Sequence) ~ "Panicum_hallii",
        grepl("^OEL", Sequence) ~ "Dicanthelium",
        grepl("^Chala", Sequence) ~ "Chasmanthum_laxum",
        grepl("^Bradi", Sequence) ~ "Brachypodium_distachyon",
        grepl("^LOC", Sequence) ~ "Oryza_sativa",
        TRUE ~ "Unknown"
      )
    )
  
  # --- species presence ---
  motif_species_presence <- all_motifs %>%
    group_by(Gene, Motif_ID, Species, Group) %>%
    summarise(hit = TRUE, .groups = "drop")
  
  # --- novel motifs ---
  novel_motifs <- motif_species_presence %>%
    group_by(Gene, Motif_ID) %>%
    summarise(
      n_c4_species = length(unique(Species[Group == "C4_PACMAD"])),
      n_c3_species = length(unique(Species[Group == "C3_PACMAD"])),
      .groups = "drop"
    ) %>%
    filter(n_c4_species >= min_c4_species & n_c3_species == 0)
  
  # --- architecture ---
  arch_df <- all_motifs %>%
    anti_join(novel_motifs, by = c("Gene", "Motif_ID")) %>%
    arrange(Gene, Sequence, Start) %>%
    group_by(Gene, Sequence, Species, Group) %>%
    summarise(motif_order = list(Motif_ID), .groups = "drop")
  
  # --- position extraction ---
  motif_positions <- do.call(rbind, lapply(1:nrow(arch_df), function(i) {
    motifs <- arch_df$motif_order[[i]]
    
    do.call(rbind, lapply(unique(motifs), function(m) {
      data.frame(
        Gene = arch_df$Gene[i],
        Sequence = arch_df$Sequence[i],
        Species = arch_df$Species[i],
        Group = arch_df$Group[i],
        Motif_ID = m,
        Position = which(motifs == m)[1]
      )
    }))
  }))
  
  # --- classification ---
  classify_motifs <- function(df) {
    
    get_consensus <- function(pos_vec) {
      if (length(pos_vec) == 0) return(list(pos = NA, support = 0))
      tab <- table(pos_vec)
      list(
        pos = as.numeric(names(tab)[which.max(tab)]),
        support = max(tab) / sum(tab)
      )
    }
    
    results <- data.frame()
    
    for (gene in unique(df$Gene)) {
      gene_df <- df[df$Gene == gene, ]
      
      for (motif in unique(gene_df$Motif_ID)) {
        
        sub <- gene_df[gene_df$Motif_ID == motif, ]
        
        c4  <- sub[sub$Group == "C4_PACMAD", ]
        c3p <- sub[sub$Group == "C3_PACMAD", ]
        c3b <- sub[sub$Group == "C3_BEP", ]
        
        c4_pos <- c4 %>% group_by(Species) %>% summarise(pos = min(Position), .groups="drop")
        c3_pos <- c3p %>% group_by(Species) %>% summarise(pos = min(Position), .groups="drop")
        c3b_pos <- c3b %>% group_by(Species) %>% summarise(pos = min(Position), .groups="drop")
        
        n_c4 <- nrow(c4_pos)
        n_c3p <- nrow(c3_pos)
        n_c3b <- nrow(c3b_pos)
        
        c4_cons <- get_consensus(c4_pos$pos)
        c3_cons <- if (n_c3p > 0) get_consensus(c3_pos$pos) else NULL
        
        c4_ok <- n_c4 >= min_c4_species && c4_cons$support >= support_threshold
        c3_ok <- n_c3p >= min_c3_species && !is.null(c3_cons)
        
        presence_class <- case_when(
          n_c4 >= min_c4_species & n_c3p == 0 ~ "Exclusive_C4",
          n_c4 > 0 & n_c3p > 0 & n_c3b > 0 ~ "pacmad_and_bep",
          n_c4 > 0 & n_c3p > 0 & n_c3b == 0 ~ "pacmad",
          TRUE ~ "partial_C4"
        )
        
        position_class <- "not_syntenic"
        
        if (presence_class == "pacmad_and_bep" && c4_ok && c3_ok) {
          if (abs(c4_cons$pos - c3_cons$pos) <= position_tolerance) {
            position_class <- "syntenic_in_all"
          } else {
            position_class <- "pacmad_and_bep"
          }
        } else if (presence_class == "pacmad" && c4_ok && c3_ok) {
          if (abs(c4_cons$pos - c3_cons$pos) <= position_tolerance) {
            position_class <- "syntenic"
          } else {
            position_class <- "c4_shifted"
          }
        }
        
        results <- rbind(results, data.frame(
          Gene = gene,
          Motif_ID = motif,
          Presence_Class = presence_class,
          Position_Class = position_class
        ))
      }
    }
    
    results
  }
  
  motif_classes <- classify_motifs(motif_positions)
  
  final <- bind_rows(
    motif_classes,
    novel_motifs %>%
      mutate(Presence_Class = "Exclusive_C4", Position_Class = "Exclusive_C4") %>%
      select(Gene, Motif_ID, Presence_Class, Position_Class)
  )
  
  final$Dataset <- dataset_name
  
  return(final)
}

test_results <- process_dataset(
  "~/angeo_C4/cismotifs/meme_res", "Test"
)

control_results <- process_dataset(
  "~/angeo_C4/cismotifs/control", "Control"
)
combined <- bind_rows(test_results, control_results)


plot_df <- combined %>%
  count(Dataset, Presence_Class, Position_Class) %>%
  group_by(Dataset) %>%
  mutate(prop = n / sum(n))

plot_df$Presence_Class <- factor(
  plot_df$Presence_Class,
  levels = c(
    "partial_C4",
    "Exclusive_C4",
    "pacmad",
    "pacmad_and_bep"
  )
)
ggplot(plot_df, aes(x = Presence_Class, y = prop, fill = Position_Class)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Dataset) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Motif classification (Test vs Control)",
    x = "Presence Class",
    y = "Proportion"
  ) +
  theme(
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.major.x = element_blank(),  # remove vertical lines
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

process_dataset <- function(base_path, dataset_name,
                            min_c4_species = 4,
                            position_tolerance = 1,
                            support_threshold = 0.5) {
  
  library(dplyr)
  library(readr)
  
  setwd(base_path)
  
  # ---------------------------
  # Load
  # ---------------------------
  meme_df <- read_tsv("meme_filtered.tsv")
  fimo_df <- read.delim("fimo_filtered.tsv", stringsAsFactors = FALSE)
  
  all_motifs <- bind_rows(
    meme_df %>% select(Gene, Sequence, Motif_ID, Start, Strand),
    fimo_df %>% select(Gene, Sequence, Motif_ID, Start, Strand)
  )
  
  # ---------------------------
  # Annotate
  # ---------------------------
  all_motifs <- all_motifs %>%
    mutate(
      Group = case_when(
        grepl("^(Sevir|Urofu|Pavag|ELECO|GRMZM|Pahal)", Sequence) ~ "C4_PACMAD",
        grepl("^(OEL|Chala)", Sequence) ~ "C3_PACMAD",
        grepl("^(Bradi|LOC)", Sequence) ~ "C3_BEP",
        TRUE ~ "Unknown"
      ),
      Species = case_when(
        grepl("^GRMZM", Sequence) ~ "Maize",
        grepl("^Sevir", Sequence) ~ "Setaria_viridis",
        grepl("^Urofu", Sequence) ~ "Urochloa_fusca",
        grepl("^Pavag", Sequence) ~ "Paspalum_vaginatum",
        grepl("^ELECO", Sequence) ~ "Eleusine_coracana",
        grepl("^Pahal", Sequence) ~ "Panicum_hallii",
        grepl("^OEL", Sequence) ~ "Dicanthelium",
        grepl("^Chala", Sequence) ~ "Chasmanthium_latifolium",
        grepl("^Bradi", Sequence) ~ "Brachypodium_distachyon",
        grepl("^LOC", Sequence) ~ "Oryza_sativa",
        TRUE ~ "Unknown"
      )
    )
  
  # ---------------------------
  # Presence classification (FIXED)
  # ---------------------------
  motif_presence <- all_motifs %>%
    group_by(Gene, Motif_ID) %>%
    summarise(
      n_c4  = n_distinct(Species[Group == "C4_PACMAD"]),
      n_c3p = n_distinct(Species[Group == "C3_PACMAD"]),
      n_c3b = n_distinct(Species[Group == "C3_BEP"]),
      .groups = "drop"
    ) %>%
    mutate(
      Presence_Class = case_when(
        n_c4 >= min_c4_species & n_c3p == 0 & n_c3b == 0 ~ "C4_only",
        n_c4 > 0 & n_c3p > 0 & n_c3b == 0 ~ "PACMAD",
        n_c4 > 0 & n_c3p > 0 & n_c3b > 0 ~ "PACMAD_BEP",
        TRUE ~ "Other"
      )
    )
  # ---------------------------
  # 🔥 Capture dropped motifs
  # ---------------------------
  dropped_motifs <- motif_presence %>%
    filter(Presence_Class == "Other")
  # clean join (no duplicates)
  all_motifs <- all_motifs %>%
    left_join(
      motif_presence %>% select(Gene, Motif_ID, Presence_Class),
      by = c("Gene", "Motif_ID")
    )
  
  # ---------------------------
  # Position extraction (per class)
  # ---------------------------
  get_positions <- function(df_group) {
    
    arch_df <- df_group %>%
      arrange(Gene, Sequence, Start) %>%
      group_by(Gene, Sequence, Species, Group) %>%
      summarise(motif_order = list(Motif_ID), .groups = "drop")
    
    do.call(rbind, lapply(1:nrow(arch_df), function(i) {
      motifs <- arch_df$motif_order[[i]]
      
      do.call(rbind, lapply(unique(motifs), function(m) {
        data.frame(
          Gene = arch_df$Gene[i],
          Sequence = arch_df$Sequence[i],
          Species = arch_df$Species[i],
          Group = arch_df$Group[i],
          Motif_ID = m,
          Position = which(motifs == m)[1]
        )
      }))
    }))
  }
  
  # ---------------------------
  # Synteny (PACMAD + BEP FIXED)
  # ---------------------------
  compute_synteny <- function(df_positions) {
    
    results <- data.frame()
    
    get_consensus <- function(pos_vec) {
      if (length(pos_vec) == 0) return(list(pos = NA, support = 0))
      tab <- table(pos_vec)
      list(
        pos = as.numeric(names(tab)[which.max(tab)]),
        support = max(tab)/sum(tab)
      )
    }
    
    for (gene in unique(df_positions$Gene)) {
      
      gene_df <- df_positions[df_positions$Gene == gene, ]
      motifs <- unique(gene_df$Motif_ID)
      
      if (length(motifs) <= 2) {
        results <- rbind(results, data.frame(
          Gene = gene,
          Motif_ID = motifs,
          Position_Class = "not_enough_motifs"
        ))
        next
      }
      
      for (motif in motifs) {
        
        sub <- gene_df[gene_df$Motif_ID == motif, ]
        
        c4  <- sub[sub$Group == "C4_PACMAD", ]
        c3p <- sub[sub$Group == "C3_PACMAD", ]
        c3b <- sub[sub$Group == "C3_BEP", ]
        
        c3_all <- rbind(c3p, c3b)
        
        c4_pos <- c4 %>% group_by(Species) %>% summarise(pos = min(Position), .groups="drop")
        c3_pos <- c3_all %>% group_by(Species) %>% summarise(pos = min(Position), .groups="drop")
        
        c4_cons <- get_consensus(c4_pos$pos)
        c3_cons <- if (nrow(c3_pos) > 0) get_consensus(c3_pos$pos) else NULL
        
        position_class <- "not_syntenic"
        
        if (!is.null(c3_cons) &&
            c4_cons$support >= support_threshold &&
            c3_cons$support >= support_threshold) {
          
          if (abs(c4_cons$pos - c3_cons$pos) <= position_tolerance) {
            position_class <- "syntenic"
          } else {
            position_class <- "c4_shifted"
          }
        }
        
        results <- rbind(results, data.frame(
          Gene = gene,
          Motif_ID = motif,
          Position_Class = position_class
        ))
      }
    }
    
    return(results)
  }
  
  # ---------------------------
  # Apply per class
  # ---------------------------
  final_list <- list()
  
  for (grp in c("C4_only", "PACMAD", "PACMAD_BEP")) {
    
    df_sub <- all_motifs %>% filter(Presence_Class == grp)
    if (nrow(df_sub) == 0) next
    
    pos_df <- get_positions(df_sub)
    
    if (grp == "C4_only") {
      
      tmp <- data.frame()
      
      for (gene in unique(pos_df$Gene)) {
        
        gene_df <- pos_df[pos_df$Gene == gene, ]
        motifs <- unique(gene_df$Motif_ID)
        
        if (length(motifs) <= 2) {
          tmp <- rbind(tmp, data.frame(
            Gene = gene,
            Motif_ID = motifs,
            Position_Class = "not_enough_motifs"
          ))
          next
        }
        
        for (motif in motifs) {
          
          sub <- gene_df[gene_df$Motif_ID == motif, ]
          
          c4_pos <- sub %>%
            group_by(Species) %>%
            summarise(pos = min(Position), .groups = "drop")
          
          tab <- table(c4_pos$pos)
          support <- max(tab) / sum(tab)
          
          position_class <- ifelse(
            support >= support_threshold,
            "syntenic",
            "not_syntenic"
          )
          
          tmp <- rbind(tmp, data.frame(
            Gene = gene,
            Motif_ID = motif,
            Position_Class = position_class
          ))
        }
      }
      
    } else {
      tmp <- compute_synteny(pos_df)
    }
    
    tmp$Presence_Class <- grp
    final_list[[grp]] <- tmp
  }
  
  final <- bind_rows(final_list)
  # ---------------------------
  # 🔥 Capture dropped genes
  # ---------------------------
  all_genes <- unique(all_motifs$Gene)
  final_genes <- unique(final$Gene)
  
  dropped_genes <- setdiff(all_genes, final_genes)
  # ---------------------------
  # 🔥 Add counts + position summary
  # ---------------------------
  pos_summary <- all_motifs %>%
    group_by(Gene, Motif_ID) %>%
    summarise(
      mean_position = mean(Start),
      n_sequences = n_distinct(Sequence),
      .groups = "drop"
    )
  
  final <- final %>%
    left_join(
      motif_presence %>% select(Gene, Motif_ID, n_c4, n_c3p, n_c3b),
      by = c("Gene", "Motif_ID")
    ) %>%
    left_join(pos_summary, by = c("Gene", "Motif_ID"))
  
  final$Dataset <- dataset_name
  
  return(list(
    final = final,
    dropped_motifs = dropped_motifs,
    dropped_genes = dropped_genes
  ))
}

#hopefully fixed function
process_dataset <- function(base_path, dataset_name,
                            min_c4_species = 4,
                            position_tolerance = 1,
                            support_threshold = 0.5) {
  
  library(dplyr)
  library(readr)
  
  setwd(base_path)
  
  # ---------------------------
  # Load
  # ---------------------------
  meme_df <- read_tsv("meme_filtered.tsv")
  fimo_df <- read.delim("fimo_filtered.tsv", stringsAsFactors = FALSE)
  
  all_motifs <- bind_rows(
    meme_df %>% select(Gene, Sequence, Motif_ID, Start, Strand),
    fimo_df %>% select(Gene, Sequence, Motif_ID, Start, Strand)
  )
  all_motifs <- all_motifs %>%
    mutate(
      Seq_clean = sub("\\|.*", "", Sequence)
    )
  
  # ---------------------------
  # Annotate
  # ---------------------------
  all_motifs <- all_motifs %>%
    mutate(
      Group = case_when(
        grepl("Sevir|Urofu|Pavag|ELECO|GRMZM|Pahal", Seq_clean) ~ "C4_PACMAD",
        grepl("OEL|Chala", Seq_clean) ~ "C3_PACMAD",
        grepl("Bradi|LOC", Seq_clean) ~ "C3_BEP",
        TRUE ~ "Unknown"
      ),
      Species = case_when(
        grepl("GRMZM", Seq_clean) ~ "Maize",
        grepl("Sevir", Seq_clean) ~ "Setaria_viridis",
        grepl("Urofu", Seq_clean) ~ "Urochloa_fusca",
        grepl("Pavag", Seq_clean) ~ "Paspalum_vaginatum",
        grepl("ELECO", Seq_clean) ~ "Eleusine_coracana",
        grepl("Pahal", Seq_clean) ~ "Panicum_hallii",
        grepl("OEL", Seq_clean) ~ "Dicanthelium",
        grepl("Chala", Seq_clean) ~ "Chasmanthium_latifolium",
        grepl("Bradi", Seq_clean) ~ "Brachypodium_distachyon",
        grepl("LOC", Seq_clean) ~ "Oryza_sativa",
        TRUE ~ "Unknown"
      )
    )
  
  # ---------------------------
  # Presence classification (FIXED)
  # ---------------------------
  motif_presence <- all_motifs %>%
    group_by(Gene, Motif_ID) %>%
    summarise(
      n_c4  = n_distinct(Species[Group == "C4_PACMAD"]),
      n_c3p = n_distinct(Species[Group == "C3_PACMAD"]),
      n_c3b = n_distinct(Species[Group == "C3_BEP"]),
      .groups = "drop"
    ) %>%
    mutate(
      Presence_Class = case_when(
        n_c4 >= min_c4_species & n_c3p == 0 & n_c3b == 0 ~ "C4_only",
        n_c4 > 0 & n_c3p > 0 & n_c3b == 0 ~ "PACMAD",
        n_c4 > 0 & n_c3p > 0 & n_c3b > 0 ~ "PACMAD_BEP",
        TRUE ~ "Other"
      )
    )
  # ---------------------------
  # 🔥 Capture dropped motifs
  # ---------------------------
  dropped_motifs <- motif_presence %>%
    filter(Presence_Class == "Other")
  # clean join (no duplicates)
  all_motifs <- all_motifs %>%
    left_join(
      motif_presence %>% select(Gene, Motif_ID, Presence_Class),
      by = c("Gene", "Motif_ID")
    )
  
  # ---------------------------
  # Position extraction (per class)
  # ---------------------------
  get_positions <- function(df_group) {
    
    arch_df <- df_group %>%
      arrange(Gene, Sequence, Start) %>%
      group_by(Gene, Sequence, Species, Group) %>%
      summarise(motif_order = list(Motif_ID), .groups = "drop")
    
    do.call(rbind, lapply(1:nrow(arch_df), function(i) {
      motifs <- arch_df$motif_order[[i]]
      
      do.call(rbind, lapply(unique(motifs), function(m) {
        data.frame(
          Gene = arch_df$Gene[i],
          Sequence = arch_df$Sequence[i],
          Species = arch_df$Species[i],
          Group = arch_df$Group[i],
          Motif_ID = m,
          Position = which(motifs == m)[1]
        )
      }))
    }))
  }
  
  # ---------------------------
  # Synteny (PACMAD + BEP FIXED)
  # ---------------------------
  compute_synteny <- function(df_positions) {
    
    results <- data.frame()
    
    get_consensus <- function(pos_vec) {
      if (length(pos_vec) == 0) return(list(pos = NA, support = 0))
      tab <- table(pos_vec)
      list(
        pos = as.numeric(names(tab)[which.max(tab)]),
        support = max(tab)/sum(tab)
      )
    }
    
    for (gene in unique(df_positions$Gene)) {
      
      gene_df <- df_positions[df_positions$Gene == gene, ]
      motifs <- unique(gene_df$Motif_ID)
      
      if (length(motifs) <= 2) {
        results <- rbind(results, data.frame(
          Gene = gene,
          Motif_ID = motifs,
          Position_Class = "not_enough_motifs"
        ))
        next
      }
      
      for (motif in motifs) {
        
        sub <- gene_df[gene_df$Motif_ID == motif, ]
        
        c4  <- sub[sub$Group == "C4_PACMAD", ]
        c3p <- sub[sub$Group == "C3_PACMAD", ]
        c3b <- sub[sub$Group == "C3_BEP", ]
        
        c3_all <- rbind(c3p, c3b)
        
        c4_pos <- c4 %>% group_by(Species) %>% summarise(pos = min(Position), .groups="drop")
        c3_pos <- c3_all %>% group_by(Species) %>% summarise(pos = min(Position), .groups="drop")
        
        c4_cons <- get_consensus(c4_pos$pos)
        c3_cons <- if (nrow(c3_pos) > 0) get_consensus(c3_pos$pos) else NULL
        
        position_class <- "not_syntenic"
        
        if (!is.null(c3_cons) &&
            c4_cons$support >= support_threshold &&
            c3_cons$support >= support_threshold) {
          
          if (abs(c4_cons$pos - c3_cons$pos) <= position_tolerance) {
            position_class <- "syntenic"
          } else {
            position_class <- "c4_shifted"
          }
        }
        
        results <- rbind(results, data.frame(
          Gene = gene,
          Motif_ID = motif,
          Position_Class = position_class
        ))
      }
    }
    
    return(results)
  }
  
  # ---------------------------
  # Apply per class
  # ---------------------------
  final_list <- list()
  
  for (grp in c("C4_only", "PACMAD", "PACMAD_BEP")) {
    
    df_sub <- all_motifs %>% filter(Presence_Class == grp)
    if (nrow(df_sub) == 0) next
    
    pos_df <- get_positions(df_sub)
    
    if (grp == "C4_only") {
      
      tmp <- data.frame()
      
      for (gene in unique(pos_df$Gene)) {
        
        gene_df <- pos_df[pos_df$Gene == gene, ]
        motifs <- unique(gene_df$Motif_ID)
        
        if (length(motifs) <= 2) {
          tmp <- rbind(tmp, data.frame(
            Gene = gene,
            Motif_ID = motifs,
            Position_Class = "not_enough_motifs"
          ))
          next
        }
        
        for (motif in motifs) {
          
          sub <- gene_df[gene_df$Motif_ID == motif, ]
          
          c4_pos <- sub %>%
            group_by(Species) %>%
            summarise(pos = min(Position), .groups = "drop")
          
          tab <- table(c4_pos$pos)
          support <- max(tab) / sum(tab)
          
          position_class <- ifelse(
            support >= support_threshold,
            "syntenic",
            "not_syntenic"
          )
          
          tmp <- rbind(tmp, data.frame(
            Gene = gene,
            Motif_ID = motif,
            Position_Class = position_class
          ))
        }
      }
      
    } else {
      tmp <- compute_synteny(pos_df)
    }
    
    tmp$Presence_Class <- grp
    final_list[[grp]] <- tmp
  }
  
  final <- bind_rows(final_list)
  # ---------------------------
  # 🔥 Capture dropped genes
  # ---------------------------
  all_genes <- unique(all_motifs$Gene)
  final_genes <- unique(final$Gene)
  
  dropped_genes <- setdiff(all_genes, final_genes)
  # ---------------------------
  # 🔥 Add counts + position summary
  # ---------------------------
  pos_summary <- all_motifs %>%
    group_by(Gene, Motif_ID) %>%
    summarise(
      mean_position = mean(Start),
      n_sequences = n_distinct(Sequence),
      .groups = "drop"
    )
  
  final <- final %>%
    left_join(
      motif_presence %>% select(Gene, Motif_ID, n_c4, n_c3p, n_c3b),
      by = c("Gene", "Motif_ID")
    ) %>%
    left_join(pos_summary, by = c("Gene", "Motif_ID"))
  
  final$Dataset <- dataset_name
  
  return(list(
    final = final,
    dropped_motifs = dropped_motifs,
    dropped_genes = dropped_genes
  ))
}
test_results <- process_dataset(
  "~/angeo_C4/cismotifs/meme_res", "Test"
)

control_results <- process_dataset(
  "~/angeo_C4/cismotifs/control", "Control"
)

view(control_results$final)
view(control_results$dropped_motifs)
library(dplyr)
library(ggplot2)

c4_only= control_results$final %>%
  filter(Presence_Class == "C4_only")

# check if they actually have C3 hits in raw FIMO
fimo_raw <- read_tsv("fimo_filtered.tsv")

fimo_raw %>%
  filter(Motif_ID %in% c4_only$Motif_ID) %>%
  distinct(Motif_ID) %>%
  nrow()

combined <- bind_rows(test_results$final, control_results$final)

# ---------------------------
# Global proportions (NOT grouped by Presence)
# ---------------------------
plot_df <- combined %>%
  count(Dataset, Presence_Class, Position_Class) %>%
  ungroup() %>%
  group_by(Dataset) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ---------------------------
# Clean factor order
# ---------------------------
plot_df$Presence_Class <- factor(
  plot_df$Presence_Class,
  levels = c("C4_only", "PACMAD", "PACMAD_BEP")
)

plot_df$Position_Class <- factor(
  plot_df$Position_Class,
  levels = c("syntenic", "c4_shifted", "not_syntenic","not_enough_motifs")
)

# ---------------------------
# Plot
# ---------------------------
ggplot(plot_df, aes(x = Presence_Class, y = prop, fill = Position_Class)) +
  
  geom_bar(stat = "identity") +
  
  facet_wrap(~Dataset) +
  
  scale_fill_manual(values = c(
    "syntenic" = "#4CAF50",       # green
    "c4_shifted" = "#2196F3",     # blue
    "not_syntenic" = "#9E9E9E",
    "not_enough_motifs"= "orange"# grey
  )) +
  
  theme_minimal(base_size = 16) +
  
  labs(
    title = "Motif Synteny Classification (Test vs Control)",
    x = "Presence Class",
    y = "Proportion",
    fill = "Position Class"
  ) +
  
  theme(
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

combined$Presence_Class
not_syntenic_df <- combined %>%
  filter(Position_Class == "not_syntenic",
         Presence_Class != "C4_only")

not_syntenic_test <- not_syntenic_df %>%
  filter(Dataset == "Test")

not_syntenic_control <- not_syntenic_df %>%
  filter(Dataset == "Control")

write_csv(not_syntenic_test, "not_syntenic_test.csv")
write_csv(not_syntenic_control, "not_syntenic_control.csv")

C4_only_df <- combined %>%
  filter(Presence_Class == "C4_only")

C4_only_test <- C4_only_df %>%
  filter(Dataset == "Test")

C4_only_control <- C4_only_df %>%
  filter(Dataset == "Control")

not_syntenic_test <- not_syntenic_test %>%
  mutate(Motif_UID = paste(Gene, Motif_ID, sep = "_"))

not_syntenic_control <- not_syntenic_control %>%
  mutate(Motif_UID = paste(Gene, Motif_ID, sep = "_"))

test_motifs <- unique(not_syntenic_test$Motif_UID)
control_motifs <- unique(not_syntenic_control$Motif_UID)

write_lines(test_motifs, "test_motif_uid.txt")
write_lines(control_motifs, "control_motif_uid.txt")


C4_shifted_df <- combined %>%
  filter(Position_Class == "c4_shifted",
         Presence_Class != "C4_only")

C4_shifted_test <- C4_shifted_df %>%
  filter(Dataset == "Test")

C4_shifted_control <- C4_shifted_df %>%
  filter(Dataset == "Control")

write_csv(not_syntenic_test, "C4_shifted_test.csv")
write_csv(not_syntenic_control, "C4_shifted_control.csv")

C4_shifted_test <- C4_shifted_test %>%
  mutate(Motif_UID = paste(Gene, Motif_ID, sep = "_"))

C4_shifted_control <- C4_shifted_control %>%
  mutate(Motif_UID = paste(Gene, Motif_ID, sep = "_"))

C4_shifted_test_motifs <- unique(C4_shifted_test$Motif_UID)
C4_shifted_control_motifs <- unique(C4_shifted_control$Motif_UID)

write_lines(C4_shifted_test_motifs, "test_C4_shifted_motif_uid.txt")
write_lines(C4_shifted_control_motifs, "control_C4_shifted_motif_uid.txt")

setwd("~/angeo_C4/cismotifs/meme_res")
meme_df <- read_tsv("meme_positions.csv")
fimo_df <- read.delim("fimo_c3_positions.tsv", stringsAsFactors = FALSE)
all_motifs <- bind_rows(
  meme_df %>% select(Gene, Sequence, Motif_ID, Start, Strand),
  fimo_df %>% select(Gene, Sequence, Motif_ID, Start, Strand)
)

unique(all_motifs$Gene)
all_motifs %>% 
  select(Gene, Motif_ID) %>% 
  n_distinct() #1427

all_motifs <- all_motifs %>%
  mutate(
    Group = case_when(
      grepl("^(Sevir|Urofu|Pavag|ELECO|GRMZM|Pahal)", Sequence) ~ "C4_PACMAD",
      grepl("^(OEL|Chala)", Sequence) ~ "C3_PACMAD",
      grepl("^(Bradi|LOC)", Sequence) ~ "C3_BEP",
      TRUE ~ "Unknown"
    ),
    Species = case_when(
      grepl("^GRMZM", Sequence) ~ "Maize",
      grepl("^Sevir", Sequence) ~ "Setaria_viridis",
      grepl("^Urofu", Sequence) ~ "Urochloa_fusca",
      grepl("^Pavag", Sequence) ~ "Paspalum_vaginatum",
      grepl("^ELECO", Sequence) ~ "Eleusine_coracana",
      grepl("^Pahal", Sequence) ~ "Panicum_hallii",
      grepl("^OEL", Sequence) ~ "Dicanthelium",
      grepl("^Chala", Sequence) ~ "Chasmanthum_laxum",
      grepl("^Bradi", Sequence) ~ "Brachypodium_distachyon",
      grepl("^LOC", Sequence) ~ "Oryza_sativa",
      TRUE ~ "Unknown"
    )
  )

table(all_motifs$Group)
table(all_motifs$Species)

motif_presence <- all_motifs %>%
  group_by(Gene, Motif_ID) %>%
  summarise(
    n_c4  = n_distinct(Species[Group == "C4_PACMAD"]),
    n_c3p = n_distinct(Species[Group == "C3_PACMAD"]),
    n_c3b = n_distinct(Species[Group == "C3_BEP"]),
    .groups = "drop"
  ) %>%
  mutate(
    Presence_Class = case_when(
      n_c4 >= 4 & n_c3p == 0 ~ "C4_only",
      n_c4 > 0 & n_c3p > 0 & n_c3b == 0 ~ "PACMAD",
      n_c4 > 0 & n_c3p > 0 & n_c3b > 0 ~ "PACMAD_BEP",
      TRUE ~ "onlyinlessthas3C4"
    )
  )

table(motif_presence$Presence_Class)

all_motifs <- all_motifs %>%
  left_join(motif_presence, by = c("Gene", "Motif_ID"))

all_motifs %>% count(Presence_Class)

table(all_motifs$Presence_Class)
arch_df <- df_sub %>%
  arrange(Gene, Sequence, Start) %>%
  group_by(Gene, Sequence) %>%
  summarise(motif_order = list(Motif_ID), .groups = "drop")

lost_genes <- setdiff(
  unique(all_motifs$Gene),
  unique(arch_df$Gene)
)

length(lost_genes)
head(lost_genes)

all_motifs %>%
  filter(Gene %in% lost_genes) %>%
  count(Presence_Class)

pos_df <- all_motifs %>%
  arrange(Gene, Sequence, Start) %>%
  group_by(Gene, Sequence) %>%
  mutate(Position = row_number()) %>%
  ungroup()

pos_df <- pos_df %>%
  left_join(motif_presence, by = c("Gene", "Motif_ID"))

pos_df <- pos_df %>%
  select(-n_c4.x, -n_c3p.x, -n_c3b.x, -Presence_Class.x) %>%
  rename(
    n_c4 = n_c4.y,
    n_c3p = n_c3p.y,
    n_c3b = n_c3b.y,
    Presence_Class = Presence_Class.y
  )
length(unique(pos_df$Gene))

df_pacmad <- pos_df %>% filter(Presence_Class == "PACMAD")
df_c4     <- pos_df %>% filter(Presence_Class == "C4_only")
df_bep    <- pos_df %>% filter(Presence_Class == "PACMAD_BEP")

compute_synteny_general <- function(df, support_threshold = 0.6) {
  
  get_consensus <- function(pos_vec) {
    if (length(pos_vec) == 0) return(list(pos = NA, support = 0))
    
    tab <- table(pos_vec)
    best <- as.numeric(names(tab)[which.max(tab)])
    support <- max(tab) / sum(tab)
    
    list(pos = best, support = support)
  }
  
  df %>%
    group_by(Gene, Motif_ID) %>%
    summarise(
      c4_pos = list(Position[Group == "C4_PACMAD"]),
      c3_pos = list(Position[Group == "C3_PACMAD"]),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      c4_cons = list(get_consensus(c4_pos)),
      c3_cons = list(get_consensus(c3_pos)),
      
      c4_pos_val = c4_cons$pos,
      c4_support = c4_cons$support,
      
      c3_pos_val = c3_cons$pos,
      c3_support = c3_cons$support,
      
      Position_Class = case_when(
        
        # C4 vs C3 shift
        length(c3_pos) > 0 &&
          c4_support >= support_threshold &&
          c3_support >= support_threshold &&
          abs(c4_pos_val - c3_pos_val) > 1 ~ "c4_shifted",
        
        # Syntenic (within C4 or general)
        c4_support >= support_threshold ~ "syntenic",
        
        TRUE ~ "not_syntenic"
      )
    ) %>%
    ungroup()
}
compute_c4_only <- function(df, support_threshold = 0.6) {
  
  get_consensus <- function(pos_vec) {
    tab <- table(pos_vec)
    best <- as.numeric(names(tab)[which.max(tab)])
    support <- max(tab) / sum(tab)
    list(pos = best, support = support)
  }
  
  df %>%
    group_by(Gene, Motif_ID) %>%
    summarise(
      pos = list(Position),
      n_instances = length(Position),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      cons = list(get_consensus(pos)),
      support = cons$support,
      
      Position_Class = case_when(
        
        # 🔥 YOUR RULE
        n_instances <= 2 ~ "not_enough_motifs",
        
        support >= support_threshold ~ "syntenic",
        
        TRUE ~ "not_syntenic"
      )
    ) %>%
    ungroup() %>%
    mutate(Presence_Class = "C4_only")
}

compute_pacmad <- function(df, support_threshold = 0.6) {
  
  get_consensus <- function(pos_vec) {
    if (length(pos_vec) == 0) return(list(pos = NA, support = 0))
    tab <- table(pos_vec)
    best <- as.numeric(names(tab)[which.max(tab)])
    support <- max(tab) / sum(tab)
    list(pos = best, support = support)
  }
  
  df %>%
    group_by(Gene, Motif_ID) %>%
    summarise(
      c4_pos = list(Position[Group == "C4_PACMAD"]),
      c3_pos = list(Position[Group == "C3_PACMAD"]),
      n_instances = length(Position),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      c4_cons = list(get_consensus(c4_pos)),
      c3_cons = list(get_consensus(c3_pos)),
      
      c4_pos_val = c4_cons$pos,
      c3_pos_val = c3_cons$pos,
      c4_sup = c4_cons$support,
      c3_sup = c3_cons$support,
      
      Position_Class = case_when(
        
        # 🔥 YOUR RULE FIRST
        n_instances <= 2 ~ "not_enough_motifs",
        
        # shift
        c4_sup >= support_threshold & c3_sup >= support_threshold &
          abs(c4_pos_val - c3_pos_val) > 1 ~ "c4_shifted",
        
        # conserved
        c4_sup >= support_threshold & c3_sup >= support_threshold ~ "syntenic",
        
        TRUE ~ "not_syntenic"
      )
    ) %>%
    ungroup() %>%
    mutate(Presence_Class = "PACMAD")
}

compute_bep <- function(df, support_threshold = 0.6) {
  compute_pacmad(df, support_threshold) %>%
    mutate(Presence_Class = "PACMAD_BEP")
}

res_c4     <- compute_c4_only(df_c4)
res_pacmad <- compute_pacmad(df_pacmad)
res_bep    <- compute_bep(df_bep)

res_c4     <- compute_c4_only(df_c4) %>%
  mutate(Presence_Class = "C4_only")

res_pacmad <- compute_pacmad(df_pacmad) %>%
  mutate(Presence_Class = "PACMAD")

res_bep    <- compute_pacmad(df_bep) %>%
  mutate(Presence_Class = "PACMAD_BEP")

final <- bind_rows(res_c4, res_pacmad, res_bep)
final$Dataset = "Test"

plot_df <- final %>%
  count(Dataset, Presence_Class, Position_Class) %>%
  group_by(Dataset) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

plot_df$Presence_Class <- factor(
  plot_df$Presence_Class,
  levels = c("C4_only", "PACMAD", "PACMAD_BEP")
)

plot_df$Position_Class <- factor(
  plot_df$Position_Class,
  levels = c("syntenic", "c4_shifted", "not_syntenic", "not_enough_motifs")
)

setwd("~/angeo_C4/cismotifs/control")
meme_df <- read_tsv("meme_positions.csv")
fimo_df <- read.delim("fimo_c3_positions.tsv", stringsAsFactors = FALSE)

all_motifs <- bind_rows(
  meme_df %>% select(Gene, Sequence, Motif_ID, Start, Strand),
  fimo_df %>% select(Gene, Sequence, Motif_ID, Start, Strand)
)
all_motifs <- all_motifs %>%
  mutate(
    Group = case_when(
      grepl("^(Sevir|Urofu|Pavag|ELECO|GRMZM|Pahal)", Sequence) ~ "C4_PACMAD",
      grepl("^(OEL|Chala)", Sequence) ~ "C3_PACMAD",
      grepl("^(Bradi|LOC)", Sequence) ~ "C3_BEP",
      TRUE ~ "Unknown"
    ),
    Species = case_when(
      grepl("^GRMZM", Sequence) ~ "Maize",
      grepl("^Sevir", Sequence) ~ "Setaria_viridis",
      grepl("^Urofu", Sequence) ~ "Urochloa_fusca",
      grepl("^Pavag", Sequence) ~ "Paspalum_vaginatum",
      grepl("^ELECO", Sequence) ~ "Eleusine_coracana",
      grepl("^Pahal", Sequence) ~ "Panicum_hallii",
      grepl("^OEL", Sequence) ~ "Dicanthelium",
      grepl("^Chala", Sequence) ~ "Chasmanthium_latifolium",
      grepl("^Bradi", Sequence) ~ "Brachypodium_distachyon",
      grepl("^LOC", Sequence) ~ "Oryza_sativa",
      TRUE ~ "Unknown"
    )
  )
motif_presence <- all_motifs %>%
  group_by(Gene, Motif_ID) %>%
  summarise(
    n_c4  = n_distinct(Species[Group == "C4_PACMAD" & !is.na(Species)]),
    n_c3p = n_distinct(Species[Group == "C3_PACMAD" & !is.na(Species)]),
    n_c3b = n_distinct(Species[Group == "C3_BEP" & !is.na(Species)]),
    
    # 🔥 NEW: direct presence flags (much safer)
    has_c4  = any(Group == "C4_PACMAD"),
    has_c3p = any(Group == "C3_PACMAD"),
    has_c3b = any(Group == "C3_BEP"),
    
    .groups = "drop"
  ) %>%
  mutate(
    Presence_Class = case_when(
      
      # 🔥 STRICT definition
      n_c4 >= 4 & !has_c3p & !has_c3b ~ "C4_only",
      
      has_c4 & has_c3p & !has_c3b ~ "PACMAD",
      
      has_c4 & has_c3p & has_c3b ~ "PACMAD_BEP",
      
      TRUE ~ "onlyinlessthas3C4"
    )
  )

# attach class FIRST
all_motifs <- all_motifs %>%
  left_join(
    motif_presence %>%
      select(Gene, Motif_ID, Presence_Class),
    by = c("Gene", "Motif_ID")
  )

# split
df_c4     <- all_motifs %>% filter(Presence_Class == "C4_only")
df_pacmad <- all_motifs %>% filter(Presence_Class == "PACMAD")
df_bep    <- all_motifs %>% filter(Presence_Class == "PACMAD_BEP")
table(df_c4$Group)
assign_positions <- function(df) {
  df %>%
    arrange(Gene, Sequence, Start) %>%
    group_by(Gene, Sequence) %>%
    mutate(Position = row_number()) %>%
    ungroup()
}

df_c4     <- assign_positions(df_c4)
df_pacmad <- assign_positions(df_pacmad)
df_bep    <- assign_positions(df_bep)

res_c4     <- compute_c4_only(df_c4) %>%
  mutate(Presence_Class = "C4_only")

res_pacmad <- compute_pacmad(df_pacmad) %>%
  mutate(Presence_Class = "PACMAD")

res_bep    <- compute_bep(df_bep) %>%
  mutate(Presence_Class = "PACMAD_BEP")

final_control <- bind_rows(res_c4, res_pacmad, res_bep) %>%
  mutate(Dataset = "Control")
combined <- bind_rows(final, final_control)

plot_df <- combined %>%
  count(Dataset, Presence_Class, Position_Class) %>%
  group_by(Dataset) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
ggplot(plot_df, aes(x = Presence_Class, y = prop, fill = Position_Class)) +
  
  geom_bar(stat = "identity") +
  
  facet_wrap(~Dataset) +
  
  scale_fill_manual(values = c(
    "syntenic" = "#4CAF50",
    "c4_shifted" = "#2196F3",
    "not_syntenic" = "#9E9E9E",
    "not_enough_motifs" = "#FFA500"
  ))+
  
  theme_minimal(base_size = 16) +
  
  labs(
    title = "Motif Synteny Classification (Test vs Control)",
    x = "Presence Class",
    y = "Proportion",
    fill = "Position Class"
  ) +
  
  theme(
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
final %>% 
  select(Gene, Motif_ID) %>% 
  n_distinct()
