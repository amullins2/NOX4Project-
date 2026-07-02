## NOX4 islet analysis confocal vs slide scanner, NOX4 vs NPC 

setwd("~/Desktop")

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(forcats)

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

paths <- list(
  Confocal_NOX4      = "16b NOX4_IdentifySecondaryObjects.xlsx",
  Confocal_NPC       = "16b NPC_IdentifySecondaryObjects.xlsx",
  SlideScanner_NOX4  = "Slide Scanner NOX4_IdentifySecondaryObjects.xlsx",
  SlideScanner_NPC   = "Slide Scanner NPC_IdentifySecondaryObjects.xlsx"
)

# donor thresholds (Average row from Thresholding.xlsx, lower limit only)
thresholds <- tibble::tribble(
  ~modality,        ~insulin_lower, ~glucagon_lower,
  "Confocal",              0.2661,          0.13356,
  "Slide Scanner",         0.03308,         0.02786
)

#channel order differs by modality
channel_map <- list(
  "Confocal"      = c(target = "C1", insulin = "C2", dapi = "C3", glucagon = "C4"),
  "Slide Scanner" = c(dapi = "C1", target = "C2", glucagon = "C3", insulin = "C4")
)

load_cp_file <- function(path, modality, stain, donor_id = "Donor1") {
  raw <- read_excel(path)
  cm <- channel_map[[modality]]
  raw %>%
    transmute(
      donor_id, modality, stain,
      ImageNumber, ObjectNumber,
      FileName     = FileName_C1,
      target_intensity   = .data[[paste0("Intensity_MeanIntensity_", cm[["target"]])]],
      insulin_intensity  = .data[[paste0("Intensity_MeanIntensity_", cm[["insulin"]])]],
      glucagon_intensity = .data[[paste0("Intensity_MeanIntensity_", cm[["glucagon"]])]]
    ) %>%
    mutate(islet_number = str_extract(FileName, "Islet\\s*\\d+") %>% str_extract("\\d+") %>% as.numeric())
}

df <- bind_rows(
  load_cp_file(paths$Confocal_NOX4,     "Confocal",      "NOX4"),
  load_cp_file(paths$Confocal_NPC,      "Confocal",      "NPC"),
  load_cp_file(paths$SlideScanner_NOX4, "Slide Scanner", "NOX4"),
  load_cp_file(paths$SlideScanner_NPC,  "Slide Scanner", "NPC")
)

stopifnot(all(!is.na(df$islet_number)))

#classify cells 
df <- df %>%
  left_join(thresholds, by = "modality") %>%
  mutate(
    insulin_pos  = insulin_intensity  >= insulin_lower,
    glucagon_pos = glucagon_intensity >= glucagon_lower,
    cell_type = case_when(
      insulin_pos & glucagon_pos  ~ "Bi-hormonal",
      insulin_pos                 ~ "Beta",
      glucagon_pos                ~ "Alpha",
      TRUE                        ~ "Other"
    ),
    cell_type = factor(cell_type, levels = c("Beta", "Alpha", "Bi-hormonal", "Other")),
    modality  = factor(modality, levels = c("Confocal", "Slide Scanner")),
    stain     = factor(stain, levels = c("NOX4", "NPC")),
    islet_id  = paste(modality, stain, islet_number, sep = "_")
  )

write.csv(df, "tables/cell_level_data.csv", row.names = FALSE)

## background correction

MIN_BG_CELLS <- 3

background_insulin_global <- df %>%
  filter(stain == "NPC", cell_type %in% c("Beta", "Bi-hormonal")) %>%
  group_by(modality) %>%
  summarise(bg_insulin_global = mean(target_intensity), n_bg_global = n(), .groups = "drop")

background_glucagon_global <- df %>%
  filter(stain == "NPC", cell_type %in% c("Alpha", "Bi-hormonal")) %>%
  group_by(modality) %>%
  summarise(bg_glucagon_global = mean(target_intensity), n_bg_global = n(), .groups = "drop")

npc_islet_bg <- df %>%
  filter(stain == "NPC") %>%
  group_by(modality, islet_number) %>%
  summarise(
    bg_insulin_local  = mean(target_intensity[cell_type %in% c("Beta", "Bi-hormonal")]),
    n_insulin_local   = sum(cell_type %in% c("Beta", "Bi-hormonal")),
    bg_glucagon_local = mean(target_intensity[cell_type %in% c("Alpha", "Bi-hormonal")]),
    n_glucagon_local  = sum(cell_type %in% c("Alpha", "Bi-hormonal")),
    .groups = "drop"
  )

nearest_npc_islet <- function(mod, islet_num, candidates) {
  pool <- candidates %>% filter(modality == mod)
  if (nrow(pool) == 0) return(NA_real_)
  pool$islet_number[which.min(abs(pool$islet_number - islet_num))]
}

nox4_insulin_pos <- df %>%
  filter(stain == "NOX4", cell_type %in% c("Beta", "Bi-hormonal")) %>%
  rowwise() %>%
  mutate(nearest_islet = nearest_npc_islet(modality, islet_number, npc_islet_bg)) %>%
  ungroup() %>%
  left_join(npc_islet_bg %>% select(modality, islet_number, bg_insulin_local, n_insulin_local),
            by = c("modality" = "modality", "nearest_islet" = "islet_number")) %>%
  left_join(background_insulin_global, by = "modality") %>%
  mutate(
    bg_used    = if_else(!is.na(n_insulin_local) & n_insulin_local >= MIN_BG_CELLS, bg_insulin_local, bg_insulin_global),
    bg_source  = if_else(!is.na(n_insulin_local) & n_insulin_local >= MIN_BG_CELLS, "local", "global fallback"),
    corrected_intensity = pmax(target_intensity - bg_used, 0)
  )

nox4_glucagon_pos <- df %>%
  filter(stain == "NOX4", cell_type %in% c("Alpha", "Bi-hormonal")) %>%
  rowwise() %>%
  mutate(nearest_islet = nearest_npc_islet(modality, islet_number, npc_islet_bg)) %>%
  ungroup() %>%
  left_join(npc_islet_bg %>% select(modality, islet_number, bg_glucagon_local, n_glucagon_local),
            by = c("modality" = "modality", "nearest_islet" = "islet_number")) %>%
  left_join(background_glucagon_global, by = "modality") %>%
  mutate(
    bg_used    = if_else(!is.na(n_glucagon_local) & n_glucagon_local >= MIN_BG_CELLS, bg_glucagon_local, bg_glucagon_global),
    bg_source  = if_else(!is.na(n_glucagon_local) & n_glucagon_local >= MIN_BG_CELLS, "local", "global fallback"),
    corrected_intensity = pmax(target_intensity - bg_used, 0)
  )

message("Insulin-positive NOX4 cells -- background source:")
print(table(nox4_insulin_pos$bg_source))
message("Glucagon-positive NOX4 cells -- background source:")
print(table(nox4_glucagon_pos$bg_source))

## per-islet composition graph 4
islet_composition <- df %>%
  group_by(donor_id, modality, stain, islet_id, islet_number) %>%
  summarise(
    n_cells        = n(),
    beta_cells     = sum(cell_type == "Beta"),
    alpha_cells    = sum(cell_type == "Alpha"),
    bihormonal_cells = sum(cell_type == "Bi-hormonal"),
    other_cells    = sum(cell_type == "Other"),
    .groups = "drop"
  ) %>%
  mutate(beta_alpha_ratio = if_else(alpha_cells > 0, beta_cells / alpha_cells, NA_real_))

n_undefined_ratio <- sum(is.na(islet_composition$beta_alpha_ratio))
message(sprintf("%d of %d islets have zero alpha cells -- excluded from Graph 4.", n_undefined_ratio, nrow(islet_composition)))

message("Zero-alpha islets, by modality x stain batch:")
print(
  islet_composition %>%
    group_by(modality, stain) %>%
    summarise(n_islets = n(), n_zero_alpha = sum(alpha_cells == 0),
              pct_zero_alpha = round(100 * mean(alpha_cells == 0), 1), .groups = "drop")
)

write.csv(islet_composition, "tables/islet_composition.csv", row.names = FALSE)

## per-islet mean corrected NOX4, by cell type graph 3 
islet_expression <- bind_rows(
    nox4_insulin_pos %>% filter(cell_type == "Beta")  %>% mutate(pop = "Beta"),
    nox4_glucagon_pos %>% filter(cell_type == "Alpha") %>% mutate(pop = "Alpha")
  ) %>%
  group_by(donor_id, modality, islet_id, islet_number, cell_type = pop) %>%
  summarise(mean_corrected_intensity = mean(corrected_intensity), n_cells = n(), .groups = "drop")

write.csv(islet_expression, "tables/islet_expression_by_celltype.csv", row.names = FALSE)

## graph 1: NOX4 in insulin-positive cells
p1 <- ggplot(nox4_insulin_pos, aes(x = 1, y = corrected_intensity)) +
  geom_violin(fill = "grey85", trim = FALSE, alpha = 0.5) +
  geom_quasirandom(size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
  facet_wrap(~ modality, scales = "free_y") +
  labs(y = "NOX4 intensity\n(background-subtracted)", x = NULL,
       title = "NOX4 expression in insulin-positive cells",
       caption = "Panels have different y-axis scales -- do not compare magnitude across modalities") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.caption = element_text(hjust = 0, size = 9, colour = "grey40"))

ggsave("figures/Graph1_NOX4_insulin_positive_corrected.png", p1, width = 9, height = 6, dpi = 300)

## graph 2: NOX4 in glucagon-positive cells 
p2 <- ggplot(nox4_glucagon_pos, aes(x = 1, y = corrected_intensity)) +
  geom_violin(fill = "grey85", trim = FALSE, alpha = 0.5) +
  geom_quasirandom(size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
  facet_wrap(~ modality, scales = "free_y") +
  labs(y = "NOX4 intensity\n(background-subtracted)", x = NULL,
       title = "NOX4 expression in glucagon-positive cells",
       caption = "Panels have different y-axis scales -- do not compare magnitude across modalities") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.caption = element_text(hjust = 0, size = 9, colour = "grey40"))

ggsave("figures/Graph2_NOX4_glucagon_positive_corrected.png", p2, width = 9, height = 6, dpi = 300)

## graph 3: per-islet NOX4 by cell type 
p3 <- ggplot(islet_expression, aes(cell_type, mean_corrected_intensity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ modality, scales = "free_y") +
  labs(y = "Mean NOX4 intensity per islet\n(background-subtracted)", x = "Cell type",
       title = "Per-islet NOX4 expression by cell type") +
  theme_classic(base_size = 13)

ggsave("figures/Graph3_per_islet_expression_by_celltype.png", p3, width = 9, height = 6, dpi = 300)

## graph 4: beta:alpha ratio per islet
p4 <- ggplot(islet_composition %>% filter(!is.na(beta_alpha_ratio)), aes(modality, beta_alpha_ratio)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  scale_y_log10() +
  labs(y = "Beta:Alpha ratio (log scale)", x = NULL,
       title = "Islet composition: Beta:Alpha ratio per islet",
       caption = paste0(n_undefined_ratio, " islet(s) with zero alpha cells excluded -- see console for breakdown")) +
  theme_classic(base_size = 13) +
  theme(plot.caption = element_text(hjust = 0, size = 9, colour = "grey40"))

ggsave("figures/Graph4_BetaAlpha_ratio_by_modality.png", p4, width = 8, height = 6, dpi = 300)

## outputs: tables/cell_level_data.csv, tables/islet_composition.csv,
## tables/islet_expression_by_celltype.csv, figures/Graph1-4_*.png
