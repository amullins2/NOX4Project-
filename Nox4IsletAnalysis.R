## NOX4 islet analysis: confocal vs slide scanner, NOX4 vs. NPC 

setwd("~/Desktop")

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(forcats)
library(lme4)
library(lmerTest)   
library(emmeans)

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

# if read_excel fails, resave the file in Excel first (strict OOXML issue)
paths <- list(
  Confocal_NOX4      = "16b NOX4_IdentifySecondaryObjects.xlsx",
  Confocal_NPC       = "16b NPC_IdentifySecondaryObjects.xlsx",
  SlideScanner_NOX4  = "Slide Scanner NOX4_IdentifySecondaryObjects.xlsx",
  SlideScanner_NPC   = "Slide Scanner NPC_IdentifySecondaryObjects.xlsx"
)

# donor thresholds (lower limit only)
thresholds <- tibble::tribble(
  ~modality,        ~insulin_lower, ~glucagon_lower,
  "Confocal",              0.2661,          0.13356,
  "Slide Scanner",         0.01648,         0.01524
)

threshold_overrides <- tibble::tribble(
  ~modality,        ~stain, ~islet_number, ~insulin_lower_ov, ~glucagon_lower_ov,
  "Slide Scanner",  "NPC",  5,             0.0168,            0.0136,
  "Confocal",       "NOX4", 1,             0.145,             0.0464,
  "Confocal",       "NOX4", 24,            0.178,             0.0325,
  "Confocal",       "NOX4", 26,            0.2759,            0.063,
  "Confocal",       "NPC",  5,             0.1571,            0.074,
  "Slide Scanner",  "NOX4", 16,            0.0369,            0.0161,
  "Slide Scanner",  "NOX4", 17,            0.034,             0.0169
)

# channel order differs by modality
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

# classify cells
df <- df %>%
  left_join(thresholds, by = "modality") %>%
  left_join(threshold_overrides, by = c("modality", "stain", "islet_number")) %>%
  mutate(
    used_override  = !is.na(insulin_lower_ov),
    insulin_lower  = coalesce(insulin_lower_ov, insulin_lower),
    glucagon_lower = coalesce(glucagon_lower_ov, glucagon_lower),
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
  ) %>%
  select(-insulin_lower_ov, -glucagon_lower_ov)

message(sprintf("%d cells used a per-islet override threshold.", sum(df$used_override)))

# excluded
excluded_islets <- tibble::tribble(
  ~modality,       ~stain, ~islet_number,
  "Slide Scanner", "NOX4", 16,
  "Slide Scanner", "NOX4", 17
)
n_before <- nrow(df)
df <- df %>% anti_join(excluded_islets, by = c("modality", "stain", "islet_number"))
message(sprintf("Excluded %d cells from Islets 16/17 (Slide Scanner NOX4).", n_before - nrow(df)))

write.csv(df, "tables/cell_level_data.csv", row.names = FALSE)

#background correction

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

#per-islet composition 
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

#per-islet mean corrected NOX4, by cell type 
islet_expression <- bind_rows(
    nox4_insulin_pos %>% filter(cell_type == "Beta")  %>% mutate(pop = "Beta"),
    nox4_glucagon_pos %>% filter(cell_type == "Alpha") %>% mutate(pop = "Alpha")
  ) %>%
  group_by(donor_id, modality, islet_id, islet_number, cell_type = pop) %>%
  summarise(mean_corrected_intensity = mean(corrected_intensity), n_cells = n(), .groups = "drop")

write.csv(islet_expression, "tables/islet_expression_by_celltype.csv", row.names = FALSE)

lmm_data <- bind_rows(
  nox4_insulin_pos  %>% filter(cell_type == "Beta")  %>% mutate(pop = "Beta"),
  nox4_glucagon_pos %>% filter(cell_type == "Alpha") %>% mutate(pop = "Alpha")
) %>%
  mutate(pop = factor(pop, levels = c("Beta", "Alpha")))

run_lmm <- function(mod_name) {
  d <- lmm_data %>% filter(modality == mod_name)
  m <- lmer(corrected_intensity ~ pop + (1 | islet_id), data = d)
  em <- emmeans(m, pairwise ~ pop, adjust = "none")
  list(model = m, emmeans = em, p_raw = summary(em$contrasts)$p.value)
}

lmm_confocal      <- run_lmm("Confocal")
lmm_slidescanner  <- run_lmm("Slide Scanner")

p_raw <- c(Confocal = lmm_confocal$p_raw, `Slide Scanner` = lmm_slidescanner$p_raw)
p_bonferroni <- p.adjust(p_raw, method = "bonferroni")

stats_summary <- tibble::tibble(
  modality      = names(p_raw),
  comparison    = "Alpha vs Beta (NOX4 corrected intensity)",
  p_raw         = p_raw,
  p_bonferroni  = p_bonferroni
)

message("Alpha vs Beta NOX4 expression, LMM with islet as random effect:")
print(stats_summary)
write.csv(stats_summary, "tables/stats_alpha_vs_beta.csv", row.names = FALSE)


comp_data <- islet_composition %>%
  filter(!is.na(beta_alpha_ratio), beta_alpha_ratio > 0)

n_excluded_comp <- nrow(islet_composition) - nrow(comp_data)
message(sprintf("%d islets excluded from composition model (zero alpha or beta cells).", n_excluded_comp))

m_comp <- lmer(log(beta_alpha_ratio) ~ modality + (1 | stain), data = comp_data)
comp_summary <- summary(m_comp)$coefficients
p_comp <- comp_summary["modalitySlide Scanner", "Pr(>|t|)"]
message("Beta:Alpha ratio, Confocal vs Slide Scanner (log scale):")
print(comp_summary)
write.csv(as.data.frame(comp_summary), "tables/stats_composition.csv", row.names = TRUE)


confocal_islets      <- sort(unique(df$islet_id[df$modality == "Confocal"]))
slidescanner_islets  <- sort(unique(df$islet_id[df$modality == "Slide Scanner"]))
confocal_pal     <- colorRampPalette(c("#08306B", "#9ECAE1"))(length(confocal_islets))
slidescanner_pal <- colorRampPalette(c("#67000D", "#FC9272"))(length(slidescanner_islets))
names(confocal_pal)     <- confocal_islets
names(slidescanner_pal) <- slidescanner_islets
islet_colours <- c(confocal_pal, slidescanner_pal)


p1 <- ggplot(nox4_insulin_pos, aes(x = 1, y = corrected_intensity, colour = islet_id, shape = cell_type)) +
  geom_violin(aes(shape = NULL), fill = "grey85", colour = "grey40", trim = FALSE, alpha = 0.5) +
  geom_quasirandom(size = 1, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  stat_summary(aes(shape = NULL), fun = mean, geom = "point", shape = 18, size = 3, colour = "black") +
  facet_wrap(~ modality, scales = "free_y") +
  scale_colour_manual(values = islet_colours, guide = "none") +
  scale_shape_manual(values = c(Beta = 16, "Bi-hormonal" = 17), drop = TRUE) +
  labs(y = "NOX4 intensity", x = NULL, shape = "Cell type",
       title = "NOX4 expression in insulin-positive cells") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("figures/Graph1_NOX4_insulin_positive_corrected.png", p1, width = 9, height = 6, dpi = 300)


p2 <- ggplot(nox4_glucagon_pos, aes(x = 1, y = corrected_intensity, colour = islet_id, shape = cell_type)) +
  geom_violin(aes(shape = NULL), fill = "grey85", colour = "grey40", trim = FALSE, alpha = 0.5) +
  geom_quasirandom(size = 1, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  stat_summary(aes(shape = NULL), fun = mean, geom = "point", shape = 18, size = 3, colour = "black") +
  facet_wrap(~ modality, scales = "free_y") +
  scale_colour_manual(values = islet_colours, guide = "none") +
  scale_shape_manual(values = c(Alpha = 16, "Bi-hormonal" = 17), drop = TRUE) +
  labs(y = "NOX4 intensity", x = NULL, shape = "Cell type",
       title = "NOX4 expression in glucagon-positive cells") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("figures/Graph2_NOX4_glucagon_positive_corrected.png", p2, width = 9, height = 6, dpi = 300)


sig_stars <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else NA_character_
}

p3_annot <- tibble::tibble(
  modality = factor(c("Confocal", "Slide Scanner"), levels = levels(islet_expression$modality)),
  label = c(sig_stars(p_bonferroni["Confocal"]), sig_stars(p_bonferroni["Slide Scanner"]))
) %>%
  filter(!is.na(label))

p3 <- ggplot(islet_expression, aes(cell_type, mean_corrected_intensity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(aes(colour = islet_id, shape = cell_type), width = 0.15, size = 2.2, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ modality, scales = "free_y") +
  scale_shape_manual(values = c(Beta = 16, Alpha = 17), guide = "none") +
  scale_colour_manual(values = islet_colours, guide = "none") +
  labs(y = "Mean NOX4 intensity per islet", x = "Cell type",
       title = "Per-islet NOX4 expression by cell type") +
  theme_classic(base_size = 13)

if (nrow(p3_annot) > 0) {
  p3 <- p3 + geom_text(data = p3_annot, aes(x = 1.5, y = Inf, label = label),
                        inherit.aes = FALSE, vjust = 1.3, size = 6)
}

ggsave("figures/Graph3_per_islet_expression_by_celltype.png", p3, width = 9, height = 6, dpi = 300)


comp_stars <- sig_stars(p_comp)

p4_data <- islet_composition %>% filter(!is.na(beta_alpha_ratio))
p4 <- ggplot(p4_data, aes(modality, beta_alpha_ratio, colour = islet_id)) +
  geom_boxplot(aes(colour = NULL), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.9, size = 2.5) +
  scale_y_log10() +
  scale_colour_manual(values = islet_colours, guide = "none") +
  labs(y = "Beta:Alpha ratio (log scale)", x = NULL,
       title = "Islet composition: Beta:Alpha ratio per islet") +
  theme_classic(base_size = 13)

if (!is.na(comp_stars)) {
  p4 <- p4 + annotate("text", x = 1.5, y = max(p4_data$beta_alpha_ratio, na.rm = TRUE), label = comp_stars, size = 6)
}

ggsave("figures/Graph4_BetaAlpha_ratio_by_modality.png", p4, width = 8, height = 6, dpi = 300)


