# =============================================================================
# Program  : 02_create_visualizations.R
# Purpose  : Create two AE visualizations using {ggplot2}:
#            Plot 1 – AE severity distribution by treatment group (bar chart)
#            Plot 2 – Top 10 most frequent AEs with 95% confidence intervals
# Input    : pharmaverseadam::adae, pharmaverseadam::adsl
# Output   : question_3_tlg/plot1_ae_severity.png
#            question_3_tlg/plot2_top10_ae.png
#            question_3_tlg/viz_log.txt
# =============================================================================

library(pharmaverseadam)  # CDISC Pilot ADaM datasets
library(dplyr)            # Data manipulation
library(tidyr)            # Reshaping
library(ggplot2)          # Visualizations
library(forcats)          # Factor reordering
library(scales)           # Axis formatting helpers
library(stringr)          # String wrapping

# ---- 1. Read Input Data & Define Population ----------------------------------
adae_raw <- pharmaverseadam::adae
adsl_raw <- pharmaverseadam::adsl

# Treated subjects: ADSL subjects who received treatment (have a TRTSDT)
adsl_treated <- adsl_raw %>%
  filter(!is.na(TRTSDT)) %>%
  select(USUBJID, ACTARM) %>%
  mutate(
    ACTARM = factor(
      ACTARM,
      levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"),
      labels = c("Placebo\n(N=86)", "Xanomeline\nLow Dose\n(N=96)",
                 "Xanomeline\nHigh Dose\n(N=72)")
    )
  )

# Per-arm N for denominators
arm_N <- c(
  "Placebo\n(N=86)"              = 86L,
  "Xanomeline\nLow Dose\n(N=96)" = 96L,
  "Xanomeline\nHigh Dose\n(N=72)" = 72L
)
total_N <- sum(arm_N)  # 254

# TEAE records: TRTEMFL == "Y".
# ADAE already contains ACTARM; just filter and apply factor ordering.
teae <- adae_raw %>%
  filter(TRTEMFL == "Y",
         ACTARM %in% c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")) %>%
  mutate(
    ACTARM = factor(
      ACTARM,
      levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"),
      labels = c("Placebo\n(N=86)", "Xanomeline\nLow Dose\n(N=96)",
                 "Xanomeline\nHigh Dose\n(N=72)")
    )
  )

cat("TEAE records:", nrow(teae), "\n")
cat("TEAE subjects:", n_distinct(teae$USUBJID), "\n\n")

# ============================================================================
# PLOT 1: AE Severity Distribution by Treatment Group (Stacked Bar Chart)
# ============================================================================
# Each bar = treatment group; segments = MILD / MODERATE / SEVERE
# Y axis = count of unique SUBJECT × AETERM × AESEV combinations
# (subject-level: deduplicate by USUBJID × AETERM × AESEV)
# ============================================================================

# Deduplicate to unique subject × AETERM × AESEV within each arm
sev_data <- teae %>%
  filter(!is.na(AESEV)) %>%
  distinct(USUBJID, ACTARM, AETERM, AESEV) %>%
  mutate(
    AESEV = factor(AESEV, levels = c("MILD", "MODERATE", "SEVERE"))
  ) %>%
  count(ACTARM, AESEV, name = "n_events") %>%
  # Add percentage within each arm (denominator = arm N above)
  mutate(
    N_arm = arm_N[as.character(ACTARM)],
    pct   = round(n_events / N_arm * 100, 1)
  )

cat("Severity distribution data:\n")
print(sev_data)

# Colour palette: clinical traffic-light scheme
sev_colours <- c(
  "MILD"     = "#52B788",   # green
  "MODERATE" = "#F4A261",   # amber
  "SEVERE"   = "#E63946"    # red
)

plot1 <- ggplot(sev_data, aes(x = ACTARM, y = n_events, fill = AESEV)) +
  geom_bar(
    stat     = "identity",
    position = "stack",
    width    = 0.6,
    colour   = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = paste0(n_events, "\n(", pct, "%)")),
    position = position_stack(vjust = 0.5),
    size     = 3.2,
    colour   = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = sev_colours,
    name   = "AE Severity",
    labels = c("Mild", "Moderate", "Severe")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    labels = comma_format()
  ) +
  labs(
    title    = "Distribution of Treatment-Emergent AEs by Severity and Treatment Group",
    subtitle = paste0(
      "Subject-level unique AE-severity pairs; TEAE (TRTEMFL='Y') | ",
      "Total treated subjects (N=", total_N, ")"
    ),
    x        = "Treatment Group",
    y        = "Number of Unique Subject-AE Severity Observations",
    caption  = "Source: pharmaverseadam::adae\nPercentages based on treated subjects per arm"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle    = element_text(size = 10, colour = "grey40", hjust = 0),
    plot.caption     = element_text(size = 9, colour = "grey50"),
    axis.title       = element_text(face = "bold"),
    axis.text.x      = element_text(size = 10, lineheight = 1.2),
    legend.title     = element_text(face = "bold"),
    legend.position  = "right",
    panel.grid.major.y = element_line(colour = "grey90"),
    plot.margin      = margin(12, 12, 12, 12)
  )

ggsave(
  "question_3_tlg/plot1_ae_severity.png",
  plot   = plot1,
  width  = 10,
  height = 7,
  dpi    = 300,
  bg     = "white"
)
cat("Plot 1 saved: question_3_tlg/plot1_ae_severity.png\n")

# ============================================================================
# PLOT 2: Top 10 Most Frequent AEs with 95% Confidence Intervals
# ============================================================================
# Incidence proportion = n_subjects_with_AE / N_total_treated
# 95% CI using Wilson (score) method via prop.test():
#   more accurate than normal approximation for small proportions
# Forest-plot style: horizontal bars with CI error bars
# ============================================================================

# Deduplicate: one record per subject per AETERM (any severity, any arm)
ae_subject <- teae %>%
  distinct(USUBJID, AETERM)

# Count subjects per AETERM (across all treatment arms = overall incidence)
ae_counts <- ae_subject %>%
  count(AETERM, name = "n_subjects") %>%
  arrange(desc(n_subjects)) %>%
  slice_head(n = 10)

cat("\nTop 10 AEs by subject count:\n")
print(ae_counts)

# Calculate incidence proportion and 95% Wilson CI for each top-10 AETERM
calc_wilson_ci <- function(x, n) {
  # Wilson score CI for proportions
  res <- prop.test(x, n, conf.level = 0.95, correct = FALSE)
  data.frame(
    prop  = res$estimate,
    lower = res$conf.int[1],
    upper = res$conf.int[2]
  )
}

top10_ci <- ae_counts %>%
  rowwise() %>%
  mutate(
    ci = list(calc_wilson_ci(n_subjects, total_N))
  ) %>%
  unnest(ci) %>%
  ungroup() %>%
  mutate(
    # Convert to % for display
    pct        = prop * 100,
    lower_pct  = lower * 100,
    upper_pct  = upper * 100,
    # Wrap long AE term names for readability
    AETERM_wrap = str_wrap(AETERM, width = 35),
    # Order factor by frequency (lowest at top for horizontal flip readability)
    AETERM_wrap = fct_reorder(AETERM_wrap, pct)
  )

cat("\nTop 10 AEs with 95% CI (%):\n")
print(top10_ci %>% select(AETERM, n_subjects, pct, lower_pct, upper_pct))

plot2 <- ggplot(
  top10_ci,
  aes(x = pct, y = AETERM_wrap)
) +
  # CI band (ribbon)
  geom_segment(
    aes(x = lower_pct, xend = upper_pct, yend = AETERM_wrap),
    colour    = "#2196F3",
    linewidth = 1.5,
    alpha     = 0.35
  ) +
  # CI whiskers (using geom_errorbar with orientation = "y" for horizontal bars)
  geom_errorbar(
    aes(xmin = lower_pct, xmax = upper_pct),
    width       = 0.35,
    orientation = "y",
    colour      = "#1565C0",
    linewidth   = 0.9
  ) +
  # Point estimate
  geom_point(
    aes(size = n_subjects),
    colour = "#0D47A1",
    alpha  = 0.9
  ) +
  # Label: n (xx%)
  geom_text(
    aes(
      x     = upper_pct + 0.8,
      label = paste0("n=", n_subjects, " (", round(pct, 1), "%)")
    ),
    hjust  = 0,
    size   = 3.4,
    colour = "grey25"
  ) +
  scale_x_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, max(top10_ci$upper_pct) + 12),
    expand = expansion(mult = c(0.01, 0))
  ) +
  scale_size_continuous(
    name   = "Subject\nCount (n)",
    range  = c(3, 8),
    breaks = c(10, 30, 50, 70, 90)
  ) +
  labs(
    title    = "Top 10 Most Frequent Treatment-Emergent Adverse Events",
    subtitle = paste0(
      "Incidence proportion (%) with 95% Wilson Confidence Intervals | ",
      "Overall treated population (N=", total_N, ")"
    ),
    x        = "Incidence (%)",
    y        = NULL,
    caption  = paste0(
      "Source: pharmaverseadam::adae | TEAE (TRTEMFL='Y') | ",
      "One count per subject per AE term | 95% CI by Wilson score method"
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle    = element_text(size = 10, colour = "grey40", hjust = 0),
    plot.caption     = element_text(size = 9, colour = "grey50"),
    axis.title.x     = element_text(face = "bold"),
    axis.text.y      = element_text(size = 10, lineheight = 1.1),
    axis.text.x      = element_text(size = 10),
    legend.title     = element_text(face = "bold", size = 10),
    legend.position  = "right",
    panel.grid.major.x = element_line(colour = "grey90"),
    panel.grid.minor.x = element_line(colour = "grey95"),
    plot.margin      = margin(12, 12, 12, 12)
  )

ggsave(
  "question_3_tlg/plot2_top10_ae.png",
  plot   = plot2,
  width  = 11,
  height = 7,
  dpi    = 300,
  bg     = "white"
)
cat("Plot 2 saved: question_3_tlg/plot2_top10_ae.png\n")

cat("\n=== Both visualizations created successfully ===\n")
