# =============================================================================
# Program  : 01_create_ae_summary_table.R
# Purpose  : Create a summary table of Treatment-Emergent Adverse Events (TEAEs)
#            using {gtsummary}. Style: FDA Adverse Events Summary (Table 10).
# Reference: https://pharmaverse.github.io/cardinal/quarto/index-catalog.html
#            https://www.danieldsjoberg.com/gtsummary/
# Input    : pharmaverseadam::adae, pharmaverseadam::adsl
# Output   : question_3_tlg/ae_summary_table.html
#            question_3_tlg/ae_summary_log.txt
# Table    : Rows   = AETERM (top AEs by SOC/PT), sorted descending frequency
#            Columns = Placebo | Xanomeline Low Dose | Xanomeline High Dose | Total
#            Cells   = n (%) subjects with ≥1 occurrence
#            TEAE   = TRTEMFL == "Y"
# =============================================================================

library(pharmaverseadam)  # CDISC Pilot ADaM datasets
library(dplyr)            # Data manipulation
library(tidyr)            # Pivoting
library(gtsummary)        # Clinical summary tables
library(gt)               # HTML table styling
library(stringr)          # String helpers

# ---- 1. Read Input Data ------------------------------------------------------
adae_raw <- pharmaverseadam::adae
adsl_raw <- pharmaverseadam::adsl

# ---- 2. Define Analysis Populations ------------------------------------------
# Treated subjects: exclude Screen Failure (no TRTSDT = never dosed)
adsl_treated <- adsl_raw %>%
  filter(!is.na(TRTSDT)) %>%
  select(USUBJID, ACTARM) %>%
  mutate(
    ACTARM = factor(
      ACTARM,
      levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
    )
  )

# Treatment arm N counts (denominators for percentages)
arm_n <- adsl_treated %>%
  count(ACTARM, name = "N_arm") %>%
  mutate(arm_label = paste0(ACTARM, "\n(N=", N_arm, ")"))

total_n <- nrow(adsl_treated)

cat("Treatment arm denominators:\n")
print(arm_n)
cat("Total treated subjects (denominator):", total_n, "\n\n")

# ---- 3. Select TEAE Records --------------------------------------------------
# Treatment-emergent AEs: TRTEMFL == "Y"
# Deduplicate to one record per subject × AETERM (subject-level occurrence flag)
teae <- adae_raw %>%
  filter(TRTEMFL == "Y") %>%
  distinct(USUBJID, AESOC, AETERM)

# ---- 4. Determine Top AEs by Overall Frequency (for table row selection) -----
# Rank AE terms by number of unique subjects experiencing each term
aeterm_rank <- teae %>%
  count(AETERM, name = "n_total") %>%
  arrange(desc(n_total))

cat("Total unique AETERM values in TEAE:", nrow(aeterm_rank), "\n")
cat("Top 10 AEs by subject count:\n")
print(head(aeterm_rank, 10))

# ---- 5. Build Wide Binary Flag Dataset for gtsummary -------------------------
# Approach: For each AETERM in the top AEs, create a binary (Y/N) flag column.
# Start from ADSL_treated so every subject appears (non-event = N = FALSE).
# Show all AEs sorted by descending overall frequency.

# Get ordered AETERM names (all terms, will be ordered by frequency)
ae_terms_ordered <- aeterm_rank$AETERM

# Create subject × AETERM occurrence indicator
teae_flags <- teae %>%
  mutate(occurred = TRUE) %>%
  pivot_wider(
    id_cols     = USUBJID,
    names_from  = AETERM,
    values_from = occurred,
    values_fill = FALSE
  )

# Merge to ADSL (all treated subjects); missing AE flags become FALSE
adsl_ae_wide <- adsl_treated %>%
  left_join(teae_flags, by = "USUBJID") %>%
  mutate(across(
    .cols = all_of(ae_terms_ordered),
    .fns  = ~ if_else(is.na(.), FALSE, .)
  ))

# ---- 6. Build ALSO: Any TEAE flag -------------------------------------------
adsl_ae_wide <- adsl_ae_wide %>%
  mutate(
    ANY_TEAE = rowSums(across(all_of(ae_terms_ordered))) > 0
  )

# ---- 7. Create gtsummary Table -----------------------------------------------
# Columns to include: "Any TEAE" header row + all AE terms ordered by frequency
cols_for_table <- c("ANY_TEAE", ae_terms_ordered)

# Build table: one row per AE flag; by treatment arm; dichotomous (show TRUE)
ae_table <- adsl_ae_wide %>%
  select(ACTARM, all_of(cols_for_table)) %>%
  tbl_summary(
    by          = ACTARM,
    type        = everything() ~ "dichotomous",
    value       = everything() ~ TRUE,
    label       = c(
      ANY_TEAE ~ "Subjects with any TEAE"
    ),
    statistic   = all_dichotomous() ~ "{n} ({p}%)",
    digits      = list(all_dichotomous() ~ c(0, 1)),
    missing     = "no"
  ) %>%
  add_overall(
    last    = TRUE,
    col_label = paste0("**Total**\n**(N=", total_n, ")**")
  ) %>%
  modify_header(
    label           ~ "**Adverse Event Term**",
    all_stat_cols() ~ "**{level}**\n**(N={n})**"
  ) %>%
  modify_caption(
    "**Table 1. Summary of Treatment-Emergent Adverse Events (TEAEs)**\n
    Subjects with ≥1 occurrence; n (%). Sorted by descending overall frequency.\n
    TEAE defined as TRTEMFL='Y' in ADAE."
  ) %>%
  bold_labels() %>%
  italicize_levels()

cat("\nAE Summary Table preview (first 10 rows):\n")
print(ae_table)

# ---- 8. Export to HTML -------------------------------------------------------
ae_table %>%
  as_gt() %>%
  tab_options(
    table.font.size      = px(13),
    data_row.padding     = px(4),
    table.width          = pct(100),
    column_labels.font.weight = "bold",
    heading.title.font.size   = px(15),
    heading.subtitle.font.size = px(12),
    source_notes.font.size     = px(11)
  ) %>%
  tab_source_note(
    source_note = md(
      "*Source: pharmaverseadam::adae, pharmaverseadam::adsl | TEAE: TRTEMFL='Y' |
      Denominator: treated subjects (N=254) | Percentages based on column N*"
    )
  ) %>%
  gtsave("question_3_tlg/ae_summary_table.html")

cat("\n=== AE Summary Table saved: question_3_tlg/ae_summary_table.html ===\n")
