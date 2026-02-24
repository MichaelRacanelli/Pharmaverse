# =============================================================================
# Program  : 01_create_ds_domain.R
# Purpose  : Create SDTM Disposition (DS) domain using {sdtm.oak}
# Reference: https://pharmaverse.github.io/examples/sdtm/ae.html
#            https://pharmaverse.github.io/sdtm.oak/
#            SDTMIG v3.4 - DS Domain
# Input    : pharmaverseraw::ds_raw, question_1_sdtm/sdtm_ct.csv,
#            pharmaversesdtm::dm
# Output   : question_1_sdtm/ds.csv  (CSV)
#            question_1_sdtm/ds.xpt  (SAS transport)
#            question_1_sdtm/ds_log.txt (execution log)
# Variables: STUDYID, DOMAIN, USUBJID, DSSEQ, DSTERM, DSDECOD,
#            DSCAT, VISITNUM, VISIT, DSDTC, DSSTDTC, DSSTDY
# =============================================================================

library(sdtm.oak)       # SDTM mapping algorithms
library(pharmaverseraw) # Source raw datasets
library(pharmaversesdtm)# SDTM reference domains (DM)
library(dplyr)          # Data manipulation
library(haven)          # Write SAS transport (.xpt)

# ---- 1. Read in Raw Data -----------------------------------------------------
# ds_raw: raw disposition eCRF data (pharmaverseraw package)
# dm    : Demographics SDTM domain (provides reference start date for study day)
ds_raw <- pharmaverseraw::ds_raw
dm     <- pharmaversesdtm::dm

# ---- 2. Create oak_id_vars ---------------------------------------------------
# oak_id_vars is the crucial link between the raw dataset and the mapped SDTM
# domain. It consists of: oak_id (row identifier), raw_source (source dataset
# name), and patient_number (subject identifier from pat_var).
ds_raw <- ds_raw %>%
  generate_oak_id_vars(
    pat_var = "PATNUM",
    raw_src = "ds_raw"
  )

# ---- 3. Read in Controlled Terminology ---------------------------------------
# study_ct contains codelist mappings: collected_value → term_value (CDISC CT)
study_ct <- read.csv("question_1_sdtm/sdtm_ct.csv")

# ---- 4. Map Topic Variable ---------------------------------------------------
# DSTERM is the topic variable for the DS domain — the verbatim disposition
# term as collected on the eCRF. Mapped from raw column IT.DSTERM using
# assign_no_ct() since no controlled terminology lookup is required here.
ds <- assign_no_ct(
  raw_dat = ds_raw,
  raw_var = "IT.DSTERM",
  tgt_var = "DSTERM",
  id_vars = oak_id_vars()
)

# ---- 5. Map Qualifier, Identifier, and Timing Variables ----------------------
# Starting from the topic-variable dataset (ds), pipe through assign_*
# functions to add each additional SDTM variable. raw_dat is always ds_raw
# (the source); the piped value is the growing target dataset (tgt_dat).
ds <- ds %>%

  # DSDECOD: Standardized disposition decode.
  # Mapped from IT.DSDECOD using C66727 (Disposition Reason codelist).
  # collect_value → term_value lookup (e.g., "Complete" → "COMPLETED").
  assign_ct(
    raw_dat = ds_raw,
    raw_var = "IT.DSDECOD",
    tgt_var = "DSDECOD",
    ct_spec = study_ct,
    ct_clst = "C66727",
    id_vars = oak_id_vars()
  ) %>%

  # VISIT: Visit label as collected on the eCRF (raw column: INSTANCE).
  assign_no_ct(
    raw_dat = ds_raw,
    raw_var = "INSTANCE",
    tgt_var = "VISIT",
    id_vars = oak_id_vars()
  ) %>%

  # DSDTC: Date of the disposition assessment (raw column: DSDTCOL).
  # Format "m-d-y" corresponds to "01-02-2014" = Jan 2, 2014 (ISO: 2014-01-02).
  assign_datetime(
    raw_dat = ds_raw,
    raw_var = "DSDTCOL",
    tgt_var = "DSDTC",
    raw_fmt = c("m-d-y"),
    id_vars = oak_id_vars()
  ) %>%

  # DSSTDTC: Date of study discontinuation / milestone (raw column: IT.DSSTDAT).
  assign_datetime(
    raw_dat = ds_raw,
    raw_var = "IT.DSSTDAT",
    tgt_var = "DSSTDTC",
    raw_fmt = c("m-d-y"),
    id_vars = oak_id_vars()
  )

# ---- 6. Create SDTM Derived Variables ----------------------------------------
ds <- ds %>%
  dplyr::mutate(

    # Identifier variables derived from raw source columns
    STUDYID = ds_raw$STUDY,
    DOMAIN  = "DS",
    USUBJID = paste0("01-", ds_raw$PATNUM),

    # DSCAT: Categorises each record per SDTMIG DS guidance.
    #   PROTOCOL MILESTONE  – events that are part of the planned study flow
    #                         (Randomization, study Completion, Screen Failure).
    #   STUDY DISCONTINUATION – events documenting premature early termination.
    DSCAT = dplyr::case_when(
      toupper(ds_raw$IT.DSDECOD) %in% c("RANDOMIZED", "COMPLETED",
                                         "SCREEN FAILURE")
                          ~ "PROTOCOL MILESTONE",
      !is.na(ds_raw$IT.DSDECOD) ~ "STUDY DISCONTINUATION",
      TRUE ~ NA_character_
    ),

    # VISITNUM: Numeric visit identifier derived from the visit label (INSTANCE).
    # Assigned per the planned visit schedule for CDISCPILOT01.
    VISITNUM = dplyr::case_when(
      ds_raw$INSTANCE == "Screening 1"       ~  1,
      ds_raw$INSTANCE == "Baseline"          ~  2,
      ds_raw$INSTANCE == "Week 2"            ~  3,
      ds_raw$INSTANCE == "Week 4"            ~  4,
      ds_raw$INSTANCE == "Week 6"            ~  5,
      ds_raw$INSTANCE == "Week 8"            ~  6,
      ds_raw$INSTANCE == "Week 12"           ~  7,
      ds_raw$INSTANCE == "Week 16"           ~  8,
      ds_raw$INSTANCE == "Week 20"           ~  9,
      ds_raw$INSTANCE == "Week 24"           ~ 10,
      ds_raw$INSTANCE == "Week 26"           ~ 11,
      ds_raw$INSTANCE == "Retrieval"         ~ 12,
      ds_raw$INSTANCE == "Ambul Ecg Removal" ~ 13,
      TRUE ~ NA_real_
    )

  ) %>%

  # DSSEQ: Sequence number within each subject, ordered by DSSTDTC.
  derive_seq(
    tgt_var  = "DSSEQ",
    rec_vars = c("USUBJID", "DSSTDTC")
  ) %>%

  # DSSTDY: Study day of the disposition event.
  # Calculated as: DSSTDTC - RFXSTDTC (reference first dose date from DM) + 1.
  derive_study_day(
    sdtm_in       = .,
    dm_domain     = dm,
    tgdt          = "DSSTDTC",
    refdt         = "RFXSTDTC",
    study_day_var = "DSSTDY"
  ) %>%

  # Select and order final variables per SDTMIG DS domain specification
  dplyr::select(
    STUDYID, DOMAIN, USUBJID, DSSEQ,
    DSTERM, DSDECOD, DSCAT,
    VISITNUM, VISIT,
    DSDTC, DSSTDTC, DSSTDY
  )

# ---- 7. Save Output ----------------------------------------------------------
# CSV (human-readable)
write.csv(ds, "question_1_sdtm/ds.csv", row.names = FALSE, na = "")

# SAS transport file (.xpt) — standard SDTM submission format
haven::write_xpt(ds, "question_1_sdtm/ds.xpt")

# Summary to console
cat("\n=== DS domain created successfully ===\n")
cat("Rows    :", nrow(ds), "\n")
cat("Columns :", ncol(ds), "\n\n")
print(head(ds, 10))