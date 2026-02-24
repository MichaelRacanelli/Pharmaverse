# =============================================================================
# Program  : create_adsl.R
# Purpose  : Create ADaM Subject-Level Analysis Dataset (ADSL) using {admiral}
# Reference: https://pharmaverse.github.io/admiral/cran-release/articles/adsl.html
#            https://pharmaverse.github.io/examples/adam/adsl
# Input    : pharmaversesdtm::dm, ::vs, ::ex, ::ds, ::ae
# Output   : question_2_adam/adsl.csv  (CSV)
#            question_2_adam/adsl.xpt  (SAS transport)
#            question_2_adam/adsl_log.txt (execution log)
# Custom variables derived (per spec):
#   AGEGR9 / AGEGR9N  – Age groups: "<18" (1), "18 - 50" (2), ">50" (3)
#   TRTSDTM / TRTSTMF – First exposure datetime (impute missing time to 00:00:00;
#                        suppress flag when only seconds imputed)
#   ITTFL             – "Y" if DM.ARM populated, else "N"
#   LSTAVLDT          – Last known alive date from VS, AE, DS, TRTEDTM
# =============================================================================

library(admiral)          # ADaM derivation functions
library(pharmaversesdtm)  # CDISC Pilot SDTM source data
library(dplyr)            # Data manipulation
library(lubridate)        # Date/datetime helpers
library(stringr)          # String functions (str_detect)
library(tidyr)            # Reshaping helpers
library(haven)            # Write SAS transport (.xpt)

# ---- 1. Read in SDTM Source Data ---------------------------------------------
# All character blanks in SAS XPT-origin datasets are converted to NA
dm  <- convert_blanks_to_na(pharmaversesdtm::dm)
vs  <- convert_blanks_to_na(pharmaversesdtm::vs)
ex  <- convert_blanks_to_na(pharmaversesdtm::ex)
ds  <- convert_blanks_to_na(pharmaversesdtm::ds)
ae  <- convert_blanks_to_na(pharmaversesdtm::ae)

# ---- 2. Build ADSL Basis from DM ---------------------------------------------
# The DM domain forms the one-record-per-subject skeleton of ADSL.
# DOMAIN is removed as it is not an ADSL variable.
adsl <- dm %>%
  select(-DOMAIN)

# ---- 3. Derive Treatment Variables (TRT01P / TRT01A) -------------------------
# TRT01P = Planned treatment arm; TRT01A = Actual treatment arm
adsl <- adsl %>%
  mutate(
    TRT01P = ARM,
    TRT01A = ACTARM
  )

# ---- 4. Derive Treatment Start / End Datetime (TRTSDTM, TRTEDTM) ------------
# Pre-process EX: convert EXSTDTC / EXENDTC to datetime, imputing missing time.
#
# Per spec for TRTSDTM / TRTSTMF:
#   - Use first exposure record per subject where valid dose and complete date.
#   - Impute completely missing time to 00:00:00 (time_imputation = "first").
#   - Impute partially missing time: 00 for missing hours, 00 for minutes.
#   - If ONLY seconds were missing, do NOT populate TRTSTMF (post-process below).
#
# For TRTEDTM / TRTETMF:
#   - Use last exposure record per subject; impute time to 23:59:59 (per ADSL example).
ex_ext <- ex %>%
  derive_vars_dtm(
    dtc              = EXSTDTC,
    new_vars_prefix  = "EXST",
    time_imputation  = "first",    # impute missing time to 00:00:00
    flag_imputation  = "time"      # generate EXSTTMF flag
  ) %>%
  derive_vars_dtm(
    dtc              = EXENDTC,
    new_vars_prefix  = "EXEN",
    time_imputation  = "last",     # impute missing end time to 23:59:59
    flag_imputation  = "time"      # generate EXENTMF flag
  )

# A valid dose: EXDOSE > 0 OR (EXDOSE == 0 AND EXTRT contains "PLACEBO")
adsl <- adsl %>%
  # Treatment start: first valid dosing record with complete date
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add  = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) &
                  !is.na(EXSTDTM),
    new_vars    = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order       = exprs(EXSTDTM, EXSEQ),
    mode        = "first",
    by_vars     = exprs(STUDYID, USUBJID)
  ) %>%
  # Per spec: suppress imputation flag when only seconds were imputed
  mutate(TRTSTMF = if_else(TRTSTMF == "S", NA_character_, TRTSTMF)) %>%

  # Treatment end: last valid dosing record with complete date
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add  = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) &
                  !is.na(EXENDTM),
    new_vars    = exprs(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
    order       = exprs(EXENDTM, EXSEQ),
    mode        = "last",
    by_vars     = exprs(STUDYID, USUBJID)
  )

# Derive date-only versions of treatment start/end (TRTSDT / TRTEDT)
adsl <- adsl %>%
  derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM))

# Derive treatment duration (TRTDURD = TRTEDT - TRTSDT + 1)
adsl <- adsl %>%
  derive_var_trtdurd()

# ---- 5. Derive Disposition Variables -----------------------------------------

# 5a. Extend DS with a numeric date version of DSSTDTC (no imputation)
ds_ext <- derive_vars_dt(
  ds,
  dtc             = DSSTDTC,
  new_vars_prefix = "DSST"
)

# 5b. End of Study date (EOSDT): last disposition event excluding screen failures
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars     = exprs(STUDYID, USUBJID),
    new_vars    = exprs(EOSDT = DSSTDT),
    filter_add  = DSCAT == "DISPOSITION EVENT" & DSDECOD != "SCREEN FAILURE"
  )

# 5c. End of Study status (EOSSTT)
format_eosstt <- function(x) {
  case_when(
    x %in% "COMPLETED"     ~ "COMPLETED",
    x %in% "SCREEN FAILURE" ~ NA_character_,
    TRUE                   ~ "DISCONTINUED"
  )
}

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add    = ds,
    by_vars        = exprs(STUDYID, USUBJID),
    filter_add     = DSCAT == "DISPOSITION EVENT",
    new_vars       = exprs(EOSSTT = format_eosstt(DSDECOD)),
    missing_values = exprs(EOSSTT = "ONGOING")
  )

# 5d. Discontinuation reason (DCSREAS) — populated only for discontinued subjects
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars     = exprs(STUDYID, USUBJID),
    new_vars    = exprs(DCSREAS = DSDECOD, DCSREASP = DSTERM),
    filter_add  = DSCAT == "DISPOSITION EVENT" &
                  !(DSDECOD %in% c("SCREEN FAILURE", "COMPLETED", NA))
  )

# 5e. Randomization date (RANDDT)
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    filter_add  = DSDECOD == "RANDOMIZED",
    by_vars     = exprs(STUDYID, USUBJID),
    new_vars    = exprs(RANDDT = DSSTDT)
  )

# ---- 6. Derive Birth Date and Analysis Age -----------------------------------
# BRTHDT: numeric birth date from BRTHDTC
adsl <- adsl %>%
  derive_vars_dt(
    new_vars_prefix = "BRTH",
    dtc             = BRTHDTC
  )

# AAGE / AAGEU: analysis age in years relative to randomization date
adsl <- adsl %>%
  derive_vars_aage(
    start_date = BRTHDT,
    end_date   = RANDDT
  )

# ---- 7. Derive Death Variables -----------------------------------------------
# DTHDT: numeric death date from DTHDTC (no imputation)
adsl <- adsl %>%
  derive_vars_dt(
    new_vars_prefix = "DTH",
    dtc             = DTHDTC
  )

# ---- 8. Derive Last Known Alive Date (LSTAVLDT) ------------------------------
# Per spec – max of:
#   (1) Last complete VS date where valid test result exists
#   (2) Last complete AE onset date (AESTDTC)
#   (3) Last complete DS disposition date (DSSTDTC)
#   (4) Last treatment date (datepart of ADSL.TRTEDTM)
adsl <- adsl %>%
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events  = list(
      # (1) Vital Signs: complete VSDTC and at least one valid result
      event(
        dataset_name  = "vs",
        order         = exprs(VSDTC, VSSEQ),
        condition     = !is.na(VSDTC) &
                        !(is.na(VSSTRESN) & is.na(VSSTRESC)),
        set_values_to = exprs(
          LSTAVLDT = convert_dtc_to_dt(VSDTC),
          seq      = VSSEQ
        )
      ),
      # (2) Adverse Events: complete AE onset date
      event(
        dataset_name  = "ae",
        order         = exprs(AESTDTC, AESEQ),
        condition     = !is.na(AESTDTC),
        set_values_to = exprs(
          LSTAVLDT = convert_dtc_to_dt(AESTDTC),
          seq      = AESEQ
        )
      ),
      # (3) Disposition: complete disposition date
      event(
        dataset_name  = "ds",
        order         = exprs(DSSTDTC, DSSEQ),
        condition     = !is.na(DSSTDTC),
        set_values_to = exprs(
          LSTAVLDT = convert_dtc_to_dt(DSSTDTC),
          seq      = DSSEQ
        )
      ),
      # (4) Treatment: datepart of last exposure end date
      event(
        dataset_name  = "adsl",
        condition     = !is.na(TRTEDT),
        set_values_to = exprs(LSTAVLDT = TRTEDT, seq = 0L)
      )
    ),
    source_datasets = list(vs = vs, ae = ae, ds = ds, adsl = adsl),
    tmp_event_nr_var = event_nr,
    order   = exprs(LSTAVLDT, seq, event_nr),
    mode    = "last",
    new_vars = exprs(LSTAVLDT)
  )

# ---- 9. Derive AGEGR9 / AGEGR9N (Custom Age Groups) -------------------------
# Categories per spec: "<18", "18 - 50", ">50"
# AGEGR9N numeric codes:  1,        2,         3
# Based on DM.AGE (age at screening/baseline)
agegr9_lookup <- exprs(
  ~condition,           ~AGEGR9,    ~AGEGR9N,
  AGE < 18,             "<18",       1L,
  between(AGE, 18, 50), "18 - 50",   2L,
  AGE > 50,             ">50",       3L,
  is.na(AGE),           "Missing",   NA_integer_
)

adsl <- adsl %>%
  derive_vars_cat(definition = agegr9_lookup)

# ---- 10. Derive ITTFL (Intent-to-Treat Flag) ----------------------------------
# ITT population: all randomized subjects (ARM is populated in DM)
adsl <- adsl %>%
  mutate(ITTFL = if_else(!is.na(ARM) & ARM != "", "Y", "N"))

# ---- 11. Derive Safety Flag (SAFFL) ------------------------------------------
# Safety population: subjects who received at least one valid dose
adsl <- adsl %>%
  derive_var_merged_exist_flag(
    dataset_add   = ex,
    by_vars       = exprs(STUDYID, USUBJID),
    new_var       = SAFFL,
    false_value   = "N",
    missing_value = "N",
    condition     = EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))
  )

# ---- 12. Final Variable Selection and Ordering --------------------------------
adsl <- adsl %>%
  select(
    # Identifiers
    STUDYID, USUBJID, SUBJID, SITEID, COUNTRY,
    # Demographics
    AGE, AGEU, SEX, RACE, ETHNIC,
    BRTHDT, BRTHDTC,
    # Age groups (custom)
    AGEGR9, AGEGR9N,
    # Analysis age
    AAGE, AAGEU,
    # Treatment assignment
    ARM, ARMCD, ACTARM, ACTARMCD,
    TRT01P, TRT01A,
    # Reference dates from DM
    RFSTDTC, RFENDTC, RFXSTDTC, RFXENDTC, RFICDTC, RFPENDTC,
    # Treatment dates/times
    TRTSDTM, TRTSTMF, TRTSDT,
    TRTEDTM, TRTETMF, TRTEDT,
    TRTDURD,
    # Randomization
    RANDDT,
    # Disposition
    EOSDT, EOSSTT, DCSREAS, DCSREASP,
    # Death
    DTHDTC, DTHDT, DTHFL,
    # Last known alive date (custom)
    LSTAVLDT,
    # Population flags
    ITTFL, SAFFL,
    # DM dates
    DMDTC, DMDY
  )

# ---- 13. Save Output ---------------------------------------------------------
# CSV (human-readable)
write.csv(adsl, "question_2_adam/adsl.csv", row.names = FALSE, na = "")

# SAS transport file (.xpt) — regulatory submission standard
haven::write_xpt(adsl, "question_2_adam/adsl.xpt")

# Summary to console / log
cat("\n=== ADSL dataset created successfully ===\n")
cat("Rows    :", nrow(adsl), "\n")
cat("Columns :", ncol(adsl), "\n\n")

cat("--- Custom variable summary ---\n")
cat("AGEGR9 / AGEGR9N counts:\n")
print(table(adsl$AGEGR9, adsl$AGEGR9N, useNA = "ifany"))

cat("\nTRTSDTM sample (first 5):\n")
print(adsl[1:5, c("USUBJID", "TRTSDTM", "TRTSTMF", "TRTSDT")])

cat("\nITTFL counts:\n")
print(table(adsl$ITTFL, useNA = "ifany"))

cat("\nLSTAVLDT sample (first 5):\n")
print(adsl[1:5, c("USUBJID", "TRTEDTM", "TRTEDT", "LSTAVLDT")])
