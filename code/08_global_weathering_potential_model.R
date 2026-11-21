
# Global weathering potential model
# Function-only microbial suitability + V18 fixed-domain substrate scaling and AI+PET climate
# Evidence score uses two tiers only: 0.5 or 1.0

# Data root shared with code/00_setup.R. Defined here only when that file
# has not already been sourced, so this script does not load its full package list.
if (!exists("data_dir", inherits = TRUE)) {
  data_dir <- path.expand(Sys.getenv("DESERT_DATA_DIR", unset = "~"))
}

library(tidyverse)
library(readxl)
library(janitor)
library(mgcv)


# 1. Paths


# Curated input directory. Expected files:
#   sand-matrix.xlsx
#   climate_table_updated.xlsx
#   soilgrids_0_5cm_summary.xlsx
# Required for the full model run:
#   microbe_tolerance_table.csv
# Optional for separate 84-day anchor sensitivity:
#   Alxa daily climate POWER csv
input_dir  <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "v16")
parent_input_dir <- dirname(input_dir)
legacy_update_dir <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "file_input", "update")
legacy_file_input_dir <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "file_input")

output_dir <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "v18_model_outputs")
fig_output_dir <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "v18_figure_outputs")
out_plot_dir <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_plot_dir, showWarnings = FALSE, recursive = TRUE)

substrate_file <- file.path(input_dir, "sand-matrix.xlsx")
climate_file   <- file.path(input_dir, "climate_table_updated.xlsx")
soilgrids_file <- file.path(input_dir, "soilgrids_0_5cm_summary.xlsx")

required_input_files <- c(substrate_file, climate_file, soilgrids_file)
missing_input_files <- required_input_files[!file.exists(path.expand(required_input_files))]
if (length(missing_input_files) > 0) {
  stop(
    "Missing required V18 input file(s):\n",
    paste(missing_input_files, collapse = "\n")
  )
}

microbe_candidates <- c(
  
  file.path(input_dir, "microbe_tolerance_table-be.csv")
  
)

microbe_candidates <- microbe_candidates[file.exists(path.expand(microbe_candidates))]

if (length(microbe_candidates) == 0) {
  stop(
    "Cannot find microbe_tolerance_table.csv. Put it in input_dir, parent_input_dir, legacy update dir, or legacy file_input dir."
  )
}

microbe_file <- microbe_candidates[1]

alxa_field_candidates <- c(
  list.files(
    input_dir,
    pattern = "Alxa.*daily.*climate.*POWER.*\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  ),
  list.files(
    parent_input_dir,
    pattern = "Alxa.*daily.*climate.*POWER.*\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  ),
  list.files(
    legacy_update_dir,
    pattern = "Alxa.*daily.*climate.*POWER.*\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  ),
  list.files(
    legacy_file_input_dir,
    pattern = "Alxa.*daily.*climate.*POWER.*\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
)

if (length(alxa_field_candidates) == 0) {
  message(
    "No Alxa 84-day daily climate POWER csv found. ",
    "This is optional in V18 because the main climate anchor uses long-term Tengger climate values."
  )
  alxa_field_file <- NA_character_
} else {
  alxa_field_file <- alxa_field_candidates[1]
}

# 2. Helper functions


clamp01 <- function(x, lower = 0.001, upper = 1) {
  pmin(pmax(as.numeric(x), lower), upper)
}

scale_ai <- function(x) {
  x <- as.numeric(x)
  case_when(
    is.na(x) ~ NA_real_,
    x > 10 ~ x / 10000,
    TRUE ~ x
  )
}

scale01 <- function(x, min_score = 0.05, max_score = 1) {
  x <- as.numeric(x)
  
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  
  rng <- range(x, na.rm = TRUE, finite = TRUE)
  
  if (!is.finite(rng[1]) || !is.finite(rng[2])) {
    return(rep(NA_real_, length(x)))
  }
  
  if (rng[1] == rng[2]) {
    return(rep((min_score + max_score) / 2, length(x)))
  }
  
  out <- (x - rng[1]) / (rng[2] - rng[1])
  out * (max_score - min_score) + min_score
}

scale_to_prediction_domain <- function(
    x,
    site_role,
    full_model_evaluated,
    prediction_role = "global_prediction",
    min_score = 0.05,
    max_score = 1
) {
  idx <- site_role == prediction_role &
    full_model_evaluated &
    !is.na(x)
  
  out <- rep(NA_real_, length(x))
  
  if (sum(idx) < 2) return(out)
  
  rng <- range(x[idx], na.rm = TRUE, finite = TRUE)
  
  if (!is.finite(rng[1]) || !is.finite(rng[2])) return(out)
  
  if (rng[1] == rng[2]) {
    out[!is.na(x)] <- (min_score + max_score) / 2
    return(out)
  }
  
  out[!is.na(x)] <- (x[!is.na(x)] - rng[1]) / (rng[2] - rng[1])
  out[!is.na(x)] <- out[!is.na(x)] * (max_score - min_score) + min_score
  
  # Important for Alxa measured/background and other non-prediction rows:
  # they are projected onto the prediction-domain scale and then bounded.
  out[!is.na(out)] <- pmin(pmax(out[!is.na(out)], min_score), max_score)
  
  out
}

geo_mean <- function(...) {
  x <- cbind(...)
  
  apply(x, 1, function(z) {
    if (any(is.na(z))) return(NA_real_)
    exp(mean(log(pmax(z, 1e-6))))
  })
}

weighted_geo_mean <- function(mat, weights) {
  mat <- as.matrix(mat)
  weights <- weights / sum(weights)
  
  apply(mat, 1, function(z) {
    if (any(is.na(z))) return(NA_real_)
    exp(sum(weights * log(pmax(z, 1e-6))))
  })
}

sample_random_weights <- function(n) {
  w <- rexp(n, rate = 1)
  w / sum(w)
}

pick_col <- function(df, col_x, col_y, col_plain = NULL) {
  if (is.null(col_plain)) {
    col_plain <- gsub("\\.x$", "", col_x)
  }
  
  if (col_x %in% names(df) && col_y %in% names(df)) {
    return(dplyr::coalesce(df[[col_x]], df[[col_y]]))
  }
  
  if (col_x %in% names(df)) return(df[[col_x]])
  if (col_y %in% names(df)) return(df[[col_y]])
  if (col_plain %in% names(df)) return(df[[col_plain]])
  
  rep(NA, nrow(df))
}

get_numeric_col <- function(df, possible_names) {
  nms <- names(df)
  hit <- intersect(possible_names, nms)
  
  if (length(hit) == 0) {
    return(rep(NA_real_, nrow(df)))
  }
  
  as.numeric(df[[hit[1]]])
}

safe_date_parse <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.Date(x))
  
  if (all(is.na(out))) {
    out <- suppressWarnings(as.Date(x, format = "%Y%m%d"))
  }
  
  if (all(is.na(out))) {
    out <- suppressWarnings(as.Date(x, format = "%m/%d/%Y"))
  }
  
  out
}

# Estimate daily reference evapotranspiration when POWER files do not contain PET.
# This is used only for the 84-day Tengger field-climate-anchor sensitivity.
# Hargreaves-Samani ET0 depends on Tmin/Tmax/Tmean, latitude, and day of year.
estimate_hargreaves_pet <- function(date, tmean, tmax, tmin, latitude_deg = 38.85) {
  date <- as.Date(date)
  tmean <- as.numeric(tmean)
  tmax <- as.numeric(tmax)
  tmin <- as.numeric(tmin)
  latitude_rad <- latitude_deg * pi / 180
  j <- as.numeric(format(date, "%j"))
  gsc <- 0.0820
  dr <- 1 + 0.033 * cos(2 * pi * j / 365)
  delta <- 0.409 * sin(2 * pi * j / 365 - 1.39)
  ws <- suppressWarnings(acos(pmax(pmin(-tan(latitude_rad) * tan(delta), 1), -1)))
  ra <- (24 * 60 / pi) * gsc * dr *
    (ws * sin(latitude_rad) * sin(delta) +
       cos(latitude_rad) * cos(delta) * sin(ws))
  td <- pmax(tmax - tmin, 0)
  et0 <- 0.0023 * (tmean + 17.8) * sqrt(td) * ra
  et0[!is.finite(et0) | is.na(date) | is.na(tmean) | is.na(tmax) | is.na(tmin)] <- NA_real_
  pmax(et0, 0)
}

perturb_lognormal <- function(x, log_sd) {
  x <- as.numeric(x)
  out <- x * exp(rnorm(length(x), mean = 0, sd = log_sd))
  out[is.na(x)] <- NA_real_
  out
}

confidence_text_to_log_sd <- function(x) {
  x <- str_to_lower(as.character(x))
  
  case_when(
    str_detect(x, "high") ~ 0.05,
    str_detect(x, "medium|moderate") ~ 0.12,
    str_detect(x, "low|proxy|lower") ~ 0.25,
    TRUE ~ 0.12
  )
}

anchor_lower_ratio <- function(x, anchor_value, min_score = 0.05) {
  # For variables where too-low value is limiting.
  # Example: aridity index.
  x <- as.numeric(x)
  score <- pmin(x / anchor_value, 1)
  score[is.na(score)] <- NA_real_
  pmax(score, min_score)
}

anchor_upper_ratio <- function(x, anchor_value, min_score = 0.05) {
  # For variables where too-high value is limiting.
  # Example: PET or Tmax.
  x <- as.numeric(x)
  score <- pmin(anchor_value / x, 1)
  score[is.na(score)] <- NA_real_
  pmax(score, min_score)
}




# V18 helper functions: fixed-domain scaling, aggregation, grids

MODEL_FLOOR_MAIN <- 0.05

scale_fixed_domain <- function(x, bounds_tbl, predictor, min_score = MODEL_FLOOR_MAIN, max_score = 1) {
  # V18.1 robust fixed-domain scaling.
  # This function intentionally avoids dplyr filtering, because both the bounds
  # table and the calling environment contain objects named `predictor`.
  # It also has a conservative order-based fallback for older/intermediate
  # bounds objects that were created before the predictor labels were finalized.
  x <- as.numeric(x)

  expected_predictors <- c(
    "feldspar_proxy",
    "inverse_quartz_proxy",
    "size_surface_proxy",
    "silt_proxy",
    "cec_proxy"
  )

  predictor_key <- trimws(as.character(predictor)[1])

  if (is.null(bounds_tbl) || !is.data.frame(bounds_tbl) || nrow(bounds_tbl) == 0) {
    stop("Fixed-domain bounds table is empty, NULL, or not a data frame.")
  }

  # Normalize common column names from possible older versions.
  names(bounds_tbl) <- trimws(as.character(names(bounds_tbl)))
  if (!"fixed_min" %in% names(bounds_tbl) && "xmin" %in% names(bounds_tbl)) {
    bounds_tbl$fixed_min <- bounds_tbl$xmin
  }
  if (!"fixed_max" %in% names(bounds_tbl) && "xmax" %in% names(bounds_tbl)) {
    bounds_tbl$fixed_max <- bounds_tbl$xmax
  }
  if (!all(c("fixed_min", "fixed_max") %in% names(bounds_tbl))) {
    stop(
      "Fixed-domain bounds table is missing fixed_min/fixed_max. Current columns: ",
      paste(names(bounds_tbl), collapse = ", ")
    )
  }

  # Try exact predictor-label matching first.
  b <- bounds_tbl[FALSE, , drop = FALSE]
  available_predictors <- character(0)
  if ("predictor" %in% names(bounds_tbl)) {
    available_predictors <- trimws(as.character(bounds_tbl[["predictor"]]))
    b <- bounds_tbl[available_predictors == predictor_key, , drop = FALSE]
  }

  # Try row names as a second option.
  if (nrow(b) != 1) {
    rn <- trimws(as.character(rownames(bounds_tbl)))
    if (length(rn) == nrow(bounds_tbl)) {
      b2 <- bounds_tbl[rn == predictor_key, , drop = FALSE]
      if (nrow(b2) == 1) b <- b2
    }
  }

  # Conservative fallback: if the table has exactly the expected five bounds,
  # use the canonical order. This prevents older bounds objects from crashing
  # when predictor labels are missing or were not saved correctly.
  if (nrow(b) != 1 && nrow(bounds_tbl) >= length(expected_predictors)) {
    pos <- match(predictor_key, expected_predictors)
    if (!is.na(pos) && pos <= nrow(bounds_tbl)) {
      b <- bounds_tbl[pos, , drop = FALSE]
      message(
        "V18 fixed-domain bounds: using canonical-order fallback for ",
        predictor_key,
        ". Please check v18_substrate_fixed_domain_scaling_bounds.csv."
      )
    }
  }

  if (nrow(b) != 1) {
    stop(
      "Cannot find exactly one fixed-domain bound for predictor: ", predictor_key,
      ". Matches found: ", nrow(b),
      ". Available predictors: ",
      ifelse(length(available_predictors) == 0, "<none>", paste(unique(available_predictors), collapse = ", ")),
      ". Bounds columns: ", paste(names(bounds_tbl), collapse = ", ")
    )
  }

  xmin <- suppressWarnings(as.numeric(b[["fixed_min"]][1]))
  xmax <- suppressWarnings(as.numeric(b[["fixed_max"]][1]))
  if (!is.finite(xmin) || !is.finite(xmax) || xmin == xmax) {
    stop(
      "Invalid fixed-domain bounds for predictor: ", predictor_key,
      ". fixed_min=", xmin, ", fixed_max=", xmax
    )
  }

  out <- min_score + (max_score - min_score) * (x - xmin) / (xmax - xmin)
  out[!is.na(out)] <- pmin(pmax(out[!is.na(out)], min_score), max_score)
  out
}
aggregate_scores <- function(mat, weights = NULL, method = "geometric", min_score = MODEL_FLOOR_MAIN) {
  method <- match.arg(method, c("arithmetic", "geometric", "harmonic", "minimum"))
  mat <- as.matrix(mat)
  if (is.null(weights)) weights <- rep(1 / ncol(mat), ncol(mat))
  weights <- as.numeric(weights) / sum(weights)
  apply(mat, 1, function(z) {
    if (any(is.na(z))) return(NA_real_)
    z <- pmax(as.numeric(z), min_score)
    if (method == "arithmetic") return(sum(weights * z))
    if (method == "geometric")  return(exp(sum(weights * log(z))))
    if (method == "harmonic")   return(1 / sum(weights / z))
    if (method == "minimum")    return(min(z))
    NA_real_
  })
}

make_simplex_grid <- function(step = 0.1, include_zero = FALSE) {
  vals <- seq(0, 1, by = step)
  grid <- expand.grid(w1 = vals, w2 = vals) %>%
    mutate(w3 = 1 - w1 - w2) %>%
    filter(w3 >= -1e-9, w3 <= 1 + 1e-9) %>%
    mutate(w3 = pmax(w3, 0))
  if (!include_zero) {
    grid <- grid %>% filter(w1 > 0, w2 > 0, w3 > 0)
  }
  grid %>% arrange(w1, w2, w3)
}

perturb_additive_normal <- function(x, sd) {
  x <- as.numeric(x)
  out <- x + rnorm(length(x), mean = 0, sd = sd)
  out[is.na(x)] <- NA_real_
  out
}

soilgrids_ph_sd_from_evidence <- function(evidence_score) {
  case_when(
    is.na(evidence_score) ~ 0.20,
    evidence_score >= 1 ~ 0.10,
    TRUE ~ 0.25
  )
}

soilgrids_log_sd_from_evidence <- function(evidence_score) {
  case_when(
    is.na(evidence_score) ~ 0.12,
    evidence_score >= 1 ~ 0.05,
    TRUE ~ 0.18
  )
}


# 3. Read V18 input files


# V18 separates sand-matrix, climate, and SoilGrids covariates.
substrate_raw <- read_excel(substrate_file, sheet = 1) %>%
  clean_names() %>%
  mutate(
    site = str_trim(as.character(site))
  )

climate_raw <- read_excel(climate_file, sheet = "Climate_final_14") %>%
  clean_names() %>%
  mutate(
    site = str_trim(as.character(site))
  )

soilgrids_raw <- read_excel(soilgrids_file, sheet = "SoilGrids_wide_model") %>%
  clean_names() %>%
  mutate(
    site = str_trim(as.character(site)),
    # The SoilGrids extraction file used the short/legacy site names for two rows.
    site = case_when(
      site == "Alxa_Left_Banner_project_base" ~ "Alxa_Left_Banner_project_base_measured",
      site == "Argentina_MedanosGra" ~ "Argentina_MedanosGrandes_NE",
      TRUE ~ site
    )
  )

# Join local SoilGrids covariates into the sand-matrix table.
substrate_raw <- substrate_raw %>%
  left_join(soilgrids_raw, by = "site")

missing_soilgrids_sites <- substrate_raw %>%
  filter(is.na(median_silt) | is.na(median_cec) | is.na(median_ph)) %>%
  pull(site)

if (length(missing_soilgrids_sites) > 0) {
  warning(
    "The following V18 sites lack complete SoilGrids covariates after joining: ",
    paste(missing_soilgrids_sites, collapse = ", ")
  )
}

microbe_raw <- read.csv(
  microbe_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  clean_names()

# 4. Optional: read Alxa 84-day field climate for later sensitivity analysis

# V18 main model does NOT use the 84-day field climate as the primary
# climate anchor. All three climate-anchor terms (Tmax, AI, PET) are
# taken from the long-term Tengger/Alxa anchor values in the climate table.
# This optional block only parses the 84-day file if present, so that
# a separate sensitivity analysis can be run downstream if needed.

read_alxa_field_climate <- function(file) {
  
  raw0 <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (nrow(raw0) == 1 && any(str_detect(names(raw0), "\\.\\d{8}$"))) {
    
    out <- raw0 %>%
      pivot_longer(
        cols = everything(),
        names_to = "key",
        values_to = "value"
      ) %>%
      filter(str_detect(key, "\\.\\d{8}$")) %>%
      mutate(
        variable = str_replace(key, "\\.\\d{8}$", ""),
        date = str_extract(key, "\\d{8}$"),
        date = as.Date(date, format = "%Y%m%d"),
        value = as.numeric(value)
      ) %>%
      select(date, variable, value) %>%
      pivot_wider(names_from = variable, values_from = value) %>%
      clean_names()
    
    return(out)
  }
  
  raw <- raw0 %>% clean_names()
  
  if (all(c("year", "mo", "dy") %in% names(raw))) {
    raw <- raw %>%
      mutate(date = as.Date(sprintf("%04d-%02d-%02d", year, mo, dy)))
  } else if ("date" %in% names(raw)) {
    raw <- raw %>%
      mutate(date = safe_date_parse(date))
  } else if ("yyyymmdd" %in% names(raw)) {
    raw <- raw %>%
      mutate(date = as.Date(as.character(yyyymmdd), format = "%Y%m%d"))
  } else {
    stop("Cannot identify date column in Alxa field climate file.")
  }
  
  raw
}

if (!is.na(alxa_field_file) && file.exists(path.expand(alxa_field_file))) {

  alxa_field_daily <- read_alxa_field_climate(alxa_field_file) %>%
    clean_names() %>%
    mutate(
      t2m = get_numeric_col(
        .,
        c("t2m", "temperature", "temp", "mean_temperature")
      ),
      t2m_max = get_numeric_col(
        .,
        c("t2m_max", "tmax", "t_max", "maximum_temperature", "max_temperature")
      ),
      t2m_min = get_numeric_col(
        .,
        c("t2m_min", "tmin", "t_min", "minimum_temperature", "min_temperature")
      ),
      prectotcorr = get_numeric_col(
        .,
        c("prectotcorr", "precipitation", "precip", "ppt", "rain")
      ),
      rh2m = get_numeric_col(
        .,
        c("rh2m", "relative_humidity", "rh", "humidity")
      ),
      allsky_sfc_sw_dwn = get_numeric_col(
        .,
        c("allsky_sfc_sw_dwn", "solar_radiation", "radiation", "shortwave_radiation")
      ),
      pet = get_numeric_col(
        .,
        c(
          "pet", "pet_mm", "pet_daily", "potential_evapotranspiration",
          "potential_evapotranspiration_mm", "evpot", "evapotranspiration"
        )
      ),
      # Most POWER daily files used here do not include PET. In that case,
      # estimate daily ET0 from Tmin/Tmax/Tmean using Hargreaves-Samani for
      # the field-period climate-anchor sensitivity only.
      pet_estimated_hargreaves = estimate_hargreaves_pet(
        date = date,
        tmean = t2m,
        tmax = t2m_max,
        tmin = t2m_min,
        latitude_deg = 38.85
      ),
      pet_for_anchor = dplyr::coalesce(pet, pet_estimated_hargreaves)
    )

  alxa_field_summary <- alxa_field_daily %>%
    summarise(
      n_days = n(),
      field_tmax_mean = mean(t2m_max, na.rm = TRUE),
      field_tmax_q90  = quantile(t2m_max, 0.90, na.rm = TRUE),
      field_tmax_q95  = quantile(t2m_max, 0.95, na.rm = TRUE),
      field_tmax_q98  = quantile(t2m_max, 0.98, na.rm = TRUE),
      field_tmax_max  = max(t2m_max, na.rm = TRUE),
      field_tmin_mean = mean(t2m_min, na.rm = TRUE),
      field_tmin_q05  = quantile(t2m_min, 0.05, na.rm = TRUE),
      field_t2m_mean = mean(t2m, na.rm = TRUE),
      field_precip_total = sum(prectotcorr, na.rm = TRUE),
      field_precip_days_1mm = sum(prectotcorr >= 1, na.rm = TRUE),
      field_precip_max = max(prectotcorr, na.rm = TRUE),
      field_rh_mean = mean(rh2m, na.rm = TRUE),
      field_rh_q25 = quantile(rh2m, 0.25, na.rm = TRUE),
      field_rad_mean = mean(allsky_sfc_sw_dwn, na.rm = TRUE),
      field_pet_total = ifelse(any(is.finite(pet_for_anchor)), sum(pet_for_anchor, na.rm = TRUE), NA_real_),
      field_pet_mean = ifelse(any(is.finite(pet_for_anchor)), mean(pet_for_anchor, na.rm = TRUE), NA_real_),
      field_pet_annualized = ifelse(
        any(is.finite(pet_for_anchor)) & n() > 0,
        sum(pet_for_anchor, na.rm = TRUE) / n() * 365,
        NA_real_
      ),
      field_ai_p_over_pet = ifelse(
        any(is.finite(pet_for_anchor)) & sum(pet_for_anchor, na.rm = TRUE) > 0,
        sum(prectotcorr, na.rm = TRUE) / sum(pet_for_anchor, na.rm = TRUE),
        NA_real_
      ),
      field_pet_source = ifelse(any(is.finite(pet)), "daily PET column", "Hargreaves-Samani ET0 from Tmean/Tmax/Tmin")
    )

  write.csv(
    alxa_field_daily,
    file.path(output_dir, "v18_alxa_84day_field_climate_daily_parsed.csv"),
    row.names = FALSE
  )

  write.csv(
    alxa_field_summary,
    file.path(output_dir, "v18_alxa_84day_field_climate_summary.csv"),
    row.names = FALSE
  )

} else {

  alxa_field_daily <- tibble()
  alxa_field_summary <- tibble()
}


# 5. Prepare V18 model-ready sand-matrix table


# Ensure optional columns used downstream exist.
optional_cols <- c(
  "site",
  "site_label",
  "region",
  "country_or_region",
  "longitude",
  "latitude",
  "coordinate_basis",
  "quartz",
  "feldspar_total",
  "quartz_to_feldspar",
  "representative_mean_um",
  "size_surface_proxy",
  "median_sand",
  "median_silt",
  "median_clay",
  "median_coarse_frag",
  "median_cec",
  "median_ph",
  "median_p_h",
  "mean_ph",
  "mean_p_h",
  "evidence_score",
  "confidence_tier",
  "mineral_source",
  "psd_source",
  "mineral_confidence",
  "psd_confidence",
  "mineral_source_url",
  "psd_source_url",
  "main_source_urls",
  "soilgrids_update_status",
  "update_status",
  "replaced_from",
  "update_note",
  "update_note_v18",
  "model_use",
  "keep_in_final_analysis",
  "use_in_main_model"
)

for (nm in optional_cols) {
  if (!nm %in% names(substrate_raw)) {
    substrate_raw[[nm]] <- NA
  }
}

substrate_raw <- substrate_raw %>%
  mutate(
    site = str_trim(as.character(site)),
    site_original = site,
    mineral_source = as.character(mineral_source),
    psd_source = as.character(psd_source),
    mineral_confidence = as.character(mineral_confidence),
    psd_confidence = as.character(psd_confidence),
    evidence_score = as.numeric(evidence_score),
    # V18 evidence score uses two tiers only.
    evidence_score = case_when(
      is.na(evidence_score) ~ 0.50,
      evidence_score >= 0.75 ~ 1.00,
      TRUE ~ 0.50
    ),
    confidence_tier = case_when(
      evidence_score >= 1 ~ "High confidence",
      TRUE ~ "Lower confidence / proxy"
    ),
    mineral_confidence = case_when(
      !is.na(mineral_confidence) & mineral_confidence != "" & mineral_confidence != "NA" ~ mineral_confidence,
      evidence_score >= 1 ~ "High confidence",
      TRUE ~ "Lower confidence / proxy"
    ),
    psd_confidence = case_when(
      !is.na(psd_confidence) & psd_confidence != "" & psd_confidence != "NA" ~ psd_confidence,
      evidence_score >= 1 ~ "High confidence",
      TRUE ~ "Lower confidence / proxy"
    )
  ) %>%
  filter(!is.na(site), site != "NA") %>%
  filter(!site %in% c("Etosha", "Qaidam", "Mongolian_gobi", "Tengger"))


# 6. Define site roles and final inclusion


substrate_raw <- substrate_raw %>%
  mutate(
    site_role = case_when(
      site == "Alxa_Left_Banner_project_base_measured" ~ "source_anchor_measured",
      site == "Alxa_Left_Banner_project_base_downloaded" ~ "source_background_downloaded",
      TRUE ~ "global_prediction"
    ),
    climate_site = case_when(
      site %in% c(
        "Alxa_Left_Banner_project_base_measured",
        "Alxa_Left_Banner_project_base_downloaded"
      ) ~ "Alxa_Left_Banner_project_base_measured",
      TRUE ~ site
    ),
    include_in_main_model = case_when(
      site == "Alxa_Left_Banner_project_base_downloaded" ~ FALSE,
      !is.na(keep_in_final_analysis) &
        str_detect(str_to_lower(as.character(keep_in_final_analysis)), "^false|^no|0|exclude") ~ FALSE,
      !is.na(use_in_main_model) &
        str_detect(str_to_lower(as.character(use_in_main_model)), "^false|^no|0|exclude") ~ FALSE,
      TRUE ~ TRUE
    )
  )


# 7. Merge climate table


# The V18 climate table was re-extracted for the final 14 sites.
# It is joined by the final site names.
climate_for_join <- climate_raw %>%
  mutate(
    site = as.character(site),
    climate_site = case_when(
      site == "Alxa_Left_Banner_project_base_measured" ~ "Alxa_Left_Banner_project_base_measured",
      TRUE ~ site
    )
  ) %>%
  select(-any_of(c("longitude_x", "latitude_x", "longitude_y", "latitude_y")))

site_master <- substrate_raw %>%
  left_join(climate_for_join, by = "climate_site", suffix = c("", ".climate"))

# Prefer coordinates from the sand-matrix table, but fill from climate table if needed.
site_master$longitude <- as.numeric(pick_col(site_master, "longitude", "longitude.climate", "longitude"))
site_master$latitude  <- as.numeric(pick_col(site_master, "latitude", "latitude.climate", "latitude"))

site_master <- site_master %>%
  select(
    -any_of(c(
      "longitude.climate",
      "latitude.climate",
      "site.climate"
    ))
  )


# 8. Numeric conversion and derived variables


num_cols <- c(
  "quartz",
  "feldspar_total",
  "quartz_to_feldspar",
  "representative_mean_um",
  "median_sand",
  "median_silt",
  "median_clay",
  "median_coarse_frag",
  "median_cec",
  "median_p_h",
  "mean_p_h",
  "median_ph",
  "mean_ph",
  "mat",
  "tmax_warmest_month",
  "tmin_coldest_month",
  "annual_precip",
  "precip_driest_month",
  "precip_seasonality",
  "aridity_index",
  "pet_annual",
  "ai_min",
  "ai_max",
  "ai_mean",
  "ai_sd",
  "ai_cv",
  "pet_min",
  "pet_max",
  "pet_mean",
  "pet_sd",
  "pet_cv",
  "aridity_p_over_pet",
  "pet_for_model"
)

for (nm in num_cols) {
  if (!nm %in% names(site_master)) {
    site_master[[nm]] <- NA_real_
  }
}

site_master <- site_master %>%
  mutate(
    across(all_of(num_cols), as.numeric),
    quartz_to_feldspar = case_when(
      is.na(quartz_to_feldspar) &
        !is.na(quartz) &
        !is.na(feldspar_total) &
        feldspar_total > 0 ~ quartz / feldspar_total,
      TRUE ~ quartz_to_feldspar
    ),
    site_p_h = coalesce(
      median_p_h,
      mean_p_h,
      median_ph,
      mean_ph
    ),
    aridity_index_scaled = scale_ai(aridity_index),
    ai_min_scaled  = scale_ai(ai_min),
    ai_max_scaled  = scale_ai(ai_max),
    ai_mean_scaled = scale_ai(ai_mean),
    aridity_p_over_pet = case_when(
      !is.na(aridity_p_over_pet) ~ scale_ai(aridity_p_over_pet),
      !is.na(aridity_index_scaled) ~ aridity_index_scaled,
      !is.na(ai_mean_scaled) ~ ai_mean_scaled,
      !is.na(annual_precip) & !is.na(pet_annual) & pet_annual > 0 ~ annual_precip / pet_annual,
      TRUE ~ NA_real_
    ),
    pet_for_model = coalesce(pet_for_model, pet_annual, pet_mean),
    size_surface_proxy = 1 / representative_mean_um
  )


# 9. Confidence, domain, completeness


site_master <- site_master %>%
  mutate(
    mineral_complete = !is.na(quartz) &
      !is.na(feldspar_total) &
      !is.na(quartz_to_feldspar),
    # V18 completeness logic:
    # sand and clay are retained as metadata but are no longer required by the model.
    # The substrate module uses mineralogy, representative particle size, silt, and CEC.
    size_complete = !is.na(representative_mean_um),
    surface_exchange_complete = !is.na(median_silt) &
      !is.na(median_cec),
    texture_complete = surface_exchange_complete,
    psd_complete = size_complete,
    substrate_evaluated = mineral_complete & size_complete & surface_exchange_complete,
    substrate_domain = case_when(
      site_role != "global_prediction" ~ "source_or_reference",
      site %in% c("Badain_Jaran", "Taklamakan") ~ "gobi_transition",
      site %in% c(
        "Sonoran_Altar",
        "Argentina_MedanosGrandes_NE",
        "Mu_Us",
        "Hobq",
        "Mojave",
        "Lower_Draa",
        "Namib",
        "SE_Qatar",
        "Kalahari_Bagatelle",
        "Karakum_AshkAbad",
        "Thal_Mankera"
      ) ~ "primary_or_arid_sandy_desert",
      TRUE ~ "other_or_unclassified"
    ),
    # Evidence score is read from the curated V18 input table and restricted
    # to two tiers only: 1.0 or 0.5.
    evidence_score = case_when(
      site == "Alxa_Left_Banner_project_base_measured" ~ 1.00,
      is.na(evidence_score) ~ 0.50,
      evidence_score >= 0.75 ~ 1.00,
      TRUE ~ 0.50
    ),
    confidence_tier = case_when(
      evidence_score >= 1.00 ~ "High confidence",
      TRUE ~ "Lower confidence / proxy"
    ),
    data_confidence = case_when(
      evidence_score >= 1.00 ~ "high",
      TRUE ~ "lower_proxy"
    ),
    plot_domain = case_when(
      site %in% c(
        "Alxa_Left_Banner_project_base_measured",
        "Alxa_Left_Banner_project_base_downloaded"
      ) ~ "Source/reference",
      substrate_domain == "gobi_transition" ~ "Gobi transition",
      substrate_domain == "primary_or_arid_sandy_desert" ~ "Sandy desert",
      TRUE ~ "Other arid domain"
    ),
    full_model_evaluated =
      substrate_evaluated &
      !is.na(site_p_h) &
      !is.na(mat) &
      !is.na(tmax_warmest_month) &
      !is.na(aridity_p_over_pet) &
      !is.na(pet_for_model)
  )

# V18 diagnostic: Alxa measured must have representative_mean_um, because it
# remains the denominator for relative_to_alxa_measured.
cat("\nCheck Alxa representative_mean_um after V18 input read:\n")
print(
  site_master %>%
    filter(site %in% c(
      "Alxa_Left_Banner_project_base_measured",
      "Alxa_Left_Banner_project_base_downloaded"
    )) %>%
    select(
      site,
      representative_mean_um,
      size_complete,
      surface_exchange_complete,
      full_model_evaluated,
      evidence_score
    )
)


# 10. Microbial response curves: pH + temperature only

# V18 uses only CFU-normalized CAS activity as the microbial-function
# suitability term. Growth and population-level function are not used in
# the global projection, to avoid mixing cell abundance with per-viable-cell
# siderophore production.

required_microbe_cols <- c(
  "condition_type",
  "condition_value",
  "function_score"
)

missing_microbe_cols <- setdiff(required_microbe_cols, names(microbe_raw))

if (length(missing_microbe_cols) > 0) {
  stop(
    "Missing microbe tolerance columns: ",
    paste(missing_microbe_cols, collapse = ", ")
  )
}

microbe_df <- microbe_raw %>%
  mutate(
    condition_type = as.character(condition_type),
    condition_value = as.numeric(condition_value),
    function_score = as.numeric(function_score),
    condition_type_std = case_when(
      str_detect(str_to_lower(condition_type), "ph") ~ "pH",
      str_detect(str_to_lower(condition_type), "temp") ~ "temperature",
      str_detect(str_to_lower(condition_type), "nacl|salt|salinity") ~ "NaCl",
      TRUE ~ condition_type
    )
  )

microbe_df_v18 <- microbe_df %>%
  filter(condition_type_std %in% c("pH", "temperature"))

# Exclude extrapolative edge conditions used only as tolerance checks.
microbe_function_df_v18 <- microbe_df_v18 %>%
  filter(
    !(condition_type_std == "pH" & condition_value >= 10),
    !(condition_type_std == "temperature" & condition_value >= 45)
  )

fit_response_curve <- function(dat, yvar) {

  dat <- dat %>%
    filter(!is.na(condition_value), !is.na(.data[[yvar]]))

  n_unique <- length(unique(dat$condition_value))

  if (nrow(dat) < 4 || n_unique < 4) {
    model <- loess(
      as.formula(paste(yvar, "~ condition_value")),
      data = dat,
      span = 0.8,
      degree = 2,
      control = loess.control(surface = "direct")
    )
    model_type <- "loess"
  } else {
    model <- gam(
      as.formula(paste(yvar, "~ s(condition_value, k = min(5, n_unique))")),
      data = dat,
      method = "REML"
    )
    model_type <- "gam"
  }

  list(
    model = model,
    model_type = model_type,
    yvar = yvar,
    x_min = min(dat$condition_value, na.rm = TRUE),
    x_max = max(dat$condition_value, na.rm = TRUE)
  )
}

predict_response <- function(model_obj, new_x, clamp = TRUE) {

  original_x <- as.numeric(new_x)
  pred_x <- original_x

  if (clamp) {
    pred_x <- pmin(pmax(pred_x, model_obj$x_min), model_obj$x_max)
  }

  new_df <- data.frame(condition_value = pred_x)
  pred <- predict(model_obj$model, newdata = new_df)
  pred <- as.numeric(pred)
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1

  pred
}

response_models_v18 <- list()

for (cond in c("pH", "temperature")) {

  dat_func <- microbe_function_df_v18 %>%
    filter(condition_type_std == cond)

  if (nrow(dat_func) == 0) {
    stop("No CFU-normalized function data found for condition: ", cond)
  }

  response_models_v18[[cond]] <- list(
    function_score = fit_response_curve(dat_func, "function_score")
  )
}

response_models <- response_models_v18
microbe_function_df <- microbe_function_df_v18


# 11. Anchor values for climate module only; V18 does not use Alxa size anchor

# V18 update:
# All climate-anchor terms are taken from the long-term Tengger/Alxa
# anchor row in the climate table. The main model no longer uses the
# 84-day field-period Tmax as the heat anchor.

alxa_tmax_ref <- site_master %>%
  filter(site == "Alxa_Left_Banner_project_base_measured") %>%
  pull(tmax_warmest_month) %>%
  .[1] %>%
  as.numeric()

alxa_ai_ref <- site_master %>%
  filter(site == "Alxa_Left_Banner_project_base_measured") %>%
  pull(aridity_p_over_pet) %>%
  .[1] %>%
  as.numeric()

alxa_pet_ref <- site_master %>%
  filter(site == "Alxa_Left_Banner_project_base_measured") %>%
  pull(pet_for_model) %>%
  .[1] %>%
  as.numeric()

if (is.na(alxa_tmax_ref) || alxa_tmax_ref <= 0) {
  stop("Alxa long-term Tmax reference is missing or invalid.")
}

if (is.na(alxa_ai_ref) || alxa_ai_ref <= 0) {
  stop("Alxa long-term aridity reference is missing or invalid.")
}

if (is.na(alxa_pet_ref) || alxa_pet_ref <= 0) {
  stop("Alxa long-term PET reference is missing or invalid.")
}






# 12. V18 fixed-domain scaling bounds for sand-matrix predictors

# V18 keeps min-max scaling for sand-matrix variables, but the min/max
# values are calculated once from the final candidate-site dataset plus the
# Tengger experimental anchor and then held fixed across the main model,
# Monte Carlo runs, and all sensitivity analyses.

make_substrate_fixed_bounds_v18 <- function(df, buffer_fraction = 0) {
  bound_df <- df %>%
    filter(
      include_in_main_model,
      full_model_evaluated,
      site != "Alxa_Left_Banner_project_base_downloaded"
    ) %>%
    mutate(
      feldspar_proxy = feldspar_total,
      inverse_quartz_proxy = 1 / quartz_to_feldspar,
      size_surface_proxy = 1 / representative_mean_um,
      silt_proxy = median_silt,
      cec_proxy = median_cec
    )
  
  predictors <- c(
    "feldspar_proxy",
    "inverse_quartz_proxy",
    "size_surface_proxy",
    "silt_proxy",
    "cec_proxy"
  )
  
  map_dfr(predictors, function(p) {
    rng <- range(bound_df[[p]], na.rm = TRUE, finite = TRUE)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) {
      stop("Invalid fixed-domain range for ", p)
    }
    width <- rng[2] - rng[1]
    tibble(
      predictor = p,
      fixed_min = rng[1] - buffer_fraction * width,
      fixed_max = rng[2] + buffer_fraction * width,
      empirical_min = rng[1],
      empirical_max = rng[2],
      buffer_fraction = buffer_fraction,
      n_sites_for_bounds = nrow(bound_df)
    )
  })
}

validate_substrate_fixed_bounds_v18 <- function(bounds_tbl) {
  required_predictors <- c(
    "feldspar_proxy",
    "inverse_quartz_proxy",
    "size_surface_proxy",
    "silt_proxy",
    "cec_proxy"
  )
  if (is.null(bounds_tbl) || nrow(bounds_tbl) == 0) {
    stop("V18 fixed-domain bounds table was not created.")
  }
  if (!all(c("predictor", "fixed_min", "fixed_max") %in% names(bounds_tbl))) {
    stop("V18 fixed-domain bounds table lacks required columns.")
  }
  bounds_tbl <- bounds_tbl %>%
    mutate(
      predictor = trimws(as.character(predictor)),
      fixed_min = as.numeric(fixed_min),
      fixed_max = as.numeric(fixed_max)
    )
  missing_predictors <- setdiff(required_predictors, bounds_tbl$predictor)
  duplicated_predictors <- bounds_tbl$predictor[duplicated(bounds_tbl$predictor)]
  invalid_predictors <- bounds_tbl %>%
    filter(!is.finite(fixed_min) | !is.finite(fixed_max) | fixed_min == fixed_max) %>%
    pull(predictor)
  if (length(missing_predictors) > 0 || length(duplicated_predictors) > 0 || length(invalid_predictors) > 0) {
    stop(
      "Invalid V18 fixed-domain bounds. Missing: ", paste(missing_predictors, collapse = ", "),
      "; duplicated: ", paste(unique(duplicated_predictors), collapse = ", "),
      "; invalid ranges: ", paste(unique(invalid_predictors), collapse = ", ")
    )
  }
  bounds_tbl
}

substrate_fixed_bounds_v18 <- validate_substrate_fixed_bounds_v18(
  make_substrate_fixed_bounds_v18(
    site_master,
    buffer_fraction = 0
  )
)

write.csv(
  substrate_fixed_bounds_v18,
  file.path(output_dir, "v18_substrate_fixed_domain_scaling_bounds.csv"),
  row.names = FALSE
)

cat("\nV18 fixed-domain substrate bounds created with predictors:\n")
print(substrate_fixed_bounds_v18)


# 13. Core model computation function


compute_model_once <- function(
    df,
    substrate_axis_weights = rep(1 / 3, 3),
    scm_weights = rep(1 / 3, 3),
    fixed_bounds = substrate_fixed_bounds_v18,
    min_score = MODEL_FLOOR_MAIN,
    aggregation_method = "geometric",
    climate_mode = c("ai_pet", "ai_only", "ai_pet_min", "ai_pet_tmax", "conservative_min"),
    alxa_tmax_anchor = alxa_tmax_ref,
    alxa_ai_anchor = alxa_ai_ref,
    alxa_pet_anchor = alxa_pet_ref,
    perturb_raw_substrate = FALSE,
    perturb_soilgrids = FALSE
) {
  
  aggregation_method <- match.arg(
    aggregation_method,
    c("arithmetic", "geometric", "harmonic", "minimum")
  )
  climate_mode <- match.arg(climate_mode)
  
  df2 <- df
  
  if (perturb_raw_substrate || perturb_soilgrids) {
    
    mineral_sd <- confidence_text_to_log_sd(df2$mineral_confidence)
    psd_sd <- confidence_text_to_log_sd(df2$psd_confidence)
    soil_log_sd <- soilgrids_log_sd_from_evidence(df2$evidence_score)
    soil_ph_sd <- soilgrids_ph_sd_from_evidence(df2$evidence_score)
    
    # The Tengger experimental anchor is treated as the fixed denominator and
    # is not perturbed in the main Monte Carlo projection.
    anchor_idx <- df2$site == "Alxa_Left_Banner_project_base_measured"
    mineral_sd[anchor_idx] <- 0
    psd_sd[anchor_idx] <- 0
    soil_log_sd[anchor_idx] <- 0
    soil_ph_sd[anchor_idx] <- 0
    
    if (perturb_raw_substrate) {
      df2 <- df2 %>%
        mutate(
          quartz = perturb_lognormal(quartz, mineral_sd),
          feldspar_total = perturb_lognormal(feldspar_total, mineral_sd),
          representative_mean_um = perturb_lognormal(representative_mean_um, psd_sd),
          median_sand = perturb_lognormal(median_sand, psd_sd),
          median_clay = perturb_lognormal(median_clay, psd_sd)
        )
    }
    
    if (perturb_soilgrids) {
      df2 <- df2 %>%
        mutate(
          median_silt = perturb_lognormal(median_silt, soil_log_sd),
          median_cec = perturb_lognormal(median_cec, soil_log_sd),
          site_p_h = perturb_additive_normal(site_p_h, soil_ph_sd),
          site_p_h = pmin(pmax(site_p_h, 3.5), 10.5)
        )
    }
    
    df2 <- df2 %>%
      mutate(
        texture_sum = median_sand + median_silt + median_clay,
        median_sand = ifelse(!is.na(texture_sum) & texture_sum > 0, median_sand / texture_sum * 100, median_sand),
        median_silt = ifelse(!is.na(texture_sum) & texture_sum > 0, median_silt / texture_sum * 100, median_silt),
        median_clay = ifelse(!is.na(texture_sum) & texture_sum > 0, median_clay / texture_sum * 100, median_clay),
        quartz_to_feldspar = case_when(
          !is.na(quartz) & !is.na(feldspar_total) & feldspar_total > 0 ~
            quartz / feldspar_total,
          TRUE ~ quartz_to_feldspar
        )
      )
  }
  
  df2 <- df2 %>%
    mutate(
      # 
      # V18 sand-matrix transformations
      # Each predictor is first direction-adjusted so that larger values mean
      # more favorable substrate conditions. Direction-adjusted predictors are
      # then scaled by fixed-domain min-max bounds, not by dynamic site-set
      # min/max values.
      # 
      feldspar_proxy = feldspar_total,
      inverse_quartz_proxy = 1 / quartz_to_feldspar,
      size_surface_proxy = 1 / representative_mean_um,
      silt_proxy = median_silt,
      cec_proxy = median_cec,
      
      feldspar_score = scale_fixed_domain(
        feldspar_proxy,
        fixed_bounds,
        "feldspar_proxy",
        min_score = min_score
      ),
      inverse_quartz_score = scale_fixed_domain(
        inverse_quartz_proxy,
        fixed_bounds,
        "inverse_quartz_proxy",
        min_score = min_score
      ),
      size_surface_score = scale_fixed_domain(
        size_surface_proxy,
        fixed_bounds,
        "size_surface_proxy",
        min_score = min_score
      ),
      silt_score = scale_fixed_domain(
        silt_proxy,
        fixed_bounds,
        "silt_proxy",
        min_score = min_score
      ),
      cec_score = scale_fixed_domain(
        cec_proxy,
        fixed_bounds,
        "cec_proxy",
        min_score = min_score
      )
    )
  
  mineral_reactivity_axis <- aggregate_scores(
    cbind(df2$feldspar_score, df2$inverse_quartz_score),
    weights = rep(1 / 2, 2),
    method = aggregation_method,
    min_score = min_score
  )
  
  physical_accessibility_axis <- df2$size_surface_score
  
  surface_exchange_axis <- aggregate_scores(
    cbind(df2$silt_score, df2$cec_score),
    weights = rep(1 / 2, 2),
    method = aggregation_method,
    min_score = min_score
  )
  
  substrate_axis <- aggregate_scores(
    cbind(
      mineral_reactivity_axis,
      physical_accessibility_axis,
      surface_exchange_axis
    ),
    weights = substrate_axis_weights,
    method = aggregation_method,
    min_score = min_score
  )
  
  substrate_fe_axis <- aggregate_scores(
    cbind(mineral_reactivity_axis, surface_exchange_axis),
    weights = rep(1 / 2, 2),
    method = aggregation_method,
    min_score = min_score
  )
  
  substrate_texture_axis <- aggregate_scores(
    cbind(physical_accessibility_axis, surface_exchange_axis),
    weights = rep(1 / 2, 2),
    method = aggregation_method,
    min_score = min_score
  )
  
  # 
  # V18 climate module
  # Main model uses AI + PET only. Tmax is still calculated and retained as a
  # diagnostic/sensitivity covariate, but it is not included in the main climate
  # score unless climate_mode = "ai_pet_tmax" or "conservative_min".
  # 
  tmax_suitability <- anchor_upper_ratio(
    df2$tmax_warmest_month,
    alxa_tmax_anchor,
    min_score = min_score
  )
  
  aridity_suitability <- anchor_lower_ratio(
    df2$aridity_p_over_pet,
    alxa_ai_anchor,
    min_score = min_score
  )
  
  pet_suitability <- anchor_upper_ratio(
    df2$pet_for_model,
    alxa_pet_anchor,
    min_score = min_score
  )
  
  if (climate_mode == "ai_only") {
    climate_suitability <- aridity_suitability
  } else if (climate_mode == "ai_pet") {
    climate_suitability <- aggregate_scores(
      cbind(aridity_suitability, pet_suitability),
      weights = rep(1 / 2, 2),
      method = aggregation_method,
      min_score = min_score
    )
  } else if (climate_mode == "ai_pet_min") {
    # V18 strict climate constraint without reintroducing Tmax:
    # the climate score is determined by the more limiting of AI and PET.
    climate_suitability <- pmin(aridity_suitability, pet_suitability)
  } else if (climate_mode == "ai_pet_tmax") {
    climate_suitability <- aggregate_scores(
      cbind(tmax_suitability, aridity_suitability, pet_suitability),
      weights = rep(1 / 3, 3),
      method = aggregation_method,
      min_score = min_score
    )
  } else if (climate_mode == "conservative_min") {
    climate_suitability <- pmin(tmax_suitability, aridity_suitability, pet_suitability)
  } else {
    stop("Unknown climate_mode: ", climate_mode)
  }
  
  # 
  # Microbial module from experimental response curves
  # V18 perturbs SoilGrids-derived pH as an input when requested, but the
  # coordinate-derived temperature input is held fixed in MC.
  # 
  microbe_ph_function <- predict_response(
    response_models_v18[["pH"]][["function_score"]],
    df2$site_p_h,
    clamp = TRUE
  )
  microbe_ph_function <- pmax(microbe_ph_function, min_score)
  
  microbe_mat_function <- predict_response(
    response_models_v18[["temperature"]][["function_score"]],
    df2$mat,
    clamp = TRUE
  )
  microbe_mat_function <- pmax(microbe_mat_function, min_score)
  
  microbial_suitability <- aggregate_scores(
    cbind(microbe_ph_function, microbe_mat_function),
    weights = rep(1 / 2, 2),
    method = aggregation_method,
    min_score = min_score
  )
  
  microbe_function_suitability <- microbial_suitability
  
  weathering_potential_raw <- aggregate_scores(
    cbind(substrate_axis, climate_suitability, microbial_suitability),
    weights = scm_weights,
    method = aggregation_method,
    min_score = min_score
  )
  
  fe_mobilization_potential_raw <- aggregate_scores(
    cbind(substrate_fe_axis, climate_suitability, microbial_suitability),
    weights = scm_weights,
    method = aggregation_method,
    min_score = min_score
  )
  
  texture_restructuring_potential_raw <- aggregate_scores(
    cbind(substrate_texture_axis, climate_suitability, microbial_suitability),
    weights = scm_weights,
    method = aggregation_method,
    min_score = min_score
  )
  
  alxa_idx <- which(df2$site == "Alxa_Left_Banner_project_base_measured")[1]
  
  df2 %>%
    mutate(
      feldspar_score = feldspar_score,
      inverse_quartz_score = inverse_quartz_score,
      size_surface_proxy = size_surface_proxy,
      size_surface_score = size_surface_score,
      silt_score = silt_score,
      cec_score = cec_score,
      
      mineral_reactivity_axis = mineral_reactivity_axis,
      physical_accessibility_axis = physical_accessibility_axis,
      surface_exchange_axis = surface_exchange_axis,
      substrate_axis = substrate_axis,
      substrate_fe_axis = substrate_fe_axis,
      substrate_texture_axis = substrate_texture_axis,
      
      tmax_suitability = tmax_suitability,
      aridity_suitability = aridity_suitability,
      pet_suitability = pet_suitability,
      climate_suitability = climate_suitability,
      climate_mode = climate_mode,
      
      microbe_ph_function = microbe_ph_function,
      microbe_mat_function = microbe_mat_function,
      microbe_function_suitability = microbe_function_suitability,
      microbial_suitability = microbial_suitability,
      
      weathering_potential_raw = weathering_potential_raw,
      fe_mobilization_potential_raw = fe_mobilization_potential_raw,
      texture_restructuring_potential_raw = texture_restructuring_potential_raw,
      
      relative_to_alxa_measured =
        weathering_potential_raw / weathering_potential_raw[alxa_idx],
      
      relative_fe_to_alxa_measured =
        fe_mobilization_potential_raw / fe_mobilization_potential_raw[alxa_idx],
      
      relative_texture_to_alxa_measured =
        texture_restructuring_potential_raw / texture_restructuring_potential_raw[alxa_idx],
      
      site_label_clean = case_when(
        site == "Alxa_Left_Banner_project_base_measured" ~ "Tengger experimental anchor",
        site == "Alxa_Left_Banner_project_base_downloaded" ~ "Alxa background",
        !is.na(site_label) & site_label != "" ~ as.character(site_label),
        TRUE ~ str_replace_all(site, "_", " ")
      ),
      
      use_in_main_map =
        include_in_main_model &
        !is.na(relative_to_alxa_measured) &
        !is.na(longitude) &
        !is.na(latitude) &
        site != "Alxa_Left_Banner_project_base_downloaded",
      
      use_in_main_barplot =
        include_in_main_model &
        site != "Alxa_Left_Banner_project_base_measured" &
        site != "Alxa_Left_Banner_project_base_downloaded" &
        !is.na(relative_to_alxa_measured) &
        full_model_evaluated,
      
      use_in_main_scatter =
        include_in_main_model &
        !is.na(relative_fe_to_alxa_measured) &
        !is.na(relative_texture_to_alxa_measured) &
        site != "Alxa_Left_Banner_project_base_downloaded"
    )
}

# 
# 13. Deterministic V18
# 

site_master_det <- compute_model_once(
  site_master,
  substrate_axis_weights = rep(1 / 3, 3),
  scm_weights = rep(1 / 3, 3),
  fixed_bounds = substrate_fixed_bounds_v18,
  min_score = MODEL_FLOOR_MAIN,
  aggregation_method = "geometric",
  climate_mode = "ai_pet",
  perturb_raw_substrate = FALSE,
  perturb_soilgrids = FALSE
)

write.csv(
  site_master_det,
  file.path(output_dir, "v18_deterministic_site_master_full.csv"),
  row.names = FALSE
)

cat("\nCheck denominator Alxa weathering potential:\n")
print(
  site_master_det %>%
    filter(site == "Alxa_Left_Banner_project_base_measured") %>%
    select(
      site,
      representative_mean_um,
      size_surface_proxy,
      size_surface_score,
      substrate_axis,
      climate_suitability,
      microbial_suitability,
      weathering_potential_raw,
      relative_to_alxa_measured
    )
)

deterministic_ranking_v18 <- site_master_det %>%
  filter(use_in_main_barplot) %>%
  arrange(desc(relative_to_alxa_measured)) %>%
  select(
    site,
    site_label_clean,
    relative_to_alxa_measured,
    substrate_axis,
    mineral_reactivity_axis,
    physical_accessibility_axis,
    surface_exchange_axis,
    climate_suitability,
    microbial_suitability,
    relative_fe_to_alxa_measured,
    relative_texture_to_alxa_measured,
    confidence_tier,
    evidence_score,
    include_in_main_model,
    substrate_domain,
    plot_domain
  )

write.csv(
  deterministic_ranking_v18,
  file.path(output_dir, "v18_deterministic_main_ranking_with_modules.csv"),
  row.names = FALSE
)

component_scores_v18 <- site_master_det %>%
  select(
    site,
    site_label_clean,
    region,
    longitude,
    latitude,
    include_in_main_model,
    use_in_main_barplot,
    substrate_domain,
    plot_domain,
    confidence_tier,
    evidence_score,
    
    relative_to_alxa_measured,
    substrate_axis,
    mineral_reactivity_axis,
    physical_accessibility_axis,
    surface_exchange_axis,
    climate_suitability,
    microbial_suitability,
    relative_fe_to_alxa_measured,
    relative_texture_to_alxa_measured,
    
    quartz,
    feldspar_total,
    quartz_to_feldspar,
    representative_mean_um,
    median_sand,
    median_silt,
    median_clay,
    median_cec,
    site_p_h,
    
    feldspar_score,
    inverse_quartz_score,
    size_surface_proxy,
    size_surface_score,
    silt_score,
    cec_score,
    
    mat,
    tmax_warmest_month,
    aridity_index,
    aridity_index_scaled,
    aridity_p_over_pet,
    pet_for_model,
    
    tmax_suitability,
    aridity_suitability,
    pet_suitability,
    
    microbe_ph_function,
    microbe_mat_function,
    microbe_function_suitability,
    
    mineral_source,
    psd_source,
    mineral_confidence,
    psd_confidence,
    psd_source_url
   
  ) %>%
  arrange(desc(relative_to_alxa_measured))

write.csv(
  component_scores_v18,
  file.path(output_dir, "v18_deterministic_component_scores_wide.csv"),
  row.names = FALSE
)

psd_values_used_v18 <- site_master_det %>%
  select(
    site,
    site_label_clean,
    representative_mean_um,
    median_sand,
    median_silt,
    median_clay,
    size_surface_proxy,
    size_surface_score,
    silt_score,
    psd_confidence,
    psd_source,
    psd_source_url
    
  ) %>%
  arrange(site)

write.csv(
  psd_values_used_v18,
  file.path(output_dir, "v18_psd_values_used_in_model.csv"),
  row.names = FALSE
)


# 14. Monte Carlo uncertainty model: input-data uncertainty only


set.seed(202606)

n_iter <- 5000

run_one_mc <- function(iter_id) {
  
  # V18 MC propagates input-data uncertainty only.
  # Substrate sub-axis weights and SCM weights remain fixed in the main MC
  # and are evaluated separately by deterministic simplex-grid sensitivity.
  alxa_tmax_anchor_i <- alxa_tmax_ref
  alxa_ai_anchor_i <- alxa_ai_ref
  alxa_pet_anchor_i <- alxa_pet_ref
  
  model_i <- compute_model_once(
    site_master,
    substrate_axis_weights = rep(1 / 3, 3),
    scm_weights = rep(1 / 3, 3),
    fixed_bounds = substrate_fixed_bounds_v18,
    min_score = MODEL_FLOOR_MAIN,
    aggregation_method = "geometric",
    climate_mode = "ai_pet",
    alxa_tmax_anchor = alxa_tmax_anchor_i,
    alxa_ai_anchor = alxa_ai_anchor_i,
    alxa_pet_anchor = alxa_pet_anchor_i,
    perturb_raw_substrate = TRUE,
    perturb_soilgrids = TRUE
  )
  
  main_idx <- model_i$include_in_main_model &
    model_i$site != "Alxa_Left_Banner_project_base_measured" &
    model_i$site != "Alxa_Left_Banner_project_base_downloaded" &
    model_i$full_model_evaluated
  
  rank_vec <- rep(NA_real_, nrow(model_i))
  rank_vec[main_idx] <- rank(
    -model_i$relative_to_alxa_measured[main_idx],
    ties.method = "average"
  )
  
  model_i %>%
    transmute(
      iter = iter_id,
      site,
      site_label_clean,
      include_in_main_model,
      use_in_main_barplot,
      confidence_tier,
      evidence_score,
      
      relative_to_alxa = relative_to_alxa_measured,
      relative_fe_to_alxa = relative_fe_to_alxa_measured,
      relative_texture_to_alxa = relative_texture_to_alxa_measured,
      rank_main = rank_vec,
      
      substrate_axis,
      mineral_reactivity_axis,
      physical_accessibility_axis,
      surface_exchange_axis,
      climate_suitability,
      microbial_suitability,
      
      feldspar_score,
      inverse_quartz_score,
      size_surface_proxy,
      size_surface_score,
      silt_score,
      cec_score,
      
      tmax_suitability,
      aridity_suitability,
      pet_suitability,
      
      alxa_tmax_anchor = alxa_tmax_anchor_i,
      alxa_ai_anchor = alxa_ai_anchor_i,
      alxa_pet_anchor = alxa_pet_anchor_i
    )
}

mc_results_v18 <- map_dfr(
  seq_len(n_iter),
  run_one_mc
)

write.csv(
  mc_results_v18,
  file.path(output_dir, "v18_monte_carlo_all_iterations_long.csv"),
  row.names = FALSE
)


# 15. Monte Carlo summaries


mc_summary_v18 <- mc_results_v18 %>%
  group_by(
    site,
    site_label_clean,
    include_in_main_model,
    use_in_main_barplot,
    confidence_tier,
    evidence_score
  ) %>%
  summarise(
    n_iter = n(),
    
    relative_median = median(relative_to_alxa, na.rm = TRUE),
    relative_iqr_low = quantile(relative_to_alxa, 0.25, na.rm = TRUE),
    relative_iqr_high = quantile(relative_to_alxa, 0.75, na.rm = TRUE),
    relative_ci95_low = quantile(relative_to_alxa, 0.025, na.rm = TRUE),
    relative_ci95_high = quantile(relative_to_alxa, 0.975, na.rm = TRUE),
    
    substrate_median = median(substrate_axis, na.rm = TRUE),
    substrate_iqr_low = quantile(substrate_axis, 0.25, na.rm = TRUE),
    substrate_iqr_high = quantile(substrate_axis, 0.75, na.rm = TRUE),
    
    mineral_median = median(mineral_reactivity_axis, na.rm = TRUE),
    mineral_iqr_low = quantile(mineral_reactivity_axis, 0.25, na.rm = TRUE),
    mineral_iqr_high = quantile(mineral_reactivity_axis, 0.75, na.rm = TRUE),
    
    physical_median = median(physical_accessibility_axis, na.rm = TRUE),
    physical_iqr_low = quantile(physical_accessibility_axis, 0.25, na.rm = TRUE),
    physical_iqr_high = quantile(physical_accessibility_axis, 0.75, na.rm = TRUE),
    
    surface_median = median(surface_exchange_axis, na.rm = TRUE),
    surface_iqr_low = quantile(surface_exchange_axis, 0.25, na.rm = TRUE),
    surface_iqr_high = quantile(surface_exchange_axis, 0.75, na.rm = TRUE),
    
    climate_median = median(climate_suitability, na.rm = TRUE),
    climate_iqr_low = quantile(climate_suitability, 0.25, na.rm = TRUE),
    climate_iqr_high = quantile(climate_suitability, 0.75, na.rm = TRUE),
    
    microbial_median = median(microbial_suitability, na.rm = TRUE),
    microbial_iqr_low = quantile(microbial_suitability, 0.25, na.rm = TRUE),
    microbial_iqr_high = quantile(microbial_suitability, 0.75, na.rm = TRUE),
    
    relative_fe_median = median(relative_fe_to_alxa, na.rm = TRUE),
    relative_fe_iqr_low = quantile(relative_fe_to_alxa, 0.25, na.rm = TRUE),
    relative_fe_iqr_high = quantile(relative_fe_to_alxa, 0.75, na.rm = TRUE),
    
    relative_texture_median = median(relative_texture_to_alxa, na.rm = TRUE),
    relative_texture_iqr_low = quantile(relative_texture_to_alxa, 0.25, na.rm = TRUE),
    relative_texture_iqr_high = quantile(relative_texture_to_alxa, 0.75, na.rm = TRUE),
    
    rank_median = median(rank_main, na.rm = TRUE),
    rank_iqr_low = quantile(rank_main, 0.25, na.rm = TRUE),
    rank_iqr_high = quantile(rank_main, 0.75, na.rm = TRUE),
    rank_ci95_low = quantile(rank_main, 0.025, na.rm = TRUE),
    rank_ci95_high = quantile(rank_main, 0.975, na.rm = TRUE),
    
    top3_probability = mean(rank_main <= 3, na.rm = TRUE),
    top5_probability = mean(rank_main <= 5, na.rm = TRUE),
    above_alxa_probability = mean(relative_to_alxa > 1, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(desc(relative_median))

write.csv(
  mc_summary_v18,
  file.path(output_dir, "v18_monte_carlo_summary_median_IQR_rank_stability.csv"),
  row.names = FALSE
)

mc_main_summary_v18 <- mc_summary_v18 %>%
  filter(use_in_main_barplot) %>%
  arrange(desc(relative_median))

write.csv(
  mc_main_summary_v18,
  file.path(output_dir, "v18_monte_carlo_main_ranking_summary.csv"),
  row.names = FALSE
)





# 15b. Sensitivity analyses

# V18 sensitivity analyses are generated together in the final section of this script.
# The final section writes table-only outputs and does not create sensitivity figures.


# 15b. Collinearity and model-stability diagnostics, V18


library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)
library(grid)

diagnostic_dir <- file.path(output_dir, "v18_diagnostics")
dir.create(diagnostic_dir, showWarnings = FALSE, recursive = TRUE)


# 1. Helper functions


safe_spearman_cor <- function(df, vars) {
  
  vars <- vars[vars %in% names(df)]
  
  if (length(vars) < 2) {
    stop("Fewer than two variables found for correlation analysis.")
  }
  
  mat <- df %>%
    dplyr::select(all_of(vars)) %>%
    dplyr::mutate(across(everything(), as.numeric))
  
  cor_mat <- cor(
    mat,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  
  cor_mat
}

cor_matrix_to_long <- function(cor_mat) {
  
  as.data.frame(cor_mat) %>%
    tibble::rownames_to_column("var1") %>%
    tidyr::pivot_longer(
      cols = -var1,
      names_to = "var2",
      values_to = "rho"
    ) %>%
    dplyr::mutate(
      abs_rho = abs(rho)
    )
}

make_site_label_diag <- function(site) {
  dplyr::case_when(
    site == "Hobq" ~ "Kubuqi",
    site == "Sonoran_Altar" ~ "Gran Desierto",
    site == "Argentina_MedanosGrandes_NE" ~ "Medanos Grandes",
    site == "Argentina_MedanosGra" ~ "Medanos Grandes",
    site == "Kalahari_Bagatelle" ~ "Kalahari",
    site == "Karakum_AshkAbad" ~ "Karakum",
    site == "Thal_Mankera" ~ "Thal",
    site == "Alxa_Left_Banner_project_base_measured" ~ "Tengger anchor",
    TRUE ~ stringr::str_replace_all(site, "_", " ")
  )
}


# 2. Build diagnostic input table

# Use only global prediction sites included in the main model.
# The Tengger anchor is excluded from correlation diagnostics because
# it is the denominator for the relative score.

diag_input_v18 <- site_master_det %>%
  dplyr::filter(
    use_in_main_barplot,
    full_model_evaluated,
    site != "Alxa_Left_Banner_project_base_measured",
    site != "Alxa_Left_Banner_project_base_downloaded"
  ) %>%
  dplyr::mutate(
    site_label_diag = make_site_label_diag(site)
  )

write.csv(
  diag_input_v18,
  file.path(diagnostic_dir, "v18_diagnostic_input_table.csv"),
  row.names = FALSE
)


# 3. Raw / near-raw input collinearity

# These are the actual environmental and substrate inputs before
# final module aggregation.
#
# Sand-matrix:
#   quartz, feldspar_total, representative_mean_um, median_silt, median_cec
# Climate:
#   tmax_warmest_month, aridity_p_over_pet, pet_for_model
# Microbial environmental inputs:
#   site_p_h, mat
# Microbial response scores:
#   microbe_ph_function, microbe_mat_function

raw_input_vars_v18 <- c(
  "quartz",
  "feldspar_total",
  "quartz_to_feldspar",
  "representative_mean_um",
  "median_silt",
  "median_cec",
  "site_p_h",
  "mat",
  "tmax_warmest_month",
  "aridity_p_over_pet",
  "pet_for_model",
  "microbe_ph_function",
  "microbe_mat_function"
)

raw_input_labels_v18 <- c(
  quartz = "Quartz",
  feldspar_total = "Feldspar",
  quartz_to_feldspar = "Q/F ratio",
  representative_mean_um = "Grain size",
  median_silt = "Silt",
  median_cec = "CEC",
  site_p_h = "pH",
  mat = "MAT",
  tmax_warmest_month = "Tmax",
  aridity_p_over_pet = "AI",
  pet_for_model = "PET",
  microbe_ph_function = "F8 pH response",
  microbe_mat_function = "F8 temperature response"
)

raw_cor_mat_v18 <- safe_spearman_cor(
  diag_input_v18,
  raw_input_vars_v18
)

raw_cor_long_v18 <- cor_matrix_to_long(raw_cor_mat_v18) %>%
  dplyr::mutate(
    var1_label = dplyr::recode(var1, !!!raw_input_labels_v18),
    var2_label = dplyr::recode(var2, !!!raw_input_labels_v18),
    var1_label = factor(var1_label, levels = raw_input_labels_v18[rownames(raw_cor_mat_v18)]),
    var2_label = factor(var2_label, levels = raw_input_labels_v18[colnames(raw_cor_mat_v18)])
  )

write.csv(
  raw_cor_long_v18,
  file.path(diagnostic_dir, "v18_raw_input_spearman_correlation_long.csv"),
  row.names = FALSE
)

raw_high_cor_pairs_v18 <- raw_cor_long_v18 %>%
  dplyr::filter(as.character(var1) < as.character(var2)) %>%
  dplyr::arrange(desc(abs_rho)) %>%
  dplyr::mutate(
    high_correlation_flag = dplyr::case_when(
      abs_rho >= 0.90 ~ "very high |rho| >= 0.90",
      abs_rho >= 0.75 ~ "high |rho| >= 0.75",
      abs_rho >= 0.60 ~ "moderate |rho| >= 0.60",
      TRUE ~ "low"
    )
  )

write.csv(
  raw_high_cor_pairs_v18,
  file.path(diagnostic_dir, "v18_raw_input_high_correlation_pairs.csv"),
  row.names = FALSE
)


# 4. Module-score collinearity

# These are the scores actually entering the final model structure.

module_vars_v18 <- c(
  "mineral_reactivity_axis",
  "physical_accessibility_axis",
  "surface_exchange_axis",
  "tmax_suitability",
  "aridity_suitability",
  "pet_suitability",
  "microbe_ph_function",
  "microbe_mat_function",
  "substrate_axis",
  "climate_suitability",
  "microbial_suitability"
)

module_labels_v18 <- c(
  mineral_reactivity_axis = "Mineral reactivity",
  physical_accessibility_axis = "Physical accessibility",
  surface_exchange_axis = "Surface exchange",
  tmax_suitability = "Heat suitability",
  aridity_suitability = "Aridity suitability",
  pet_suitability = "PET suitability",
  microbe_ph_function = "F8 pH response",
  microbe_mat_function = "F8 temperature response",
  substrate_axis = "Sand-matrix reactivity",
  climate_suitability = "Climate suitability",
  microbial_suitability = "Microbial suitability"
)

module_cor_mat_v18 <- safe_spearman_cor(
  diag_input_v18,
  module_vars_v18
)

module_cor_long_v18 <- cor_matrix_to_long(module_cor_mat_v18) %>%
  dplyr::mutate(
    var1_label = dplyr::recode(var1, !!!module_labels_v18),
    var2_label = dplyr::recode(var2, !!!module_labels_v18),
    var1_label = factor(var1_label, levels = module_labels_v18[rownames(module_cor_mat_v18)]),
    var2_label = factor(var2_label, levels = module_labels_v18[colnames(module_cor_mat_v18)])
  )

write.csv(
  module_cor_long_v18,
  file.path(diagnostic_dir, "v18_module_score_spearman_correlation_long.csv"),
  row.names = FALSE
)

module_high_cor_pairs_v18 <- module_cor_long_v18 %>%
  dplyr::filter(as.character(var1) < as.character(var2)) %>%
  dplyr::arrange(desc(abs_rho)) %>%
  dplyr::mutate(
    high_correlation_flag = dplyr::case_when(
      abs_rho >= 0.90 ~ "very high |rho| >= 0.90",
      abs_rho >= 0.75 ~ "high |rho| >= 0.75",
      abs_rho >= 0.60 ~ "moderate |rho| >= 0.60",
      TRUE ~ "low"
    )
  )

write.csv(
  module_high_cor_pairs_v18,
  file.path(diagnostic_dir, "v18_module_score_high_correlation_pairs.csv"),
  row.names = FALSE
)


# 5. Correlation heatmaps


theme_diag_heatmap <- theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "#222222",
      size = 9
    ),
    axis.text.y = element_text(
      colour = "#222222",
      size = 9
    ),
    axis.title = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 14,
      colour = "#222222"
    ),
    plot.subtitle = element_text(
      size = 10,
      colour = "#555555"
    ),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(8, 12, 8, 8)
  )

p_v18_raw_collinearity <- ggplot(
  raw_cor_long_v18,
  aes(x = var1_label, y = var2_label, fill = rho)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(
    aes(label = sprintf("%.2f", rho)),
    size = 2.5,
    colour = "black"
  ) +
  scale_fill_gradient2(
    low = "#2D7A78",
    mid = "white",
    high = "#A77A18",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho"
  ) +
  coord_equal() +
  labs(
    title = "Input-variable correlation structure",
    subtitle = "Spearman correlations among raw and derived model inputs"
  ) +
  theme_diag_heatmap

print(p_v18_raw_collinearity)

ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "p_v18_raw_collinearity.pdf"),
  plot = p_v18_raw_collinearity,
  width = 8,
  height = 6.5,
  bg = "white")



p_v18_module_collinearity <- ggplot(
  module_cor_long_v18,
  aes(x = var1_label, y = var2_label, fill = rho)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(
    aes(label = sprintf("%.2f", rho)),
    size = 2.5,
    colour = "black"
  ) +
  scale_fill_gradient2(
    low = "#2D7A78",
    mid = "white",
    high = "#A77A18",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho"
  ) +
  coord_equal() +
  labs(
    title = "Module-score correlation structure",
    subtitle = "Spearman correlations among model component scores"
  ) +
  theme_diag_heatmap

print(p_v18_module_collinearity)


# V18 module-score collinearity heatmap
# Remove Heat suitability
# Reorder to match Fig. 7C main order


module_order_v18 <- c(
  "Mineral reactivity",
  "Physical accessibility",
  "Surface exchange",
  "Aridity suitability",
  "PET suitability",
  "F8 pH response",
  "F8 temperature response",
  "Sand-matrix reactivity",
  "Climate suitability",
  "Microbial suitability"
)

module_cor_long_v18_plot <- module_cor_long_v18 %>%
  filter(
    !var1_label %in% c("Heat suitability", "Heat constraint (Tmax)", "Heat (Tmax)"),
    !var2_label %in% c("Heat suitability", "Heat constraint (Tmax)", "Heat (Tmax)")
  ) %>%
  mutate(
    var1_label = factor(var1_label, levels = module_order_v18),
    var2_label = factor(var2_label, levels = rev(module_order_v18))
  )

p_v18_module_collinearity <- ggplot(
  module_cor_long_v18_plot,
  aes(x = var1_label, y = var2_label, fill = rho)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(
    aes(label = sprintf("%.2f", rho)),
    size = 2.5,
    colour = "black"
  ) +
  scale_fill_gradient2(
    low = "#2D7A78",
    mid = "white",
    high = "#A77A18",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho"
  ) +
  coord_equal() +
  labs(
    title = "Module-score correlation structure",
    subtitle = "Spearman correlations among model component scores"
  ) +
  theme_diag_heatmap +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE)

print(p_v18_module_collinearity)


ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "p_v18_module_collinearity.pdf"),
  plot = p_v18_module_collinearity,
  width = 8,
  height = 6.5,
  bg = "white")


ggsave(
  file.path(diagnostic_dir, "v18_module_score_collinearity_heatmap.pdf"),
  p_v18_module_collinearity,
  width = 7.6,
  height = 6.6,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  file.path(diagnostic_dir, "v18_module_score_collinearity_heatmap.png"),
  p_v18_module_collinearity,
  width = 7.6,
  height = 6.6,
  dpi = 600,
  bg = "white"
)


# 6. Monte Carlo rank-stability diagnostics


det_rank_v18 <- deterministic_ranking_v18 %>%
  dplyr::mutate(
    deterministic_rank = dplyr::min_rank(dplyr::desc(relative_to_alxa_measured))
  ) %>%
  dplyr::select(
    site,
    site_label_clean,
    deterministic_relative = relative_to_alxa_measured,
    deterministic_rank
  )

mc_rank_diag_v18 <- mc_results_v18 %>%
  dplyr::filter(
    use_in_main_barplot,
    !is.na(relative_to_alxa),
    !is.na(rank_main)
  ) %>%
  dplyr::left_join(det_rank_v18, by = c("site", "site_label_clean")) %>%
  dplyr::mutate(
    rank_shift_from_deterministic = rank_main - deterministic_rank,
    relative_shift_from_deterministic =
      relative_to_alxa - deterministic_relative
  )

write.csv(
  mc_rank_diag_v18,
  file.path(diagnostic_dir, "v18_mc_rank_stability_long.csv"),
  row.names = FALSE
)

# Per-iteration agreement with deterministic model
mc_iteration_stability_v18 <- mc_rank_diag_v18 %>%
  dplyr::group_by(iter) %>%
  dplyr::summarise(
    n_sites = dplyr::n(),
    spearman_rank_vs_deterministic = suppressWarnings(
      cor(
        deterministic_rank,
        rank_main,
        method = "spearman",
        use = "complete.obs"
      )
    ),
    spearman_score_vs_deterministic = suppressWarnings(
      cor(
        deterministic_relative,
        relative_to_alxa,
        method = "spearman",
        use = "complete.obs"
      )
    ),
    mean_abs_rank_shift = mean(abs(rank_shift_from_deterministic), na.rm = TRUE),
    max_abs_rank_shift = max(abs(rank_shift_from_deterministic), na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  mc_iteration_stability_v18,
  file.path(diagnostic_dir, "v18_mc_iteration_stability_summary.csv"),
  row.names = FALSE
)

mc_stability_summary_v18 <- mc_iteration_stability_v18 %>%
  dplyr::summarise(
    n_iter = dplyr::n(),
    median_spearman_rank = median(spearman_rank_vs_deterministic, na.rm = TRUE),
    iqr_low_spearman_rank = quantile(spearman_rank_vs_deterministic, 0.25, na.rm = TRUE),
    iqr_high_spearman_rank = quantile(spearman_rank_vs_deterministic, 0.75, na.rm = TRUE),
    ci95_low_spearman_rank = quantile(spearman_rank_vs_deterministic, 0.025, na.rm = TRUE),
    ci95_high_spearman_rank = quantile(spearman_rank_vs_deterministic, 0.975, na.rm = TRUE),
    median_mean_abs_rank_shift = median(mean_abs_rank_shift, na.rm = TRUE),
    median_max_abs_rank_shift = median(max_abs_rank_shift, na.rm = TRUE)
  )

print(mc_stability_summary_v18)

write.csv(
  mc_stability_summary_v18,
  file.path(diagnostic_dir, "v18_mc_rank_stability_overall_summary.csv"),
  row.names = FALSE
)

# Site-level rank stability
site_rank_stability_v18 <- mc_rank_diag_v18 %>%
  dplyr::group_by(
    site,
    site_label_clean,
    deterministic_rank,
    deterministic_relative,
    evidence_score
  ) %>%
  dplyr::summarise(
    rank_median = median(rank_main, na.rm = TRUE),
    rank_iqr_low = quantile(rank_main, 0.25, na.rm = TRUE),
    rank_iqr_high = quantile(rank_main, 0.75, na.rm = TRUE),
    rank_ci95_low = quantile(rank_main, 0.025, na.rm = TRUE),
    rank_ci95_high = quantile(rank_main, 0.975, na.rm = TRUE),
    rank_shift_median = median(rank_shift_from_deterministic, na.rm = TRUE),
    rank_shift_abs_median = median(abs(rank_shift_from_deterministic), na.rm = TRUE),
    relative_median = median(relative_to_alxa, na.rm = TRUE),
    relative_iqr_low = quantile(relative_to_alxa, 0.25, na.rm = TRUE),
    relative_iqr_high = quantile(relative_to_alxa, 0.75, na.rm = TRUE),
    top3_probability = mean(rank_main <= 3, na.rm = TRUE),
    top5_probability = mean(rank_main <= 5, na.rm = TRUE),
    above_alxa_probability = mean(relative_to_alxa > 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(deterministic_rank)

write.csv(
  site_rank_stability_v18,
  file.path(diagnostic_dir, "v18_site_level_rank_stability_summary.csv"),
  row.names = FALSE
)


# 7. Rank-stability plots


theme_diag_stability <- theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(colour = "#222222"),
    axis.title = element_text(colour = "#222222"),
    plot.title = element_text(
      face = "bold",
      size = 15,
      colour = "#222222"
    ),
    plot.subtitle = element_text(
      size = 13,
      colour = "#555555"
    ),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(8, 12, 8, 8)
  )

p_v18_mc_rank_agreement <- ggplot(
  mc_iteration_stability_v18,
  aes(x = spearman_rank_vs_deterministic)
) +
  geom_histogram(
    bins = 40,
    fill = "#5F987B",
    colour = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = median(
      mc_iteration_stability_v18$spearman_rank_vs_deterministic,
      na.rm = TRUE
    ),
    linetype = "dashed",
    colour = "#B23A2F",
    linewidth = 0.55
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    title = "Monte Carlo rank agreement with deterministic model",
    subtitle = "Each value represents one Monte Carlo iteration",
    x = "Spearman rank correlation with deterministic ranking",
    y = "Number of iterations"
  ) +
  theme_diag_stability +
  theme(axis.line = element_line(colour = "#333333", linewidth = 0.45))

print(p_v18_mc_rank_agreement)

ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "p_v18_mc_rank_agreement.pdf"),
  plot = p_v18_mc_rank_agreement,
  width = 6.5,
  height = 4.5,
  bg = "white")



# Rank interval plot: lower rank number = higher potential

# Site-level rank stability plot, V18
# Domain shape added:
#   Sandy desert    = circle
#   Gobi transition = square


library(dplyr)
library(ggplot2)
library(scales)
library(grid)


# 1. Prepare plotting table


site_rank_stability_plot_v18 <- site_rank_stability_v18 %>%
  dplyr::mutate(
    site_label_plot = make_site_label_diag(site),
    
    # Keep domain style consistent with other sensitivity panels
    plot_domain = dplyr::case_when(
      site %in% c("Badain_Jaran", "Taklamakan") ~ "Gobi transition",
      TRUE ~ "Sandy desert"
    ),
    plot_domain = factor(
      plot_domain,
      levels = c("Sandy desert", "Gobi transition")
    ),
    
    site_label_plot = factor(
      site_label_plot,
      levels = rev(make_site_label_diag(site_rank_stability_v18$site))
    ),
    
    evidence_score = dplyr::case_when(
      is.na(evidence_score) ~ 0.50,
      evidence_score >= 0.75 ~ 1.00,
      TRUE ~ 0.50
    )
  )


# 2. Domain shape setting


domain_shapes_rank <- c(
  "Sandy desert" = 21,
  "Gobi transition" = 22
)

# If you already defined domain_shapes elsewhere and want to reuse it:
# domain_shapes_rank <- domain_shapes[c("Sandy desert", "Gobi transition")]


# 3. Plot


p_v18_site_rank_stability <- ggplot(
  site_rank_stability_plot_v18,
  aes(y = site_label_plot)
) +
  

  
  # Median rank point
  geom_point(
    aes(
      x = rank_median,
      fill = relative_median,
      size = evidence_score,
      shape = plot_domain
    ),
    colour = "#3A3A3A",
    stroke = 0.65
  ) +
  
  # Lower rank number = higher predicted potential
  scale_x_reverse(
    breaks = seq(
      1,
      nrow(site_rank_stability_plot_v18),
      by = 2
    ),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.03, 0.04))
  ) +
  
  scale_fill_gradientn(
    colours = potential_palette,
    name = "MC median\nrelative score"
  ) +
  # 95% interval
  geom_linerange(
    aes(
      xmin = rank_ci95_low,
      xmax = rank_ci95_high
    ),
    colour = "#BDBDBD",
    linewidth = 0.45
  ) +
  
  # IQR
  geom_linerange(
    aes(
      xmin = rank_iqr_low,
      xmax = rank_iqr_high
    ),
    colour = "#4F4F4F",
    linewidth = 1.40
  ) +
  
  scale_size_continuous(
    name = "Evidence\nscore",
    limits = c(0.50, 1.00),
    breaks = c(0.50, 1.00),
    labels = c("0.50", "1.00"),
    range = c(2.8, 5.2),
    guide = guide_legend(
      order = 2,
      title.position = "top",
      nrow = 1,
      override.aes = list(
        shape = 21,
        fill = "white",
        colour = "#3A3A3A",
        stroke = 0.65
      )
    )
  ) +
  
  scale_shape_manual(
    name = "Domain",
    values = domain_shapes_rank,
    breaks = c("Sandy desert", "Gobi transition"),
    labels = c("Sandy desert", "Gobi transition"),
    drop = FALSE,
    guide = guide_legend(
      order = 1,
      title.position = "top",
      nrow = 1,
      override.aes = list(
        size = 3.8,
        fill = "white",
        colour = "#3A3A3A",
        stroke = 0.75
      )
    )
  ) +
  
  labs(
    title = "Site-level rank stability",
    subtitle = "Points show median ranks; thick bars indicate IQRs; thin bars indicate 95% intervals.",
    x = "Rank in Monte Carlo simulations",
    y = NULL
  ) +
  
  theme_diag_stability +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    axis.text.y = element_text(size = 12, colour = "#222222"),
    axis.text.x = element_text(size = 12, colour = "#222222"),
    axis.title.x = element_text(size = 12, colour = "#222222"),
    plot.title = element_text(
      size = 14,
      face = "bold",
      colour = "#222222"
    ),
    plot.subtitle = element_text(
      size = 10,
      colour = "#555555"
    )
  )

print(p_v18_site_rank_stability)


# 4. Save




ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "p_v18_site_rank_stability.pdf"),
  plot = p_v18_site_rank_stability,
  width = 6,
  height = 5,
  bg = "white")


ggsave(
  file.path(diagnostic_dir, "v18_site_level_rank_stability.pdf"),
  p_v18_site_rank_stability,
  width = 7.2,
  height = 5.8,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  file.path(diagnostic_dir, "v18_site_level_rank_stability.png"),
  p_v18_site_rank_stability,
  width = 7.2,
  height = 5.8,
  dpi = 600,
  bg = "white"
)


# 8. Console diagnostic summary


cat("\n================ V18 diagnostics finished ================\n")
cat("Diagnostic output directory:\n")
cat(diagnostic_dir, "\n\n")

cat("Top raw-input correlations:\n")
print(
  raw_high_cor_pairs_v18 %>%
    dplyr::select(var1_label, var2_label, rho, abs_rho, high_correlation_flag) %>%
    dplyr::slice_head(n = 10)
)

cat("\nTop module-score correlations:\n")
print(
  module_high_cor_pairs_v18 %>%
    dplyr::select(var1_label, var2_label, rho, abs_rho, high_correlation_flag) %>%
    dplyr::slice_head(n = 10)
)

cat("\nMonte Carlo rank-stability summary:\n")
print(mc_stability_summary_v18)

cat("\nKey diagnostic files:\n")
cat("1. v18_raw_input_collinearity_heatmap.pdf/png\n")
cat("2. v18_module_score_collinearity_heatmap.pdf/png\n")
cat("3. v18_raw_input_high_correlation_pairs.csv\n")
cat("4. v18_module_score_high_correlation_pairs.csv\n")
cat("5. v18_mc_rank_agreement_histogram.pdf/png\n")
cat("6. v18_site_level_rank_stability.pdf/png\n")
cat("7. v18_mc_rank_stability_overall_summary.csv\n")
cat("8. v18_site_level_rank_stability_summary.csv\n")




# 16. Lookup tables


temp_lookup_v18 <- tibble(
  temperature = seq(
    response_models_v18[["temperature"]][["function_score"]]$x_min,
    response_models_v18[["temperature"]][["function_score"]]$x_max,
    by = 1
  )
) %>%
  mutate(
    function_score = predict_response(
      response_models_v18[["temperature"]][["function_score"]],
      temperature,
      clamp = TRUE
    )
  )

write.csv(
  temp_lookup_v18,
  file.path(output_dir, "v18_temperature_function_response_lookup_table.csv"),
  row.names = FALSE
)

ph_lookup_v18 <- tibble(
  pH = seq(
    response_models_v18[["pH"]][["function_score"]]$x_min,
    response_models_v18[["pH"]][["function_score"]]$x_max,
    by = 0.1
  )
) %>%
  mutate(
    function_score = predict_response(
      response_models_v18[["pH"]][["function_score"]],
      pH,
      clamp = TRUE
    )
  )

write.csv(
  ph_lookup_v18,
  file.path(output_dir, "v18_pH_function_response_lookup_table.csv"),
  row.names = FALSE
)


# 17. Console report



cat("\n================ V18 model finished ================\n")
cat("Output directory:\n")
cat(output_dir, "\n\n")

cat("Alxa long-term climate anchor values:\n")
cat("Alxa long-term Tmax warmest month:", alxa_tmax_ref, "\n")
cat("Alxa AI:", alxa_ai_ref, "\n")
cat("Alxa PET:", alxa_pet_ref, "\n")
cat("Particle-size treatment: inverse-diameter surface-area proxy, scaled to the prediction domain.\n\n")

cat("Deterministic main ranking:\n")
print(
  deterministic_ranking_v18 %>%
    select(
      site,
      relative_to_alxa_measured,
      substrate_axis,
      mineral_reactivity_axis,
      physical_accessibility_axis,
      surface_exchange_axis,
      climate_suitability,
      microbial_suitability,
      confidence_tier
    )
)

cat("\nMonte Carlo main ranking summary:\n")
print(
  mc_main_summary_v18 %>%
    select(
      site,
      relative_median,
      relative_iqr_low,
      relative_iqr_high,
      top3_probability,
      top5_probability,
      above_alxa_probability,
      rank_median
    )
)

cat("\nKey output files:\n")
cat("1. v18_deterministic_main_ranking_with_modules.csv\n")
cat("2. v18_deterministic_component_scores_wide.csv\n")
cat("3. v18_psd_values_used_in_model.csv\n")
cat("4. v18_monte_carlo_main_ranking_summary.csv\n")
cat("5. v18_monte_carlo_summary_median_IQR_rank_stability.csv\n")
cat("6. v18_monte_carlo_all_iterations_long.csv\n")
cat("7. v18_temperature_response_lookup_table.csv\n")
cat("8. v18_pH_response_lookup_table.csv\n")



#### 7a ####

# Figure 7A world map, V18
# Global projection of DFOE-mediated sand weathering potential
#
# Required objects from V18 model:
#   site_master_det
#   mc_main_summary_v18
#
# Visual encoding:
#   background = arid core + semi-arid transition
#   point fill colour = MC median relative potential
#   point border colour = model stability, P(> anchor)
#   point size = evidence score
#   point shape = dryland domain
#   red point = Tengger experimental anchor

# 1. Packages

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(terra)
library(maps)
library(scales)
library(grid)

if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}
library(ggnewscale)

# 2. Output directory

fig_output_dir <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "v18_figure_outputs")
dir.create(fig_output_dir, showWarnings = FALSE, recursive = TRUE)

# 3. Map extent

xlim_use <- c(-125, 155)
ylim_use <- c(-56, 61.5)

# 4. Aridity background

aridity_bg_path <- file.path(data_dir, "Desktop", "replot", "global_weathering", "background_layers", "global_aridity_classes_AI_lt_050_v31_lowres_df_land_only_latcrop.rds")

if (file.exists(aridity_bg_path)) {
  
  aridity_bg_df <- readRDS(aridity_bg_path) %>%
    filter(
      longitude >= xlim_use[1],
      longitude <= xlim_use[2],
      latitude  >= ylim_use[1],
      latitude  <= ylim_use[2]
    ) %>%
    mutate(
      aridity_class = factor(
        as.character(aridity_class),
        levels = c(
          "Arid core (AI < 0.20)",
          "Semi-arid transition (0.20–0.50)"
        )
      )
    )
  
} else {
  
  ai_path <- file.path(data_dir, "Desktop", "replot", "global_weathering", "aridity_pet", "7504448", "Global-AI_ET0__annual_v3_1", "ai_v31_yr.tif")
  
  ai <- terra::rast(ai_path)
  ai_range <- terra::global(ai, range, na.rm = TRUE)
  print(ai_range)
  
  if (ai_range[1, 2] > 10) {
    ai <- ai / 10000
  }
  
  ai_crop <- terra::crop(
    ai,
    terra::ext(xlim_use[1], xlim_use[2], ylim_use[1], ylim_use[2])
  )
  
  ai_class <- terra::classify(
    ai_crop,
    rcl = matrix(
      c(
        -Inf, 0.20, 1,
        0.20, 0.50, 2,
        0.50, Inf, NA
      ),
      ncol = 3,
      byrow = TRUE
    )
  )
  
  ai_class_ag <- terra::aggregate(
    ai_class,
    fact = 2,
    fun = "modal",
    na.rm = TRUE
  )
  
  aridity_bg_df <- as.data.frame(ai_class_ag, xy = TRUE, na.rm = TRUE)
  colnames(aridity_bg_df) <- c("longitude", "latitude", "aridity_class")
  
  aridity_bg_df$land_region <- maps::map.where(
    database = "world",
    x = aridity_bg_df$longitude,
    y = aridity_bg_df$latitude
  )
  
  aridity_bg_df <- aridity_bg_df %>%
    filter(!is.na(land_region)) %>%
    select(-land_region) %>%
    mutate(
      aridity_class = case_when(
        aridity_class == 1 ~ "Arid core (AI < 0.20)",
        aridity_class == 2 ~ "Semi-arid transition (0.20–0.50)",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(aridity_class)) %>%
    mutate(
      aridity_class = factor(
        aridity_class,
        levels = c(
          "Arid core (AI < 0.20)",
          "Semi-arid transition (0.20–0.50)"
        )
      )
    )
}

tile_width <- median(diff(sort(unique(aridity_bg_df$longitude))), na.rm = TRUE)
tile_height <- median(diff(sort(unique(aridity_bg_df$latitude))), na.rm = TRUE)

if (!is.finite(tile_width) || is.na(tile_width)) {
  tile_width <- 30 / 3600 * 10
}

if (!is.finite(tile_height) || is.na(tile_height)) {
  tile_height <- 30 / 3600 * 10
}

# 5. Clean site table for V18


if (!exists("site_master_det")) {
  stop("Object `site_master_det` not found. Please run the V18 model first.")
}

if (!exists("mc_main_summary_v18")) {
  stop("Object `mc_main_summary_v18` not found. Please run the V18 Monte Carlo summary first.")
}

score_from_text <- function(x) {
  x <- str_to_lower(as.character(x))
  case_when(
    str_detect(x, "high") ~ 1.00,
    str_detect(x, "medium|moderate") ~ 0.75,
    str_detect(x, "proxy|low|lower") ~ 0.50,
    TRUE ~ NA_real_
  )
}

main_sites_v18 <- c(
  "Sonoran_Altar",
  "Argentina_MedanosGrandes_NE",
  "Mu_Us",
  "Hobq",
  "Badain_Jaran",
  "Taklamakan",
  "Mojave",
  "Lower_Draa",
  "Namib",
  "SE_Qatar",
  "Kalahari_Bagatelle",
  "Karakum_AshkAbad",
  "Thal_Mankera"
)

anchor_site <- "Alxa_Left_Banner_project_base_measured"

optional_cols <- c(
  "confidence_tier",
  "data_confidence",
  "mineral_confidence",
  "psd_confidence"
)

for (nm in optional_cols) {
  if (!nm %in% names(site_master_det)) {
    site_master_det[[nm]] <- NA_character_
  }
}

plot_df <- site_master_det %>%
  filter(site %in% c(main_sites_v18, anchor_site)) %>%
  select(
    site,
    site_label,
    country_or_region,
    longitude,
    latitude,
    confidence_tier,
    data_confidence,
    mineral_confidence,
    psd_confidence,
    evidence_score
  ) %>%
  left_join(
    mc_main_summary_v18 %>%
      select(
        site,
        relative_median,
        relative_iqr_low,
        relative_iqr_high,
        above_alxa_probability
      ),
    by = "site"
  ) %>%
  mutate(
    # 
    # Country labels
    # 
    country = case_when(
      site == anchor_site ~ "China",
      !is.na(country_or_region) & country_or_region != "" ~ as.character(country_or_region),
      TRUE ~ NA_character_
    ),
    
    site_label_base = case_when(
      site == anchor_site ~ "Tengger experimental anchor",
      site == "Hobq" ~ "Kubuqi",
      !is.na(site_label) & site_label != "" ~ as.character(site_label),
      TRUE ~ str_replace_all(site, "_", " ")
    ),
    
    site_label_clean = case_when(
      !is.na(country) ~ paste0(site_label_base, ", ", country),
      TRUE ~ site_label_base
    ),
    
    # 
    # Evidence score
    # 
    mineral_score    = score_from_text(mineral_confidence),
    psd_score        = score_from_text(psd_confidence),
    confidence_score = score_from_text(confidence_tier),
    data_score       = score_from_text(data_confidence),
    
    evidence_score = case_when(
      !is.na(mineral_score) | !is.na(psd_score) ~ rowMeans(
        cbind(mineral_score, psd_score),
        na.rm = TRUE
      ),
      !is.na(confidence_score) ~ confidence_score,
      !is.na(data_score) ~ data_score,
      TRUE ~ 0.50
    ),
    
    evidence_score = pmin(pmax(evidence_score, 0.50), 1.00),
    
    # 
    # Domain
    # Basin/Gobi transition is merged into Gobi transition
    # 
    plot_domain = case_when(
      site %in% c("Badain_Jaran", "Taklamakan") ~ "Gobi transition",
      site == anchor_site ~ "Source/reference",
      TRUE ~ "Sandy desert"
    ),
    
    plot_domain = factor(
      plot_domain,
      levels = c(
        "Gobi transition",
        "Sandy desert",
        "Source/reference"
      )
    ),
    
    # 
    # Model stability shown by point border
    # 
    stability_class = case_when(
      site == anchor_site ~ "Anchor",
      above_alxa_probability >= 0.95 ~ "Stable > anchor",
      above_alxa_probability >= 0.50 ~ "Uncertain / near anchor",
      above_alxa_probability < 0.50 ~ "Mostly <= anchor",
      TRUE ~ "Not estimated"
    ),
    
    stability_class = factor(
      stability_class,
      levels = c(
        "Stable > anchor",
        "Uncertain / near anchor",
        "Mostly <= anchor",
        "Not estimated",
        "Anchor"
      )
    )
  ) %>%
  filter(!is.na(longitude), !is.na(latitude))

regular_sites <- plot_df %>%
  filter(site %in% main_sites_v18) %>%
  filter(!is.na(relative_median))

anchor_df <- plot_df %>%
  filter(site == anchor_site)

if (nrow(anchor_df) == 0) {
  anchor_df <- tibble(
    site = anchor_site,
    longitude = 105.7,
    latitude = 38.8,
    relative_median = 1,
    relative_iqr_low = 1,
    relative_iqr_high = 1,
    above_alxa_probability = 1,
    country = "China",
    site_label_base = "Tengger experimental anchor",
    site_label_clean = "Tengger experimental anchor, China",
    evidence_score = 1.00,
    plot_domain = factor("Source/reference", levels = levels(plot_df$plot_domain)),
    stability_class = factor("Anchor", levels = levels(plot_df$stability_class))
  )
}

# 6. Label offsets

label_offsets <- tibble::tribble(
  ~site, ~dx, ~dy,
  "Mojave",                         5.8,  2.5,
  "Sonoran_Altar",                   5.0, -2.2,
  "Argentina_MedanosGrandes_NE",     4.2, -2.2,
  "Lower_Draa",                     -6.2,  2.4,
  "SE_Qatar",                        3.6, -2.5,
  "Namib",                          -7.2, -2.0,
  "Kalahari_Bagatelle",              4.0, 4.2,
  "Karakum_AshkAbad",                6.8,  2.0,
  "Thal_Mankera",                    3.8, -2.0,
  "Taklamakan",                     -5.2,  3.0,
  "Badain_Jaran",                   -6.2,  2.0,
  "Hobq",                            4.3,  1.6,
  "Mu_Us",                           4.3, -1.8,
  "Alxa_Left_Banner_project_base_measured", 4.7, -4.0
)

regular_sites <- regular_sites %>%
  left_join(label_offsets, by = "site") %>%
  mutate(
    dx = replace_na(dx, 0),
    dy = replace_na(dy, 0)
  )

anchor_df <- anchor_df %>%
  left_join(label_offsets, by = "site") %>%
  mutate(
    dx = replace_na(dx, 0),
    dy = replace_na(dy, 0)
  )


# 7. Palettes and scales

world_map <- ggplot2::map_data("world")

aridity_bg_colors <- c(
  "Arid core (AI < 0.20)" = "#E2CF80",
  "Semi-arid transition (0.20-0.50)" = "#F4EFD8"
)

# Make names match the factor level generated above.
levels(aridity_bg_df$aridity_class) <- c(
  "Arid core (AI < 0.20)",
  "Semi-arid transition (0.20-0.50)"
)

potential_palette <- c(
  "#A77A18",
  "#C39A37",
  "#8FAA68",
  "#5F987B",
  "#2D7A78",
  "#123F4B"
)

domain_shapes <- c(
  "Gobi transition" = 22,
  "Sandy desert" = 21,
  "Source/reference" = 21
)

stability_border_colors <- c(
  "Stable > anchor" = "#1B1B1B",
  "Uncertain / near anchor" = "#777777",
  "Mostly <= anchor" = "#C9C9C9",
  "Not estimated" = "#C9C9C9",
  "Anchor" = "#8A1C16"
)

stability_border_colors <- c(
  "Stable > anchor" = "#1B1B1B",
  "Uncertain / near anchor" = "#1B1B1B",
  "Mostly <= anchor" = "#1B1B1B",
  "Not estimated" = "#1B1B1B",
  "Anchor" = "#8A1C16"
)


rel_limits <- range(
  c(regular_sites$relative_median, 1),
  na.rm = TRUE,
  finite = TRUE
)

rel_limits <- c(
  min(0.75, rel_limits[1]),
  max(1.60, rel_limits[2])
)

rel_breaks <- c(0.8, 1.0, 1.2, 1.5)
rel_breaks <- rel_breaks[
  rel_breaks >= rel_limits[1] &
    rel_breaks <= rel_limits[2]
]

# 8. Theme

theme_fig7a <- theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    
    plot.title = element_text(
      size = 22,
      face = "bold",
      hjust = 0,
      colour = "#222222"
    ),
    plot.subtitle = element_text(
      size = 12,
      hjust = 0,
      colour = "grey30",
      margin = margin(b = 8)
    ),
    plot.title.position = "plot",
    plot.margin = margin(10, 18, 10, 10),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.x = unit(3.5, "mm"),
    legend.spacing.y = unit(1.5, "mm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA)
  )


# 9. Build Figure 7A

fig7a_world_v18 <- ggplot() +
  
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "white",
    colour = NA
  ) +
  
  geom_tile(
    data = aridity_bg_df,
    aes(
      x = longitude,
      y = latitude,
      fill = aridity_class
    ),
    width = tile_width,
    height = tile_height,
    alpha = 0.96
  ) +
  
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = NA,
    colour = "#D0CCC2",
    linewidth = 0.30
  ) +
  
  scale_fill_manual(
    name = "Aridity\nbackground",
    values = aridity_bg_colors,
    drop = FALSE,
    guide = "legend"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  geom_point(
    data = regular_sites,
    aes(
      x = longitude,
      y = latitude,
      fill = relative_median,
      size = evidence_score,
      shape = plot_domain,
      colour = stability_class
    ),
    stroke = 0.95,
    alpha = 1
  ) +
  
  geom_point(
    data = anchor_df,
    aes(x = longitude, y = latitude),
    shape = 21,
    size = 4.5,
    fill = "#D73A31",
    colour = "#8A1C16",
    stroke = 0.8,
    show.legend = FALSE
  ) +
  
   # geom_segment(
   # data = regular_sites,
   # aes(
   #    x = longitude,
   #    y = latitude,
   #    xend = longitude + dx * 0.52,
   #    yend = latitude + dy * 0.52
   # ),
   #   colour = alpha("#666666", 0.60),
   #  linewidth = 0.22,
   #   show.legend = FALSE
   # ) +
  
 #   geom_text(
 #     data = regular_sites,
 #    aes(
 #       x = longitude + dx,
 #       y = latitude + dy,
 #      label = site_label_clean
 #    ),
 #     size = 3.10,
 #    colour = "#2F2F2F",
 #    show.legend = FALSE
 #   ) +
 # 
 #   geom_segment(
 #    data = anchor_df,
 #    aes(
 #      x = longitude,
 #       y = latitude,
 #     xend = longitude + dx * 0.52,
 #     yend = latitude + dy * 0.52
 #   ),
 #    colour = alpha("#8A1C16", 0.75),
 #  linewidth = 0.25,
 #    show.legend = FALSE
 # ) +
 # 
 #  geom_text(
 #   data = anchor_df,
 #     aes(
 #    x = longitude + dx,
 #     y = latitude + dy,
 #    label = site_label_clean
 #     ),
 #     size = 3.2,
 #    colour = "#8A1C16",
 #     show.legend = FALSE
 #   ) +
  
  scale_fill_gradientn(
    colours = potential_palette,
    name = "Relative\npotential\n(MC median)",
    limits = rel_limits,
    breaks = rel_breaks,
    labels = paste0(rel_breaks, "x"),
    oob = scales::squish,
    guide = "colourbar"
  ) +
  
  scale_colour_manual(
    name = "P(> anchor)",
    values = stability_border_colors,
    breaks = c(
      "Stable > anchor",
      "Uncertain / near anchor",
      "Mostly <= anchor"
    ),
    labels = c(
      "≥0.95",
      "0.50–0.95",
      "<0.50"
    ),
    guide = guide_legend(
      order = 2,
      title.position = "top",
      nrow = 1,
      override.aes = list(
        shape = 21,
        fill = "white",
        size = 4.5,
        stroke = 1.2
      )
    )
  ) +
  
  scale_size_continuous(
    name = "Evidence\nscore",
    limits = c(0.50, 1.00),
    breaks = c(0.50, 1.00),
    labels = c("0.50", "1.00"),
    range = c(4, 6),
    guide = guide_legend(
      order = 3,
      title.position = "top",
      nrow = 1,
      override.aes = list(
        shape = 21,
        fill = "white",
        colour = "grey25",
        stroke = 1.1
      )
    )
  ) +
  
  scale_shape_manual(
    name = "Domain",
    values = domain_shapes,
    guide = "legend"
  ) +
  
  coord_quickmap(
    xlim = xlim_use,
    ylim = ylim_use,
    expand = FALSE
  ) +
  
  theme_fig7a

print(fig7a_world_v18)




ggsave(
  file.path(out_plot_dir, "fig7a-no-label.pdf"),
  plot = fig7a_world_v18,
  width = 18,
  height = 10,
  bg = "white"
)

ggsave(
  filename = file.path(out_plot_dir, "fig7a_no_label.pdf"),
  plot = fig7a_world_v18,
  width = 18,
  height = 10,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(out_plot_dir, "fig7a_no_label.eps"),
  plot = fig7a_world_v18,
  width = 18,
  height = 10,
  device = cairo_ps,
  bg = "white"
)
ggsave(
  filename = file.path(out_plot_dir, "fig7a_no_label.svg"),
  plot = fig7a_world_v18,
  width = 18,
  height = 10,
  bg = "white"
)


#### figure 7b ####



# Figure 7B ranked point-range plot, V18
# Scheme B:
#   median = colored point
#   IQR    = thick grey line
#   no 95% interval
#   point border = P(> anchor)
#   right text   = P(>A)
#   highest relative_median shown at top
#
# Required object:
#   mc_main_summary_v18
#
# Expected columns in mc_main_summary_v18:
#   site
#   relative_median
#   relative_iqr_low
#   relative_iqr_high
#   above_alxa_probability
#   rank_median

# 1. Packages

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(grid)


# Figure 7B | Ranked MC median relative potential
# Point colour = MC median
# Point size   = evidence score
# Point shape  = dryland domain
# Right text   = P(> anchor)


fig_output_dir <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "v18_figure_outputs")
dir.create(fig_output_dir, showWarnings = FALSE, recursive = TRUE)


# 3. Input check


if (!exists("mc_main_summary_v18")) {
  stop("Object `mc_main_summary_v18` not found. Please run the V18 Monte Carlo summary first.")
}

required_cols <- c(
  "site",
  "relative_median",
  "relative_iqr_low",
  "relative_iqr_high",
  "above_alxa_probability",
  "rank_median",
  "evidence_score"
)

missing_cols <- setdiff(required_cols, names(mc_main_summary_v18))
if (length(missing_cols) > 0) {
  stop(
    "The following required columns are missing in `mc_main_summary_v18`: ",
    paste(missing_cols, collapse = ", ")
  )
}


# 4. Clean table

site_label_map <- c(
  "Sonoran_Altar" = "Gran Desierto",
  "Argentina_MedanosGrandes_NE" = "Medanos Grandes",
  "Mu_Us" = "Mu Us",
  "Hobq" = "Kubuqi",
  "Badain_Jaran" = "Badain Jaran",
  "Taklamakan" = "Taklamakan",
  "Mojave" = "Mojave",
  "Lower_Draa" = "Lower Draa",
  "Namib" = "Namib",
  "SE_Qatar" = "SE Qatar",
  "Kalahari_Bagatelle" = "Kalahari",
  "Karakum_AshkAbad" = "Karakum",
  "Thal_Mankera" = "Thal",
  "Alxa_Left_Banner_project_base_measured" = "Tengger anchor"
)

plot_df_b <- mc_main_summary_v18 %>%
  mutate(
    site_label = case_when(
      site %in% names(site_label_map) ~ unname(site_label_map[site]),
      TRUE ~ str_replace_all(site, "_", " ")
    ),
    
    plot_domain = case_when(
      site %in% c("Badain_Jaran", "Taklamakan") ~ "Gobi transition",
      TRUE ~ "Sandy desert"
    ),
    plot_domain = factor(
      plot_domain,
      levels = c("Sandy desert", "Gobi transition")
    ),
    
    evidence_score = as.numeric(evidence_score),
    p_anchor_label = sprintf("%.3f", above_alxa_probability)
  ) %>%
  arrange(desc(relative_median), rank_median) %>%
  mutate(
    # highest site should appear at the top
    y_id = rev(seq_len(n()))
  )


# 5. Palette / scales

potential_palette <- c(
  "#A77A18",
  "#C39A37",
  "#8FAA68",
  "#5F987B",
  "#2D7A78",
  "#123F4B"
)

x_data_min <- min(plot_df_b$relative_iqr_low, na.rm = TRUE)
x_data_max <- max(plot_df_b$relative_iqr_high, na.rm = TRUE)

# Reserve space for right-side P(> anchor)
x_left  <- min(0.65, x_data_min - 0.05)
x_pcol  <- x_data_max + 0.12
x_right <- x_data_max + 0.26

rel_limits <- c(x_left, x_right)

x_breaks <- c(0.8, 1.0, 1.2, 1.5)
x_breaks <- x_breaks[
  x_breaks >= rel_limits[1] &
    x_breaks <= rel_limits[2]
]

fill_limits <- range(plot_df_b$relative_median, na.rm = TRUE)
fill_breaks <- c(0.8, 1.0, 1.2, 1.5)
fill_breaks <- fill_breaks[
  fill_breaks >= fill_limits[1] &
    fill_breaks <= fill_limits[2]
]


# 6. Theme

theme_fig7b <- theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(
      colour = "#E6E6E6",
      linewidth = 0.35
    ),
    
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 11, colour = "#222222"),
    axis.text.y  = element_text(size = 11, colour = "#222222"),
    axis.text.x  = element_text(size = 10, colour = "#333333"),
    axis.ticks.y = element_blank(),
    
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = margin(10, 22, 10, 10),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.x = unit(3.5, "mm"),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA)
  )


# 7. Build Figure 7B

fig7b_ranked_v18 <- ggplot(plot_df_b, aes(x = relative_median, y = y_id)) +
  
  # Reference line: Tengger experimental anchor
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "#B23A2F"
  ) +
  
  # IQR

  
  # Median points
  geom_point(
    aes(
      fill = relative_median,
      size = evidence_score,
      shape = plot_domain
    ),
    colour = "#1B1B1B",
    stroke = 0.60
  ) +
  
  # Right-side text: P(> anchor)
  geom_text(
    data = plot_df_b,
    aes(
      x = x_pcol,
      y = y_id,
      label = p_anchor_label
    ),
    inherit.aes = FALSE,
    hjust = 0.9,
    size = 3.85,
    colour = "#333333"
  ) +
  
  # Right-side header
  annotate(
    "text",
    x = x_pcol,
    y = max(plot_df_b$y_id) + 0.72,
    label = "P(> anchor)",
    hjust = 0.9,
    vjust = 0,
    size = 4,
    colour = "#222222"
  ) +
  
  # Anchor label
  annotate(
    "text",
    x = 1,
    y = 0.30,
    label = "Tengger experimental anchor",
    colour = "#8A1C16",
    size = 4.25,
    hjust = 0.5,
    vjust = 0.4
  ) +
  
  scale_fill_gradientn(
    colours = potential_palette,
    name = "Relative\npotential\n(MC median)",
    limits = fill_limits,
    breaks = fill_breaks,
    labels = paste0(fill_breaks, "x"),
    oob = scales::squish,
    guide = guide_colourbar(
      title.position = "top",
      label.position = "bottom",
      barheight = unit(4, "mm"),
      barwidth  = unit(34, "mm")
    )
  ) +
  
  scale_size_continuous(
    name = "Evidence\nscore",
    limits = c(0.5, 1.0),
    breaks = c(0.5,  1.0),
    labels = c("0.50", "1.00"),
    range = c(4.0, 6),
    guide = guide_legend(
      title.position = "top",
      nrow = 1,
      override.aes = list(
        shape = 21,
        fill = "white",
        colour = "#1B1B1B",
        stroke = 0.60
      )
    )
  ) +
  
  scale_shape_manual(
    name = "Domain",
    values = c(
      "Sandy desert" = 21,
      "Gobi transition" = 22
    ),
    guide = guide_legend(
      title.position = "top",
      nrow = 1,
      override.aes = list(
        fill = "white",
        colour = "#1B1B1B",
        size = 4.2,
        stroke = 0.60
      )
    )
  ) +
  geom_segment(
    data = plot_df_b,
    aes(
      x = relative_iqr_low,
      xend = relative_iqr_high,
      y = y_id,
      yend = y_id
    ),
    linewidth = 1.25,
    colour = "#37474f",
    lineend = "round"
  ) +
  
  scale_x_continuous(
    limits = rel_limits,
    breaks = x_breaks,
    labels = paste0(x_breaks, "x"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  scale_y_continuous(
    breaks = plot_df_b$y_id,
    labels = plot_df_b$site_label,
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  
  labs(
    x = "Relative potential to Tengger anchor"
  ) +
  
  coord_cartesian(clip = "off") +
  theme_fig7b

print(fig7b_ranked_v18)




#8. Save



ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "fig7b.pdf"),
  plot = fig7b_ranked_v18,
  width = 6.5,
  height = 5
 
)

####figure 7c ####


# Figure 7C grouped heatmap, V18
# Compact component-score heatmap
#
# Updates:
# - Anchor row label: "Tengger experimental anchor"
# - Climate group renamed to "Climate constraint"
# - Climate components renamed to:
#     Heat constraint (Tmax)
#     Aridity constraint (AI)
#     PET constraint
#
# Required objects:
#   site_master_det
#   mc_main_summary_v18


library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(grid)


# 1. Output directory


fig_output_dir <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "v18_figure_outputs")
dir.create(fig_output_dir, showWarnings = FALSE, recursive = TRUE)


# 2. Input checks

if (!exists("site_master_det")) {
  stop("Object `site_master_det` not found.")
}

if (!exists("mc_main_summary_v18")) {
  stop("Object `mc_main_summary_v18` not found.")
}

required_component_cols <- c(
  "site",
  "mineral_reactivity_axis",
  "physical_accessibility_axis",
  "surface_exchange_axis",
  "tmax_suitability",
  "aridity_suitability",
  "pet_suitability",
  "microbe_ph_function",
  "microbe_mat_function"
)

missing_component_cols <- setdiff(required_component_cols, names(site_master_det))

if (length(missing_component_cols) > 0) {
  stop(
    "Missing required columns in `site_master_det`: ",
    paste(missing_component_cols, collapse = ", ")
  )
}


# 3. Helper


geo_mean_vec <- function(a, b) {
  out <- exp(rowMeans(log(pmax(cbind(a, b), 1e-6)), na.rm = FALSE))
  out[is.na(a) | is.na(b)] <- NA_real_
  out
}


# 4. Site labels and ranking order from 7B


site_label_map <- c(
  "Sonoran_Altar" = "Gran Desierto",
  "Argentina_MedanosGrandes_NE" = "Medanos Grandes",
  "Mu_Us" = "Mu Us",
  "Hobq" = "Hobq",
  "Badain_Jaran" = "Badain Jaran",
  "Taklamakan" = "Taklamakan",
  "Mojave" = "Mojave",
  "Lower_Draa" = "Lower Draa",
  "Namib" = "Namib",
  "SE_Qatar" = "SE Qatar",
  "Kalahari_Bagatelle" = "Kalahari",
  "Karakum_AshkAbad" = "Karakum",
  "Thal_Mankera" = "Thal",
  "Alxa_Left_Banner_project_base_measured" = "Tengger anchor"
)

ranked_sites_df <- mc_main_summary_v18 %>%
  mutate(
    site_label = case_when(
      site %in% names(site_label_map) ~ unname(site_label_map[site]),
      TRUE ~ str_replace_all(site, "_", " ")
    )
  ) %>%
  arrange(desc(relative_median), rank_median)

ranked_site_order  <- ranked_sites_df$site
ranked_site_labels <- ranked_sites_df$site_label

anchor_site_name <- "Alxa_Left_Banner_project_base_measured"

if (!(anchor_site_name %in% site_master_det$site)) {
  stop("Anchor site `Alxa_Left_Banner_project_base_measured` not found in `site_master_det`.")
}

# Final order: ranked global sites + Tengger experimental anchor
site_order_final  <- c(ranked_site_order, anchor_site_name)
site_label_final  <- c(ranked_site_labels, "Tengger experimental anchor")


# 5. Build component table


heat_wide <- site_master_det %>%
  filter(site %in% site_order_final) %>%
  mutate(
    ph_component_score = microbe_ph_function,
    mat_component_score = microbe_mat_function,
    site_label = case_when(
      site %in% names(site_label_map) ~ unname(site_label_map[site]),
      TRUE ~ str_replace_all(site, "_", " ")
    ),
    site_label = ifelse(
      site == anchor_site_name,
      "Tengger experimental anchor",
      site_label
    )
  ) %>%
  select(
    site,
    site_label,
    mineral_reactivity_axis,
    physical_accessibility_axis,
    surface_exchange_axis,
    tmax_suitability,
    aridity_suitability,
    pet_suitability,
    ph_component_score,
    mat_component_score
  )

component_vars <- c(
  "mineral_reactivity_axis",
  "physical_accessibility_axis",
  "surface_exchange_axis",
  "tmax_suitability",
  "aridity_suitability",
  "pet_suitability",
  "ph_component_score",
  "mat_component_score"
)

component_labels <- c(
  "Mineral\nreactivity",
  "Physical\naccessibility",
  "Surface\nexchange",
  "Heat constraint\n(Tmax)",
  "Aridity constraint\n(AI)",
  "PET constraint",
  "pH\nsuitability",
  "Thermal\nsuitability"
)

heat_df <- heat_wide %>%
  mutate(
    site = factor(site, levels = site_order_final),
    site_label = factor(site_label, levels = rev(site_label_final))
  ) %>%
  pivot_longer(
    cols = all_of(component_vars),
    names_to = "component",
    values_to = "score"
  ) %>%
  mutate(
    component_group = case_when(
      component %in% c(
        "mineral_reactivity_axis",
        "physical_accessibility_axis",
        "surface_exchange_axis"
      ) ~ "Substrate",
      component %in% c(
        "tmax_suitability",
        "aridity_suitability",
        "pet_suitability"
      ) ~ "Climate constraint",
      component %in% c(
        "ph_component_score",
        "mat_component_score"
      ) ~ "Microbial",
      TRUE ~ "Other"
    ),
    
    component_label = case_when(
      component == "mineral_reactivity_axis"     ~ "Mineral\nreactivity",
      component == "physical_accessibility_axis" ~ "Physical\naccessibility",
      component == "surface_exchange_axis"       ~ "Surface\nexchange",
      component == "tmax_suitability"            ~ "Heat constraint\n(Tmax)",
      component == "aridity_suitability"         ~ "Aridity constraint\n(AI)",
      component == "pet_suitability"             ~ "PET constraint",
      component == "ph_component_score"          ~ "pH\nsuitability",
      component == "mat_component_score"         ~ "Thermal\nsuitability",
      TRUE ~ component
    ),
    
    component = factor(component, levels = component_vars),
    component_label = factor(component_label, levels = component_labels),
    x_id = match(component, component_vars)
  )


# 6. Palette


component_palette <- c(
  "#A77A18",
  "#C39A37",
  "#8FAA68",
  "#5F987B",
  "#2D7A78",
  "#123F4B"
)

fill_limits <- c(0, 1)
fill_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)


# 7. Theme

theme_fig7c_compact <- theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    
    axis.title = element_blank(),
    axis.text.y = element_text(size = 8.8, colour = "#222222"),
    axis.text.x = element_text(
      size = 8.0,
      colour = "#222222",
      angle = 35,
      hjust = 1,
      vjust = 1,
      lineheight = 0.92
    ),
    axis.ticks = element_blank(),
    
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = margin(8, 6, 4, 6),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.x = unit(2.5, "mm"),
    legend.title = element_text(size = 8.8),
    legend.text  = element_text(size = 8.0),
    legend.key.height = unit(3.2, "mm"),
    legend.key.width  = unit(7, "mm"),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA)
  )

# 8. Build Figure 7C

n_sites_total <- length(site_label_final)

fig7c_heatmap_v18 <- ggplot(
  heat_df,
  aes(x = x_id, y = site_label, fill = score)
) +
  
  geom_tile(
    width = 0.94,
    height = 0.92,
    colour = "white",
    linewidth = 0.55
  ) +
  
  # vertical separators between groups
  geom_vline(
    xintercept = c(3.5, 6.5),
    colour = "#D0D0D0",
    linewidth = 0.7
  ) +
  
  # horizontal separator above Tengger anchor row
  geom_hline(
    yintercept = 1.5,
    colour = "#D0D0D0",
    linewidth = 0.8
  ) +
  
  # group labels
  annotate(
    "text",
    x = 2.0,
    y = n_sites_total + 0.92,
    label = "Substrate",
    fontface = "bold",
    size = 3.2,
    colour = "#222222"
  ) +
  annotate(
    "text",
    x = 5.0,
    y = n_sites_total + 0.92,
    label = "Climate constraint",
    fontface = "bold",
    size = 3.2,
    colour = "#222222"
  ) +
  annotate(
    "text",
    x = 7.5,
    y = n_sites_total + 0.92,
    label = "Microbial",
    fontface = "bold",
    size = 3.2,
    colour = "#222222"
  ) +
  
  # group underlines
  annotate(
    "segment",
    x = 0.65,
    xend = 3.35,
    y = n_sites_total + 0.62,
    yend = n_sites_total + 0.62,
    linewidth = 0.55,
    colour = "#8C8C8C"
  ) +
  annotate(
    "segment",
    x = 3.65,
    xend = 6.35,
    y = n_sites_total + 0.62,
    yend = n_sites_total + 0.62,
    linewidth = 0.55,
    colour = "#8C8C8C"
  ) +
  annotate(
    "segment",
    x = 6.65,
    xend = 8.35,
    y = n_sites_total + 0.62,
    yend = n_sites_total + 0.62,
    linewidth = 0.55,
    colour = "#8C8C8C"
  ) +
  
  scale_x_continuous(
    breaks = 1:8,
    labels = component_labels,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  scale_y_discrete(
    labels = rev(site_label_final)
  ) +
  
  scale_fill_gradientn(
    colours = component_palette,
    limits = fill_limits,
    breaks = fill_breaks,
    labels = c("0", "0.25", "0.50", "0.75", "1.00"),
    oob = scales::squish,
    name = "Component\nscore",
    guide = guide_colourbar(
      title.position = "top",
      label.position = "bottom",
      barheight = unit(3.6, "mm"),
      barwidth  = unit(28, "mm")
    )
  ) +
  
  coord_cartesian(
    ylim = c(0.5, n_sites_total + 1.1),
    clip = "off"
  ) +
  
  theme_fig7c_compact

print(fig7c_heatmap_v18)

ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "fig7c-1.pdf"),
  plot = fig7c_heatmap_v18,
  width = 6,
  height = 7
  
)

# 9. Save

ggsave(
  filename = file.path(fig_output_dir, "figure_7C_component_heatmap_v18_climate_constraint_anchor.pdf"),
  plot = fig7c_heatmap_v18,
  width = 5.2,
  height = 4.4,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(fig_output_dir, "figure_7C_component_heatmap_v18_climate_constraint_anchor.png"),
  plot = fig7c_heatmap_v18,
  width = 5.2,
  height = 4.4,
  units = "in",
  dpi = 600
)


####figure 7c-1 ####

# Figure 7C grouped heatmap, V18
# Portrait / rotated version
#
# Layout:
#   x-axis = sites
#   y-axis = components
#   width < height
#
# Updates included:
# - Anchor label: "Tengger experimental anchor"
# - Climate group: "Climate constraint"
# - Climate sublabels:
#     Heat constraint (Tmax)
#     Aridity constraint (AI)
#     PET constraint

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(grid)

# 1. Output directory

fig_output_dir <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "v18_figure_outputs")
dir.create(fig_output_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Input checks

if (!exists("site_master_det")) {
  stop("Object `site_master_det` not found.")
}

if (!exists("mc_main_summary_v18")) {
  stop("Object `mc_main_summary_v18` not found.")
}

required_component_cols <- c(
  "site",
  "mineral_reactivity_axis",
  "physical_accessibility_axis",
  "surface_exchange_axis",
  "tmax_suitability",
  "aridity_suitability",
  "pet_suitability",
  "microbe_ph_function",
  "microbe_mat_function"
)

missing_component_cols <- setdiff(required_component_cols, names(site_master_det))

if (length(missing_component_cols) > 0) {
  stop(
    "Missing required columns in `site_master_det`: ",
    paste(missing_component_cols, collapse = ", ")
  )
}

# 3. Helper

geo_mean_vec <- function(a, b) {
  out <- exp(rowMeans(log(pmax(cbind(a, b), 1e-6)), na.rm = FALSE))
  out[is.na(a) | is.na(b)] <- NA_real_
  out
}

# 4. Site labels and order (same ranking as 7B)

site_label_map <- c(
  "Sonoran_Altar" = "Gran Desierto",
  "Argentina_MedanosGrandes_NE" = "Medanos Grandes",
  "Mu_Us" = "Mu Us",
  "Hobq" = "Hobq",
  "Badain_Jaran" = "Badain Jaran",
  "Taklamakan" = "Taklamakan",
  "Mojave" = "Mojave",
  "Lower_Draa" = "Lower Draa",
  "Namib" = "Namib",
  "SE_Qatar" = "SE Qatar",
  "Kalahari_Bagatelle" = "Kalahari",
  "Karakum_AshkAbad" = "Karakum",
  "Thal_Mankera" = "Thal",
  "Alxa_Left_Banner_project_base_measured" = "Tengger anchor"
)

ranked_sites_df <- mc_main_summary_v18 %>%
  mutate(
    site_label = case_when(
      site %in% names(site_label_map) ~ unname(site_label_map[site]),
      TRUE ~ str_replace_all(site, "_", " ")
    )
  ) %>%
  arrange(desc(relative_median), rank_median)

ranked_site_order  <- ranked_sites_df$site
ranked_site_labels <- ranked_sites_df$site_label

anchor_site_name <- "Alxa_Left_Banner_project_base_measured"

if (!(anchor_site_name %in% site_master_det$site)) {
  stop("Anchor site `Alxa_Left_Banner_project_base_measured` not found in `site_master_det`.")
}

# Final order: ranked sites + anchor at the end
site_order_final <- c(ranked_site_order, anchor_site_name)
site_label_final <- c(ranked_site_labels, "Tengger experimental anchor")

# 5. Build component table

heat_wide <- site_master_det %>%
  filter(site %in% site_order_final) %>%
  mutate(
    ph_component_score = microbe_ph_function,
    mat_component_score = microbe_mat_function,
    site_label = case_when(
      site %in% names(site_label_map) ~ unname(site_label_map[site]),
      TRUE ~ str_replace_all(site, "_", " ")
    ),
    site_label = ifelse(
      site == anchor_site_name,
      "Tengger experimental anchor",
      site_label
    )
  ) %>%
  select(
    site,
    site_label,
    mineral_reactivity_axis,
    physical_accessibility_axis,
    surface_exchange_axis,
    tmax_suitability,
    aridity_suitability,
    pet_suitability,
    ph_component_score,
    mat_component_score
  )

component_vars <- c(
  "mineral_reactivity_axis",
  "physical_accessibility_axis",
  "surface_exchange_axis",
  "tmax_suitability",
  "aridity_suitability",
  "pet_suitability",
  "ph_component_score",
  "mat_component_score"
)

component_label_order <- c(
  "Mineral reactivity",
  "Physical accessibility",
  "Surface exchange",
  "Heat constraint (Tmax)",
  "Aridity constraint (AI)",
  "PET constraint",
  "pH suitability",
  "Thermal suitability"
)

heat_df <- heat_wide %>%
  mutate(
    site = factor(site, levels = site_order_final),
    site_label = factor(site_label, levels = site_label_final)
  ) %>%
  pivot_longer(
    cols = all_of(component_vars),
    names_to = "component",
    values_to = "score"
  ) %>%
  mutate(
    component_group = case_when(
      component %in% c(
        "mineral_reactivity_axis",
        "physical_accessibility_axis",
        "surface_exchange_axis"
      ) ~ "Substrate",
      component %in% c(
        "tmax_suitability",
        "aridity_suitability",
        "pet_suitability"
      ) ~ "Climate constraint",
      component %in% c(
        "ph_component_score",
        "mat_component_score"
      ) ~ "Microbial",
      TRUE ~ "Other"
    ),
    
    component_label = case_when(
      component == "mineral_reactivity_axis"     ~ "Mineral reactivity",
      component == "physical_accessibility_axis" ~ "Physical accessibility",
      component == "surface_exchange_axis"       ~ "Surface exchange",
      component == "tmax_suitability"            ~ "Heat constraint (Tmax)",
      component == "aridity_suitability"         ~ "Aridity constraint (AI)",
      component == "pet_suitability"             ~ "PET constraint",
      component == "ph_component_score"          ~ "pH suitability",
      component == "mat_component_score"         ~ "Thermal suitability",
      TRUE ~ component
    ),
    
    component_group = factor(
      component_group,
      levels = c("Substrate", "Climate constraint", "Microbial")
    ),
    
    component_label = factor(
      component_label,
      levels = rev(component_label_order)
    )
  )

# anchor separator location
anchor_x_pos <- length(site_label_final) - 0.5

# 6. Palette

component_palette <- c(
  "#A77A18",
  "#C39A37",
  "#8FAA68",
  "#5F987B",
  "#2D7A78",
  "#123F4B"
)

fill_limits <- c(0, 1)
fill_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)

# 7. Theme (portrait)

theme_fig7c_portrait <- theme_minimal(base_size = 10.5) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    
    axis.title = element_blank(),
    
    axis.text.x = element_text(
      size = 7.8,
      colour = "#222222",
      angle = 50,
      hjust = 1,
      vjust = 1,
      lineheight = 0.92
    ),
    axis.text.y = element_text(
      size = 8.2,
      colour = "#222222",
      lineheight = 0.94
    ),
    axis.ticks = element_blank(),
    
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      size = 9.0,
      face = "bold",
      angle = 0,
      colour = "#222222"
    ),
    
    panel.spacing.y = unit(1.6, "mm"),
    
    plot.margin = margin(8, 8, 4, 6),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.x = unit(2.5, "mm"),
    legend.title = element_text(size = 8.8),
    legend.text  = element_text(size = 8.0),
    legend.key.height = unit(3.2, "mm"),
    legend.key.width  = unit(7, "mm"),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA)
  )

# 8. Build portrait / rotated heatmap

fig7c_heatmap_v18_portrait <- ggplot(
  heat_df,
  aes(x = site_label, y = component_label, fill = score)
) +
  
  geom_tile(
    width = 0.94,
    height = 0.90,
    colour = "white",
    linewidth = 0.55
  ) +
  
  # separator before anchor column
  geom_vline(
    xintercept = anchor_x_pos,
    colour = "#D0D0D0",
    linewidth = 0.85
  ) +
  
  facet_grid(
    rows = vars(component_group),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  scale_fill_gradientn(
    colours = component_palette,
    limits = fill_limits,
    breaks = fill_breaks,
    labels = c("0", "0.25", "0.50", "0.75", "1.00"),
    oob = scales::squish,
    name = "Component\nscore",
    guide = guide_colourbar(
      title.position = "top",
      label.position = "bottom",
      barheight = unit(3.6, "mm"),
      barwidth  = unit(28, "mm")
    )
  ) +
  
  theme_fig7c_portrait

print(fig7c_heatmap_v18_portrait)


#remove PET

# 只去掉 Climate constraint 里面的 Heat constraint (Tmax) 这一行
# 其他绘图参数、顺序、配色、theme 全部不变

heat_df_v18 <- heat_df %>%
  filter(
    !(component_group == "Climate constraint" &
        component_label == "Heat constraint (Tmax)")
  ) %>%
  droplevels()

fig7c_heatmap_v18_portrait <- ggplot(
  heat_df_v18,
  aes(x = site_label, y = component_label, fill = score)
) +
  
  geom_tile(
    width = 0.94,
    height = 0.90,
    colour = "white",
    linewidth = 0.55
  ) +
  
  # separator before anchor column
  geom_vline(
    xintercept = anchor_x_pos,
    colour = "#D0D0D0",
    linewidth = 0.85
  ) +
  
  facet_grid(
    rows = vars(component_group),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  scale_fill_gradientn(
    colours = component_palette,
    limits = fill_limits,
    breaks = fill_breaks,
    labels = c("0", "0.25", "0.50", "0.75", "1.00"),
    oob = scales::squish,
    name = "Component\nscore",
    guide = guide_colourbar(
      title.position = "top",
      label.position = "bottom",
      barheight = unit(3.6, "mm"),
      barwidth  = unit(28, "mm")
    )
  ) +
  
  theme_fig7c_portrait



ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "fig7c.pdf"),
  plot = fig7c_heatmap_v18_portrait,
  width = 8,
  height = 5
  
)



#### S1 model construction diagnostics ####

# Fig. S15. Model construction diagnostics, V18
#
# Purpose:
#   A. Substrate PCA diagnostic:
#      show that candidate sites occupy different substrate space.
#
#   B. Climate PCA diagnostic:
#      show that candidate sites occupy different climate-stress space.
#
#   C. Climate one-sided constraint diagnostic:
#      show how Heat / Aridity / PET constraint scores were constructed
#      relative to the Tengger experimental anchor.
#
#   D. Microbial response curves:
#      show experimental pH and temperature response curves used for
#      microbial suitability projection.
#
# Required objects:
#   site_master_det
#   mc_main_summary_v18
#
# Optional objects for panel D:
#   response_models_v18
#   predict_response()
#   response_points_v18 / response_data_v18 / microbial_response_points_v18
# 1. Packages

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(grid)

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
library(patchwork)

# 2. Output directory

fig_output_dir <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "v18_figure_outputs")
dir.create(fig_output_dir, showWarnings = FALSE, recursive = TRUE)


# 3. Input checks

if (!exists("site_master_det")) {
  stop("Object `site_master_det` not found. Please run the V18 model first.")
}

if (!exists("mc_main_summary_v18")) {
  stop("Object `mc_main_summary_v18` not found. Please run the V18 Monte Carlo summary first.")
}

# 4. Helper functions

pick_first_existing <- function(candidates, nm) {
  hit <- candidates[candidates %in% nm][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  hit
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

rescale_01 <- function(x, min_score = 0, max_score = 1) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep((min_score + max_score) / 2, length(x)))
  }
  out <- (x - rng[1]) / diff(rng)
  out <- out * (max_score - min_score) + min_score
  out
}

geo_mean_vec <- function(...) {
  mat <- cbind(...)
  out <- exp(rowMeans(log(pmax(mat, 1e-6)), na.rm = FALSE))
  out[rowSums(is.na(mat)) > 0] <- NA_real_
  out
}

clamp_score <- function(x, lower = 0.05, upper = 1) {
  pmin(pmax(x, lower), upper)
}

anchor_upper_constraint <- function(x, anchor, lower = 0.05) {
  out <- ifelse(
    is.na(x) | is.na(anchor) | anchor <= 0,
    NA_real_,
    ifelse(x <= anchor, 1, anchor / x)
  )
  clamp_score(out, lower = lower, upper = 1)
}

anchor_lower_constraint <- function(x, anchor, lower = 0.05) {
  out <- ifelse(
    is.na(x) | is.na(anchor) | anchor <= 0,
    NA_real_,
    ifelse(x >= anchor, 1, x / anchor)
  )
  clamp_score(out, lower = lower, upper = 1)
}

make_text_panel <- function(message) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = message,
      hjust = 0.5,
      vjust = 0.5,
      size = 3.4,
      colour = "grey30"
    ) +
    theme_void()
}

# 5. Site order, labels, domain, evidence score

anchor_site <- "Alxa_Left_Banner_project_base_measured"

site_label_map <- c(
  "Sonoran_Altar" = "Gran Desierto",
  "Argentina_MedanosGrandes_NE" = "Medanos Grandes",
  "Mu_Us" = "Mu Us",
  "Hobq" = "Hobq",
  "Badain_Jaran" = "Badain Jaran",
  "Taklamakan" = "Taklamakan",
  "Mojave" = "Mojave",
  "Lower_Draa" = "Lower Draa",
  "Namib" = "Namib",
  "SE_Qatar" = "SE Qatar",
  "Kalahari_Bagatelle" = "Kalahari",
  "Karakum_AshkAbad" = "Karakum",
  "Thal_Mankera" = "Thal",
  "Alxa_Left_Banner_project_base_measured" = "Tengger anchor"
)

site_order_df <- mc_main_summary_v18 %>%
  mutate(
    site_label = case_when(
      site %in% names(site_label_map) ~ unname(site_label_map[site]),
      TRUE ~ str_replace_all(site, "_", " ")
    )
  ) %>%
  arrange(desc(relative_median), rank_median)

site_order <- site_order_df$site
site_label_order <- site_order_df$site_label

if (anchor_site %in% site_master_det$site) {
  site_order_final <- c(site_order, anchor_site)
  site_label_final <- c(site_label_order, "Tengger experimental anchor")
} else {
  site_order_final <- site_order
  site_label_final <- site_label_order
}

score_from_text <- function(x) {
  x <- str_to_lower(as.character(x))
  case_when(
    str_detect(x, "high") ~ 1.00,
    str_detect(x, "medium|moderate") ~ 0.75,
    str_detect(x, "proxy|low|lower") ~ 0.50,
    TRUE ~ NA_real_
  )
}

optional_cols <- c(
  "confidence_tier",
  "data_confidence",
  "mineral_confidence",
  "psd_confidence",
  "domain",
  "dryland_domain",
  "site_domain",
  "plot_domain"
)

for (cc in optional_cols) {
  if (!cc %in% names(site_master_det)) {
    site_master_det[[cc]] <- NA
  }
}

site_meta <- site_master_det %>%
  filter(site %in% site_order_final) %>%
  mutate(
    site_label = case_when(
      site %in% names(site_label_map) ~ unname(site_label_map[site]),
      TRUE ~ str_replace_all(site, "_", " ")
    ),
    site_label = ifelse(site == anchor_site, "Tengger experimental anchor", site_label),
    
    mineral_score = score_from_text(mineral_confidence),
    psd_score = score_from_text(psd_confidence),
    confidence_score = score_from_text(confidence_tier),
    data_score = score_from_text(data_confidence),
    
    evidence_score = case_when(
      !is.na(mineral_score) | !is.na(psd_score) ~ rowMeans(
        cbind(mineral_score, psd_score),
        na.rm = TRUE
      ),
      !is.na(confidence_score) ~ confidence_score,
      !is.na(data_score) ~ data_score,
      TRUE ~ 0.50
    ),
    evidence_score = pmin(pmax(evidence_score, 0.50), 1.00),
    
    raw_domain = coalesce(
      as.character(plot_domain),
      as.character(domain),
      as.character(dryland_domain),
      as.character(site_domain)
    ),
    plot_domain = case_when(
      site == anchor_site ~ "Source/reference",
      site %in% c("Badain_Jaran", "Taklamakan") ~ "Gobi transition",
      str_detect(str_to_lower(raw_domain), "gobi") ~ "Gobi transition",
      TRUE ~ "Sandy desert"
    ),
    plot_domain = factor(
      plot_domain,
      levels = c("Sandy desert", "Gobi transition", "Source/reference")
    )
  ) %>%
  select(site, site_label, evidence_score, plot_domain)

# 6. Shared palettes and themes

potential_palette <- c(
  "#A77A18",
  "#C39A37",
  "#8FAA68",
  "#5F987B",
  "#2D7A78",
  "#123F4B"
)

domain_shapes <- c(
  "Sandy desert" = 21,
  "Gobi transition" = 22,
  "Source/reference" = 21
)


theme_pca_s15 <- theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "#E8E8E8", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(size = 13, colour = "#222222"),
    axis.text  = element_text(size = 12, colour = "#222222"),
    
    plot.title = element_text(size = 12.5, face = "bold", colour = "#222222"),
    plot.subtitle = element_text(size = 8.5, colour = "grey35", margin = margin(b = 4)),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.margin = margin(6, 8, 6, 6)
  )

theme_dot_s15 <- theme_minimal(base_size = 10.5) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#E8E8E8", linewidth = 0.35),
    
    axis.title = element_text(size = 9.5, colour = "#222222"),
    axis.text = element_text(size = 8.5, colour = "#222222"),
    
    plot.title = element_text(size = 12.5, face = "bold", colour = "#222222"),
    plot.subtitle = element_text(size = 8.5, colour = "grey35", margin = margin(b = 4)),
    
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.0),
    plot.margin = margin(6, 8, 6, 6)
  )

theme_curve_s15 <- theme_minimal(base_size = 10.5) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "#E8E8E8", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(size = 9.5, colour = "#222222"),
    axis.text = element_text(size = 8.5, colour = "#222222"),
    
    plot.title = element_text(size = 12.5, face = "bold", colour = "#222222"),
    plot.subtitle = element_text(size = 8.5, colour = "grey35", margin = margin(b = 4)),
    
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.0),
    plot.margin = margin(6, 8, 6, 6)
  )

# 7. Generic PCA biplot builder

build_pca_biplot <- function(
    data_wide,
    variable_cols,
    variable_labels,
    colour_col,
    colour_name,
    title,
    subtitle,
    arrow_label_size = 3.0
) {
  
  pca_input <- data_wide %>%
    select(site, site_label, evidence_score, plot_domain, all_of(variable_cols), all_of(colour_col)) %>%
    filter(if_all(all_of(variable_cols), ~ !is.na(.x)))
  
  if (nrow(pca_input) < 4) {
    return(make_text_panel(paste0("Not enough complete rows for PCA: ", title)))
  }
  
  pca_mat <- pca_input %>%
    select(all_of(variable_cols)) %>%
    as.data.frame()
  
  pca_fit <- prcomp(pca_mat, center = TRUE, scale. = TRUE)
  
  var_exp <- summary(pca_fit)$importance[2, 1:2] * 100
  
  score_df <- as.data.frame(pca_fit$x[, 1:2]) %>%
    rename(PC1 = PC1, PC2 = PC2) %>%
    bind_cols(
      pca_input %>%
        select(site, site_label, evidence_score, plot_domain, all_of(colour_col))
    ) %>%
    mutate(
      is_anchor = site == anchor_site
    )
  
  ## 关键修改：普通点和 anchor 分开
  score_regular <- score_df %>%
    filter(!is_anchor)
  
  score_anchor <- score_df %>%
    filter(is_anchor)
  
  loading_df <- as.data.frame(pca_fit$rotation[, 1:2]) %>%
    tibble::rownames_to_column("variable") %>%
    mutate(
      variable_label = variable_labels[variable],
      variable_label = ifelse(is.na(variable_label), variable, variable_label)
    )
  
  score_range <- max(
    diff(range(score_df$PC1, na.rm = TRUE)),
    diff(range(score_df$PC2, na.rm = TRUE))
  )
  
  loading_range <- max(abs(c(loading_df$PC1, loading_df$PC2)), na.rm = TRUE)
  
  arrow_scale <- ifelse(
    is.finite(loading_range) && loading_range > 0,
    score_range * 0.35 / loading_range,
    1
  )
  
  loading_df <- loading_df %>%
    mutate(
      xend = PC1 * arrow_scale,
      yend = PC2 * arrow_scale
    )
  
  ggplot(score_df, aes(x = PC1, y = PC2)) +
    geom_hline(yintercept = 0, colour = "#D7D7D7", linewidth = 0.35) +
    geom_vline(xintercept = 0, colour = "#D7D7D7", linewidth = 0.35) +
    
    geom_segment(
      data = loading_df,
      aes(x = 0, y = 0, xend = xend, yend = yend),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(2.0, "mm")),
      colour = "#6E6E6E",
      linewidth = 0.45
    ) +
    
    ## 普通点：不再包含 anchor
    geom_point(
      data = score_regular,
      aes(
        x = PC1,
        y = PC2,
        fill = .data[[colour_col]],
        size = evidence_score,
        shape = plot_domain
      ),
      colour = "#4C4C4C",
      stroke = 0.65,
      alpha = 1,
      inherit.aes = FALSE
    ) +
    
    ## anchor：固定样式，不受 colour_col 影响
    geom_point(
      data = score_anchor,
      aes(x = PC1, y = PC2),
      shape = 21,
      fill = "#D73A31",
      colour = "#8A1C16",
      size = 4.3,
      stroke = 0.80,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    
    scale_fill_gradientn(
      colours = potential_palette,
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      oob = scales::squish,
      name = colour_name,
      guide = "colourbar"
    ) +
    
    scale_size_continuous(
      name = "Evidence\nscore",
      limits = c(0.50, 1.00),
      breaks = c(0.50, 1.00),
      range = c(4, 6),
      guide = guide_legend(
        override.aes = list(
          shape = 21,
          fill = "white",
          colour = "grey30",
          stroke = 0.8
        )
      )
    ) +
    
    scale_shape_manual(
      name = "Domain",
      values = domain_shapes,
      drop = FALSE,
      guide = guide_legend(
        override.aes = list(
          size = 3,
          fill = "white",
          colour = "grey30",
          stroke = 0.9
        )
      )
    ) +
    
    labs(
      title = title,
      subtitle = subtitle,
      x = paste0("PC1 (", sprintf("%.1f", var_exp[1]), "%)"),
      y = paste0("PC2 (", sprintf("%.1f", var_exp[2]), "%)")
    ) +
    theme_pca_s15
}

# 8. Panel A: substrate PCA diagnostic

nm <- names(site_master_det)

feldspar_col <- pick_first_existing(
  c("feldspar_total", "feldspar", "feldspar_pct", "feldspar_percent"),
  nm
)

quartz_col <- pick_first_existing(
  c("quartz", "quartz_pct", "quartz_percent", "quartz_total"),
  nm
)

size_col <- pick_first_existing(
  c("representative_mean_um", "recommended_representative_mean_um", "current_representative_mean_um"),
  nm
)

silt_col <- pick_first_existing(
  c("median_silt", "silt", "silt_percent", "silt_pct"),
  nm
)

cec_col <- pick_first_existing(
  c("median_cec", "cec", "CEC"),
  nm
)

required_substrate_raw <- c(feldspar_col, quartz_col, size_col, silt_col, cec_col)

if (any(is.na(required_substrate_raw))) {
  
  missing_names <- c(
    Feldspar = feldspar_col,
    Quartz = quartz_col,
    Particle_size = size_col,
    Silt = silt_col,
    CEC = cec_col
  )
  missing_names <- names(missing_names)[is.na(missing_names)]
  
  p_s15a <- make_text_panel(
    paste0(
      "Substrate PCA could not be calculated.\nMissing raw variables: ",
      paste(missing_names, collapse = ", ")
    )
  )
  
} else {
  
  substrate_pca_df <- site_master_det %>%
    filter(site %in% site_order_final) %>%
    transmute(
      site,
      feldspar_raw = safe_num(.data[[feldspar_col]]),
      quartz_raw   = safe_num(.data[[quartz_col]]),
      size_um      = safe_num(.data[[size_col]]),
      silt_raw     = safe_num(.data[[silt_col]]),
      cec_raw      = safe_num(.data[[cec_col]]),
      substrate_axis = if ("substrate_axis" %in% nm) safe_num(substrate_axis) else NA_real_
    ) %>%
    left_join(site_meta, by = "site") %>%
    mutate(
      # These are diagnostic variables for PCA only.
      # Direction is standardized so higher values indicate higher presumed
      # weatherability / accessibility.
      feldspar_score_diag = rescale_01(feldspar_raw),
      inverse_quartz_diag = rescale_01(-quartz_raw),
      particle_accessibility_diag = rescale_01(1 / pmax(size_um, 1e-6)),
      silt_score_diag = rescale_01(silt_raw),
      cec_score_diag = rescale_01(cec_raw),
      
      substrate_axis_diag = ifelse(
        is.na(substrate_axis),
        geo_mean_vec(
          feldspar_score_diag,
          inverse_quartz_diag,
          particle_accessibility_diag,
          silt_score_diag,
          cec_score_diag
        ),
        substrate_axis
      )
    )
  
  p_s15a <- build_pca_biplot(
    data_wide = substrate_pca_df,
    variable_cols = c(
      "feldspar_score_diag",
      "inverse_quartz_diag",
      "particle_accessibility_diag",
      "silt_score_diag",
      "cec_score_diag"
    ),
    variable_labels = c(
      feldspar_score_diag = "Feldspar",
      inverse_quartz_diag = "Inverse quartz",
      particle_accessibility_diag = "Particle\naccessibility",
      silt_score_diag = "Silt",
      cec_score_diag = "CEC"
    ),
    colour_col = "substrate_axis_diag",
    colour_name = "Substrate\nscore",
    title = "Substrate PCA diagnostic",
    subtitle = "PCA of V18 substrate predictors showing separation among candidate sites."
  )
} 

p <- p_s15a + theme(axis.line = element_line(colour = "#333333", linewidth = 0.45))

ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "s1.pdf"),
  plot = p,
  width = 6,
  height = 5,
  bg = "white"
)

# 9. Climate scores recalculated for Panels B-C

required_climate_raw <- c(
  "tmax_warmest_month",
  "aridity_p_over_pet",
  "pet_for_model"
)

missing_climate_raw <- setdiff(required_climate_raw, names(site_master_det))

if (length(missing_climate_raw) > 0) {
  climate_diag_df <- NULL
} else {
  
  # Prefer V18 anchor objects if available.
  # Otherwise fall back to the anchor site's values.
  if (!exists("alxa_tmax_ref")) {
    alxa_tmax_ref <- site_master_det %>%
      filter(site == anchor_site) %>%
      pull(tmax_warmest_month) %>%
      .[1] %>%
      as.numeric()
  }
  
  if (!exists("alxa_ai_ref")) {
    alxa_ai_ref <- site_master_det %>%
      filter(site == anchor_site) %>%
      pull(aridity_p_over_pet) %>%
      .[1] %>%
      as.numeric()
  }
  
  if (!exists("alxa_pet_ref")) {
    alxa_pet_ref <- site_master_det %>%
      filter(site == anchor_site) %>%
      pull(pet_for_model) %>%
      .[1] %>%
      as.numeric()
  }
  
  climate_diag_df <- site_master_det %>%
    filter(site %in% site_order_final) %>%
    transmute(
      site,
      tmax_warmest_month = safe_num(tmax_warmest_month),
      aridity_p_over_pet = safe_num(aridity_p_over_pet),
      pet_for_model = safe_num(pet_for_model)
    ) %>%
    left_join(site_meta, by = "site") %>%
    mutate(
      heat_constraint_score = anchor_upper_constraint(
        tmax_warmest_month,
        alxa_tmax_ref
      ),
      aridity_constraint_score = anchor_lower_constraint(
        aridity_p_over_pet,
        alxa_ai_ref
      ),
      pet_constraint_score = anchor_upper_constraint(
        pet_for_model,
        alxa_pet_ref
      ),
      climate_constraint_score = geo_mean_vec(
        heat_constraint_score,
        aridity_constraint_score,
        pet_constraint_score
      ),
      
      # Diagnostic PCA variables. Higher values = stronger climate stress.
      heat_stress_diag = rescale_01(tmax_warmest_month),
      aridity_stress_diag = rescale_01(1 / pmax(aridity_p_over_pet, 1e-6)),
      pet_stress_diag = rescale_01(pet_for_model)
    )
}

# 10. Panel B: climate PCA diagnostic

if (is.null(climate_diag_df)) {
  
  p_s15b <- make_text_panel(
    paste0(
      "Climate PCA could not be calculated.\nMissing variables: ",
      paste(missing_climate_raw, collapse = ", ")
    )
  )
  
} else {
  
  p_s15b <- build_pca_biplot(
    data_wide = climate_diag_df,
    variable_cols = c(
      "heat_stress_diag",
      "aridity_stress_diag",
      "pet_stress_diag"
    ),
    variable_labels = c(
      heat_stress_diag = "Heat\n(Tmax)",
      aridity_stress_diag = "Aridity\nstress",
      pet_stress_diag = "PET"
    ),
    colour_col = "climate_constraint_score",
    colour_name = "Climate\nconstraint",
    title = "Climate PCA diagnostic",
    subtitle = "PCA of heat, aridity-stress and evaporative-demand dimensions."
  )
}
p2 <- p_s15b + theme(axis.line = element_line(colour = "#333333", linewidth = 0.45))
p2
ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "s2.pdf"),
  plot = p2,
  width = 6,
  height = 5,
  bg = "white"
)

# 11. Panel C: one-sided climate constraint diagnostic

if (is.null(climate_diag_df)) {
  
  p_s15c <- make_text_panel("Climate constraint diagnostic could not be calculated.")
  
} else {
  
  climate_dot_df <- climate_diag_df %>%
    select(
      site,
      site_label,
      heat_constraint_score,
      aridity_constraint_score,
      pet_constraint_score,
      climate_constraint_score
    ) %>%
    pivot_longer(
      cols = c(
        heat_constraint_score,
        aridity_constraint_score,
        pet_constraint_score,
        climate_constraint_score
      ),
      names_to = "constraint",
      values_to = "score"
    ) %>%
    mutate(
      constraint_label = case_when(
        constraint == "heat_constraint_score" ~ "Heat (Tmax)",
        constraint == "aridity_constraint_score" ~ "Aridity (AI)",
        constraint == "pet_constraint_score" ~ "PET",
        constraint == "climate_constraint_score" ~ "Final climate constraint",
        TRUE ~ constraint
      ),
      constraint_label = factor(
        constraint_label,
        levels = c(
          "Heat (Tmax)",
          "Aridity (AI)",
          "PET",
          "Final climate constraint"
        )
      ),
      site_label = factor(site_label, levels = rev(site_label_final))
    )
  
  climate_constraint_palette <- c(
    "Heat (Tmax)" = "#123F4B",
    "Aridity (AI)" = "#C39A37",
    "PET" = "#5F987B",
    "Final climate constraint" = "#333333"
  )
  
  p_s15c <- ggplot(
    climate_dot_df,
    aes(x = score, y = site_label, colour = constraint_label)
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "#B23A2F",
      linewidth = 0.45
    ) +
    geom_point(
      size = 2.3,
      alpha = 0.95,
      position = position_dodge(width = 0.55)
    ) +
    scale_colour_manual(
      values = climate_constraint_palette,
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(0, 1.04),
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      title = "One-sided climate constraint scores",
      subtitle = "Higher values indicate weaker limitation relative to the Tengger experimental anchor.",
      x = "Constraint score",
      y = NULL
    ) +
    theme_dot_s15
}

p_s15c 

# 11. Panel C: one-sided climate constraint diagnostic
#    Alternative visualization: faceted horizontal lollipop plot

if (is.null(climate_diag_df)) {
  
  p_s15c <- make_text_panel("Climate constraint diagnostic could not be calculated.")
  
} else {
  
  climate_dot_df <- climate_diag_df %>%
    select(
      site,
      site_label,
      heat_constraint_score,
      aridity_constraint_score,
      pet_constraint_score,
      climate_constraint_score
    ) %>%
    pivot_longer(
      cols = c(
        heat_constraint_score,
        aridity_constraint_score,
        pet_constraint_score,
        climate_constraint_score
      ),
      names_to = "constraint",
      values_to = "score"
    ) %>%
    mutate(
      constraint_label = case_when(
        constraint == "heat_constraint_score" ~ "Heat (Tmax)",
        constraint == "aridity_constraint_score" ~ "Aridity (AI)",
        constraint == "pet_constraint_score" ~ "PET",
        constraint == "climate_constraint_score" ~ "Final climate\nconstraint",
        TRUE ~ constraint
      ),
      constraint_label = factor(
        constraint_label,
        levels = c(
          "Heat (Tmax)",
          "Aridity (AI)",
          "PET",
          "Final climate\nconstraint"
        )
      ),
      site_label = factor(site_label, levels = rev(site_label_final))
    )
  
  climate_constraint_palette <- c(
    "Heat (Tmax)" = "#123F4B",
    "Aridity (AI)" = "#C39A37",
    "PET" = "#5F987B",
    "Final climate\nconstraint" = "#333333"
  )
  
  p_s15c <- ggplot(
    climate_dot_df,
    aes(x = score, y = site_label)
  ) +
    geom_segment(
      aes(x = 0, xend = score, y = site_label, yend = site_label),
      colour = "#CFCFCF",
      linewidth = 0.55
    ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "#B23A2F",
      linewidth = 0.45
    ) +
    geom_point(
      aes(fill = constraint_label),
      shape = 21,
      size = 2.8,
      stroke = 0.35,
      colour = "#333333"
    ) +
    facet_wrap(
      ~ constraint_label,
      nrow = 1
    ) +
    scale_fill_manual(
      values = climate_constraint_palette,
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(0, 1.04),
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      title = "One-sided climate constraint scores",
      subtitle = "Higher values indicate weaker limitation relative to the Tengger experimental anchor.",
      x = "Constraint score",
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#E8E8E8", linewidth = 0.35),
      
      axis.title = element_text(size = 12, colour = "#222222"),
      axis.text.x = element_text(size = 10, colour = "#222222"),
      axis.text.y = element_text(size = 12, colour = "#222222"),
      
      strip.background = element_rect(fill = "white", colour = NA),
      strip.text = element_text(size = 12,  colour = "#222222"),
      
      plot.title = element_text(size = 12, colour = "#222222"),
      plot.subtitle = element_text(size = 12, colour = "grey35", margin = margin(b = 4)),
      
      legend.position = "none",
      plot.margin = margin(6, 8, 6, 6)
    )
}


# remove Tmax

# Climate constraint diagnostic plot
# V18: remove Heat (Tmax), keep AI + PET + final climate constraint

if (is.null(climate_diag_df)) {
  
  p_s15c <- make_text_panel("Climate constraint diagnostic could not be calculated.")
  
} else {
  
  climate_dot_df <- climate_diag_df %>%
    select(
      site,
      site_label,
      aridity_constraint_score,
      pet_constraint_score,
      climate_constraint_score
    ) %>%
    pivot_longer(
      cols = c(
        aridity_constraint_score,
        pet_constraint_score,
        climate_constraint_score
      ),
      names_to = "constraint",
      values_to = "score"
    ) %>%
    mutate(
      constraint_label = case_when(
        constraint == "aridity_constraint_score" ~ "Aridity (AI)",
        constraint == "pet_constraint_score" ~ "PET",
        constraint == "climate_constraint_score" ~ "Final climate\nconstraint",
        TRUE ~ constraint
      ),
      constraint_label = factor(
        constraint_label,
        levels = c(
          "Aridity (AI)",
          "PET",
          "Final climate\nconstraint"
        )
      ),
      site_label = factor(site_label, levels = rev(site_label_final))
    )
  
  climate_constraint_palette <- c(
    "Aridity (AI)" = "#C39A37",
    "PET" = "#5F987B",
    "Final climate\nconstraint" = "#333333"
  )
  
  p_s15c <- ggplot(
    climate_dot_df,
    aes(x = score, y = site_label)
  ) +
    geom_segment(
      aes(x = 0, xend = score, y = site_label, yend = site_label),
      colour = "#CFCFCF",
      linewidth = 0.55
    ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "#B23A2F",
      linewidth = 0.45
    ) +
    geom_point(
      aes(fill = constraint_label),
      shape = 21,
      size = 2.8,
      stroke = 0.35,
      colour = "#333333"
    ) +
    facet_wrap(
      ~ constraint_label,
      nrow = 1
    ) +
    scale_fill_manual(
      values = climate_constraint_palette,
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(0, 1.04),
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      title = "One-sided climate constraint scores",
      subtitle = "Higher values indicate weaker limitation relative to the Tengger experimental anchor.",
      x = "Constraint score",
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#E8E8E8", linewidth = 0.35),
      
      axis.title = element_text(size = 12, colour = "#222222"),
      axis.text.x = element_text(size = 10, colour = "#222222"),
      axis.text.y = element_text(size = 12, colour = "#222222"),
      
      strip.background = element_rect(fill = "white", colour = NA),
      strip.text = element_text(size = 12, colour = "#222222"),
      
      plot.title = element_text(size = 12, colour = "#222222"),
      plot.subtitle = element_text(size = 12, colour = "grey35", margin = margin(b = 4)),
      
      legend.position = "none",
      plot.margin = margin(6, 8, 6, 6)
    )
}

ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "s3.pdf"),
  plot = p_s15c,
  width = 8,
  height = 5,
  bg = "white"
)
p_s15c

# 12. Panel D: microbial function response curves
#     V18:
#     - microbial suitability uses CFU-normalized CAS function only
#     - one fitted response line is shown for each condition
#     - raw experimental points are shown

# Helper for S15D microbial response curves

# 12. Panel D: microbial function response curves
# Function-only V18 version

# This panel visualizes the microbial suitability term used in the model.
# Microbial suitability is based only on CFU-normalized CAS function score.
# No growth score or population-level function is shown.

# 12.1 Locate response-model object

if (exists("response_models_v13")) {
  response_models_use <- response_models_v13
} else if (exists("response_models")) {
  response_models_use <- response_models
} else {
  stop("No response model object found. Expected `response_models_v13` or `response_models`.")
}

# 12.2 Helper functions

get_response_axis <- function(axis_candidates) {
  
  hit <- intersect(axis_candidates, names(response_models_use))
  
  if (length(hit) == 0) {
    stop(
      "None of the candidate response axes were found in the response model object: ",
      paste(axis_candidates, collapse = ", "),
      "\nAvailable axes are: ",
      paste(names(response_models_use), collapse = ", ")
    )
  }
  
  hit[1]
}


get_function_model <- function(axis_name) {
  
  model_list <- response_models_use[[axis_name]]
  
  possible_function_names <- c(
    "function_score",
    "function",
    "func_per_cell",
    "per_cell_function_score"
  )
  
  hit <- intersect(possible_function_names, names(model_list))
  
  if (length(hit) == 0) {
    stop(
      "No function model found for axis `", axis_name, "`. ",
      "Expected one of: ",
      paste(possible_function_names, collapse = ", "),
      ". Available model names are: ",
      paste(names(model_list), collapse = ", ")
    )
  }
  
  model_list[[hit[1]]]
}

predict_curve_safe <- function(model_obj, new_x) {
  
  if (exists("predict_response")) {
    return(predict_response(model_obj, new_x, clamp = TRUE))
  }
  
  original_x <- as.numeric(new_x)
  pred_x <- original_x
  
  if (!is.null(model_obj$x_min) && !is.null(model_obj$x_max)) {
    pred_x <- pmin(pmax(pred_x, model_obj$x_min), model_obj$x_max)
  }
  
  new_df <- data.frame(condition_value = pred_x)
  pred <- predict(model_obj$model, newdata = new_df)
  pred <- as.numeric(pred)
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  pred
}

make_response_curve_data <- function(
    axis_candidates,
    x_range,
    condition_label
) {
  
  axis_name <- get_response_axis(axis_candidates)
  model_function <- get_function_model(axis_name)
  
  x_grid <- seq(x_range[1], x_range[2], length.out = 250)
  
  tibble(
    condition = condition_label,
    x = x_grid,
    response_score = predict_curve_safe(model_function, x_grid),
    response_type = "CFU-normalized CAS function"
  ) %>%
    filter(!is.na(response_score))
}

# 12.3 Fitted function-response curves

curve_ph <- make_response_curve_data(
  axis_candidates = c("pH", "ph"),
  x_range = c(6, 10),
  condition_label = "pH"
)

curve_temp <- make_response_curve_data(
  axis_candidates = c("temperature", "Temperature", "temp"),
  x_range = c(15, 45),
  condition_label = "Temperature"
)

curve_df_s15d <- bind_rows(curve_ph, curve_temp) %>%
  mutate(
    condition = factor(condition, levels = c("pH", "Temperature")),
    response_type = factor(
      response_type,
      levels = "CFU-normalized CAS function"
    )
  )

# 12.4 Raw experimental points

if (exists("microbe_function_df")) {
  microbe_raw_use <- microbe_function_df
} else if (exists("microbe_function_df_v13")) {
  microbe_raw_use <- microbe_function_df_v13
} else if (exists("microbe_df_v13")) {
  microbe_raw_use <- microbe_df_v13
} else if (exists("microbe_df")) {
  microbe_raw_use <- microbe_df
} else {
  stop("No microbial raw data object found. Expected `microbe_function_df`, `microbe_function_df_v13`, `microbe_df_v13`, or `microbe_df`.")
}

point_df_s15d <- microbe_raw_use %>%
  mutate(
    condition_type_std = case_when(
      "condition_type_std" %in% names(.) ~ as.character(condition_type_std),
      str_detect(str_to_lower(as.character(condition_type)), "ph") ~ "pH",
      str_detect(str_to_lower(as.character(condition_type)), "temp") ~ "temperature",
      TRUE ~ as.character(condition_type)
    ),
    condition = case_when(
      condition_type_std %in% c("pH", "ph") ~ "pH",
      condition_type_std %in% c("temperature", "Temperature", "temp") ~ "Temperature",
      TRUE ~ NA_character_
    ),
    x = as.numeric(condition_value),
    response_score = as.numeric(function_score),
    response_type = "CFU-normalized CAS function"
  ) %>%
  filter(
    condition %in% c("pH", "Temperature"),
    !is.na(x),
    !is.na(response_score)
  ) %>%
  mutate(
    condition = factor(condition, levels = c("pH", "Temperature")),
    response_type = factor(
      response_type,
      levels = "CFU-normalized CAS function"
    )
  )

# 12.5 Theme

theme_s15d <- theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    
    panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    
    strip.text = element_text(size = 11, face = "bold", colour = "#222222"),
    strip.background = element_blank(),
    
    axis.title = element_text(size = 10, colour = "#222222"),
    axis.text  = element_text(size = 9, colour = "#333333"),
    
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key = element_rect(fill = "white", colour = NA),
    
    plot.margin = margin(5, 8, 5, 5)
  )

# 12.6 Build Panel D

p_s15d <- ggplot() +
  
  geom_point(
    data = point_df_s15d,
    aes(
      x = x,
      y = response_score
    ),
    shape = 21,
    fill = "#5F987B",
    colour = "#2F4F4F",
    size = 2.2,
    stroke = 0.45,
    alpha = 0.80
  ) +
  
  geom_line(
    data = curve_df_s15d,
    aes(
      x = x,
      y = response_score
    ),
    colour = "#2D7A78",
    linewidth = 0.85
  ) +
  
  facet_wrap(
    ~ condition,
    ncol = 1,
    scales = "free_x"
  ) +
  
  scale_y_continuous(
    limits = c(0, 1.02),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  
  labs(
    x = "Condition value",
    y = "Function score"
  ) +
  
  theme_s15d

print(p_s15d)




p_s15d

ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "s4.pdf"),
  plot = p_s15d,
  width = 4,
  height = 7,
  bg = "white"
)


#### figuer s16 ####

# Fig. S16. Monte Carlo uncertainty and rank stability, V18
#
# Panels:
#   A. Monte Carlo relative-potential distributions
#   B. Median, IQR and 95% interval
#   C. Rank-stability heatmap
#   D. Probability summary
#
# Required objects:
#   mc_results_v18
#   mc_main_summary_v18

# 1. Packages


# S16A. Monte Carlo relative-potential distributions
# Violin + boxplot + MC median point
#
# Domain shape:
#   Sandy desert     = circle
#   Gobi transition  = square
#
# Required objects:
#   site_master_det
#   mc_results_v18 and mc_main_summary_v18
#   OR mc_results_v15 and mc_main_summary_v15

# S16A. Monte Carlo relative-potential distributions
# Violin + boxplot + MC median point
# Point border fixed to match previous figures
#   Sandy desert    = circle
#   Gobi transition = square

library(dplyr)
library(ggplot2)
library(stringr)
library(scales)
library(grid)

# 0. Object-name compatibility

if (!exists("mc_results_v18") && exists("mc_results_v15")) {
  mc_results_v18 <- mc_results_v15
}
if (!exists("mc_main_summary_v18") && exists("mc_main_summary_v15")) {
  mc_main_summary_v18 <- mc_main_summary_v15
}

if (!exists("mc_results_v18")) stop("Object `mc_results_v18` not found.")
if (!exists("mc_main_summary_v18")) stop("Object `mc_main_summary_v18` not found.")
if (!exists("site_master_det")) stop("Object `site_master_det` not found.")

# 1. Lookup tables

domain_lookup_s16 <- site_master_det %>%
  dplyr::select(site, dplyr::any_of("plot_domain")) %>%
  dplyr::distinct()

if (!"plot_domain" %in% names(domain_lookup_s16)) {
  domain_lookup_s16 <- site_master_det %>%
    dplyr::select(site) %>%
    dplyr::distinct() %>%
    dplyr::mutate(plot_domain = NA_character_)
}

label_lookup_s16 <- site_master_det %>%
  dplyr::select(site, dplyr::any_of(c("site_label", "site_label_clean"))) %>%
  dplyr::distinct()

if (!"site_label" %in% names(label_lookup_s16)) {
  label_lookup_s16$site_label <- NA_character_
}
if (!"site_label_clean" %in% names(label_lookup_s16)) {
  label_lookup_s16$site_label_clean <- NA_character_
}

label_lookup_s16 <- label_lookup_s16 %>%
  dplyr::mutate(
    site_label_clean = dplyr::case_when(
      site == "Hobq" ~ "Kubuqi",
      site == "Sonoran_Altar" ~ "Gran Desierto",
      site == "Argentina_MedanosGrandes_NE" ~ "Medanos Grandes",
      site == "Argentina_MedanosGra" ~ "Medanos Grandes",
      site == "Kalahari_Bagatelle" ~ "Kalahari",
      site == "Karakum_AshkAbad" ~ "Karakum",
      site == "Thal_Mankera" ~ "Thal",
      site == "Alxa_Left_Banner_project_base_measured" ~ "Tengger experimental anchor",
      !is.na(site_label_clean) ~ as.character(site_label_clean),
      !is.na(site_label) ~ as.character(site_label),
      TRUE ~ stringr::str_replace_all(site, "_", " ")
    )
  ) %>%
  dplyr::select(site, site_label_clean) %>%
  dplyr::distinct()

# 2. Build MC plot tables

mc_plot_s16 <- mc_results_v18 %>%
  dplyr::filter(
    include_in_main_model,
    site != "Alxa_Left_Banner_project_base_measured",
    site != "Alxa_Left_Banner_project_base_downloaded",
    !is.na(relative_to_alxa)
  ) %>%
  dplyr::select(-dplyr::any_of(c("plot_domain", "site_label_clean"))) %>%
  dplyr::left_join(domain_lookup_s16, by = "site") %>%
  dplyr::left_join(label_lookup_s16, by = "site") %>%
  dplyr::mutate(
    plot_domain = dplyr::case_when(
      site %in% c("Badain_Jaran", "Taklamakan") ~ "Gobi transition",
      TRUE ~ "Sandy desert"
    ),
    plot_domain = factor(
      plot_domain,
      levels = c("Sandy desert", "Gobi transition")
    )
  )

summary_s16 <- mc_main_summary_v18 %>%
  dplyr::filter(
    include_in_main_model,
    site != "Alxa_Left_Banner_project_base_measured",
    site != "Alxa_Left_Banner_project_base_downloaded",
    !is.na(relative_median)
  ) %>%
  dplyr::select(-dplyr::any_of(c("plot_domain", "site_label_clean"))) %>%
  dplyr::left_join(domain_lookup_s16, by = "site") %>%
  dplyr::left_join(label_lookup_s16, by = "site") %>%
  dplyr::mutate(
    plot_domain = dplyr::case_when(
      site %in% c("Badain_Jaran", "Taklamakan") ~ "Gobi transition",
      TRUE ~ "Sandy desert"
    ),
    plot_domain = factor(
      plot_domain,
      levels = c("Sandy desert", "Gobi transition")
    ),
    evidence_score = dplyr::case_when(
      is.na(evidence_score) ~ 0.50,
      evidence_score >= 0.75 ~ 1.00,
      TRUE ~ 0.50
    )
  )

# 3. Order sites by MC median

site_order_s16 <- summary_s16 %>%
  dplyr::arrange(relative_median) %>%
  dplyr::pull(site_label_clean)

mc_plot_s16 <- mc_plot_s16 %>%
  dplyr::mutate(
    site_label = factor(site_label_clean, levels = site_order_s16)
  )

summary_s16 <- summary_s16 %>%
  dplyr::mutate(
    site_label = factor(site_label_clean, levels = site_order_s16)
  )

# 4. Axis range and breaks

relative_xlim <- range(mc_plot_s16$relative_to_alxa, na.rm = TRUE)
relative_xlim <- c(
  min(0.25, relative_xlim[1]),
  max(2.00, relative_xlim[2])
)

relative_breaks <- c(0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0)
relative_breaks <- relative_breaks[
  relative_breaks >= relative_xlim[1] &
    relative_breaks <= relative_xlim[2]
]

# 5. Palettes and shapes

if (!exists("potential_palette")) {
  potential_palette <- c(
    "#A77A18",
    "#C39A37",
    "#8FAA68",
    "#5F987B",
    "#2D7A78",
    "#123F4B"
  )
}

domain_shapes_s16 <- c(
  "Sandy desert" = 21,
  "Gobi transition" = 22
)

# 6. Theme

theme_s16_base <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#E2E2E2", linewidth = 0.35),
    
    axis.text.y = element_text(size = 12, colour = "#222222"),
    axis.text.x = element_text(size = 11, colour = "#222222"),
    axis.title.x = element_text(size = 13, colour = "#222222", margin = margin(t = 8)),
    
    plot.title = element_text(size = 16, face = "bold", colour = "#222222"),
    plot.subtitle = element_text(size = 11.5, colour = "#555555", margin = margin(b = 10)),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.spacing.x = unit(3.0, "mm"),
    legend.key = element_rect(fill = "white", colour = NA),
    
    plot.margin = margin(8, 16, 8, 8)
  )

# 7. Plot

p_s16a <- ggplot(
  mc_plot_s16,
  aes(x = relative_to_alxa, y = site_label)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "#B23A2F",
    linewidth = 0.50
  ) +
  geom_violin(
    fill = "#F1F1F1",
    colour = "#A7A7A7",
    linewidth = 0.35,
    width = 0.82,
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = "white",
    colour = "#666666",
    linewidth = 0.38
  ) +
  geom_point(
    data = summary_s16,
    aes(
      x = relative_median,
      y = site_label,
      fill = relative_median,
      size = evidence_score,
      shape = plot_domain
    ),
    inherit.aes = FALSE,
    colour = "#4C4C4C",
    stroke = 0.65
  ) +
  scale_fill_gradientn(
    colours = potential_palette,
    limits = range(summary_s16$relative_median, na.rm = TRUE),
    oob = scales::squish,
    name = "MC median"
  ) +
  scale_shape_manual(
    name = "Domain",
    values = domain_shapes_s16,
    breaks = c("Sandy desert", "Gobi transition"),
    labels = c("Sandy desert", "Gobi transition"),
    drop = FALSE,
    guide = guide_legend(
      title.position = "top",
      nrow = 1,
      override.aes = list(
        fill = "white",
        colour = "#4C4C4C",
        size = 3.6,
        stroke = 0.65
      )
    )
  ) +
  scale_size_continuous(
    name = "Evidence\nscore",
    limits = c(0.50, 1.00),
    breaks = c(0.50, 1.00),
    labels = c("0.50", "1.00"),
    range = c(2.0, 4),
    guide = guide_legend(
      title.position = "top",
      nrow = 1,
      override.aes = list(
        shape = 21,
        fill = "white",
        colour = "#4C4C4C",
        stroke = 0.65
      )
    )
  ) +
  scale_x_continuous(
    limits = relative_xlim,
    breaks = relative_breaks,
    labels = paste0(relative_breaks, "x"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "Monte Carlo relative-potential distributions",
    subtitle = "Each distribution represents 5000 Monte Carlo estimates relative to the Tengger experimental anchor.",
    x = "Relative potential to Tengger anchor",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_s16_base

print(p_s16a)

ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "s16a.pdf"),
  plot = p_s16a,
  width = 6,
  height = 5,
  bg = "white"
)

# 7. Panel B: Median, IQR and 95% interval

# S16B. Median, IQR and 95% interval
# Consistent with Panel A
# Point border fixed; no stability mapping

library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(stringr)

# 0. Check required objects

if (!exists("summary_s16")) {
  stop("Object `summary_s16` not found. Please run the Panel A preparation code first.")
}

if (!exists("potential_palette")) {
  potential_palette <- c(
    "#A77A18",
    "#C39A37",
    "#8FAA68",
    "#5F987B",
    "#2D7A78",
    "#123F4B"
  )
}

domain_shapes_s16 <- c(
  "Sandy desert" = 21,
  "Gobi transition" = 22
)

# 1. Standardize interval columns

get_col_or_na <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    return(rep(NA_real_, nrow(df)))
  }
  as.numeric(df[[hit[1]]])
}

summary_s16_b <- summary_s16 %>%
  dplyr::mutate(
    relative_iqr_low_std = get_col_or_na(
      .,
      c("relative_iqr_low", "relative_q25", "q25")
    ),
    relative_iqr_high_std = get_col_or_na(
      .,
      c("relative_iqr_high", "relative_q75", "q75")
    ),
    relative_ci95_low_std = get_col_or_na(
      .,
      c("relative_ci95_low", "relative_q025", "q025", "relative_2_5")
    ),
    relative_ci95_high_std = get_col_or_na(
      .,
      c("relative_ci95_high", "relative_q975", "q975", "relative_97_5")
    ),
    evidence_score = dplyr::case_when(
      is.na(evidence_score) ~ 0.50,
      evidence_score >= 0.75 ~ 1.00,
      TRUE ~ 0.50
    ),
    plot_domain = dplyr::case_when(
      site %in% c("Badain_Jaran", "Taklamakan") ~ "Gobi transition",
      TRUE ~ "Sandy desert"
    ),
    plot_domain = factor(
      plot_domain,
      levels = c("Sandy desert", "Gobi transition")
    )
  )

# 2. Make sure site labels exist and match Panel A

if (!"site_label" %in% names(summary_s16_b)) {
  if ("site_label_clean" %in% names(summary_s16_b)) {
    summary_s16_b <- summary_s16_b %>%
      dplyr::mutate(site_label = site_label_clean)
  } else {
    summary_s16_b <- summary_s16_b %>%
      dplyr::mutate(
        site_label = dplyr::case_when(
          site == "Hobq" ~ "Kubuqi",
          site == "Sonoran_Altar" ~ "Gran Desierto",
          site == "Argentina_MedanosGrandes_NE" ~ "Medanos Grandes",
          site == "Argentina_MedanosGra" ~ "Medanos Grandes",
          site == "Kalahari_Bagatelle" ~ "Kalahari",
          site == "Karakum_AshkAbad" ~ "Karakum",
          site == "Thal_Mankera" ~ "Thal",
          TRUE ~ stringr::str_replace_all(site, "_", " ")
        )
      )
  }
}

if (!inherits(summary_s16_b$site_label, "factor")) {
  site_order_s16_b <- summary_s16_b %>%
    dplyr::arrange(relative_median) %>%
    dplyr::pull(site_label)
  
  summary_s16_b <- summary_s16_b %>%
    dplyr::mutate(
      site_label = factor(site_label, levels = site_order_s16_b)
    )
}

# 3. Axis range and breaks

relative_xlim_b <- range(
  c(
    summary_s16_b$relative_ci95_low_std,
    summary_s16_b$relative_ci95_high_std,
    summary_s16_b$relative_median
  ),
  na.rm = TRUE
)

relative_xlim_b <- c(
  min(0.25, relative_xlim_b[1]),
  max(2.00, relative_xlim_b[2])
)

relative_breaks_b <- c(0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0)
relative_breaks_b <- relative_breaks_b[
  relative_breaks_b >= relative_xlim_b[1] &
    relative_breaks_b <= relative_xlim_b[2]
]

# 4. Theme

if (!exists("theme_s16_base")) {
  theme_s16_base <- theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#E2E2E2", linewidth = 0.35),
      axis.text.y = element_text(size = 12, colour = "#222222"),
      axis.text.x = element_text(size = 11, colour = "#222222"),
      axis.title.x = element_text(size = 13, colour = "#222222", margin = margin(t = 8)),
      plot.title = element_text(size = 16, face = "bold", colour = "#222222"),
      plot.subtitle = element_text(size = 11.5, colour = "#555555", margin = margin(b = 10)),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.spacing.x = unit(3.0, "mm"),
      legend.key = element_rect(fill = "white", colour = NA),
      plot.margin = margin(8, 16, 8, 8)
    )
}

# 5. Panel B plot

p_s16b <- ggplot(
  summary_s16_b,
  aes(y = site_label)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "#B23A2F",
    linewidth = 0.50
  ) +
 
  geom_point(
    aes(
      x = relative_median,
      fill = relative_median,
      size = evidence_score,
      shape = plot_domain
    ),
    colour = "#4C4C4C",
    stroke = 0.65
  ) +
  scale_fill_gradientn(
    colours = potential_palette,
    limits = range(summary_s16_b$relative_median, na.rm = TRUE),
    oob = scales::squish,
    name = "MC median"
  ) +
  geom_linerange(
    aes(
      xmin = relative_ci95_low_std,
      xmax = relative_ci95_high_std
    ),
    linewidth = 0.45,
    colour = "#8C8C8C"
  ) +
  geom_linerange(
    aes(
      xmin = relative_iqr_low_std,
      xmax = relative_iqr_high_std
    ),
    linewidth = 1.45,
    colour = "#4F4F4F"
  ) +
  scale_shape_manual(
    name = "Domain",
    values = domain_shapes_s16,
    breaks = c("Sandy desert", "Gobi transition"),
    labels = c("Sandy desert", "Gobi transition"),
    drop = FALSE,
    guide = guide_legend(
      title.position = "top",
      nrow = 1,
      override.aes = list(
        fill = "white",
        colour = "#4C4C4C",
        size = 3.6,
        stroke = 0.65
      )
    )
  ) +
  scale_size_continuous(
    name = "Evidence\nscore",
    limits = c(0.50, 1.00),
    breaks = c(0.50, 1.00),
    labels = c("0.50", "1.00"),
    range = c(2.0, 4),
    guide = guide_legend(
      title.position = "top",
      nrow = 1,
      override.aes = list(
        shape = 21,
        fill = "white",
        colour = "#4C4C4C",
        stroke = 0.65
      )
    )
  ) +
  scale_x_continuous(
    limits = relative_xlim_b,
    breaks = relative_breaks_b,
    labels = paste0(relative_breaks_b, "x"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "Monte Carlo relative-potential summary",
    subtitle = "Points show medians; thick bars indicate IQRs; thin bars indicate 95% intervals.",
    x = "Relative potential to Tengger anchor",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_s16_base

print(p_s16b)


p_s16b
ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "s16b.pdf"),
  plot = p_s16b,
  width = 6,
  height = 5,
  bg = "white"
)

# 8. Panel C: Rank-stability heatmap

rank_prob_s16 <- mc_results_v18 %>%
  filter(
    site %in% site_order,
    use_in_main_barplot,
    !is.na(rank_main)
  ) %>%
  mutate(
    rank_bin = as.integer(round(rank_main)),
    site_label = case_when(
      site %in% names(site_label_map) ~ unname(site_label_map[site]),
      TRUE ~ str_replace_all(site, "_", " ")
    )
  ) %>%
  count(site, site_label, rank_bin, name = "n") %>%
  group_by(site) %>%
  mutate(rank_probability = n / sum(n)) %>%
  ungroup() %>%
  complete(
    site = site_order,
    rank_bin = seq_len(length(site_order)),
    fill = list(n = 0, rank_probability = 0)
  ) %>%
  mutate(
    site_label = case_when(
      site %in% names(site_label_map) ~ unname(site_label_map[site]),
      TRUE ~ str_replace_all(site, "_", " ")
    ),
    site_label = factor(site_label, levels = rev(site_label_order))
  )

prob_palette <- c(
  "#F2E8C9",
  "#C39A37",
  "#8FAE7A",
  "#4F8F7A",
  "#123F4B"
)

p_s16c <- ggplot(
  rank_prob_s16,
  aes(x = rank_bin, y = site_label, fill = rank_probability)
) +
  geom_tile(
    colour = "white",
    linewidth = 0.35,
    width = 0.92,
    height = 0.92
  ) +
  scale_fill_gradientn(
    colours = prob_palette,
    limits = c(0, max(rank_prob_s16$rank_probability, na.rm = TRUE)),
    oob = scales::squish,
    name = "Rank\nprobability"
  ) +
  scale_x_continuous(
    breaks = seq_len(length(site_order)),
    labels = seq_len(length(site_order)),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "Rank-stability heatmap",
    subtitle = "Colour indicates the probability of occupying each rank across Monte Carlo iterations.",
    x = "Monte Carlo rank",
    y = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    
    axis.title = element_text(size = 14, colour = "#222222"),
    axis.text = element_text(size = 14, colour = "#222222"),
    axis.ticks = element_blank(),
    
    plot.title = element_text(size = 13, face = "bold", colour = "#222222"),
    plot.subtitle = element_text(size = 9, colour = "grey35", margin = margin(b = 4)),
    
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.key.height = unit(3.5, "mm"),
    legend.key.width = unit(8, "mm"),
    
    plot.margin = margin(6, 8, 6, 6)
  )

p_s16c

ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "s15g.pdf"),
  plot = p_s16c,
  width = 7,
  height = 5,
  bg = "white"
)

# 9. Panel D: Probability summary

prob_summary_s16 <- summary_s16 %>%
  as.data.frame() %>%
  select(
    site,
    site_label,
    above_alxa_probability,
    top3_probability,
    top5_probability
  ) %>%
  pivot_longer(
    cols = c(
      above_alxa_probability,
      top3_probability,
      top5_probability
    ),
    names_to = "metric",
    values_to = "probability"
  ) %>%
  mutate(
    metric_label = case_when(
      metric == "above_alxa_probability" ~ "P(> anchor)",
      metric == "top3_probability" ~ "P(top 3)",
      metric == "top5_probability" ~ "P(top 5)",
      TRUE ~ metric
    ),
    metric_label = factor(
      metric_label,
      levels = c("P(> anchor)", "P(top 3)", "P(top 5)")
    ),
    site_label = factor(as.character(site_label), levels = rev(site_label_order)),
    probability_label = sprintf("%.2f", probability)
  )

# Rank-probability heatmap palette
prob_palette <- c(
  "#F4F1E6",
  "#D9C58A",
  "#B88A1D",
  "#5F987B",
  "#123F4B"
)

p_s16d <- ggplot(
  prob_summary_s16,
  aes(x = metric_label, y = site_label, fill = probability)
) +
  geom_tile(
    colour = "white",
    linewidth = 0.45,
    width = 0.92,
    height = 0.92
  ) +
  geom_text(
    aes(label = probability_label),
    size = 4.0,
    colour = "#222222"
  ) +
  scale_fill_gradientn(
    colours = prob_palette,
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "0.25", "0.50", "0.75", "1.00"),
    oob = scales::squish,
    name = "Probability"
  ) +
  labs(
    title = "Probability summary",
    subtitle = "Stability metrics derived from Monte Carlo ranks and relative potential distributions.",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      size = 14,
      colour = "#222222",
      angle = 0,
      hjust = 0.5,
      vjust = 1
    ),
    axis.text.y = element_text(size = 14, colour = "#222222"),
    axis.ticks = element_blank(),
    
    plot.title = element_text(size = 13, face = "bold", colour = "#222222"),
    plot.subtitle = element_text(size = 9, colour = "grey35", margin = margin(b = 4)),
    
    legend.position = "bottom",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.key.height = unit(3.5, "mm"),
    legend.key.width = unit(8, "mm"),
    
    plot.margin = margin(6, 8, 6, 6)
  )

p_s16d

ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "s16c.pdf"),
  plot = p_s16d,
  width = 6,
  height = 5,
  bg = "white"
)

# FINAL SECTION. V18 rank intervals and table-only sensitivity analyses

# This section is the authoritative sensitivity-analysis output for V18.
# It recalculates rank 95% simulation intervals and writes all sensitivity
# analyses as tables. No sensitivity visualizations are generated here.

# F1. Recalculate Monte Carlo summaries with rank 95% intervals

mc_summary_v18 <- mc_results_v18 %>%
  group_by(
    site,
    site_label_clean,
    include_in_main_model,
    use_in_main_barplot,
    confidence_tier,
    evidence_score
  ) %>%
  summarise(
    n_iter = n(),

    relative_median = median(relative_to_alxa, na.rm = TRUE),
    relative_iqr_low = quantile(relative_to_alxa, 0.25, na.rm = TRUE),
    relative_iqr_high = quantile(relative_to_alxa, 0.75, na.rm = TRUE),
    relative_ci95_low = quantile(relative_to_alxa, 0.025, na.rm = TRUE),
    relative_ci95_high = quantile(relative_to_alxa, 0.975, na.rm = TRUE),

    substrate_median = median(substrate_axis, na.rm = TRUE),
    substrate_iqr_low = quantile(substrate_axis, 0.25, na.rm = TRUE),
    substrate_iqr_high = quantile(substrate_axis, 0.75, na.rm = TRUE),

    mineral_median = median(mineral_reactivity_axis, na.rm = TRUE),
    mineral_iqr_low = quantile(mineral_reactivity_axis, 0.25, na.rm = TRUE),
    mineral_iqr_high = quantile(mineral_reactivity_axis, 0.75, na.rm = TRUE),

    physical_median = median(physical_accessibility_axis, na.rm = TRUE),
    physical_iqr_low = quantile(physical_accessibility_axis, 0.25, na.rm = TRUE),
    physical_iqr_high = quantile(physical_accessibility_axis, 0.75, na.rm = TRUE),

    surface_median = median(surface_exchange_axis, na.rm = TRUE),
    surface_iqr_low = quantile(surface_exchange_axis, 0.25, na.rm = TRUE),
    surface_iqr_high = quantile(surface_exchange_axis, 0.75, na.rm = TRUE),

    climate_median = median(climate_suitability, na.rm = TRUE),
    climate_iqr_low = quantile(climate_suitability, 0.25, na.rm = TRUE),
    climate_iqr_high = quantile(climate_suitability, 0.75, na.rm = TRUE),

    microbial_median = median(microbial_suitability, na.rm = TRUE),
    microbial_iqr_low = quantile(microbial_suitability, 0.25, na.rm = TRUE),
    microbial_iqr_high = quantile(microbial_suitability, 0.75, na.rm = TRUE),

    relative_fe_median = median(relative_fe_to_alxa, na.rm = TRUE),
    relative_fe_iqr_low = quantile(relative_fe_to_alxa, 0.25, na.rm = TRUE),
    relative_fe_iqr_high = quantile(relative_fe_to_alxa, 0.75, na.rm = TRUE),

    relative_texture_median = median(relative_texture_to_alxa, na.rm = TRUE),
    relative_texture_iqr_low = quantile(relative_texture_to_alxa, 0.25, na.rm = TRUE),
    relative_texture_iqr_high = quantile(relative_texture_to_alxa, 0.75, na.rm = TRUE),

    rank_median = median(rank_main, na.rm = TRUE),
    rank_iqr_low = quantile(rank_main, 0.25, na.rm = TRUE),
    rank_iqr_high = quantile(rank_main, 0.75, na.rm = TRUE),
    rank_ci95_low = quantile(rank_main, 0.025, na.rm = TRUE),
    rank_ci95_high = quantile(rank_main, 0.975, na.rm = TRUE),

    top3_probability = mean(rank_main <= 3, na.rm = TRUE),
    top5_probability = mean(rank_main <= 5, na.rm = TRUE),
    above_alxa_probability = mean(relative_to_alxa > 1, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  arrange(desc(relative_median))

write.csv(
  mc_summary_v18,
  file.path(output_dir, "v18_monte_carlo_summary_median_IQR_rank95_stability.csv"),
  row.names = FALSE
)

# Overwrite the previous summary filename too, so downstream plotting code can
# use the same object/file name while now containing rank_ci95_low/high.
write.csv(
  mc_summary_v18,
  file.path(output_dir, "v18_monte_carlo_summary_median_IQR_rank_stability.csv"),
  row.names = FALSE
)

mc_main_summary_v18 <- mc_summary_v18 %>%
  filter(use_in_main_barplot) %>%
  arrange(desc(relative_median))

write.csv(
  mc_main_summary_v18,
  file.path(output_dir, "v18_monte_carlo_main_ranking_summary.csv"),
  row.names = FALSE
)

# Optional rank-uncertainty plot data only. Plotting remains outside sensitivity.
rank_uncertainty_table_v18 <- mc_main_summary_v18 %>%
  select(
    site,
    site_label_clean,
    evidence_score,
    relative_median,
    rank_median,
    rank_iqr_low,
    rank_iqr_high,
    rank_ci95_low,
    rank_ci95_high,
    top5_probability,
    above_alxa_probability
  ) %>%
  arrange(rank_median)

write.csv(
  rank_uncertainty_table_v18,
  file.path(output_dir, "v18_rank_uncertainty_MC_median_IQR_95interval_table.csv"),
  row.names = FALSE
)

# F2. Shared helpers for deterministic sensitivity scenarios

main_rank_v18 <- deterministic_ranking_v18 %>%
  arrange(desc(relative_to_alxa_measured)) %>%
  mutate(main_rank = row_number()) %>%
  select(site, main_relative = relative_to_alxa_measured, main_rank)

get_scenario_rank_v18 <- function(model_df, scenario_type, scenario, extra_cols = tibble()) {
  out <- model_df %>%
    filter(use_in_main_barplot) %>%
    arrange(desc(relative_to_alxa_measured)) %>%
    mutate(
      scenario_type = scenario_type,
      scenario = scenario,
      scenario_rank = row_number(),
      scenario_relative = relative_to_alxa_measured
    ) %>%
    select(scenario_type, scenario, site, site_label_clean, scenario_rank, scenario_relative)

  if (nrow(extra_cols) > 0) {
    out <- bind_cols(out, extra_cols[rep(1, nrow(out)), , drop = FALSE])
  }
  out
}

summarize_scenario_vs_main_v18 <- function(rank_df, scenario_type, scenario) {
  tab <- rank_df %>%
    left_join(main_rank_v18, by = "site") %>%
    mutate(rank_shift = scenario_rank - main_rank)

  main_top5 <- main_rank_v18 %>% filter(main_rank <= 5) %>% pull(site)
  scen_top5 <- tab %>% filter(scenario_rank <= 5) %>% pull(site)

  tibble(
    scenario_type = scenario_type,
    scenario = scenario,
    n_sites = nrow(tab),
    spearman_rho_vs_main = suppressWarnings(cor(tab$main_rank, tab$scenario_rank, method = "spearman", use = "complete.obs")),
    top5_overlap = length(intersect(main_top5, scen_top5)),
    median_abs_rank_shift = median(abs(tab$rank_shift), na.rm = TRUE),
    max_abs_rank_shift = max(abs(tab$rank_shift), na.rm = TRUE),
    most_shifted_sites = paste(
      tab %>%
        arrange(desc(abs(rank_shift))) %>%
        slice_head(n = 3) %>%
        mutate(label = paste0(site, " (", sprintf("%+.1f", rank_shift), ")")) %>%
        pull(label),
      collapse = "; "
    )
  )
}

run_single_scenario_v18 <- function(scenario_type, scenario, ...) {
  m <- compute_model_once(
    site_master,
    substrate_axis_weights = rep(1 / 3, 3),
    scm_weights = rep(1 / 3, 3),
    fixed_bounds = substrate_fixed_bounds_v18,
    perturb_raw_substrate = FALSE,
    perturb_soilgrids = FALSE,
    ...
  )
  get_scenario_rank_v18(m, scenario_type, scenario)
}

# Build scenario rank directly from a vector of alternative relative scores.
build_rank_from_relative_v18 <- function(df, relative_col, scenario_type, scenario, extra_cols = tibble()) {
  relative_col <- rlang::ensym(relative_col)
  out <- df %>%
    filter(use_in_main_barplot) %>%
    arrange(desc(!!relative_col)) %>%
    mutate(
      scenario_type = scenario_type,
      scenario = scenario,
      scenario_rank = row_number(),
      scenario_relative = as.numeric(!!relative_col)
    ) %>%
    select(scenario_type, scenario, site, site_label_clean, scenario_rank, scenario_relative)

  if (nrow(extra_cols) > 0) {
    out <- bind_cols(out, extra_cols[rep(1, nrow(out)), , drop = FALSE])
  }
  out
}

# F3. Weight-grid sensitivities

sand_weight_grid <- make_simplex_grid(step = 0.1, include_zero = FALSE) %>%
  mutate(
    scenario_id = paste0("sand_w_", row_number()),
    w_mineral = w1,
    w_physical = w2,
    w_exchange = w3
  ) %>%
  select(scenario_id, w_mineral, w_physical, w_exchange)

sand_weight_rank_v18 <- purrr::pmap_dfr(
  sand_weight_grid,
  function(scenario_id, w_mineral, w_physical, w_exchange) {
    m <- compute_model_once(
      site_master,
      substrate_axis_weights = c(w_mineral, w_physical, w_exchange),
      scm_weights = rep(1 / 3, 3),
      fixed_bounds = substrate_fixed_bounds_v18,
      min_score = MODEL_FLOOR_MAIN,
      aggregation_method = "geometric",
      climate_mode = "ai_pet",
      perturb_raw_substrate = FALSE,
      perturb_soilgrids = FALSE
    )
    get_scenario_rank_v18(
      m,
      "sand_subaxis_weight_grid",
      scenario_id,
      tibble(w_mineral = w_mineral, w_physical = w_physical, w_exchange = w_exchange)
    )
  }
)

sand_weight_scenario_summary_v18 <- sand_weight_rank_v18 %>%
  group_split(scenario_type, scenario) %>%
  purrr::map_dfr(function(.x) {
    summarize_scenario_vs_main_v18(.x, .x$scenario_type[1], .x$scenario[1])
  }) %>%
  left_join(sand_weight_grid, by = c("scenario" = "scenario_id"))

sand_weight_site_summary_v18 <- sand_weight_rank_v18 %>%
  left_join(main_rank_v18, by = "site") %>%
  group_by(site, site_label_clean, main_rank) %>%
  summarise(
    median_rank_across_grid = median(scenario_rank, na.rm = TRUE),
    rank_min = min(scenario_rank, na.rm = TRUE),
    rank_max = max(scenario_rank, na.rm = TRUE),
    rank_range = rank_max - rank_min,
    top5_retention = mean(scenario_rank <= 5, na.rm = TRUE),
    median_relative = median(scenario_relative, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(main_rank)

scm_weight_grid <- make_simplex_grid(step = 0.1, include_zero = FALSE) %>%
  mutate(
    scenario_id = paste0("scm_w_", row_number()),
    w_sand = w1,
    w_climate = w2,
    w_microbe = w3
  ) %>%
  select(scenario_id, w_sand, w_climate, w_microbe)

scm_weight_rank_v18 <- purrr::pmap_dfr(
  scm_weight_grid,
  function(scenario_id, w_sand, w_climate, w_microbe) {
    m <- compute_model_once(
      site_master,
      substrate_axis_weights = rep(1 / 3, 3),
      scm_weights = c(w_sand, w_climate, w_microbe),
      fixed_bounds = substrate_fixed_bounds_v18,
      min_score = MODEL_FLOOR_MAIN,
      aggregation_method = "geometric",
      climate_mode = "ai_pet",
      perturb_raw_substrate = FALSE,
      perturb_soilgrids = FALSE
    )
    get_scenario_rank_v18(
      m,
      "SCM_weight_grid",
      scenario_id,
      tibble(w_sand = w_sand, w_climate = w_climate, w_microbe = w_microbe)
    )
  }
)

scm_weight_scenario_summary_v18 <- scm_weight_rank_v18 %>%
  group_split(scenario_type, scenario) %>%
  purrr::map_dfr(function(.x) {
    summarize_scenario_vs_main_v18(.x, .x$scenario_type[1], .x$scenario[1])
  }) %>%
  left_join(scm_weight_grid, by = c("scenario" = "scenario_id"))

scm_weight_site_summary_v18 <- scm_weight_rank_v18 %>%
  left_join(main_rank_v18, by = "site") %>%
  group_by(site, site_label_clean, main_rank) %>%
  summarise(
    median_rank_across_grid = median(scenario_rank, na.rm = TRUE),
    rank_min = min(scenario_rank, na.rm = TRUE),
    rank_max = max(scenario_rank, na.rm = TRUE),
    rank_range = rank_max - rank_min,
    top5_retention = mean(scenario_rank <= 5, na.rm = TRUE),
    median_relative = median(scenario_relative, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(main_rank)

# F4. Single-assumption sensitivities

# Climate-variable sensitivity: AI only and strict AI/PET minimum.
climate_variable_rank_v18 <- bind_rows(
  run_single_scenario_v18(
    "climate_variable",
    "AI_only",
    min_score = MODEL_FLOOR_MAIN,
    aggregation_method = "geometric",
    climate_mode = "ai_only"
  ),
  run_single_scenario_v18(
    "climate_variable",
    "AI_PET_min_constraint",
    min_score = MODEL_FLOOR_MAIN,
    aggregation_method = "geometric",
    climate_mode = "ai_pet_min"
  )
)

# Aggregation sensitivity.
aggregation_rank_v18 <- purrr::map_dfr(
  c("arithmetic", "geometric", "harmonic", "minimum"),
  function(method_i) {
    run_single_scenario_v18(
      "aggregation_method",
      method_i,
      min_score = MODEL_FLOOR_MAIN,
      aggregation_method = method_i,
      climate_mode = "ai_pet"
    )
  }
)

# Numerical-floor sensitivity.
floor_rank_v18 <- purrr::map_dfr(
  c(0.001, 0.01, 0.05, 0.10),
  function(floor_i) {
    run_single_scenario_v18(
      "numerical_floor",
      paste0("floor_", floor_i),
      min_score = floor_i,
      aggregation_method = "geometric",
      climate_mode = "ai_pet"
    )
  }
)

# 84-day Tengger field-period climate-anchor sensitivity.
field84_rank_v18 <- tibble()
if (exists("alxa_field_summary") && nrow(alxa_field_summary) > 0 &&
    "field_ai_p_over_pet" %in% names(alxa_field_summary) &&
    "field_pet_annualized" %in% names(alxa_field_summary) &&
    is.finite(as.numeric(alxa_field_summary$field_ai_p_over_pet[1])) &&
    is.finite(as.numeric(alxa_field_summary$field_pet_annualized[1])) &&
    as.numeric(alxa_field_summary$field_ai_p_over_pet[1]) > 0 &&
    as.numeric(alxa_field_summary$field_pet_annualized[1]) > 0) {

  alxa_ai_anchor_84day <- as.numeric(alxa_field_summary$field_ai_p_over_pet[1])
  alxa_pet_anchor_84day <- as.numeric(alxa_field_summary$field_pet_annualized[1])

  field84_rank_v18 <- run_single_scenario_v18(
    "field_84day_climate_anchor",
    "Tengger_84day_AI_PET_anchor",
    min_score = MODEL_FLOOR_MAIN,
    aggregation_method = "geometric",
    climate_mode = "ai_pet",
    alxa_ai_anchor = alxa_ai_anchor_84day,
    alxa_pet_anchor = alxa_pet_anchor_84day
  ) %>%
    mutate(
      field_ai_anchor = alxa_ai_anchor_84day,
      field_pet_anchor = alxa_pet_anchor_84day
    )
} else {
  warning(
    "Skipping field_84day_climate_anchor sensitivity because the 84-day climate anchor could not be computed."
  )
}

# Simplified-mineral alternative: mineral reactivity is represented by feldspar alone.
simplified_mineral_base_v18 <- site_master_det %>%
  mutate(
    simplified_mineral_axis = feldspar_score,
    simplified_substrate_axis = aggregate_scores(
      cbind(simplified_mineral_axis, physical_accessibility_axis, surface_exchange_axis),
      weights = rep(1 / 3, 3),
      method = "geometric",
      min_score = MODEL_FLOOR_MAIN
    ),
    simplified_weathering_raw = aggregate_scores(
      cbind(simplified_substrate_axis, climate_suitability, microbial_suitability),
      weights = rep(1 / 3, 3),
      method = "geometric",
      min_score = MODEL_FLOOR_MAIN
    )
  )

simplified_anchor_raw <- simplified_mineral_base_v18 %>%
  filter(site == "Alxa_Left_Banner_project_base_measured") %>%
  pull(simplified_weathering_raw) %>%
  .[1]

if (!is.finite(simplified_anchor_raw) || simplified_anchor_raw <= 0) {
  stop("Invalid Tengger anchor denominator for simplified-mineral sensitivity.")
}

simplified_mineral_base_v18 <- simplified_mineral_base_v18 %>%
  mutate(simplified_mineral_relative = simplified_weathering_raw / simplified_anchor_raw)

simplified_mineral_rank_v18 <- build_rank_from_relative_v18(
  simplified_mineral_base_v18,
  simplified_mineral_relative,
  "simplified_mineral_alternative",
  "feldspar_only_mineral_axis"
)

# PCA-derived substrate alternative: PC1 from the five substrate component scores.
rescale_to_training_range_v18 <- function(x, train_idx, min_score = MODEL_FLOOR_MAIN, max_score = 1) {
  x <- as.numeric(x)
  rng <- range(x[train_idx], na.rm = TRUE, finite = TRUE)
  if (!all(is.finite(rng)) || rng[1] == rng[2]) {
    stop("Invalid PCA training range for substrate sensitivity.")
  }
  out <- min_score + (max_score - min_score) * (x - rng[1]) / (rng[2] - rng[1])
  out[!is.na(out)] <- pmin(pmax(out[!is.na(out)], min_score), max_score)
  out
}

pca_substrate_vars_v18 <- c(
  "feldspar_score",
  "inverse_quartz_score",
  "size_surface_score",
  "silt_score",
  "cec_score"
)

pca_train_v18 <- site_master_det %>%
  filter(use_in_main_barplot) %>%
  filter(if_all(all_of(pca_substrate_vars_v18), ~ !is.na(.x)))

if (nrow(pca_train_v18) < 4) {
  stop("Not enough complete main sites for PCA-derived substrate sensitivity.")
}

pca_fit_v18 <- prcomp(
  pca_train_v18 %>% select(all_of(pca_substrate_vars_v18)),
  center = TRUE,
  scale. = TRUE
)

pca_scores_v18 <- predict(
  pca_fit_v18,
  newdata = site_master_det %>% select(all_of(pca_substrate_vars_v18))
)

pca_substrate_base_v18 <- site_master_det %>%
  mutate(pca_substrate_pc1_raw = as.numeric(pca_scores_v18[, 1]))

orient_cor_v18 <- suppressWarnings(cor(
  pca_substrate_base_v18$pca_substrate_pc1_raw[pca_substrate_base_v18$use_in_main_barplot],
  pca_substrate_base_v18$substrate_axis[pca_substrate_base_v18$use_in_main_barplot],
  method = "spearman",
  use = "complete.obs"
))

if (!is.na(orient_cor_v18) && orient_cor_v18 < 0) {
  pca_substrate_base_v18$pca_substrate_pc1_raw <- -pca_substrate_base_v18$pca_substrate_pc1_raw
}

train_idx_pca_v18 <- pca_substrate_base_v18$use_in_main_barplot &
  !is.na(pca_substrate_base_v18$pca_substrate_pc1_raw)

pca_substrate_base_v18 <- pca_substrate_base_v18 %>%
  mutate(
    pca_substrate_axis = rescale_to_training_range_v18(
      pca_substrate_pc1_raw,
      train_idx = train_idx_pca_v18,
      min_score = MODEL_FLOOR_MAIN,
      max_score = 1
    ),
    pca_substrate_weathering_raw = aggregate_scores(
      cbind(pca_substrate_axis, climate_suitability, microbial_suitability),
      weights = rep(1 / 3, 3),
      method = "geometric",
      min_score = MODEL_FLOOR_MAIN
    )
  )

pca_anchor_raw <- pca_substrate_base_v18 %>%
  filter(site == "Alxa_Left_Banner_project_base_measured") %>%
  pull(pca_substrate_weathering_raw) %>%
  .[1]

if (!is.finite(pca_anchor_raw) || pca_anchor_raw <= 0) {
  stop("Invalid Tengger anchor denominator for PCA-derived substrate sensitivity.")
}

pca_substrate_base_v18 <- pca_substrate_base_v18 %>%
  mutate(pca_substrate_relative = pca_substrate_weathering_raw / pca_anchor_raw)

pca_substrate_rank_v18 <- build_rank_from_relative_v18(
  pca_substrate_base_v18,
  pca_substrate_relative,
  "PCA_derived_substrate",
  "PC1_from_five_substrate_scores",
  tibble(pca_pc1_variance_explained = summary(pca_fit_v18)$importance[2, 1])
)

# F5. Consolidate and write sensitivity tables

single_scenario_rank_v18 <- bind_rows(
  climate_variable_rank_v18,
  aggregation_rank_v18,
  floor_rank_v18,
  field84_rank_v18,
  simplified_mineral_rank_v18,
  pca_substrate_rank_v18
)

single_scenario_summary_v18 <- single_scenario_rank_v18 %>%
  group_split(scenario_type, scenario) %>%
  purrr::map_dfr(function(.x) {
    summarize_scenario_vs_main_v18(.x, .x$scenario_type[1], .x$scenario[1])
  })

# Weight-grid scenario summaries already have their weights joined. Keep their
# weight columns in the detailed scenario tables, but use a shared set of columns
# for the overview table.
sensitivity_overview_table_v18 <- bind_rows(
  sand_weight_scenario_summary_v18 %>% select(any_of(names(single_scenario_summary_v18))),
  scm_weight_scenario_summary_v18 %>% select(any_of(names(single_scenario_summary_v18))),
  single_scenario_summary_v18
) %>%
  arrange(scenario_type, scenario)

all_scenario_ranks_v18 <- bind_rows(
  sand_weight_rank_v18,
  scm_weight_rank_v18,
  single_scenario_rank_v18
) %>%
  left_join(main_rank_v18, by = "site") %>%
  mutate(rank_shift = scenario_rank - main_rank) %>%
  arrange(scenario_type, scenario, scenario_rank)

site_stability_by_type_v18 <- all_scenario_ranks_v18 %>%
  group_by(scenario_type, site, site_label_clean, main_rank) %>%
  summarise(
    n_scenarios = n_distinct(scenario),
    median_rank = median(scenario_rank, na.rm = TRUE),
    rank_min = min(scenario_rank, na.rm = TRUE),
    rank_max = max(scenario_rank, na.rm = TRUE),
    rank_range = rank_max - rank_min,
    top5_retention = mean(scenario_rank <= 5, na.rm = TRUE),
    median_relative = median(scenario_relative, na.rm = TRUE),
    relative_min = min(scenario_relative, na.rm = TRUE),
    relative_max = max(scenario_relative, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(scenario_type, main_rank)

# Main output: one table for all scenario-level sensitivity results.
write.csv(
  sensitivity_overview_table_v18,
  file.path(output_dir, "v18_sensitivity_overview_table_for_methods_and_supp.csv"),
  row.names = FALSE
)

# Long detailed table: all ranks under all sensitivity scenarios.
write.csv(
  all_scenario_ranks_v18,
  file.path(output_dir, "v18_sensitivity_all_scenario_ranks_long.csv"),
  row.names = FALSE
)

# Site-level stability by sensitivity type.
write.csv(
  site_stability_by_type_v18,
  file.path(output_dir, "v18_sensitivity_site_stability_summary_by_type.csv"),
  row.names = FALSE
)

# Keep detailed grid tables for supplement because they summarize 36-scenario scans.
write.csv(
  sand_weight_rank_v18,
  file.path(output_dir, "v18_sensitivity_sand_subaxis_weight_grid_all_ranks.csv"),
  row.names = FALSE
)
write.csv(
  sand_weight_scenario_summary_v18,
  file.path(output_dir, "v18_sensitivity_sand_subaxis_weight_grid_scenario_summary.csv"),
  row.names = FALSE
)
write.csv(
  sand_weight_site_summary_v18,
  file.path(output_dir, "v18_sensitivity_sand_subaxis_weight_grid_site_summary.csv"),
  row.names = FALSE
)
write.csv(
  scm_weight_rank_v18,
  file.path(output_dir, "v18_sensitivity_SCM_weight_grid_all_ranks.csv"),
  row.names = FALSE
)
write.csv(
  scm_weight_scenario_summary_v18,
  file.path(output_dir, "v18_sensitivity_SCM_weight_grid_scenario_summary.csv"),
  row.names = FALSE
)
write.csv(
  scm_weight_site_summary_v18,
  file.path(output_dir, "v18_sensitivity_SCM_weight_grid_site_summary.csv"),
  row.names = FALSE
)
write.csv(
  single_scenario_rank_v18,
  file.path(output_dir, "v18_sensitivity_single_change_all_ranks.csv"),
  row.names = FALSE
)
write.csv(
  single_scenario_summary_v18,
  file.path(output_dir, "v18_sensitivity_single_change_summary_table.csv"),
  row.names = FALSE
)

cat("\nV18 final table-only sensitivity outputs written to:\n")
cat(output_dir, "\n")
cat("Sensitivity types included:\n")
print(sort(unique(sensitivity_overview_table_v18$scenario_type)))

# remove sensitivity analysis #

# V18 component-ablation sensitivity analysis
# Direct recombination version
#
# This version does NOT call compute_model_once().
# It uses already-computed deterministic V18 component/module scores
# from site_master_det and recombines them for ablation tests.

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(grid)

# 0. Output directories

if (!exists("output_dir")) {
  output_dir <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "v18_model_outputs")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

if (!exists("fig_output_dir")) {
  fig_output_dir <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "v18_figure_outputs")
  dir.create(fig_output_dir, showWarnings = FALSE, recursive = TRUE)
}

if (!exists("site_master_det")) {
  stop("Object `site_master_det` not found. Run the V18 deterministic model first.")
}

anchor_site <- "Alxa_Left_Banner_project_base_measured"
background_site <- "Alxa_Left_Banner_project_base_downloaded"

# 1. Helper: find columns robustly

find_col_v18 <- function(df, candidates, label) {
  hit <- candidates[candidates %in% names(df)]
  
  if (length(hit) >= 1) {
    message(label, " -> ", hit[1])
    return(hit[1])
  }
  
  message("\nCould not find column for: ", label)
  message("Tried candidates:")
  message(paste(candidates, collapse = ", "))
  message("\nAvailable columns:")
  message(paste(names(df), collapse = ", "))
  stop("Column not found: ", label)
}

make_v18_site_label <- function(site) {
  dplyr::case_when(
    site == "Hobq" ~ "Kubuqi",
    site == "Sonoran_Altar" ~ "Gran Desierto",
    site == "Argentina_MedanosGrandes_NE" ~ "Medanos Grandes",
    site == "Argentina_MedanosGra" ~ "Medanos Grandes",
    site == "Kalahari_Bagatelle" ~ "Kalahari",
    site == "Karakum_AshkAbad" ~ "Karakum",
    site == "Thal_Mankera" ~ "Thal",
    site == anchor_site ~ "Tengger experimental anchor",
    TRUE ~ stringr::str_replace_all(site, "_", " ")
  )
}

geom_mean_v18 <- function(x, floor = 0.05) {
  x <- as.numeric(x)
  x <- pmax(x, floor)
  exp(mean(log(x), na.rm = TRUE))
}

geom_mean_row_v18 <- function(mat, floor = 0.05) {
  mat <- as.matrix(mat)
  mat <- apply(mat, 2, as.numeric)
  mat <- pmax(mat, floor)
  exp(rowMeans(log(mat), na.rm = TRUE))
}

MODEL_FLOOR_MAIN <- ifelse(exists("MODEL_FLOOR_MAIN"), MODEL_FLOOR_MAIN, 0.05)

# 2. Identify required columns
#    Add more candidate names here if your script used different names.

col_mineral <- find_col_v18(
  site_master_det,
  candidates = c(
    "mineral_reactivity_axis",
    "mineral_reactivity_score",
    "mineral_reactivity",
    "score_mineral_reactivity"
  ),
  label = "Mineral reactivity"
)

col_physical <- find_col_v18(
  site_master_det,
  candidates = c(
    "physical_accessibility_axis",
    "physical_accessibility_score",
    "physical_accessibility",
    "particle_accessibility_axis",
    "particle_accessibility_score",
    "grain_size_accessibility_score"
  ),
  label = "Physical accessibility"
)

col_surface <- find_col_v18(
  site_master_det,
  candidates = c(
    "surface_exchange_axis",
    "surface_exchange_score",
    "surface_exchange",
    "surface_exchange_capacity_score"
  ),
  label = "Surface exchange"
)

col_climate <- find_col_v18(
  site_master_det,
  candidates = c(
    "climate_suitability",
    "climate_suitability_score",
    "climate_constraint_score",
    "climate_score"
  ),
  label = "Climate suitability"
)

col_microbial <- find_col_v18(
  site_master_det,
  candidates = c(
    "microbial_suitability",
    "microbial_suitability_score",
    "microbial_response_score",
    "microbial_function_score",
    "microbial_score"
  ),
  label = "Microbial suitability"
)

col_relative <- find_col_v18(
  site_master_det,
  candidates = c(
    "relative_to_alxa_measured",
    "relative_potential",
    "relative_weathering_potential"
  ),
  label = "Reference relative potential"
)

# 3. Prepare component table

required_basic_cols <- c(
  "site",
  "include_in_main_model",
  "full_model_evaluated"
)

missing_basic_cols <- setdiff(required_basic_cols, names(site_master_det))

if (length(missing_basic_cols) > 0) {
  stop(
    "Missing basic columns: ",
    paste(missing_basic_cols, collapse = ", ")
  )
}

component_df_v18 <- site_master_det %>%
  filter(
    include_in_main_model,
    full_model_evaluated,
    !site %in% c(anchor_site, background_site)
  ) %>%
  mutate(
    site_label = if ("site_label_clean" %in% names(.)) {
      ifelse(
        !is.na(site_label_clean) & site_label_clean != "",
        as.character(site_label_clean),
        make_v18_site_label(site)
      )
    } else {
      make_v18_site_label(site)
    },
    mineral_reactivity = as.numeric(.data[[col_mineral]]),
    physical_accessibility = as.numeric(.data[[col_physical]]),
    surface_exchange = as.numeric(.data[[col_surface]]),
    climate_suitability = as.numeric(.data[[col_climate]]),
    microbial_suitability = as.numeric(.data[[col_microbial]]),
    reference_relative = as.numeric(.data[[col_relative]])
  )

# Check anchor raw score from deterministic model
# For rank shift, relative values only need consistent ranking.
# For ablation relative values, we normalize each ablated raw score to the ablated Tengger anchor raw score if anchor is available.

anchor_component_v18 <- site_master_det %>%
  filter(site == anchor_site) %>%
  mutate(
    mineral_reactivity = as.numeric(.data[[col_mineral]]),
    physical_accessibility = as.numeric(.data[[col_physical]]),
    surface_exchange = as.numeric(.data[[col_surface]]),
    climate_suitability = as.numeric(.data[[col_climate]]),
    microbial_suitability = as.numeric(.data[[col_microbial]])
  )

if (nrow(anchor_component_v18) != 1) {
  stop("Could not find exactly one Tengger experimental anchor row.")
}

# 4. Reference ranking

ref_rank_v18 <- component_df_v18 %>%
  mutate(
    ref_rank = rank(-reference_relative, ties.method = "average")
  ) %>%
  select(site, site_label, ref_relative = reference_relative, ref_rank)

site_order_by_ref_v18 <- ref_rank_v18 %>%
  arrange(ref_rank) %>%
  pull(site)

site_label_order_by_ref_v18 <- ref_rank_v18 %>%
  arrange(ref_rank) %>%
  pull(site_label)

message("\nReference ranking:")
print(ref_rank_v18 %>% arrange(ref_rank))

# 5. Calculate ablated scores

calc_ablation_relative_v18 <- function(df, anchor_df, scenario) {
  
  if (scenario == "remove_mineral") {
    
    raw_site <- geom_mean_row_v18(
      df %>%
        select(
          physical_accessibility,
          surface_exchange,
          climate_suitability,
          microbial_suitability
        ),
      floor = MODEL_FLOOR_MAIN
    )
    
    raw_anchor <- geom_mean_v18(
      c(
        anchor_df$physical_accessibility,
        anchor_df$surface_exchange,
        anchor_df$climate_suitability,
        anchor_df$microbial_suitability
      ),
      floor = MODEL_FLOOR_MAIN
    )
    
  } else if (scenario == "remove_particle_size") {
    
    raw_site <- geom_mean_row_v18(
      df %>%
        select(
          mineral_reactivity,
          surface_exchange,
          climate_suitability,
          microbial_suitability
        ),
      floor = MODEL_FLOOR_MAIN
    )
    
    raw_anchor <- geom_mean_v18(
      c(
        anchor_df$mineral_reactivity,
        anchor_df$surface_exchange,
        anchor_df$climate_suitability,
        anchor_df$microbial_suitability
      ),
      floor = MODEL_FLOOR_MAIN
    )
    
  } else if (scenario == "remove_surface_exchange") {
    
    raw_site <- geom_mean_row_v18(
      df %>%
        select(
          mineral_reactivity,
          physical_accessibility,
          climate_suitability,
          microbial_suitability
        ),
      floor = MODEL_FLOOR_MAIN
    )
    
    raw_anchor <- geom_mean_v18(
      c(
        anchor_df$mineral_reactivity,
        anchor_df$physical_accessibility,
        anchor_df$climate_suitability,
        anchor_df$microbial_suitability
      ),
      floor = MODEL_FLOOR_MAIN
    )
    
  } else if (scenario == "remove_climate") {
    
    raw_site <- geom_mean_row_v18(
      df %>%
        select(
          mineral_reactivity,
          physical_accessibility,
          surface_exchange,
          microbial_suitability
        ),
      floor = MODEL_FLOOR_MAIN
    )
    
    raw_anchor <- geom_mean_v18(
      c(
        anchor_df$mineral_reactivity,
        anchor_df$physical_accessibility,
        anchor_df$surface_exchange,
        anchor_df$microbial_suitability
      ),
      floor = MODEL_FLOOR_MAIN
    )
    
  } else if (scenario == "remove_microbial_function") {
    
    raw_site <- geom_mean_row_v18(
      df %>%
        select(
          mineral_reactivity,
          physical_accessibility,
          surface_exchange,
          climate_suitability
        ),
      floor = MODEL_FLOOR_MAIN
    )
    
    raw_anchor <- geom_mean_v18(
      c(
        anchor_df$mineral_reactivity,
        anchor_df$physical_accessibility,
        anchor_df$surface_exchange,
        anchor_df$climate_suitability
      ),
      floor = MODEL_FLOOR_MAIN
    )
    
  } else {
    stop("Unknown ablation scenario: ", scenario)
  }
  
  raw_site / raw_anchor
}

ablation_scenarios_v18 <- tibble::tibble(
  scenario = c(
    "remove_mineral",
    "remove_particle_size",
    "remove_surface_exchange",
    "remove_climate",
    "remove_microbial_function"
  ),
  scenario_label = c(
    "Remove\nmineral",
    "Remove\nparticle size",
    "Remove\nsurface exchange",
    "Remove\nclimate",
    "Remove\nmicrobial function"
  )
)

ablation_rank_tables_v18 <- lapply(
  seq_len(nrow(ablation_scenarios_v18)),
  function(i) {
    
    scenario_i <- ablation_scenarios_v18$scenario[i]
    
    component_df_v18 %>%
      mutate(
        scenario = scenario_i,
        alt_relative = calc_ablation_relative_v18(
          df = .,
          anchor_df = anchor_component_v18,
          scenario = scenario_i
        ),
        alt_rank = rank(-alt_relative, ties.method = "average")
      ) %>%
      select(site, site_label, scenario, alt_relative, alt_rank)
  }
) %>%
  bind_rows()

# 6. Rank shift table

ablation_rank_shift_v18 <- ref_rank_v18 %>%
  left_join(
    ablation_rank_tables_v18,
    by = c("site", "site_label")
  ) %>%
  left_join(
    ablation_scenarios_v18,
    by = "scenario"
  ) %>%
  mutate(
    rank_shift = alt_rank - ref_rank,
    abs_rank_shift = abs(rank_shift),
    site = factor(site, levels = rev(site_order_by_ref_v18)),
    site_label = factor(site_label, levels = rev(site_label_order_by_ref_v18)),
    scenario_label = factor(
      scenario_label,
      levels = ablation_scenarios_v18$scenario_label
    ),
    rank_shift_label = case_when(
      rank_shift > 0 ~ paste0("+", rank_shift),
      rank_shift == 0 ~ "+0",
      TRUE ~ as.character(rank_shift)
    )
  )

write.csv(
  ablation_rank_shift_v18,
  file.path(output_dir, "v18_component_ablation_rank_shift_table.csv"),
  row.names = FALSE
)

# 7. Scenario-level robustness metrics

ablation_metrics_v18 <- ablation_rank_shift_v18 %>%
  group_by(scenario, scenario_label) %>%
  summarize(
    spearman_rho_vs_main = suppressWarnings(
      cor(ref_rank, alt_rank, method = "spearman", use = "complete.obs")
    ),
    top5_overlap = length(
      intersect(
        as.character(site)[ref_rank <= 5],
        as.character(site)[alt_rank <= 5]
      )
    ) / 5,
    median_abs_rank_shift = median(abs_rank_shift, na.rm = TRUE),
    max_abs_rank_shift = max(abs_rank_shift, na.rm = TRUE),
    most_shifted_sites = paste(
      as.character(site_label)[order(-abs_rank_shift)][1:min(3, n())],
      collapse = "; "
    ),
    .groups = "drop"
  )

write.csv(
  ablation_metrics_v18,
  file.path(output_dir, "v18_component_ablation_metrics.csv"),
  row.names = FALSE
)

message("\nAblation robustness metrics:")
print(ablation_metrics_v18)

# 8. Plot: rank-shift heatmap

rank_shift_limit_v18 <- max(
  abs(ablation_rank_shift_v18$rank_shift),
  na.rm = TRUE
)

rank_shift_limit_v18 <- max(rank_shift_limit_v18, 1)

rank_shift_breaks_v18 <- seq(
  -rank_shift_limit_v18,
  rank_shift_limit_v18,
  by = 1
)

theme_v18_ablation_heatmap <- theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    
    axis.title = element_blank(),
    axis.text.y = element_text(size = 9.2, colour = "#222222"),
    axis.text.x = element_text(
      size = 8.8,
      colour = "#222222",
      angle = 35,
      hjust = 1,
      vjust = 1,
      lineheight = 0.9
    ),
    axis.ticks = element_blank(),
    
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = margin(6, 6, 4, 6),
    
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    legend.key.height = unit(3.5, "mm"),
    legend.key.width  = unit(8, "mm"),
    legend.background = element_rect(fill = "white", colour = NA)
  )

p_v18_component_ablation_rank_shift <- ggplot(
  ablation_rank_shift_v18,
  aes(
    x = scenario_label,
    y = site_label,
    fill = rank_shift
  )
) +
  geom_tile(
    width = 0.94,
    height = 0.92,
    colour = "white",
    linewidth = 0.55
  ) +
  geom_text(
    aes(label = rank_shift_label),
    size = 3.0,
    colour = "#222222"
  ) +
  scale_fill_gradient2(
    low = "#2D7A78",
    mid = "white",
    high = "#A77A18",
    midpoint = 0,
    limits = c(-rank_shift_limit_v18, rank_shift_limit_v18),
    breaks = rank_shift_breaks_v18,
    name = "Rank shift",
    guide = guide_colourbar(
      title.position = "left",
      label.position = "bottom",
      barheight = unit(3.2, "mm"),
      barwidth = unit(34, "mm")
    )
  ) +
  coord_equal() +
  theme_v18_ablation_heatmap

print(p_v18_component_ablation_rank_shift)

ggsave(
  file.path(data_dir, "Desktop", "replot", "global_weathering_input", "out-v18", "p_v18_component_ablation_rank_shift.pdf"),
  plot = p_v18_component_ablation_rank_shift,
  width = 8,
  height = 5,
  bg = "white"
)

ggsave(
  file.path(fig_output_dir, "p_v18_component_ablation_rank_shift.png"),
  plot = p_v18_component_ablation_rank_shift,
  width = 5.8,
  height = 5.3,
  dpi = 600,
  bg = "white"
)

# 9. Optional combined plot with module collinearity heatmap

if (exists("p_v18_module_collinearity")) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
  }
  
  library(patchwork)
  
  p_v18_collinearity_ablation_combined <-
    p_v18_module_collinearity +
    p_v18_component_ablation_rank_shift +
    plot_layout(widths = c(1.08, 0.92)) +
    plot_annotation(tag_levels = "A")
  
  print(p_v18_collinearity_ablation_combined)
  
  ggsave(
    file.path(fig_output_dir, "p_v18_collinearity_ablation_combined.pdf"),
    plot = p_v18_collinearity_ablation_combined,
    width = 10.5,
    height = 5.5,
    bg = "white"
  )
  
  ggsave(
    file.path(fig_output_dir, "p_v18_collinearity_ablation_combined.png"),
    plot = p_v18_collinearity_ablation_combined,
    width = 10.5,
    height = 5.5,
    dpi = 600,
    bg = "white"
  )
}

message("\nFinished V18 component-ablation sensitivity analysis.")
message("Saved table: ", file.path(output_dir, "v18_component_ablation_rank_shift_table.csv"))
message("Saved plot: ", file.path(fig_output_dir, "p_v18_component_ablation_rank_shift.pdf"))

# V18 flattened component-level sensitivity analysis
#
# Purpose:
#   Test whether the hierarchical model structure affects global rankings.
#
# Main hierarchical model:
#   S = geom_mean(mineral, physical, surface)
#   C = geom_mean(AI, PET)
#   M = geom_mean(pH, temperature)
#   P = geom_mean(S, C, M)
#
# Flattened component-level alternative:
#   P_flat = geom_mean(
#     mineral,
#     physical,
#     surface,
#     AI,
#     PET,
#     pH,
#     temperature
#   )
#
# Each terminal component receives equal weight = 1/7.

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(grid)

# 0. Output directories

if (!exists("output_dir")) {
  output_dir <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "v18_model_outputs")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

if (!exists("fig_output_dir")) {
  fig_output_dir <- file.path(data_dir, "Desktop", "replot", "global_weathering_input", "v18_figure_outputs")
  dir.create(fig_output_dir, showWarnings = FALSE, recursive = TRUE)
}

if (!exists("site_master_det")) {
  stop("Object `site_master_det` not found. Run the V18 deterministic model first.")
}

if (!exists("MODEL_FLOOR_MAIN")) {
  MODEL_FLOOR_MAIN <- 0.05
}

anchor_site <- "Alxa_Left_Banner_project_base_measured"
background_site <- "Alxa_Left_Banner_project_base_downloaded"

# 1. Helper functions

find_col_v18 <- function(df, candidates, label) {
  hit <- candidates[candidates %in% names(df)]
  
  if (length(hit) >= 1) {
    message(label, " -> ", hit[1])
    return(hit[1])
  }
  
  message("\nCould not find column for: ", label)
  message("Tried candidates:")
  message(paste(candidates, collapse = ", "))
  message("\nAvailable columns:")
  message(paste(names(df), collapse = ", "))
  stop("Column not found: ", label)
}

make_v18_site_label <- function(site) {
  dplyr::case_when(
    site == "Hobq" ~ "Kubuqi",
    site == "Sonoran_Altar" ~ "Gran Desierto",
    site == "Argentina_MedanosGrandes_NE" ~ "Medanos Grandes",
    site == "Argentina_MedanosGra" ~ "Medanos Grandes",
    site == "Kalahari_Bagatelle" ~ "Kalahari",
    site == "Karakum_AshkAbad" ~ "Karakum",
    site == "Thal_Mankera" ~ "Thal",
    site == anchor_site ~ "Tengger experimental anchor",
    TRUE ~ stringr::str_replace_all(site, "_", " ")
  )
}

geom_mean_row_v18 <- function(mat, floor = 0.05) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat <- pmax(mat, floor)
  exp(rowMeans(log(mat), na.rm = TRUE))
}

geom_mean_v18 <- function(x, floor = 0.05) {
  x <- as.numeric(x)
  x <- pmax(x, floor)
  exp(mean(log(x), na.rm = TRUE))
}

# 2. Identify required terminal component columns
#    Add candidates here if your V18 object used different names.

col_mineral <- find_col_v18(
  site_master_det,
  candidates = c(
    "mineral_reactivity_axis",
    "mineral_reactivity_score",
    "mineral_reactivity",
    "score_mineral_reactivity"
  ),
  label = "Mineral reactivity"
)

col_physical <- find_col_v18(
  site_master_det,
  candidates = c(
    "physical_accessibility_axis",
    "physical_accessibility_score",
    "physical_accessibility",
    "particle_accessibility_axis",
    "particle_accessibility_score",
    "grain_size_accessibility_score"
  ),
  label = "Physical accessibility"
)

col_surface <- find_col_v18(
  site_master_det,
  candidates = c(
    "surface_exchange_axis",
    "surface_exchange_score",
    "surface_exchange",
    "surface_exchange_capacity_score"
  ),
  label = "Surface exchange"
)

col_ai <- find_col_v18(
  site_master_det,
  candidates = c(
    "aridity_constraint_score",
    "aridity_suitability",
    "aridity_suitability_score",
    "ai_constraint_score",
    "ai_suitability_score",
    "SAI",
    "s_ai"
  ),
  label = "Aridity suitability"
)

col_pet <- find_col_v18(
  site_master_det,
  candidates = c(
    "pet_constraint_score",
    "pet_suitability",
    "pet_suitability_score",
    "SPET",
    "s_pet"
  ),
  label = "PET suitability"
)

col_ph <- find_col_v18(
  site_master_det,
  candidates = c(
    "microbe_ph_function",
    "microbial_ph_score",
    "ph_response_score",
    "f8_ph_response_score",
    "f8_ph_score",
    "pH_response_score",
    "microbial_pH_score"
  ),
  label = "F8 pH response"
)

col_temp <- find_col_v18(
  site_master_det,
  candidates = c(
    "microbe_mat_function",
    "microbial_temp_score",
    "microbial_temperature_score",
    "temperature_response_score",
    "f8_temperature_response_score",
    "f8_temp_response_score",
    "f8_temperature_score",
    "mat_response_score"
  ),
  label = "F8 temperature response"
)


col_relative <- find_col_v18(
  site_master_det,
  candidates = c(
    "relative_to_alxa_measured",
    "relative_potential",
    "relative_weathering_potential"
  ),
  label = "Main reference relative potential"
)

# 3. Prepare component table

required_basic_cols <- c(
  "site",
  "include_in_main_model",
  "full_model_evaluated"
)

missing_basic_cols <- setdiff(required_basic_cols, names(site_master_det))

if (length(missing_basic_cols) > 0) {
  stop(
    "Missing basic columns in `site_master_det`: ",
    paste(missing_basic_cols, collapse = ", ")
  )
}

flat_component_all_v18 <- site_master_det %>%
  mutate(
    site_label = if ("site_label_clean" %in% names(.)) {
      ifelse(
        !is.na(site_label_clean) & site_label_clean != "",
        as.character(site_label_clean),
        make_v18_site_label(site)
      )
    } else {
      make_v18_site_label(site)
    },
    mineral_reactivity = as.numeric(.data[[col_mineral]]),
    physical_accessibility = as.numeric(.data[[col_physical]]),
    surface_exchange = as.numeric(.data[[col_surface]]),
    aridity_suitability = as.numeric(.data[[col_ai]]),
    pet_suitability = as.numeric(.data[[col_pet]]),
    f8_ph_response = as.numeric(.data[[col_ph]]),
    f8_temperature_response = as.numeric(.data[[col_temp]]),
    main_relative_potential = as.numeric(.data[[col_relative]])
  )

anchor_flat_v18 <- flat_component_all_v18 %>%
  filter(site == anchor_site)

if (nrow(anchor_flat_v18) != 1) {
  stop("Could not find exactly one Tengger experimental anchor row.")
}

flat_component_df_v18 <- flat_component_all_v18 %>%
  filter(
    include_in_main_model,
    full_model_evaluated,
    !site %in% c(anchor_site, background_site),
    !is.na(main_relative_potential)
  )

# 4. Calculate flattened component-level relative potential

flat_anchor_raw_v18 <- geom_mean_v18(
  c(
    anchor_flat_v18$mineral_reactivity,
    anchor_flat_v18$physical_accessibility,
    anchor_flat_v18$surface_exchange,
    anchor_flat_v18$aridity_suitability,
    anchor_flat_v18$pet_suitability,
    anchor_flat_v18$f8_ph_response,
    anchor_flat_v18$f8_temperature_response
  ),
  floor = MODEL_FLOOR_MAIN
)

flat_result_v18 <- flat_component_df_v18 %>%
  mutate(
    flat_raw_score = geom_mean_row_v18(
      select(
        .,
        mineral_reactivity,
        physical_accessibility,
        surface_exchange,
        aridity_suitability,
        pet_suitability,
        f8_ph_response,
        f8_temperature_response
      ),
      floor = MODEL_FLOOR_MAIN
    ),
    flat_relative_potential = flat_raw_score / flat_anchor_raw_v18,
    main_rank = rank(-main_relative_potential, ties.method = "average"),
    flat_rank = rank(-flat_relative_potential, ties.method = "average"),
    rank_shift = flat_rank - main_rank,
    abs_rank_shift = abs(rank_shift),
    rank_shift_label = case_when(
      rank_shift > 0 ~ paste0("+", rank_shift),
      rank_shift == 0 ~ "+0",
      TRUE ~ as.character(rank_shift)
    )
  ) %>%
  arrange(main_rank)

# 5. Site-level result table

flat_site_table_v18 <- flat_result_v18 %>%
  transmute(
    site,
    site_label,
    main_relative_potential,
    flat_relative_potential,
    main_rank,
    flat_rank,
    rank_shift,
    abs_rank_shift,
    mineral_reactivity,
    physical_accessibility,
    surface_exchange,
    aridity_suitability,
    pet_suitability,
    f8_ph_response,
    f8_temperature_response
  )

write.csv(
  flat_site_table_v18,
  file.path(output_dir, "v18_flattened_component_level_site_table.csv"),
  row.names = FALSE
)

print(flat_site_table_v18)

# 6. Scenario-level summary table

flat_summary_v18 <- flat_result_v18 %>%
  summarize(
    scenario_type = "hierarchical_structure",
    scenario = "flattened_component_equal",
    scenario_description = "Seven terminal component scores were combined directly using an equal-weight geometric mean.",
    spearman_rho_vs_main = suppressWarnings(
      cor(main_rank, flat_rank, method = "spearman", use = "complete.obs")
    ),
    top3_overlap = length(
      intersect(
        site[main_rank <= 3],
        site[flat_rank <= 3]
      )
    ) / 3,
    top5_overlap = length(
      intersect(
        site[main_rank <= 5],
        site[flat_rank <= 5]
      )
    ) / 5,
    median_abs_rank_shift = median(abs_rank_shift, na.rm = TRUE),
    max_abs_rank_shift = max(abs_rank_shift, na.rm = TRUE),
    most_shifted_sites = paste(
      site_label[order(-abs_rank_shift)][1:min(5, n())],
      collapse = "; "
    ),
    .groups = "drop"
  )

write.csv(
  flat_summary_v18,
  file.path(output_dir, "v18_flattened_component_level_summary.csv"),
  row.names = FALSE
)

print(flat_summary_v18)

# 7. Optional: append flattened model to Table S6 overview

overview_path_v18 <- file.path(
  output_dir,
  "v18_sensitivity_overview_table_for_methods_and_supp.csv"
)

if (file.exists(overview_path_v18)) {
  
  sensitivity_overview_v18 <- read.csv(
    overview_path_v18,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  flat_overview_row_v18 <- flat_summary_v18 %>%
    transmute(
      scenario_type,
      scenario,
      spearman_rho_vs_main,
      top5_overlap,
      median_abs_rank_shift,
      max_abs_rank_shift,
      most_shifted_sites
    )
  
  # Add missing columns if the existing overview has extra fields
  missing_in_flat <- setdiff(names(sensitivity_overview_v18), names(flat_overview_row_v18))
  for (cc in missing_in_flat) {
    flat_overview_row_v18[[cc]] <- NA
  }
  
  # Add missing columns to existing overview if needed
  missing_in_overview <- setdiff(names(flat_overview_row_v18), names(sensitivity_overview_v18))
  for (cc in missing_in_overview) {
    sensitivity_overview_v18[[cc]] <- NA
  }
  
  flat_overview_row_v18 <- flat_overview_row_v18[, names(sensitivity_overview_v18)]
  
  sensitivity_overview_v18_updated <- sensitivity_overview_v18 %>%
    filter(
      !(
        scenario_type == "hierarchical_structure" &
          scenario == "flattened_component_equal"
      )
    ) %>%
    bind_rows(flat_overview_row_v18)
  
  write.csv(
    sensitivity_overview_v18_updated,
    overview_path_v18,
    row.names = FALSE
  )
  
  message("Updated overview table: ", overview_path_v18)
}

# 8. Plot 1: main vs flattened relative potential

site_order_flat_v18 <- flat_result_v18 %>%
  arrange(main_rank) %>%
  pull(site_label)

flat_plot_long_v18 <- flat_result_v18 %>%
  select(
    site_label,
    main_relative_potential,
    flat_relative_potential
  ) %>%
  pivot_longer(
    cols = c(main_relative_potential, flat_relative_potential),
    names_to = "model",
    values_to = "relative_potential"
  ) %>%
  mutate(
    model = recode(
      model,
      main_relative_potential = "Hierarchical model",
      flat_relative_potential = "Flattened component-level model"
    ),
    model = factor(
      model,
      levels = c(
        "Hierarchical model",
        "Flattened component-level model"
      )
    ),
    site_label = factor(site_label, levels = rev(site_order_flat_v18))
  )

theme_v18_flat <- theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9.2, colour = "#222222"),
    axis.text.x = element_text(size = 9, colour = "#222222"),
    axis.title.x = element_text(size = 10, colour = "#222222"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    plot.title = element_blank(),
    plot.margin = margin(6, 6, 4, 6)
  )

p_v18_flattened_relative_comparison <- ggplot(
  flat_plot_long_v18,
  aes(
    x = relative_potential,
    y = site_label,
    shape = model,
    fill = model
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "#B22222"
  ) +
  geom_point(
    size = 2.8,
    colour = "#222222",
    alpha = 0.88,
    stroke = 0.35
  ) +
  scale_x_continuous(
    labels = function(x) paste0(round(x, 2), "x"),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  scale_shape_manual(
    values = c(
      "Hierarchical model" = 21,
      "Flattened component-level model" = 24
    )
  ) +
  scale_fill_manual(
    values = c(
      "Hierarchical model" = "#2D7A78",
      "Flattened component-level model" = "#A77A18"
    )
  ) +
  labs(
    x = "Relative potential to Tengger anchor"
  ) +
  theme_v18_flat

print(p_v18_flattened_relative_comparison)

ggsave(
  file.path(fig_output_dir, "p_v18_flattened_component_relative_comparison.pdf"),
  plot = p_v18_flattened_relative_comparison,
  width = 5.6,
  height = 4.6,
  bg = "white"
)

ggsave(
  file.path(fig_output_dir, "p_v18_flattened_component_relative_comparison.png"),
  plot = p_v18_flattened_relative_comparison,
  width = 5.6,
  height = 4.6,
  dpi = 600,
  bg = "white"
)

# 9. Plot 2: rank-shift heatmap, consistent with ablation style

rank_shift_limit_flat_v18 <- max(
  abs(flat_result_v18$rank_shift),
  na.rm = TRUE
)

rank_shift_limit_flat_v18 <- max(rank_shift_limit_flat_v18, 1)

flat_rank_shift_plot_v18 <- flat_result_v18 %>%
  mutate(
    site_label = factor(site_label, levels = rev(site_order_flat_v18)),
    scenario_label = factor(
      "Flattened\ncomponent-level",
      levels = "Flattened\ncomponent-level"
    )
  )

p_v18_flattened_component_rank_shift <- ggplot(
  flat_rank_shift_plot_v18,
  aes(
    x = scenario_label,
    y = site_label,
    fill = rank_shift
  )
) +
  geom_tile(
    width = 0.90,
    height = 0.92,
    colour = "white",
    linewidth = 0.55
  ) +
  geom_text(
    aes(label = rank_shift_label),
    size = 3.0,
    colour = "#222222"
  ) +
  scale_fill_gradient2(
    low = "#2D7A78",
    mid = "white",
    high = "#A77A18",
    midpoint = 0,
    limits = c(-rank_shift_limit_flat_v18, rank_shift_limit_flat_v18),
    breaks = seq(
      -rank_shift_limit_flat_v18,
      rank_shift_limit_flat_v18,
      by = 1
    ),
    name = "Rank shift",
    guide = guide_colourbar(
      title.position = "left",
      label.position = "bottom",
      barheight = unit(3.2, "mm"),
      barwidth = unit(28, "mm")
    )
  ) +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 9.2, colour = "#222222"),
    axis.text.x = element_text(
      size = 8.8,
      colour = "#222222",
      angle = 35,
      hjust = 1,
      vjust = 1,
      lineheight = 0.9
    ),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    legend.key.height = unit(3.5, "mm"),
    legend.key.width  = unit(8, "mm"),
    legend.background = element_rect(fill = "white", colour = NA),
    plot.title = element_blank(),
    plot.margin = margin(6, 6, 4, 6)
  )

print(p_v18_flattened_component_rank_shift)

ggsave(
  file.path(fig_output_dir, "p_v18_flattened_component_rank_shift.pdf"),
  plot = p_v18_flattened_component_rank_shift,
  width = 2.4,
  height = 5.3,
  bg = "white"
)

ggsave(
  file.path(fig_output_dir, "p_v18_flattened_component_rank_shift.png"),
  plot = p_v18_flattened_component_rank_shift,
  width = 2.4,
  height = 5.3,
  dpi = 600,
  bg = "white"
)

# 10. Optional: combine with component-ablation heatmap

if (exists("p_v18_component_ablation_rank_shift")) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
  }
  
  library(patchwork)
  
  p_v18_ablation_flattened_combined <-
    p_v18_component_ablation_rank_shift +
    p_v18_flattened_component_rank_shift +
    plot_layout(widths = c(1.0, 0.32)) +
    plot_annotation(tag_levels = "A")
  
  print(p_v18_ablation_flattened_combined)
  
  ggsave(
    file.path(fig_output_dir, "p_v18_ablation_flattened_combined.pdf"),
    plot = p_v18_ablation_flattened_combined,
    width = 7.2,
    height = 5.4,
    bg = "white"
  )
  
  ggsave(
    file.path(fig_output_dir, "p_v18_ablation_flattened_combined.png"),
    plot = p_v18_ablation_flattened_combined,
    width = 7.2,
    height = 5.4,
    dpi = 600,
    bg = "white"
  )
}

message("\nFinished V18 flattened component-level sensitivity analysis.")
message("Saved site table: ", file.path(output_dir, "v18_flattened_component_level_site_table.csv"))
message("Saved summary: ", file.path(output_dir, "v18_flattened_component_level_summary.csv"))
message("Saved plot: ", file.path(fig_output_dir, "p_v18_flattened_component_relative_comparison.pdf"))
message("Saved plot: ", file.path(fig_output_dir, "p_v18_flattened_component_rank_shift.pdf"))

