# =====================================================================
# fidelity_evaluation.R
# RATE-AF Synthetic Clinical Trials
# Functions for evaluating fidelity between real and synthetic data.
# =====================================================================

# Default variable lists and labels specific to RATE-AF
{
  vars_categorical <- c("sex", "mdfdehrascr", "nyha", "hfdiag", "typ2dm", "ethncty", "estmtdejctnfrctn")
  categorical_labels <- c(
    sex = "Gender", mdfdehrascr = "mEHRA class", nyha = "NYHA class",
    hfdiag = "Diagnosis of HF", typ2dm = "Type II Diabetes",
    ethncty = "Ethnicity", estmtdejctnfrctn = "Estimated LVEF")
  
  # Define numeric variable labels
  vars_numeric <- c("age", "creat", "ntprobnp", "ecghr", "bpsysav")
  numeric_labels <- c(
    age = "Age (years)", creat = "Creatinine (micromol/L)", ntprobnp = "NT-proBNP (pg/mL)",
    ecghr = "ECG Heart Rate (bpm)", bpsysav = "Systolic BP (mmHg)")
  
}


# Functions to evaluate fidelity in categorical columns
fidelity_evaluation_categorical <- function(var, real, dataset_list,labels) {
  output <- list()
  var_label <- labels[[var]]
  
  # Get non-NA levels from real and all synthetic datasets
  all_levels <- sort(unique(na.omit(as.character(real[[var]]))))
  for (syn in dataset_list) {
    all_levels <- union(all_levels, unique(na.omit(as.character(syn[[var]]))))
  }
  
  # First row: variable + p-values
  pval_row <- data.frame(Variable = var_label, check.names = FALSE)
  pval_row[["RATE-AF"]] <- ""
  for (syn_name in names(dataset_list)) {
    syn <- dataset_list[[syn_name]]
    combined <- bind_rows(
      real %>% mutate(dataset = "Real"),
      syn %>% mutate(dataset = "Synthetic")
    )
    pval <- tryCatch({
      kruskal.test(as.formula(paste("as.numeric(", var, ") ~ dataset")), data = combined)$p.value
    }, error = function(e) NA)
    pval_row[[syn_name]] <- if (is.na(pval)) {
      "NA"
    } else {
      formatted <- sprintf("p=%.3f", pval)
      if (pval < 0.05) paste0(formatted, "*") else formatted
    }
    
  }
  output[[1]] <- pval_row
  
  # Rows per level (excluding NA)
  for (lvl in all_levels) {
    row <- data.frame(Variable = paste0("  ", lvl), check.names = FALSE)
    
    # RATE-AF
    count_r <- sum(as.character(real[[var]]) == lvl, na.rm = TRUE)
    total_r <- nrow(real)
    row[["RATE-AF"]] <- sprintf("%d (%.1f%%)", count_r, 100 * count_r / total_r)
    
    # Synthetic datasets
    for (syn_name in names(dataset_list)) {
      syn <- dataset_list[[syn_name]]
      count_s <- sum(as.character(syn[[var]]) == lvl, na.rm = TRUE)
      total_s <- nrow(syn)
      row[[syn_name]] <- sprintf("%d (%.1f%%)", count_s, 100 * count_s / total_s)
    }
    
    output[[length(output) + 1]] <- row
  }
  
  return(bind_rows(output))
}

# Functions to evaluate fidelity in numeric columns
fidelity_evaluation_numeric <- function(var, real, dataset_list,labels) {
  output <- list()
  var_label <- labels[[var]]
  
  # First row: variable + p-values
  p_row <- data.frame(Variable = var_label, check.names = FALSE)
  p_row[["RATE-AF"]] <- ""
  for (syn_name in names(dataset_list)) {
    syn <- dataset_list[[syn_name]]
    combined <- bind_rows(
      real %>% mutate(dataset = "Real"),
      syn %>% mutate(dataset = "Synthetic")
    )
    pval <- tryCatch({
      kruskal.test(as.formula(paste(var, "~ dataset")), data = combined)$p.value
    }, error = function(e) NA)
    
    # Format with asterisk if significant
    p_row[[syn_name]] <- if (is.na(pval)) {
      "NA"
    } else {
      formatted <- sprintf("p=%.3f", pval)
      if (pval < 0.05) paste0(formatted, "*") else formatted
    }
  }
  output[[1]] <- p_row
  
  # Row: Mean [SD]
  mean_row <- data.frame(Variable = "  Mean [SD]", check.names = FALSE)
  mean_row[["RATE-AF"]] <- sprintf("%.1f [%.1f]", mean(real[[var]], na.rm = TRUE), sd(real[[var]], na.rm = TRUE))
  for (syn_name in names(dataset_list)) {
    syn <- dataset_list[[syn_name]]
    mean_row[[syn_name]] <- sprintf("%.1f [%.1f]", mean(syn[[var]], na.rm = TRUE), sd(syn[[var]], na.rm = TRUE))
  }
  output[[2]] <- mean_row
  
  # Row: Median {IQR}
  med_row <- data.frame(Variable = "  Median {IQR}", check.names = FALSE)
  med_row[["RATE-AF"]] <- {
    q <- quantile(real[[var]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    sprintf("%.0f {%.0f–%.0f}", q[2], q[1], q[3])
  }
  for (syn_name in names(dataset_list)) {
    syn <- dataset_list[[syn_name]]
    q <- quantile(syn[[var]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    med_row[[syn_name]] <- sprintf("%.0f {%.0f–%.0f}", q[2], q[1], q[3])
  }
  output[[3]] <- med_row
  
  return(bind_rows(output))
}

# Evaluate fidelity across multiple runs
fidelity_eval_replicates <- function(var, real, replicates, labels, type = "categorical") {
  results <- list()
  
  for (i in seq_along(replicates)) {
    ds <- replicates[[i]]
    # Skip if var not in dataset (safety check)
    if (!(var %in% names(ds))) next
    
    block <- if (type == "categorical") {
      fidelity_evaluation_categorical(var, real, list(Rep = ds), labels)
    } else {
      fidelity_evaluation_numeric(var, real, list(Rep = ds), labels)
    }
    
    results[[i]] <- block[1, ]  # keep only p-value row
  }
  
  df <- bind_rows(results)
  
  # Extract numeric p-values
  df$pval_mean <- suppressWarnings(as.numeric(str_extract(df$Rep, "[0-9\\.]+")))
  avg_p <- mean(df$pval_mean, na.rm = TRUE)
  ci_low  <- quantile(df$pval_mean, probs = 0.025, na.rm = TRUE) # to build ci
  ci_high <- quantile(df$pval_mean, probs = 0.975, na.rm = TRUE) # to build ci
  
  scenario <- attr(replicates, "scenario")
  
  return(data.frame(
    Variable   = labels[[var]],
    Scenario   = scenario,
    pval_mean  = avg_p,
    ci_low     = ci_low,
    ci_high    = ci_high,
    sig_flag   = ifelse(is.na(avg_p), "NA",
                        ifelse(avg_p < 0.05, "Significant", "NS")),
    stringsAsFactors = FALSE
  ))
}

# Aggregate fidelity results across all variables
get_fidelity_pvals <- function(real, dataset_list, vars_categorical, vars_numeric, cat_labels, num_labels) {
  all_tables <- list()
  
  for (scenario in names(dataset_list)) {
    if (scenario == "RATE-AF") next  # skip baseline
    replicates <- dataset_list[[scenario]]
    attr(replicates, "scenario") <- scenario
    
    # Categorical variables
    for (var in vars_categorical) {
      res <- fidelity_eval_replicates(var, real, replicates, cat_labels, type = "categorical")
      all_tables[[length(all_tables)+1]] <- res
    }
    
    # Numeric variables
    for (var in vars_numeric) {
      res <- fidelity_eval_replicates(var, real, replicates, num_labels, type = "numeric")
      all_tables[[length(all_tables)+1]] <- res
    }
  }
  
  bind_rows(all_tables)
}

# Compute baseline (RATE-AF) p-values as reference
get_fidelity_baseline <- function(real, vars_categorical, vars_numeric, cat_labels, num_labels) {
  base_tables <- list()
  
  for (var in vars_categorical) {
    block <- fidelity_evaluation_categorical(var, real, list(), cat_labels)
    base_tables[[length(base_tables)+1]] <- block[1, ]
  }
  for (var in vars_numeric) {
    block <- fidelity_evaluation_numeric(var, real, list(), num_labels)
    base_tables[[length(base_tables)+1]] <- block[1, ]
  }
  
  baseline_table <- bind_rows(base_tables) %>%
    mutate(Baseline = ifelse(grepl("\\*", `RATE-AF`), "Significant", "NS")) %>%
    select(Variable, Baseline)
  
  return(baseline_table)
}

# Helper Functions to summarise categorical results (Plus CI)
fidelity_summary_categorical_avg <- function(var, real, replicates, labels) {
  
  # ------------------------------------------
  # 1. Build per-replicate categorical summaries
  # ------------------------------------------
  block_list <- lapply(replicates, function(ds) {
    fidelity_evaluation_categorical(var, real, list(Rep = ds), labels)
  })
  
  # ------------------------------------------
  # 2. Extract p-values per replicate
  # ------------------------------------------
  pvals <- sapply(block_list, function(b) {
    as.numeric(stringr::str_extract(b[1, "Rep"], "[0-9\\.]+"))
  })
  avg_p <- mean(pvals, na.rm = TRUE)
  ci_p  <- quantile(pvals, c(0.025, 0.975), na.rm = TRUE)
  
  pval_row <- data.frame(
    Variable    = labels[[var]],
    `RATE-AF`   = "p=1.000",
    Replication = sprintf("p=%.3f", avg_p),
    CI          = sprintf("[%.3f–%.3f]", ci_p[1], ci_p[2]),
    check.names = FALSE
  )
  
  # ------------------------------------------
  # 3. Combine per-replicate counts and percentages
  # ------------------------------------------
  all_levels <- bind_rows(lapply(seq_along(block_list), function(i) {
    tbl <- block_list[[i]][-1, ]  # remove p-value row
    tbl$Iteration <- i
    tbl
  }))
  
  # remove leading spaces from Variable
  all_levels$Variable <- trimws(all_levels$Variable)
  
  # extract numeric counts and % for RATE-AF and Replication
  all_levels <- all_levels %>%
    mutate(
      count_r = as.numeric(str_extract(`RATE-AF`, "^[0-9]+")),
      perc_r  = as.numeric(str_extract(`RATE-AF`, "(?<=\\().+?(?=%)")),
      count_s = as.numeric(str_extract(Rep, "^[0-9]+")),
      perc_s  = as.numeric(str_extract(Rep, "(?<=\\().+?(?=%)"))
    )
  
  # ------------------------------------------
  # 4. Compute across-replicate summaries
  # ------------------------------------------
  summary_levels <- all_levels %>%
    group_by(Variable) %>%
    summarise(
      RATE_AF     = sprintf("%.0f (%.1f%%)", mean(count_r, na.rm = TRUE), mean(perc_r, na.rm = TRUE)),
      Replication = sprintf("%.0f (%.1f%%)", mean(count_s, na.rm = TRUE), mean(perc_s, na.rm = TRUE)),
      ci_count_low  = quantile(count_s, 0.025, na.rm = TRUE),
      ci_count_high = quantile(count_s, 0.975, na.rm = TRUE),
      ci_perc_low   = quantile(perc_s,  0.025, na.rm = TRUE),
      ci_perc_high  = quantile(perc_s,  0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      CI = sprintf("[%.0f–%.0f ; %.1f–%.1f%%]", ci_count_low, ci_count_high, ci_perc_low, ci_perc_high)
    ) %>%
    select(Variable, `RATE-AF` = RATE_AF, Replication, CI)
  
  # ------------------------------------------
  # 5. Bind p-value row + level summaries
  # ------------------------------------------
  final <- bind_rows(pval_row, summary_levels)
  return(final)
}

# Helper Functions to summarise numerical results (Plus CI)
fidelity_summary_numeric_avg <- function(var, real, replicates, labels) {
  blocks <- lapply(replicates, function(ds) {
    fidelity_evaluation_numeric(var, real, list(Rep = ds), labels)
  })
  
  # -------- Mean [SD] --------
  mean_rows <- bind_rows(lapply(blocks, function(b) b %>% filter(grepl("Mean\\s*\\[SD\\]", Variable))))
  mean_avg <- data.frame(Variable = "  Mean [SD]", check.names = FALSE)
  mean_avg[["RATE-AF"]] <- unique(mean_rows$`RATE-AF`)[1]
  
  for (col in names(mean_rows)[-1]) {
    vals_mean <- suppressWarnings(as.numeric(stringr::str_extract(mean_rows[[col]], "^[0-9\\.]+")))
    vals_sd   <- suppressWarnings(as.numeric(stringr::str_extract(mean_rows[[col]], "(?<=\\[)[0-9\\.]+")))
    avg_mean <- mean(vals_mean, na.rm = TRUE)
    avg_sd   <- mean(vals_sd,   na.rm = TRUE)
    mean_avg[[col]] <- sprintf("%.1f [%.1f]", avg_mean, avg_sd)
  }
  
  # -------- Median {IQR} --------
  med_rows <- bind_rows(lapply(blocks, function(b) b %>% filter(grepl("Median\\s*\\{IQR\\}", Variable))))
  med_avg <- data.frame(Variable = "  Median {IQR}", check.names = FALSE)
  med_avg[["RATE-AF"]] <- unique(med_rows$`RATE-AF`)[1]
  
  for (col in names(med_rows)[-1]) {
    vals_median <- suppressWarnings(as.numeric(stringr::str_extract(med_rows[[col]], "^[0-9\\.]+")))
    vals_q1     <- suppressWarnings(as.numeric(stringr::str_extract(med_rows[[col]], "(?<=\\{)[0-9\\.]+")))
    vals_q3     <- suppressWarnings(as.numeric(stringr::str_extract(med_rows[[col]], "(?<=–)[0-9\\.]+")))
    avg_median <- mean(vals_median, na.rm = TRUE)
    avg_q1     <- mean(vals_q1,     na.rm = TRUE)
    avg_q3     <- mean(vals_q3,     na.rm = TRUE)
    med_avg[[col]] <- sprintf("%.0f {%.0f–%.0f}", avg_median, avg_q1, avg_q3)
  }
  
  # -------- p-value row (compute CIs from numeric p-values) --------
  pvals <- bind_rows(lapply(blocks, function(b) b[1, ]))
  pval_avg <- data.frame(Variable = labels[[var]], check.names = FALSE)
  pval_avg[["RATE-AF"]] <- ""
  
  # store numeric p-values for later CI extraction
  pval_storage <- list()
  
  for (col in names(pvals)[-1]) {
    vals_p <- suppressWarnings(as.numeric(stringr::str_extract(pvals[[col]], "[0-9\\.]+")))
    avg_p  <- mean(vals_p, na.rm = TRUE)
    pval_storage[[col]] <- vals_p
    pval_avg[[col]] <- if (is.nan(avg_p)) "NA" else sprintf("p=%.3f", avg_p)
  }
  
  # -------- Combine all rows --------
  final <- bind_rows(pval_avg, mean_avg, med_avg)
  rep_col <- setdiff(names(final), c("Variable", "RATE-AF"))[1]
  
  # -------- Add standalone CI column --------
  final <- final %>%
    mutate(
      CI = case_when(
        grepl("Mean", Variable) ~ {
          vals_mean <- suppressWarnings(as.numeric(stringr::str_extract(.data[[rep_col]], "^[0-9\\.]+")))
          ci <- quantile(vals_mean, c(0.025, 0.975), na.rm = TRUE)
          sprintf("[%.1f–%.1f]", ci[1], ci[2])
        },
        grepl("Median", Variable) ~ {
          vals_med <- suppressWarnings(as.numeric(stringr::str_extract(.data[[rep_col]], "^[0-9\\.]+")))
          ci <- quantile(vals_med, c(0.025, 0.975), na.rm = TRUE)
          sprintf("[%.0f–%.0f]", ci[1], ci[2])
        },
        grepl("^p", .data[[rep_col]]) ~ {
          vals_p <- unlist(pval_storage)
          ci <- quantile(vals_p, c(0.025, 0.975), na.rm = TRUE)
          sprintf("[%.3f–%.3f]", ci[1], ci[2])
        },
        TRUE ~ NA_character_
      )
    )
  
  return(final)
}

# Wrapper to build summary averaged table
build_avg_fidelity_tables_replication <- function(rateAF_data,synthetic_results, vars_categorical, vars_numeric, categorical_labels, numeric_labels) {
  
  
  #  Helper cleaner 
  clean_fidelity_table <- function(tbl) {
    rate_cols <- grep("^RATE-AF", names(tbl), value = TRUE)
    if (length(rate_cols) > 1) {
      tbl <- tbl %>%
        select(-all_of(rate_cols[-1])) %>%
        rename(`RATE-AF` = all_of(rate_cols[1]))
    }
    tbl[] <- lapply(tbl, function(col) {
      if (is.character(col)) {
        col <- gsub("^NA$", "p=1.000", col)
        col <- gsub("NaN", "p=1.000", col)
        col <- ifelse(col == "", "p=1.000", col)
      }
      return(col)
    })
    return(tbl)
  }
  
  
  #  Build dataset list (Baseline + Replication)
  dataset_list <- list(
    `RATE-AF`    = rateAF_data,
    `Replication` = synthetic_results
  )
  
  
  #  Categorical variables
  categorical_table <- list()
  for (var in vars_categorical) {
    tbl <- fidelity_summary_categorical_avg(var, rateAF_data,
                                            dataset_list[["Replication"]],
                                            categorical_labels)
    names(tbl)[names(tbl) == "Rep"] <- "Replication"
    tbl <- clean_fidelity_table(tbl)
    categorical_table[[length(categorical_table) + 1]] <- tbl
  }
  fidelity_categorical_summary <- bind_rows(categorical_table)
  
  
  #  Numeric variables
  numeric_table <- list()
  for (var in vars_numeric) {
    tbl <- fidelity_summary_numeric_avg(var, rateAF_data,
                                        dataset_list[["Replication"]],
                                        numeric_labels)
    names(tbl)[names(tbl) == "Rep"] <- "Replication"
    tbl <- clean_fidelity_table(tbl)
    numeric_table[[length(numeric_table) + 1]] <- tbl
  }
  fidelity_numeric_summary <- bind_rows(numeric_table)
  
  
  #  Return final averaged tables
  return(list(
    categorical = fidelity_categorical_summary,
    numeric     = fidelity_numeric_summary
  ))
}


# =====================================================================
# Additional Fidelity Summary Functions (Controlled Boosting)
# =====================================================================

# Average categorical fidelity across replicates
fidelity_summary_categorical_avg <- function(var, real, replicates, labels) {
  blocks <- lapply(replicates, function(ds) {
    fidelity_evaluation_categorical(var, real, list(Rep = ds), labels)
  })
  
  # --- average p-values ---
  pvals <- bind_rows(lapply(blocks, function(b) b[1, ]))
  pval_avg <- data.frame(Variable = labels[[var]], check.names = FALSE)
  pval_avg[["RATE-AF"]] <- ""
  for (col in names(pvals)[-1]) {
    vals <- suppressWarnings(as.numeric(str_extract(pvals[[col]], "[0-9\\.]+")))
    avg_p <- mean(vals, na.rm = TRUE)
    label <- if (is.nan(avg_p)) "NA" else sprintf("p=%.3f%s", avg_p, ifelse(avg_p < 0.05, "*", ""))
    pval_avg[[col]] <- label
  }
  
  # --- average counts/percentages ---
  combined <- bind_rows(lapply(blocks, function(b) b[-1, , drop = FALSE]))
  level_rows <- combined %>%
    filter(grepl("^\\s", Variable)) %>%
    pivot_longer(-Variable, names_to = "Scenario", values_to = "value") %>%
    mutate(
      count = suppressWarnings(as.numeric(str_extract(value, "^[0-9]+"))),
      perc  = suppressWarnings(as.numeric(str_extract(value, "(?<=\\().+?(?=%)")))
    ) %>%
    group_by(Variable, Scenario) %>%
    summarise(
      mean_count = mean(count, na.rm = TRUE),
      mean_perc  = mean(perc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(value = sprintf("%.0f (%.1f%%)", mean_count, mean_perc)) %>%
    select(Variable, Scenario, value) %>%
    pivot_wider(names_from = Scenario, values_from = value)
  
  final <- bind_rows(pval_avg, level_rows)
  return(final)
}

# Average numeric fidelity across replicates
fidelity_summary_numeric_avg <- function(var, real, replicates, labels) {
  blocks <- lapply(replicates, function(ds) {
    fidelity_evaluation_numeric(var, real, list(Rep = ds), labels)
  })
  
  # Average "Mean [SD]" row
  mean_rows <- bind_rows(lapply(blocks, function(b) b %>% filter(grepl("Mean", Variable))))
  mean_avg <- data.frame(Variable = "  Mean [SD]", check.names = FALSE)
  mean_avg[["RATE-AF"]] <- unique(mean_rows$`RATE-AF`)[1]
  for (col in names(mean_rows)[-1]) {
    vals_m <- suppressWarnings(as.numeric(str_extract(mean_rows[[col]], "^[0-9\\.]+")))
    vals_s <- suppressWarnings(as.numeric(str_extract(mean_rows[[col]], "(?<=\\[)[0-9\\.]+")))
    avg_m <- mean(vals_m, na.rm = TRUE)
    avg_s <- mean(vals_s, na.rm = TRUE)
    mean_avg[[col]] <- sprintf("%.1f [%.1f]", avg_m, avg_s)
  }
  
  # Average "Median {IQR}" row
  med_rows <- bind_rows(lapply(blocks, function(b) b %>% filter(grepl("Median", Variable))))
  med_avg <- data.frame(Variable = "  Median {IQR}", check.names = FALSE)
  med_avg[["RATE-AF"]] <- unique(med_rows$`RATE-AF`)[1]
  for (col in names(med_rows)[-1]) {
    vals_m <- suppressWarnings(as.numeric(str_extract(med_rows[[col]], "^[0-9\\.]+")))
    vals_q1 <- suppressWarnings(as.numeric(str_extract(med_rows[[col]], "(?<=\\{)[0-9\\.]+")))
    vals_q3 <- suppressWarnings(as.numeric(str_extract(med_rows[[col]], "(?<=–)[0-9\\.]+")))
    avg_m <- mean(vals_m, na.rm = TRUE)
    avg_q1 <- mean(vals_q1, na.rm = TRUE)
    avg_q3 <- mean(vals_q3, na.rm = TRUE)
    med_avg[[col]] <- sprintf("%.0f {%.0f–%.0f}", avg_m, avg_q1, avg_q3)
  }
  
  # p-value row
  pvals <- bind_rows(lapply(blocks, function(b) b[1, ]))
  pval_avg <- data.frame(Variable = labels[[var]], check.names = FALSE)
  pval_avg[["RATE-AF"]] <- ""
  for (col in names(pvals)[-1]) {
    vals <- suppressWarnings(as.numeric(str_extract(pvals[[col]], "[0-9\\.]+")))
    avg_p <- mean(vals, na.rm = TRUE)
    label <- if (is.nan(avg_p)) "NA" else sprintf("p=%.3f%s", avg_p, ifelse(avg_p < 0.05, "*", ""))
    pval_avg[[col]] <- label
  }
  
  # Combine
  final <- bind_rows(pval_avg, mean_avg, med_avg)
  return(final)
}

# Wrapper to build full averaged tables across 30/60/90% for one ethnic group
build_avg_fidelity_tables_boosting <- function(rateAF_data, boosting_results, group_prefix,vars_categorical, vars_numeric,categorical_labels, numeric_labels) {
  
  # ------------------------------------------------------------
  #  Helper: internal table cleaner (remove duplicates, fix NA)
  # ------------------------------------------------------------
  clean_fidelity_table <- function(tbl) {
    # Remove duplicate RATE-AF columns
    rate_cols <- grep("^RATE-AF", names(tbl), value = TRUE)
    if (length(rate_cols) > 1) {
      tbl <- tbl %>%
        select(-all_of(rate_cols[-1])) %>%
        rename(`RATE-AF` = all_of(rate_cols[1]))
    }
    
    # Replace NA/NaN/blank p-values with p=1.000
    tbl[] <- lapply(tbl, function(col) {
      if (is.character(col)) {
        col <- gsub("^NA$", "p=1.000", col)
        col <- gsub("NaN", "p=1.000", col)
        col <- ifelse(col == "", "p=1.000", col)
      }
      return(col)
    })
    
    return(tbl)
  }
  
  # ------------------------------------------------------------
  #  Build dataset list (RATE-AF + three boosted scenarios)
  # ------------------------------------------------------------
  dataset_list <- list(
    rateAF_data,
    boosting_results[[paste0(group_prefix, "_30%")]],
    boosting_results[[paste0(group_prefix, "_60%")]],
    boosting_results[[paste0(group_prefix, "_80%")]]
  )
  names(dataset_list) <- c(
    "RATE-AF",
    paste0("Boost_", substr(group_prefix, 1, 1), "_30"),
    paste0("Boost_", substr(group_prefix, 1, 1), "_60"),
    paste0("Boost_", substr(group_prefix, 1, 1), "_90")
  )
  
  # ------------------------------------------------------------
  #  Categorical variables
  # ------------------------------------------------------------
  categorical_table <- list()
  for (var in vars_categorical) {
    avg_tables <- list()
    for (scen in names(dataset_list)[-1]) {  # skip RATE-AF
      tbl <- fidelity_summary_categorical_avg(var, rateAF_data, dataset_list[[scen]], categorical_labels)
      names(tbl)[names(tbl) == "Rep"] <- scen
      avg_tables[[scen]] <- tbl
    }
    merged <- Reduce(function(x, y) full_join(x, y, by = "Variable"), avg_tables)
    merged <- clean_fidelity_table(merged)
    categorical_table[[length(categorical_table) + 1]] <- merged
  }
  fidelity_categorical_summary <- bind_rows(categorical_table)
  
  # ------------------------------------------------------------
  #  Numeric variables
  # ------------------------------------------------------------
  numeric_table <- list()
  for (var in vars_numeric) {
    avg_tables <- list()
    for (scen in names(dataset_list)[-1]) {
      tbl <- fidelity_summary_numeric_avg(var, rateAF_data, dataset_list[[scen]], numeric_labels)
      names(tbl)[names(tbl) == "Rep"] <- scen
      avg_tables[[scen]] <- tbl
    }
    merged <- Reduce(function(x, y) full_join(x, y, by = "Variable"), avg_tables)
    merged <- clean_fidelity_table(merged)
    numeric_table[[length(numeric_table) + 1]] <- merged
  }
  fidelity_numeric_summary <- bind_rows(numeric_table)
  
  # ------------------------------------------------------------
  #  Return final cleaned summaries
  # ------------------------------------------------------------
  return(list(
    categorical = fidelity_categorical_summary,
    numeric = fidelity_numeric_summary
  ))
}




# =====================================================================
# Subgroup Fidelity Functions (used in Controlled Boosting)
# =====================================================================

# Numeric Summary Function
summarise_numeric_by_subgroup <- function(data, vars, subgroup_var = "ethncty") {
    
    summary_list <- list()
    
    # Mapping for ethnicity codes to labels
    ethnicity_labels <- c("1" = "White", "2" = "Asian", "3" = "Black")
    data[[subgroup_var]] <- recode(as.character(data[[subgroup_var]]), !!!ethnicity_labels)
    
    for (var in vars) {
      # Compute mean (sd) per subgroup
      mean_row <- data %>%
        group_by(!!sym(subgroup_var)) %>%
        summarise(value = sprintf("%.1f (%.1f)", mean(.data[[var]], na.rm = TRUE), sd(.data[[var]], na.rm = TRUE)), .groups = "drop") %>%
        pivot_wider(names_from = !!sym(subgroup_var), values_from = value) %>%
        mutate(Variable = "  mean (sd)")
      
      # Compute median (IQR) per subgroup
      median_row <- data %>%
        group_by(!!sym(subgroup_var)) %>%
        summarise(med = median(.data[[var]], na.rm = TRUE),
                  q1 = quantile(.data[[var]], probs = 0.25, na.rm = TRUE),
                  q3 = quantile(.data[[var]], probs = 0.75, na.rm = TRUE), .groups = "drop") %>%
        mutate(value = sprintf("%.0f (%.0f - %.0f)", med, q1, q3)) %>%
        select(!!sym(subgroup_var), value) %>%
        pivot_wider(names_from = !!sym(subgroup_var), values_from = value) %>%
        mutate(Variable = "  median (IQR)")
      
      # Add overall values
      overall_mean <- sprintf("%.1f (%.1f)", mean(data[[var]], na.rm = TRUE), sd(data[[var]], na.rm = TRUE))
      overall_median <- {
        q <- quantile(data[[var]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
        sprintf("%.0f (%.0f - %.0f)", q[2], q[1], q[3])
      }
      mean_row$Overall <- overall_mean
      median_row$Overall <- overall_median
      
      # Add Variable header row
      header_row <- tibble(Variable = var)
      for (col in setdiff(names(mean_row), "Variable")) {
        header_row[[col]] <- NA_character_
      }
      
      # Bind rows
      summary_list[[length(summary_list) + 1]] <- bind_rows(header_row, mean_row, median_row)
    }
    
    bind_rows(summary_list) %>%
      relocate(Variable)
  }
  
# Categorical Summary Function
summarise_categorical_by_subgroup <- function(data, vars, subgroup_var = "ethncty") {
    summary_list <- list()
    # Mapping for ethnicity codes to labels
    ethnicity_labels <- c("1" = "White", "2" = "Asian", "3" = "Black")
    data[[subgroup_var]] <- recode(as.character(data[[subgroup_var]]), !!!ethnicity_labels)
    
    for (var in vars) {
      cat_data <- data %>%
        group_by(across(all_of(c(subgroup_var, var)))) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(across(all_of(subgroup_var))) %>%
        mutate(pct = n / sum(n) * 100) %>%
        ungroup() %>%
        mutate(value = sprintf("%d (%.1f%%)", n, pct)) %>%
        select(group = all_of(subgroup_var), level = all_of(var), value)
      
      # Complete combinations with 0 (0%) if missing
      all_groups <- unique(data[[subgroup_var]])
      all_levels <- unique(data[[var]])
      complete_grid <- expand.grid(group = all_groups, level = all_levels, stringsAsFactors = FALSE)
      
      cat_data <- complete_grid %>%
        left_join(cat_data, by = c("group", "level")) %>%
        mutate(value = ifelse(is.na(value), "0 (0.0%)", value))
      
      wide_data <- cat_data %>%
        pivot_wider(names_from = group, values_from = value, names_repair = "unique") %>%
        mutate(Variable = paste0("  ", as.character(level))) %>%
        select(-level) %>%
        relocate(Variable)
      
      # Overall counts
      overall_counts <- data %>%
        group_by(across(all_of(var))) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(pct = n / sum(n) * 100,
               value = sprintf("%d (%.1f%%)", n, pct)) %>%
        rename(level = !!sym(var))
      
      wide_data <- wide_data %>%
        left_join(overall_counts %>% mutate(Variable = paste0("  ", level)) %>% select(Variable, Overall = value), by = "Variable")
      
      header_row <- tibble(Variable = var)
      for (col in setdiff(names(wide_data), "Variable")) {
        header_row[[col]] <- NA_character_
      }
      
      summary_list[[length(summary_list) + 1]] <- bind_rows(header_row, wide_data)
    }
    
    bind_rows(summary_list)
  }
  
# Combine Subgroup Summaries
combine_subgroup_summaries <- function(rateaf_df, syn_df,
                                         groups = c("White", "Asian", "Black", "Overall")) {
    rateaf_df <- rateaf_df %>% mutate(row_id = dplyr::row_number())
    syn_df    <- syn_df %>% mutate(row_id = dplyr::row_number())
    
    joined <- dplyr::full_join(rateaf_df, syn_df, by = "row_id", suffix = c("_real", "_syn"))
    
    out <- tibble::tibble(Variable = dplyr::coalesce(joined$Variable_real, joined$Variable_syn))
    for (g in groups) {
      out[[paste0(g, " (RATE-AF)")]]   <- joined[[paste0(g, "_real")]]
      out[[paste0(g, " (AvgBoosted)")]] <- joined[[paste0(g, "_syn")]]
    }
    out
  }
  
# Force categorical header rows to be blank (first row under each categorical var should be empty)
blank_categorical_headers <- function(df, cat_vars,
                                        groups = c("White","Asian","Black","Overall")) {
    hdr <- df$Variable %in% cat_vars                   # rows that are exactly the variable name
    cols <- unlist(lapply(groups, function(g) c(paste0(g," (RATE-AF)"),
                                                paste0(g," (AvgBoosted)"))))
    df[hdr, cols] <- NA
    df
  }
  
summarise_numeric_by_subgroup_averaged <- function(boosted_list, vars, subgroup_var = "ethncty") {
  summary_list <- list()
  ethnicity_labels <- c("1" = "White", "2" = "Asian", "3" = "Black")
  
  for (var in vars) {
    
    
    # For each dataset, compute subgroup stats
    per_rep <- lapply(boosted_list, function(ds) {
      ds[[subgroup_var]] <- recode(as.character(ds[[subgroup_var]]), !!!ethnicity_labels)
      
      ds %>%
        group_by(!!sym(subgroup_var)) %>%
        summarise(
          mean_val = mean(.data[[var]], na.rm = TRUE),
          sd_val   = sd(.data[[var]], na.rm = TRUE),
          med_val  = median(.data[[var]], na.rm = TRUE),
          q1_val   = quantile(.data[[var]], 0.25, na.rm = TRUE),
          q3_val   = quantile(.data[[var]], 0.75, na.rm = TRUE),
          .groups  = "drop"
        ) %>%
        mutate(
          Overall_mean = mean(ds[[var]], na.rm = TRUE),
          Overall_sd   = sd(ds[[var]], na.rm = TRUE),
          Overall_med  = median(ds[[var]], na.rm = TRUE),
          Overall_q1   = quantile(ds[[var]], 0.25, na.rm = TRUE),
          Overall_q3   = quantile(ds[[var]], 0.75, na.rm = TRUE)
        )
    })
    
    # Combine and average across iterations
    combined <- bind_rows(per_rep) %>%
      group_by(!!sym(subgroup_var)) %>%
      summarise(
        mean_val = mean(.data[["mean_val"]], na.rm = TRUE),
        sd_val   = mean(.data[["sd_val"]],   na.rm = TRUE),
        med_val  = mean(.data[["med_val"]],  na.rm = TRUE),
        q1_val   = mean(.data[["q1_val"]],   na.rm = TRUE),
        q3_val   = mean(.data[["q3_val"]],   na.rm = TRUE),
        Overall_mean = mean(.data[["Overall_mean"]], na.rm = TRUE),
        Overall_sd   = mean(.data[["Overall_sd"]],   na.rm = TRUE),
        Overall_med  = mean(.data[["Overall_med"]],  na.rm = TRUE),
        Overall_q1   = mean(.data[["Overall_q1"]],   na.rm = TRUE),
        Overall_q3   = mean(.data[["Overall_q3"]],   na.rm = TRUE),
        .groups = "drop"
      )
    
    # Format exactly like original function
    mean_row <- combined %>%
      select(!!sym(subgroup_var), mean_val, sd_val) %>%
      mutate(value = sprintf("%.1f (%.1f)", mean_val, sd_val)) %>%
      select(-mean_val, -sd_val) %>%
      pivot_wider(names_from = !!sym(subgroup_var), values_from = value) %>%
      mutate(
        Variable = "  mean (sd)",
        Overall  = sprintf("%.1f (%.1f)",
                           first(combined$Overall_mean),
                           first(combined$Overall_sd))
      )
    
    median_row <- combined %>%
      select(!!sym(subgroup_var), med_val, q1_val, q3_val) %>%
      mutate(value = sprintf("%.0f (%.0f - %.0f)", med_val, q1_val, q3_val)) %>%
      select(-med_val, -q1_val, -q3_val) %>%
      pivot_wider(names_from = !!sym(subgroup_var), values_from = value) %>%
      mutate(
        Variable = "  median (IQR)",
        Overall  = sprintf("%.0f (%.0f - %.0f)",
                           first(combined$Overall_med),
                           first(combined$Overall_q1),
                           first(combined$Overall_q3))
      )
    
    # Header row
    header_row <- tibble(Variable = var)
    for (col in setdiff(names(mean_row), "Variable")) {
      header_row[[col]] <- NA_character_
    }
    
    summary_list[[length(summary_list) + 1]] <- bind_rows(header_row, mean_row, median_row)
  }
  
  bind_rows(summary_list) %>% relocate(Variable)
}

summarise_categorical_by_subgroup_averaged <- function(boosted_list, vars, subgroup_var = "ethncty") {
  summary_list <- list()
  ethnicity_labels <- c("1" = "White", "2" = "Asian", "3" = "Black")
  
  for (var in vars) {
    
    # --- Preserve reference level order from RATE-AF if available ---
    if (exists("rateAF_data") && var %in% names(rateAF_data)) {
      ref_levels <- unique(na.omit(as.character(rateAF_data[[var]])))
    } else {
      ref_levels <- NULL
    }
    
    # ---- Collect subgroup summaries from all replicates ----
    per_rep <- lapply(boosted_list, function(ds) {
      ds[[subgroup_var]] <- recode(as.character(ds[[subgroup_var]]), !!!ethnicity_labels)
      
      ds %>%
        group_by(across(all_of(c(subgroup_var, var)))) %>%
        summarise(n = dplyr::n(), .groups = "drop") %>%
        group_by(across(all_of(subgroup_var))) %>%
        mutate(pct = 100 * n / sum(n, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(
          group = as.character(.data[[subgroup_var]]),
          level = as.character(.data[[var]]),
          n = as.numeric(n),
          pct = as.numeric(pct)
        ) %>%
        select(group, level, n, pct)
    })
    
    # ---- Combine and average across replicates ----
    combined <- dplyr::bind_rows(per_rep) %>%
      mutate(
        n = as.numeric(n),
        pct = as.numeric(pct)
      ) %>%
      group_by(group, level) %>%
      summarise(
        n_avg = mean(n, na.rm = TRUE),
        pct_avg = mean(pct, na.rm = TRUE),
        .groups = "drop"
      )
    
    # ---- Ensure full grid of categories with correct order ----
    all_groups <- c("White", "Asian", "Black")
    if (!is.null(ref_levels)) {
      all_levels <- ref_levels
    } else {
      all_levels <- sort(unique(combined$level))
    }
    
    grid <- expand.grid(group = all_groups, level = all_levels, stringsAsFactors = FALSE)
    
    combined <- grid %>%
      left_join(combined, by = c("group", "level")) %>%
      mutate(
        n_avg = as.numeric(ifelse(is.na(n_avg), 0, n_avg)),
        pct_avg = as.numeric(ifelse(is.na(pct_avg), 0, pct_avg))
      )
    
    # ---- Create display values ----
    combined <- combined %>%
      mutate(value = sprintf("%.0f (%.1f%%)", n_avg, pct_avg))
    
    # ---- Pivot to wide ----
    wide_data <- combined %>%
      select(group, level, value) %>%
      tidyr::pivot_wider(names_from = group, values_from = value) %>%
      mutate(Variable = paste0("  ", level)) %>%
      select(Variable, everything(), -level)
    
    # ---- Overall averages ----
    overall <- combined %>%
      group_by(level) %>%
      summarise(
        n_avg_total = sum(n_avg, na.rm = TRUE),
        .groups = "drop"
      )
    
    # recompute overall percentages relative to total sample size
    total_all <- sum(overall$n_avg_total, na.rm = TRUE)
    
    overall <- overall %>%
      mutate(
        pct_avg_total = 100 * n_avg_total / total_all,
        Variable = paste0("  ", level),
        Overall = sprintf("%.0f (%.1f%%)", n_avg_total, pct_avg_total)
      ) %>%
      select(Variable, Overall)
    
    wide_data <- left_join(wide_data, overall, by = "Variable")
    
    # ---- Header ----
    header <- tibble(Variable = var)
    for (col in setdiff(names(wide_data), "Variable")) header[[col]] <- NA_character_
    
    summary_list[[length(summary_list) + 1]] <- bind_rows(header, wide_data)
  }
  
  dplyr::bind_rows(summary_list)
}

