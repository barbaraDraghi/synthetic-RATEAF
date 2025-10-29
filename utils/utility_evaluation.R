# =====================================================================
# Utility Evaluation Functions
# =====================================================================
# These functions evaluate model-based utility (treatment effect
# preservation) between real and synthetic datasets — primarily used
# in the Replication experiment, but extendable to controlled or
# exploratory setups.
# =====================================================================

# Define outcome models used in RATE-AF
outcome_models <- list(
  #nyha12m = list(
  #  type = "linear",
  #  formula = nyha12m ~ alloc + sex + age + nyha + mdfdehrascr + estmtdejctnfrctn
  #),
  nyhaimprove = list(
    type = "logistic",
    formula = nyhaimprove ~ alloc + sex + age + estmtdejctnfrctn + nyha
  ),
  mehra2class12m = list(
    type = "logistic",
    formula = mehra2class12m ~ alloc + sex + age + estmtdejctnfrctn
  ),
  ecghr12m = list(
    type = "linear",
    formula = ecghr12m ~ alloc + sex + age + ecghr + mdfdehrascr + estmtdejctnfrctn
  ),
  ntprobnp12m_ln = list(
    type = "linear",
    formula = ntprobnp12m_ln ~ alloc + sex + age + ntprobnp_ln + mdfdehrascr + estmtdejctnfrctn,
    log_exp = TRUE
  ),
  ae_yn = list(
    type = "logistic",
    formula = ae_yn ~ alloc + sex + age + estmtdejctnfrctn
  ),
  cvadj = list(
    type = "logistic",
    formula = cvadj ~ alloc + sex + age + estmtdejctnfrctn
  )
)
  

# Core utility evaluation (single dataset or list of datasets)
utility_evaluation <- function(outcome_name, model_spec, dataset_list) {
    results <- list()
    
    for (dname in names(dataset_list)) {
      data <- dataset_list[[dname]]
      
      # --- Pre-check for variation (logistic only) ---
      if (model_spec$type == "logistic") {
        y <- data[[as.character(model_spec$formula[[2]])]]
        if (all(is.na(y)) || length(unique(y[!is.na(y)])) < 2) {
          message(sprintf("⚠️ %s has no variation in %s — skipping", outcome_name, dname))
          next
        }
      }
      
      # Check if model can be fit safely
      model <- tryCatch({
        if (model_spec$type == "logistic") {
          # Detect complete separation (any 0 cells)
          tab <- table(data[[as.character(model_spec$formula[[2]])]], data$alloc)
          if (any(tab == 0)) {
            message(sprintf("⚠️ Separation detected for %s in %s — switching to Firth regression",
                            outcome_name, dname))
            if (!requireNamespace("logistf", quietly = TRUE)) install.packages("logistf")
            logistf::logistf(model_spec$formula, data = data)
          } else {
            glm(model_spec$formula, data = data, family = binomial)
          }
        } else {
          lm(model_spec$formula, data = data)
        }
      }, error = function(e) NULL)
      
      if (!is.null(model)) {
        result <- tryCatch({
          broom::tidy(model, conf.int = TRUE)
        }, error = function(e) NULL)
        
        if (!is.null(result)) {
          result <- result %>%
            dplyr::filter(grepl("^alloc", term)) %>%
            dplyr::mutate(
              estimate  = ifelse(model_spec$type == "logistic" || isTRUE(model_spec$log_exp), exp(estimate), estimate),
              conf.low  = ifelse(model_spec$type == "logistic" || isTRUE(model_spec$log_exp), exp(conf.low), conf.low),
              conf.high = ifelse(model_spec$type == "logistic" || isTRUE(model_spec$log_exp), exp(conf.high), conf.high),
              label = sprintf("%.2f (95%% CI %.2f–%.2f; p=%s%s)",
                              estimate, conf.low, conf.high,
                              ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
                              ifelse(p.value < 0.05, " *", "")),
              stat = dplyr::case_when(
                model_spec$type == "logistic" ~ "OR",
                model_spec$type == "linear" & outcome_name == "ntprobnp12m_ln" ~ "RGM",
                TRUE ~ "AMD"
              ),
              outcome = outcome_name,
              dataset = dname,
              sig_flag = ifelse(p.value < 0.05, "p < 0.05", "NS")
            ) %>%
            dplyr::select(outcome, stat, dataset, estimate, conf.low, conf.high, p.value, sig_flag, label)
          results[[dname]] <- result
        }
      }
    }
    
    dplyr::bind_rows(results)
  }
  
  
# Average utility results across replicates
utility_eval_replicates <- function(outcome_name, model_spec, replicates) {
  
  results <- purrr::map(replicates, ~{
    tryCatch(
      utility_evaluation(outcome_name, model_spec, list(Rep = .x)),
      error = function(e) NULL
    )
  })
  
  df <- dplyr::bind_rows(results)
  
  # --- Handle completely missing outcome (no valid replicate) ---
  if (nrow(df) == 0) {
    message(sprintf("⚠️ No valid replicates for %s (%s) — inserting NA placeholder", 
                    outcome_name, attr(replicates, "scenario")))
    return(data.frame(
      outcome = outcome_name,
      stat = ifelse(model_spec$type == "logistic", "OR", "AMD"),
      dataset = attr(replicates, "scenario"),
      estimate = NA, conf.low = NA, conf.high = NA, p.value = NA,
      sig_flag = "NA",
      label = "Not computed"
    ))
  }
  
  # --- Robust filtering ---
  n_before <- nrow(df)
  df <- df %>%
    dplyr::filter(
      is.finite(estimate),
      is.finite(conf.low),
      is.finite(conf.high),
      abs(estimate) < 1e6,
      abs(conf.low) < 1e6,
      abs(conf.high) < 1e6
    )
  n_after <- nrow(df)
  if (n_after < n_before)
    message(sprintf("⚠️ Removed %d unstable rows in %s (%s)", 
                    n_before - n_after, outcome_name, attr(replicates, "scenario")))
  
  # --- Average across stable replicates ---
  df_summary <- df %>%
    dplyr::group_by(outcome, stat) %>%
    dplyr::summarise(
      estimate  = mean(estimate, na.rm = TRUE),
      conf.low  = mean(conf.low, na.rm = TRUE),
      conf.high = mean(conf.high, na.rm = TRUE),
      p.value   = mean(p.value, na.rm = TRUE),
      n_used    = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      dataset = attr(replicates, "scenario"),
      sig_flag = ifelse(p.value < 0.05, "p < 0.05", "NS"),
      label = sprintf("%.2f (95%% CI %.2f–%.2f; p=%s%s)",
                      estimate, conf.low, conf.high,
                      ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
                      ifelse(p.value < 0.05, " *", ""))
    )
  
  return(df_summary)
}
  
# Convert NYHA to numeric (helper)
convert_nyha_numeric <- function(obj) {
  if (is.data.frame(obj)) {
    if ("nyha12m" %in% names(obj))
      obj$nyha12m <- suppressWarnings(as.numeric(as.character(obj$nyha12m)))
    return(obj)
  }
  if (is.list(obj)) return(lapply(obj, convert_nyha_numeric))
  obj
}
  
# Forest plot function (replication experiment)
forest_plot <- function(data, outcome_label) {
  df <- data %>% dplyr::filter(outcome == outcome_label)
  
  if (nrow(df) == 0) {
    warning(sprintf("⚠️ No data found for outcome '%s' in forest_plot()", outcome_label))
    return(NULL)
  }
  
  # --- Ensure numeric values are available ---
  if (!all(c("estimate", "conf.low", "conf.high") %in% names(df))) {
    df <- df %>%
      mutate(
        estimate  = as.numeric(str_extract(label, "^[0-9.\\-]+")),
        conf.low  = as.numeric(str_extract(label, "(?<=95% CI )[-0-9.]+")),
        conf.high = as.numeric(str_extract(label, "(?<=–)[-0-9.]+")),
        sig_flag  = ifelse(str_detect(label, "\\*"), "p < 0.05", "NS")
      )
  }
  
  # --- Choose axis label and null effect line ---
  stat_type <- unique(df$stat)
  ref_line <- ifelse(stat_type %in% c("OR", "RGM"), 1, 0)
  x_label <- switch(
    stat_type,
    "OR"  = "Odds Ratio (Digoxin vs Beta-blocker)",
    "RGM" = "Geometric Mean Ratio (Digoxin vs Beta-blocker)",
    "AMD" = "Adjusted Mean Difference (Digoxin vs Beta-blocker)",
    "Effect Estimate"
  )
  
  # --- Define colour palette ---
  palette_sig <- c("p < 0.05" = "#1f78b4", "NS" = "black")
  
  # --- Plot ---
  ggplot(df, aes(x = estimate, y = reorder(dataset, desc(dataset)))) +
    geom_vline(aes(xintercept = ref_line), color = "red", linetype = "dashed", linewidth = 1) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = sig_flag),
                    size = 0.9, fatten = 2) +
    geom_text(aes(label = sprintf("%.2f [%.2f–%.2f]", estimate, conf.low, conf.high)),
              vjust = -1, size = 3.3, color = "black") +
    scale_color_manual(values = palette_sig) +
    labs(
      title = outcome_label,
      x = x_label,
      y = NULL,
      color = "Significance"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.line.x = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.text.y = element_text(face = "bold"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
}
  

# Function to compute utility within a subgroup
utility_by_subgroup <- function(outcome_name, model_spec, dataset_list, subgroup, group_var = "ethncty") {
  message(sprintf("\n▶ Subgroup %s — outcome: %s", subgroup, outcome_name))
  
  # Filter each dataset by subgroup before computing utility
  subgroup_datasets <- lapply(dataset_list, function(df) {
    df %>% dplyr::filter(.data[[group_var]] == subgroup)
  })
  
  # Run the standard utility_evaluation
  utility_evaluation(outcome_name, model_spec, subgroup_datasets) %>%
    dplyr::mutate(Subgroup = ethnic_labels[[as.character(subgroup)]])
}

