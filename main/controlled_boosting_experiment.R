set.seed(12345)
n_iterations <- 500

# -------------------------------------------------------
# Scripts, libraries and functions loading
# -------------------------------------------------------

setwd("H:/partTimePhD/Code/synthetic-RATEAF/utils")
source('functions.R')
source('bias_analysis.R')
source('synthetic_data_generation.R')
source('fidelity_evaluation.R')
source('utility_evaluation.R')

required_packages <- c(
  'dplyr', 'ggplot2', 'tidyr', 'bnlearn', 'purrr', 
  'parallel', 'stringr', 'mice', 'binom', 
  'ranger', 'pROC', 'PRROC', 'e1071',
  'readxl','haven',
  'fairness','caret', 'patchwork'
)
check_and_install_packages(required_packages)


# --------------------------------------------------------------
#  RATE AF - Data loading and preparation
# --------------------------------------------------------------

# Assuming having access to RATE-AF.dta data in data/ folder
rateAF_data <- rateAF_datapreparation_step1(data_path = "../data/RATE-AF.dta") 
rateAF_data_processed <- rateAF_datapreparation_step2(rateAF_data)
imputed_data <- mice(rateAF_data_processed, method = "pmm", m = 5, maxit = 50, seed = 500)
rateAF_data_processed <- complete(imputed_data)
rateAF_data_processed <- discretise_RATEAF_adaptive(rateAF_data_processed, default_quantiles = 8, max_bins = 10)

# ---------------------------------------------------------------------
#  Bayesian Network Learning - White and Black Lists Specification
# ---------------------------------------------------------------------


whitelist <- data.frame(
  from = c(
    "alloc", 
    "ntprobnp_discretised", 
    "nyha", 
    "ecghr_discretised", 
    "mehra2class12m"
  ),
  to = c(
    "ntprobnp12m_discretised", 
    "ntprobnp12m_discretised", 
    "ntprobnp12m_discretised",
    "ecghr12m_discretised",
    "cvadj"
  )
)

blacklist <- data.frame(
  from = c(
    "ntprobnp12m_discretised", "nyha12m", "mehra2class12m", "cvadj"
  ),
  to = c(
    "alloc", "alloc", "alloc", "alloc"
  )
)

# --------------------------------------------------------------
#   Synthetic Data Generation
# --------------------------------------------------------------


# Run Boosting Scenarios (30%, 60%, 90%)
drop_levels <- c(0.3, 0.6, 0.8)
groups <- c("1" = "White", "2" = "Asian", "3" = "Black")

boosting_results <- list()
scenario_summary <- list()

for (g in names(groups)) {
  for (d in drop_levels) {
    key <- paste0(groups[g], "_", d*100, "%")
    boosting_results[[key]] <- boosting_replicates(
      data = rateAF_data_processed,
      group_value = g,
      drop_frac = d,
      group_var = "ethncty",
      max_iter = 5,
      n_reps = n_iterations,
      whitelist = whitelist,
      blacklist = blacklist
    )
    
    # ---- Record scenario details ----
    n_total <- nrow(rateAF_data_processed)
    n_group <- sum(rateAF_data_processed$ethncty == g)
    n_removed <- round(n_group * d)
    n_remaining <- n_group - n_removed
    n_synth <- n_removed
    n_biased <- n_total - n_removed
    n_boosted <- n_total  # total after restoring full cohort
    
    scenario_summary[[key]] <- data.frame(
      subgroup = groups[g],
      bias_level = paste0(d * 100, "%"),
      removed = n_removed,
      remaining_real = n_remaining,
      synthetic_added = n_synth,
      total_biased = n_biased,
      total_boosted = n_boosted
    )
    
  }
}


scenario_summary_df <- do.call(rbind, scenario_summary)
print(scenario_summary_df)

# --------------------------------------------------------------
# Post-Processing
# --------------------------------------------------------------
rateAF_data <- rateAF_datapreparation_step1(data_path = "H:/partTimePhD/Code/synthetic-RATEAF/data/RATE-AF.dta")

for (k in names(boosting_results)) {
  boosting_results[[k]] <- lapply(boosting_results[[k]], function(ds) {
    post_processing(rateAF_data, ds)
  })
}

# --------------------------------------------------------------
#    Fidelity Heatmap (variable-level significance across scenarios)
# --------------------------------------------------------------

# Build dataset list: baseline + boosted
dataset_list <- list(`RATE-AF` = rateAF_data)
for (k in names(boosting_results)) {
  dataset_list[[k]] <- boosting_results[[k]]  # keep as list of 100 boosted datasets
}

# Extract p-values
fidelity_pvals <- get_fidelity_pvals(rateAF_data, dataset_list,
                                     vars_categorical, vars_numeric,
                                     categorical_labels, numeric_labels)

# Add baseline reference
baseline_ref <- get_fidelity_baseline(rateAF_data, vars_categorical, vars_numeric,
                                      categorical_labels, numeric_labels)

# Heatmap
{
  heatmap_data <- fidelity_pvals %>%
    left_join(baseline_ref, by = "Variable") %>%
    mutate(
      Preservation = case_when(
        sig_flag == "Significant" ~ "Lost",
        sig_flag == "NS"          ~ "Preserved",
        TRUE ~ "NA"
      )
    )
  
  # 4. Plot Heatmap
  g_fidelity_heatmap <- ggplot(heatmap_data, aes(x = Scenario, y = Variable, fill = Preservation)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", pval_mean)), size = 4) +#, fontface = "bold") +
    scale_fill_manual(
      values = c("Preserved" = "#33a02c", 
                 "Lost" = "#e31a1c", 
                 "Spurious" = "#ff7f00",
                 "NA" = "grey80")
    ) +
    labs(
      title = sprintf(
        "Fidelity: Preservation of Kruskal test results (Averaged over %d iterations)",
        n_iterations),
      x = "Scenario",
      y = "Variable"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"))#,
  #legend.position = "top",
  #legend.title = element_blank()
  
  plot(g_fidelity_heatmap)
  
  
}


# ----------------------------------------------------------------------------
#   Supplementary Material - Averaged Fidelity Results for each ethnic group
# ----------------------------------------------------------------------------

# Run for each ethnic group
white_tables <- build_avg_fidelity_tables_boosting(rateAF_data, boosting_results, "White",
                                            vars_categorical, vars_numeric,
                                            categorical_labels, numeric_labels)

asian_tables <- build_avg_fidelity_tables_boosting(rateAF_data, boosting_results, "Asian",
                                            vars_categorical, vars_numeric,
                                            categorical_labels, numeric_labels)
  
black_tables <- build_avg_fidelity_tables_boosting(rateAF_data, boosting_results, "Black",
                                            vars_categorical, vars_numeric,
                                            categorical_labels, numeric_labels)
  

View(white_tables$categorical)
View(white_tables$numeric)
View(asian_tables$categorical)
View(asian_tables$numeric)
View(black_tables$categorical)
View(black_tables$numeric)


# --------------------------------------------------------------
#   Utility
# --------------------------------------------------------------

# Prepare Dataset List (Baseline + Boosting Replicates)
rateAF_data <- convert_nyha_numeric(rateAF_data)

# Convert nyha in ALL boosted replicates
for (k in names(boosting_results)) {
  boosting_results[[k]] <- lapply(boosting_results[[k]], convert_nyha_numeric)
}
  
# Run Utility Models (Averaged Results)
utility_long <- purrr::map_dfr(names(outcome_models), function(outcome) {
  message(sprintf("\n▶ Running model for %s", outcome))
  results_list <- list()
  
  # --- Baseline single dataset ---
  baseline <- tryCatch(
    utility_evaluation(outcome, outcome_models[[outcome]], list(`RATE-AF` = rateAF_data)),
    error = function(e) NULL
  )
  if (!is.null(baseline) && nrow(baseline) > 0)
    results_list[["Baseline"]] <- baseline
  
  # --- Boosted scenarios averaged over replicates ---
  for (scen in names(boosting_results)) {
    reps <- boosting_results[[scen]]
    attr(reps, "scenario") <- scen
    res <- tryCatch(
      utility_eval_replicates(outcome, outcome_models[[outcome]], reps),
      error = function(e) NULL
    )
    if (!is.null(res) && nrow(res) > 0)
      results_list[[scen]] <- res
  }
  
  # --- Combine only valid outputs ---
  if (length(results_list) == 0) {
    message(sprintf("⚠️ No valid results for %s — skipping.", outcome))
    return(NULL)
  }
  
  # Harmonise sig_flag type before binding
  results_list <- lapply(results_list, function(x) {
    x$sig_flag <- as.character(x$sig_flag)
    return(x)
  })
  
  dplyr::bind_rows(results_list)
})
  
  
# --------------------------------------------------------------
#   Utility Preservation Heatmap
# --------------------------------------------------------------

#   Format for Heatmap
utility_summary <- utility_long %>%
  dplyr::mutate(
    Group = case_when(
      grepl("White", dataset, ignore.case = TRUE) ~ "White",
      grepl("Asian", dataset, ignore.case = TRUE) ~ "Asian",
      grepl("Black", dataset, ignore.case = TRUE) ~ "Black",
      dataset == "RATE-AF"                        ~ "Baseline",
      TRUE                                        ~ dataset
    ),
    DropLevel = suppressWarnings(as.numeric(str_extract(dataset, "(\\d+)(?=%|$)"))),
    DropLabel = ifelse(is.na(DropLevel), "", paste0(DropLevel, "%")),
    Scenario = ifelse(Group == "Baseline", "Baseline", paste0(Group, "_", DropLabel))
  )

# Baseline significance reference
baseline_sig <- utility_summary %>%
  dplyr::filter(Group == "Baseline") %>%
  dplyr::mutate(Baseline_Significance = ifelse(p.value < 0.05, "Significant", "NS")) %>%
  dplyr::select(outcome, stat, Baseline_Significance)

# Join baseline info and define preservation
heatmap_data <- utility_summary %>%
  dplyr::left_join(baseline_sig, by = c("outcome", "stat")) %>%
  dplyr::mutate(
    p_rounded = round(p.value, 3),
    Preservation = case_when(
      Group == "Baseline" ~ "Baseline",
      is.na(p_rounded)    ~ "Not enough data to compute",
      p_rounded < 0.05 & Baseline_Significance == "Significant" ~ "Preserved",
      p_rounded >= 0.05 & Baseline_Significance == "NS"         ~ "Preserved",
      p_rounded >= 0.05 & Baseline_Significance == "Significant"~ "Lost",
      TRUE                                                      ~ "Not enough data to compute"
    ),
    label_text = ifelse(is.na(p_rounded), "X", sprintf("%.2f", p_rounded))
  )

# Define ordered axis (ensures correct scenario order + % visible)
group_levels <- c("Asian", "Black", "White")
drop_levels  <- c("30%", "60%", "80%")
ordered_levels <- c("Baseline", as.vector(outer(group_levels, drop_levels, paste, sep = "_")))
heatmap_data$Scenario <- factor(heatmap_data$Scenario, levels = ordered_levels)


g_utility_heatmap <-ggplot(heatmap_data, aes(x = Scenario, y = outcome, fill = Preservation)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(
    aes(label = label_text),
    size = 4, lineheight = 0.9
  ) +
  scale_fill_manual(
    values = c(
      "Baseline"                  = "#1f78b4",
      "Preserved"                 = "#33a02c",
      "Lost"                      = "#e31a1c",
      "Not enough data to compute"= "grey60"
    ),
    name = "Utility"
  ) +
  labs(
    title = sprintf(
      "Utility: Preservation of Model Significance (Averaged over %d iterations)",
      n_iterations
    ),
    x = "Scenario",
    y = "Outcome"
  ) +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

# Plot
plot(g_utility_heatmap)



# --------------------------------------------------------------
#   Utility Forest Plots
# --------------------------------------------------------------

# Define mapping between variable names and titles
outcome_titles <- c(
  nyhaimprove      = "NYHA Improvement",
  #nyha12m          = "NYHA at 12 Months",
  mehra2class12m   = "mEHRA Improvement",
  ecghr12m         = "ECG Heart Rate at 12 Months",
  ntprobnp12m_ln   = "NT-proBNP (Log-Transformed)",
  ae_yn            = "Adverse Events",
  cvadj            = "Cardiovascular Events"
)

#  Prepare data for plotting
plot_data <- utility_long %>%
  dplyr::mutate(
    estimate  = as.numeric(str_extract(label, "^[0-9.\\-]+")),
    conf.low  = as.numeric(str_extract(label, "(?<=95% CI )[-0-9.]+")),
    conf.high = as.numeric(str_extract(label, "(?<=–)[-0-9.]+(?=;)")),
    pval      = as.numeric(str_extract(label, "(?<=p=)[0-9\\.]+")),
    sig_flag  = ifelse(str_detect(label, "\\*"), "p < 0.05", "NS")
  ) %>%
  dplyr::filter(!is.na(estimate)) %>%
  dplyr::mutate(
    dataset = ifelse(dataset == "RATE-AF", "RATE-AF", dataset),
    dataset = str_replace_all(dataset, "%", "%") # for consistent factor levels
  )

# Define scenario order (Baseline first, then boosting scenarios)
plot_data$dataset <- factor(
  plot_data$dataset,
  levels = c("RATE-AF",
             "White_30%","White_60%","White_80%",
             "Asian_30%","Asian_60%","Asian_80%",
             "Black_30%","Black_60%","Black_80%")
)
  
  
#  Generate and display plots for each outcome
unique(plot_data$outcome) %>%
  purrr::walk(function(outcome) {
    g <- forest_plot(plot_data, outcome)
    plot(g)
  })
  
# Optional: combine individual outcome plots into one panel (requires patchwork)
outcomes <- unique(plot_data$outcome)
  groups <- split(outcomes, ceiling(seq_along(outcomes) / 6))
  
  for (i in seq_along(groups)) {
    subset <- groups[[i]]
    forest_list <- lapply(subset, function(o) forest_plot(plot_data, o))
    
    
    g_forest_panel <- wrap_plots(forest_list, ncol = 2, guides = "collect") +
      plot_layout(guides = "collect") & 
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 10)) &
      plot_annotation(
        title = "Forest plots of treatment effects (Averaged over 500 iterations)",
        theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
      )
    
    
    plot(g_forest_panel)
  }
  
  



# --------------------------------------------------------------
#   Utility Results for supplementary materials
# --------------------------------------------------------------

# Filter for White + RATE-AF, and keep only desired columns
white_result <- utility_long %>%
  filter(grepl("White", dataset) | dataset == "RATE-AF") %>%
  select(outcome, stat, dataset, result = label)

# Similarly for Asian + RATE-AF
asian_result <- utility_long %>%
  filter(grepl("Asian", dataset) | dataset == "RATE-AF") %>%
  select(outcome, stat, dataset, result = label)

# And for Black + RATE-AF
black_result <- utility_long %>%
  filter(grepl("Black", dataset) | dataset == "RATE-AF") %>%
  select(outcome, stat, dataset, result = label)

# View in RStudio
View(white_result)
View(asian_result)
View(black_result)



# --------------------------------------------------------------
#   Utility Results Stratified by Ethnic Subgroup
# --------------------------------------------------------------

# Define subgroup variable
group_var <- "ethncty"

# Get all unique subgroup values and their labels
ethnic_groups <- levels(rateAF_data[[group_var]])
ethnic_labels <- c("1" = "White", "2" = "Asian", "3" = "Black")


#   Run subgroup utility evaluation
utility_subgroup_long <- list()
for (outcome in names(outcome_models)) {
  message(sprintf("\n▶ Running subgroup utility for outcome: %s", outcome))
  
  for (g in levels(rateAF_data$ethncty)) {
    subgroup_label <- c("1"="White","2"="Asian","3"="Black")[g]
    
    # Baseline RATE-AF (original data)
    base <- tryCatch(
      utility_evaluation(outcome, outcome_models[[outcome]],
                         list(`RATE-AF` = rateAF_data %>% dplyr::filter(ethncty == g))),
      error = function(e) NULL
    )
    if (!is.null(base) && nrow(base) > 0) base$dataset <- "RATE-AF"
    
    # Boosted scenarios
    boost_list <- list()
    for (scen in names(boosting_results)) {
      reps <- boosting_results[[scen]]
      reps_filtered <- lapply(reps, function(ds) ds %>% dplyr::filter(ethncty == g))
      attr(reps_filtered, "scenario") <- scen
      
      res <- tryCatch(
        utility_eval_replicates(outcome, outcome_models[[outcome]], reps_filtered),
        error = function(e) NULL
      )
      if (!is.null(res) && nrow(res) > 0) res$dataset <- scen
      boost_list[[scen]] <- res
    }
    
    # Combine baseline + boosted (only non-null)
    valid <- c(list(base), boost_list)
    valid <- valid[!sapply(valid, is.null)]
    
    if (length(valid) > 0) {
      # Ensure consistent column types before binding
      valid <- lapply(valid, function(x) {
        if ("sig_flag" %in% names(x)) x$sig_flag <- as.character(x$sig_flag)
        x
      })
      
      df <- dplyr::bind_rows(valid)
      df$Subgroup <- subgroup_label
      utility_subgroup_long[[paste0(outcome, "_", subgroup_label)]] <- df
    }
  }
}

#   Combine and clean subgroup results
utility_subgroup_long <- dplyr::bind_rows(utility_subgroup_long)
utility_subgroup_long_clean <- utility_subgroup_long %>%
  dplyr::mutate(
    p.value = suppressWarnings(as.numeric(p.value)),
    p.value = ifelse(is.na(p.value), 1.0, p.value),
    label = ifelse(is.na(label), "Not computed", label),
    Scenario = dataset
  )


#   Prepare summary status table
utility_subgroup_long_summary <- utility_subgroup_long_clean %>%
  dplyr::mutate(
    status = dplyr::case_when(
      label == "Not computed" ~ "Insufficient data",
      p.value < 0.05 ~ "Significant",
      TRUE ~ "NS"
    )
  ) %>%
  dplyr::select(Subgroup, Scenario, outcome, p.value, status)

#   Fill all subgroup × scenario × outcome combinations
all_combos <- expand.grid(
  Subgroup = unique(utility_subgroup_long_summary$Subgroup),
  Scenario = unique(utility_subgroup_long_summary$Scenario),
  outcome  = unique(utility_subgroup_long_summary$outcome)
)

utility_subgroup_complete <- all_combos %>%
  dplyr::left_join(utility_subgroup_long_summary,
                   by = c("Subgroup", "Scenario", "outcome")) %>%
  dplyr::mutate(
    status = dplyr::case_when(
      is.na(status) ~ "Insufficient data",
      TRUE ~ status
    ),
    p.value = ifelse(is.na(p.value), 1.0, p.value)
  )



#   Plot subgroup-specific utility heatmap (complete grid)
g_utility_heatmap_subgroup <- ggplot(
  utility_subgroup_complete,
  aes(x = Scenario, y = outcome, fill = status)
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = ifelse(status == "Insufficient data", "", sprintf("%.2f", p.value))),
    size = 3
  ) +
  facet_wrap(~Subgroup, ncol = 1, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "Significant" = "#33a02c",
      "NS" = "#e31a1c",
      "Insufficient data" = "grey80"
    )
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Subgroup-Specific Utility Validation (Averaged over 500 iterations)",
    subtitle = "Grey cells indicate insufficient or uncomputed data for model estimation.",
    x = "Scenario",
    y = "Outcome",
    fill = "Result"
  )

plot(g_utility_heatmap_subgroup)


# --------------------------------------------------------------
#   Subgroup Evaluation
# --------------------------------------------------------------

#vars_numeric     <- c("age", "creat", "ntprobnp", "ecghr", "bpsysav")
#vars_categorical <- c("sex", "mdfdehrascr", "nyha", "hfdiag", "typ2dm", "estmtdejctnfrctn")

avg_subgroup_tables <- list()

for (scenario in names(boosting_results)) {
  cat("Processing:", scenario, "\n")
  
  # 1) pool the 100 boosted datasets for this scenario
  boosted_list <- boosting_results[[scenario]]
  
  # 2) build subgroup summaries exactly like the single-shot code
  real_num  <- summarise_numeric_by_subgroup(rateAF_data, vars_numeric)
  real_cat  <- summarise_categorical_by_subgroup(rateAF_data, vars_categorical)
  real_tbl  <- dplyr::bind_rows(real_num, real_cat)
  
  syn_num <- summarise_numeric_by_subgroup_averaged(boosted_list, vars_numeric)
  syn_cat <- summarise_categorical_by_subgroup_averaged(boosted_list, vars_categorical)
  syn_tbl <- dplyr::bind_rows(syn_num, syn_cat)
  
  # 3) combine and blank categorical headers
  out_tbl <- combine_subgroup_summaries(real_tbl, syn_tbl) %>%
    dplyr::mutate(Scenario = scenario, .before = 1) %>%
    blank_categorical_headers(vars_categorical)
  
  avg_subgroup_tables[[scenario]] <- out_tbl
}

# Inspect
View(avg_subgroup_tables[["White_30%"]])
View(avg_subgroup_tables[["White_60%"]])
# ...