set.seed(12345)
seed <- 12345
n_iterations <- 500

# -------------------------------------------------------
# Scripts, libraries and functions loading
# -------------------------------------------------------


setwd("synthetic-RATEAF/utils")
source('functions.R')
source('bias_analysis.R')
source('synthetic_data_generation.R')
source('fidelity_evaluation.R')
source("utility_evaluation.R")

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

# --------------------------------------------------------------
#  Bayesian Network Model Learning Setting
# --------------------------------------------------------------


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

#  Run 100 synthetic replications
synthetic_results <- iterative_syngen(
  rateAF_data_processed,
  n_reps   = n_iterations,
  size     = nrow(rateAF_data_processed),
  whitelist = whitelist,
  blacklist = blacklist,
  resample = FALSE)
  

# --------------------------------------------------------------
#   Post-Processing
# --------------------------------------------------------------
rateAF_data <- rateAF_datapreparation_step1(data_path = "H:/partTimePhD/Projects 2025/Clinical Trials/Data/RATE-AF_RPFBrunel_09jan24.dta") 

synthetic_results <- lapply(synthetic_results, function(ds) {
  post_processing(rateAF_data, ds)
})

# Dataset mappings
dataset_list <- list(
  `RATE-AF` = rateAF_data,
  `Replication` = synthetic_results
)

# --------------------------------------------------------------
#   Fidelity 
# --------------------------------------------------------------


# Averaged Fidelity Table
replication_tables <- build_avg_fidelity_tables_replication(
    rateAF_data,
    synthetic_results,
    vars_categorical,
    vars_numeric,
    categorical_labels,
    numeric_labels)
  
View(replication_tables$categorical)
View(replication_tables$numeric)



# Plot Fidelity Results
{
  # Common colour palette and theme
  palette_fidelity <- c("Real" = "#1E88E5", "Synthetic" = "#F9A825")
  
  theme_fidelity <- theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 11)
    )
  
# NUMERIC VARIABLES
{
    rep_numeric <- replication_tables$numeric
    var_rows <- which(!grepl("Mean|Median|p=", rep_numeric$Variable))
    
    plot_data <- purrr::map_dfr(var_rows, function(i) {
      if (i + 1 <= nrow(rep_numeric) && grepl("Mean", rep_numeric$Variable[i + 1])) {
        data.frame(
          variable  = rep_numeric$Variable[i],
          mean_real = as.numeric(str_extract(rep_numeric$`RATE-AF`[i + 1], "^[0-9\\.]+")),
          sd_real   = as.numeric(str_extract(rep_numeric$`RATE-AF`[i + 1], "(?<=\\[)[0-9\\.]+")),
          mean_rep  = as.numeric(str_extract(rep_numeric$Replication[i + 1], "^[0-9\\.]+")),
          sd_rep    = as.numeric(str_extract(rep_numeric$Replication[i + 1], "(?<=\\[)[0-9\\.]+")),
          pval      = as.numeric(str_extract(rep_numeric$Replication[i], "(?<=p=)[0-9\\.]+"))
        )
      }
    }) %>%
      drop_na(mean_real, sd_real, mean_rep, sd_rep)
    
    skewed_vars <- c("NT-proBNP (pg/mL)", "Creatinine (µmol/L)")
    plot_data <- plot_data %>%
      mutate(
        mean_real = ifelse(variable %in% skewed_vars, log1p(mean_real), mean_real),
        mean_rep  = ifelse(variable %in% skewed_vars, log1p(mean_rep), mean_rep),
        variable  = ifelse(variable %in% skewed_vars,
                           paste0(variable, " [log-scale]"),
                           variable)
      )
    
    simulate_dist <- function(mean, sd, n = 10000) {
      if (is.na(mean) || is.na(sd) || sd == 0) return(rep(mean, n))
      rnorm(n, mean, sd)
    }
    
    sim_data <- plot_data %>%
      rowwise() %>%
      mutate(
        real = list(simulate_dist(mean_real, sd_real)),
        rep  = list(simulate_dist(mean_rep, sd_rep))
      ) %>%
      tidyr::pivot_longer(cols = c(real, rep), names_to = "dataset", values_to = "value") %>%
      unnest(value) %>%
      mutate(dataset = recode(dataset, real = "Real", rep = "Synthetic"))
    
    plots <- plot_data %>%
      split(.$variable) %>%
      purrr::map(~{
        v <- .x$variable
        pval <- .x$pval
        ggplot(filter(sim_data, variable == v),
               aes(x = value, fill = dataset, color = dataset)) +
          geom_density(alpha = 0.25, size = 0.7) +
          scale_fill_manual(values = palette_fidelity) +
          scale_color_manual(values = palette_fidelity) +
          labs(
            title = v,
            subtitle = paste0("Kruskal–Wallis p = ",
                              ifelse(is.na(pval), "NA",
                                     formatC(pval, digits = 3, format = "f"))),
            x = NULL, y = "Density"
          ) +
          theme_fidelity +
          theme(legend.position = "none")
      })
    
    panelA <- wrap_plots(plots, ncol = 2, guides = "collect") +
      plot_annotation(
        #title = "Fidelity: Real vs Synthetic RATE-AF Data",
        theme = theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0)  # NPJ-style
        )
      ) &
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 11)
      )
    
    # Optional: caption below legend
    panelA <- panelA + plot_annotation(
      #caption = "Each density estimated from mean ± SD values aggregated across 500 iterations.",
      theme = theme(plot.caption = element_text(size = 9, hjust = 0.5))
    )
    
    # Display or save
    plot(panelA)
    #ggsave("H:/partTimePhD/Code/Synthetic-Clinical-Trials/NPJ plots/Replication_Fidelity.png", panelA, width = 8, height = 8, dpi = 300)
    
    }
  
  
# CATEGORICAL VARIABLES
{
    rep_categorical <- replication_tables$categorical
    
    # 1️⃣ Detect variable blocks
    var_starts <- which(grepl("p=", rep_categorical$Replication))
    plot_data_cat <- list()
    
    #  Parse each variable and extract levels + percentages
    for (i in seq_along(var_starts)) {
      start_row <- var_starts[i]
      end_row <- ifelse(i < length(var_starts), var_starts[i + 1] - 1, nrow(rep_categorical))
      block <- rep_categorical[start_row:end_row, ]
      
      var_name <- block$Variable[1]
      pval <- suppressWarnings(as.numeric(str_extract(block$Replication[1], "(?<=p=)[0-9\\.]+")))
      
      levels_block <- block[-1, ] %>%
        mutate(
          level     = trimws(Variable),
          real_perc = as.numeric(str_extract(`RATE-AF`, "(?<=\\().+?(?=%)")),
          rep_perc  = as.numeric(str_extract(Replication, "(?<=\\().+?(?=%)"))
        ) %>%
        filter(!is.na(real_perc), !is.na(rep_perc))
      
      if (nrow(levels_block) == 0) next
      
      plot_data_cat[[length(plot_data_cat) + 1]] <- levels_block %>%
        select(level, real_perc, rep_perc) %>%
        pivot_longer(cols = c(real_perc, rep_perc),
                     names_to = "dataset", values_to = "percentage") %>%
        mutate(
          dataset = ifelse(dataset == "real_perc", "Real", "Synthetic"),
          variable = var_name,
          pval = pval
        )
    }
    
    plot_data_cat <- bind_rows(plot_data_cat)
    if (nrow(plot_data_cat) == 0) stop("No valid categorical data found.")
    
    #  Define consistent palette (same as numeric figure)
    palette_fidelity <- c("Real" = "#1E88E5", "Synthetic" = "#F9A825")
    
    #  Build per-variable plots
    plots_cat <- plot_data_cat %>%
      split(.$variable) %>%
      map(~{
        v <- unique(.x$variable)
        pval <- unique(.x$pval)
        
        ggplot(.x, aes(x = level, y = percentage, fill = dataset)) +
          geom_col(position = position_dodge(width = 0.7),
                   width = 0.6, color = "grey40") +
          scale_fill_manual(values = palette_fidelity) +
          labs(
            title = v,
            subtitle = paste0("Kruskal–Wallis p = ",
                              ifelse(is.na(pval), "NA",
                                     formatC(pval, digits = 3, format = "f"))),
            x = NULL, y = "Percentage (%)"
          ) +
          theme_minimal(base_size = 12) +
          theme(
            legend.position = "none",
            plot.title = element_text(size = 11, face = "bold"),
            plot.subtitle = element_text(size = 9),
            axis.text.x = element_text(angle = 45, hjust = 1)
          )
      })
    
    #  Assemble figure with shared legend and title
    panelA_cat <- wrap_plots(plots_cat, ncol = 2, guides = "collect") +
      plot_annotation(
        #title = "Distributional Fidelity Between Real and Synthetic RATE-AF Data (Categorical Variables)",
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))
      ) &
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 11)
      )
    
    #  Optional caption (same style as numeric)
    panelA_cat <- panelA_cat +
      plot_annotation(
        #caption = "Proportions derived from average percentages across 500 iterations.",
        theme = theme(plot.caption = element_text(size = 9, hjust = 0.5))
      )
    
    #  Display or save
    plot(panelA_cat)
    #ggsave("H:/partTimePhD/Code/Synthetic-Clinical-Trials/NPJ plots/Replication_Fidelity_Cat.png", panelA_cat, width = 8, height = 8, dpi = 300)
  }
  
  
# Combined
{
  
  # Combine numeric (left) and categorical (right)
  fig_combined <- (panelA | panelA_cat) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 11)
    )
  
  # Add a single overall caption
  fig_combined <- fig_combined +
    plot_annotation(
      title = "Fidelity of Real vs Synthetic RATE-AF Data Across Baseline Characteristics",
      subtitle = "Left: Numeric variables — Right: Categorical variables",
      caption = "Each density or proportion is aggregated across 500 iterations.",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # ← centered
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        plot.caption = element_text(size = 9, hjust = 0.5)
      )
    )
  
  
  # Display or save
  plot(fig_combined)
  ggsave(
    "H:/partTimePhD/Code/Synthetic-Clinical-Trials/NPJ plots/Replication_Fidelity.png",
     fig_combined, width = 15, height = 8, dpi = 300) 
}
}

  


# --------------------------------------------------------------
#   Utility
# --------------------------------------------------------------



# Apply conversion to real and synthetic datasets
rateAF_data <- convert_nyha_numeric(rateAF_data)
for (i in 1:length(dataset_list[[2]])) {
  dataset_list[[2]][[i]] <- convert_nyha_numeric(dataset_list[[2]][[i]])
}


# Run all outcome models
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
  
  # --- Replication (averaged over 100 iterations) ---
  #reps <- dataset_list[["Replication"]]
  reps <- dataset_list[[2]]
  attr(reps, "scenario") <- "Replication"
  res <- tryCatch(
    utility_eval_replicates(outcome, outcome_models[[outcome]], reps),
    error = function(e) NULL
  )
  if (!is.null(res) && nrow(res) > 0)
    results_list[["Replication"]] <- res
  
  # --- Combine only valid outputs ---
  if (length(results_list) == 0) {
    message(sprintf("⚠️ No valid results for %s — skipping.", outcome))
    return(NULL)
  }
  
  # Harmonise column types before binding
  results_list <- lapply(results_list, function(x) {
    x$sig_flag <- as.character(x$sig_flag)
    return(x)
  })
  
  dplyr::bind_rows(results_list)
})


# Plot all outcomes
outcome_titles <- c(
  nyhaimprove      = "NYHA Improvement",
  mehra2class12m   = "mEHRA Improvement",
  ecghr12m         = "ECG Heart Rate at 12 Months",
  ntprobnp12m_ln   = "NT-proBNP (Log-Transformed)",
  ae_yn            = "Adverse Events",
  cvadj            = "Cardiovascular Events"
)
    
# Prepare data for plotting
plot_data <- utility_long %>%
  mutate(
    # Extract numeric estimates if not already present
    estimate  = if (!"estimate" %in% names(.)) as.numeric(str_extract(label, "^[0-9.\\-]+")) else estimate,
    conf.low  = if (!"conf.low" %in% names(.)) as.numeric(str_extract(label, "(?<=95% CI )[-0-9.]+")) else conf.low,
    conf.high = if (!"conf.high" %in% names(.)) as.numeric(str_extract(label, "(?<=–)[-0-9.]+")) else conf.high,
    sig_flag  = ifelse(str_detect(label, "\\*"), "p < 0.05", "NS")
  ) %>%
  filter(dataset %in% c("RATE-AF", "Replication")) %>%
  mutate(
    dataset = recode(dataset,
                     "RATE-AF" = "Real",
                     "Replication" = "Synthetic")
  )

# Define dataset order
plot_data$dataset <- factor(plot_data$dataset, levels = c("Real", "Synthetic"))

# Generate one plot per outcome
forest_list <- unique(plot_data$outcome) %>%
  purrr::map(function(outcome) {
    title_label <- outcome_titles[[outcome]] %||% outcome
    g <- forest_plot(plot_data, outcome_label = outcome) +
      labs(title = title_label)
    return(g)
  })

# Combine all plots into a single 2×3 panel
g_forest_panel <- wrap_plots(forest_list, ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Forest plots of treatment effects (Averaged over 500 iterations)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  ) &
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

# Display
plot(g_forest_panel)


