set.seed(12345)

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
  'fairness','caret', 'patchwork', 'visNetwork', 'igraph'
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
#  Bayesian Network Learning
# --------------------------------------------------------------

# --- Whitelist: Force inclusion of clinically justified arcs ---
whitelist <- data.frame(from = c(
  "typ2dm", "ethncty", "ethncty", "ethncty", "ethncty", "ethncty", "ethncty", "ethncty",
  "ntprobnp_discretised", "signshf", "mihist", "mihist", "mihist", "mrthn8alchl",
  "bpsysav_discretised", "antiplat", "anticoag", "cabg",
  "creat_discretised", "egfr_discretised", "mehra2class12m",
  "mehra2class12m", "sex", "sex", "sex", "mihist", "cvadj", "yearsaf_discretised", "antihypertensive"
), to = c(
  "insln",
  "hfdiag", "mihist", "creat_discretised", "bpsysav12m_discretised",
  "hvecopdemphysma", "typ2dm", "strokehist",
  "hfdiag", "nyha",
  "antiplat", "pci", "cabg",
  "hvlvrdisease",
  "hfdiag", "hshdmjrblds", "hshdmjrblds", "hfdiag",
  "egfr_discretised", "urea_discretised", "nyha12m",
  "gp_cv_yn", "ht_discretised", "wt_discretised", "waist_discretised",
  "typ2dm", "alloc", "anticoag","acearb"
))

# --- Blacklist: Prevent inclusion of implausible or reversed arcs ---
blacklist <- data.frame(from = c(
  "alloc", "typ2dm",
  "mrthn8alchl", "hfdiag", "inhlrsarwydis", "mrthn8alchl",
  "hvasthma", "inhlrsarwydis", "antiplat", "antihypertensive" 
), to = c(
  "cvadj", "mihist",
  "hsrcvdtrtmntthyrd", "inhlrsarwydis", "hfdiag", "hsrcvdtrtmntthyrd",
  "hshdmjrblds", "hshdmjrblds", "hvasthma", "ppm"
))


# Fit the model using structural EM with constraints
bn_hc_mle <- structural.em(
  x = rateAF_data_processed,
  maximize = "hc",
  maximize.args = list(
    whitelist = whitelist,
    blacklist = blacklist,
    score = "bic"
  ),
  fit = "mle",
  fit.args = list(),
  impute = "parents",
  impute.args = list(),
  return.all = TRUE,
  start = NULL,
  max.iter = 1000,
  debug = FALSE
)



# Visualise the updated network (if you're using your helper functions)
bn_hc_mle_igraph <- convert_bn_to_igraph(bn_hc_mle)
visualise_network <- visualise_bn_with_visnetwork(bn_hc_mle_igraph, TRUE)
print(visualise_network)




# --------------------------------------------------------------
#   Exploratory Boosting
# --------------------------------------------------------------


sensitive_attributes <- 'ethncty'
bias_mode <- 'coverage'

  
# Augmentation Setting - Fully Synthetic  

# Option 2: Fully Synthetic
augmented_size <- 1600
augmentation_type <- 'Synthetic'
calculate_additional_data_method <- 'proportional_adjustment' # for fully Synthetic 

# Calculate data needed 
statistics <- calculate_coverage_metrics_intersection(rateAF_data_processed, sensitive_attributes)
additional_data_needed <- calculate_additional_data_needed(statistics, desired_coverage_ratio = 1, expected_proportions = expected_proportions, desired_data_size = augmented_size, bias_mode = bias_mode, method = calculate_additional_data_method)

# Create Evidence
evidence <- create_evidence_df(rateAF_data_processed, sensitive_attributes, desired_data_size = augmented_size)

# Format Evidence
evidence_formatted <- format_evidence(additional_data_needed, evidence, sensitive_attributes)

# Conditional Sampling
bn_fitted <- bn_hc_mle$fitted
bn <- bn_hc_mle$dag
aug_syn <- conditional_sampling(bayesian_network_name = 'bn_fitted', 
                                bayesian_network = bn, 
                                fitted_bayesian_network = bn_fitted,
                                evidence_formatted = evidence_formatted)


  





# --------------------------------------------------------------
#   Post-Processing
# --------------------------------------------------------------
# assuming rateAF_data is the original dataset, post-processing function will re-align synthetic and real data
aug_syn <- post_processing(rateAF_data, aug_syn)

# --------------------------------------------------------------
#   Fidelity 
# --------------------------------------------------------------


# Dataset mappings
dataset_list <- list(`Boost Syn` = aug_syn)
  
categorical_table <- list()
for (var in vars_categorical) {
  block <- fidelity_evaluation_categorical(var, rateAF_data, dataset_list,categorical_labels)
  categorical_table[[length(categorical_table) + 1]] <- block
}
fidelity_categorical_summary <- bind_rows(categorical_table)


numeric_table <- list()
for (var in vars_numeric) {
  block <- fidelity_evaluation_numeric(var, rateAF_data, dataset_list,numeric_labels)
  numeric_table[[length(numeric_table) + 1]] <- block
}
fidelity_numeric_summary <- bind_rows(numeric_table)


# Preview 
View(fidelity_categorical_summary)
View(fidelity_numeric_summary)


# --------------------------------------------------------------
#   Utility
# --------------------------------------------------------------


# Prepare dataset list
dataset_list <- list(`RATE-AF` = rateAF_data,`Boost Syn` = aug_syn)

# Ensure NYHA variable is numeric where needed
dataset_list <- lapply(dataset_list, convert_nyha_numeric)
  
# Run all outcome models
utility_long <- purrr::map_dfr(names(outcome_models), function(outcome) {
  message(sprintf("▶ Running model for %s", outcome))
  utility_evaluation(outcome, outcome_models[[outcome]], dataset_list)
})
  
# Format results into table view
utility_results <- utility_long %>%
  dplyr::select(outcome, stat, dataset, label) %>%
  tidyr::pivot_wider(names_from = dataset, values_from = label) %>%
  dplyr::arrange(outcome)

View(utility_results)
  
  
# Convert utility_long to numeric components
plot_data <- utility_results %>%
  pivot_longer(cols = -c(outcome, stat), names_to = "dataset", values_to = "label") %>%
  mutate(
    estimate = as.numeric(str_extract(label, "^[0-9.\\-]+")),
    conf.low = as.numeric(str_extract(label, "(?<=95% CI )[-0-9.]+")),
    conf.high = as.numeric(str_extract(label, "(?<=–)[-0-9.]+(?=;)")),
    pval = str_extract(label, "(?<=p=)[^\\s\\*\\)]+"),
    sig_flag = ifelse(str_detect(label, "\\*"), "p < 0.05", "NS")
  )

plot_data$dataset <- factor(plot_data$dataset, levels = c("RATE-AF", "Boost Syn" ))

forest_list <- unique(plot_data$outcome) %>%
  purrr::map(function(outcome) {
    g <- forest_plot(plot_data, outcome_label = outcome) +
      labs(title = outcome)
    print(g)
    return(g)
  })  

# Optionally combine plots
g_panel <- wrap_plots(forest_list, ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Forest plots of treatment effects (Exploratory Boosting)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

plot(g_panel)


# --------------------------------------------------------------
# Subgroup Fidelity Analysis Function
# --------------------------------------------------------------

# Categorical Summary Function
vars_categorical <- c("sex", "mdfdehrascr", "nyha", "hfdiag", "typ2dm", "estmtdejctnfrctn")

# Generate summary
numeric_summary <- summarise_numeric_by_subgroup(aug_syn, vars_numeric)
categorical_summary <- summarise_categorical_by_subgroup(aug_syn, vars_categorical)


numeric_summary_rateaf <- summarise_numeric_by_subgroup(rateAF_data, vars_numeric)
categorical_summary_rateaf <- summarise_categorical_by_subgroup(rateAF_data, vars_categorical)


# Combine and view
subgroup_analysis_boostsyn <- bind_rows(numeric_summary, categorical_summary)
subgroup_analysis_rateaf <- bind_rows(numeric_summary_rateaf, categorical_summary_rateaf)


desired_order <- c("White", "Asian", "Black", "Overall")
subgroup_analysis_boostsyn <- subgroup_analysis_boostsyn %>%
  select(Variable, all_of(desired_order))
subgroup_analysis_rateaf <- subgroup_analysis_rateaf %>%
  select(Variable, all_of(desired_order))

subgroup_analysis_boostsyn <- subgroup_analysis_boostsyn %>%
  filter(!(is.na(Variable) | trimws(Variable) == "NA"))

subgroup_analysis_rateaf <- subgroup_analysis_rateaf %>%
  filter(!(is.na(Variable) | trimws(Variable) == "NA"))

#View(categorical_summary)
View(subgroup_analysis_boostsyn)
View(subgroup_analysis_rateaf)


subgroup_rateaf_vs_syn <- combine_subgroup_summaries(
  rateaf_df = subgroup_analysis_rateaf,
  syn_df = subgroup_analysis_boostsyn
)

View(subgroup_rateaf_vs_syn)




