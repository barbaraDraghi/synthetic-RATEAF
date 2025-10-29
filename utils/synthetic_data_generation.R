# --------------------------------------------------------------------
# Define the Functions used for Synthetic Data Generation and Boosting
# --------------------------------------------------------------------

# Synthetic Data Generation Process
#
# This function generates synthetic data using either a bias-aware or bias-unaware approach.
# The bias-aware mode ensures that the generated data accounts for coverage bias based on the provided sensitive attributes.
#
# param data A dataframe containing the original dataset.
# param size Integer, the size of the synthetic dataset to generate.
# param sensitive_attributes A character vector of sensitive attributes used for bias-aware generation.
# param mode A character string specifying the mode of synthetic generation ('bias-aware' or 'bias-unaware').
# param bias_mode A character string specifying the bias mode (currently supports only 'coverage' for bias-aware mode).
# param calculate_additional_data_method A character string specifying the mode of computing the additional data needed when boosting ('proportional_adjustment' or 'match_to_maximum')
#
# return A dataframe containing the generated synthetic data.
#
# examples
# synthetic_data_generation(data = df, size = 1000, sensitive_attributes = c("gender"), mode = "bias-aware", bias_mode = "coverage")
#
# NB: - so far, only coverage mode is allowed
synthetic_data_generation <- function(data, size, sensitive_attributes, mode, bias_mode, expected_proportions = NULL, calculate_additional_data_method=NULL){
  
  # Input Validation
  tryCatch({
    if (!is.data.frame(data)) {
      stop("Error in 'synthetic_data_generation': 'data' must be a dataframe.")
    }
    if (!is.numeric(size) || size <= 0) {
      stop("Error in 'synthetic_data_generation': 'size' must be a positive integer.")
    }
    if (!is.character(sensitive_attributes) || length(sensitive_attributes) == 0) {
      stop("Error in 'synthetic_data_generation': 'sensitive_attributes' must be a non-empty character vector.")
    }
    if (!is.character(mode) || !mode %in% c('bias-aware', 'bias-unaware')) {
      stop("Error in 'synthetic_data_generation': 'mode' must be either 'bias-aware' or 'bias-unaware'.")
    }
    if (mode == 'bias-aware' && (!(bias_mode %in% c('coverage','expert_knowledge')))) {
      stop("Error in 'synthetic_data_generation': 'bias_mode' must be 'coverage' for bias-aware mode.")
    }
  }, error = function(e) {
    stop("Input validation failed in 'synthetic_data_generation': ", conditionMessage(e))
  })
  
  # Main Logic for Synthetic Data Generation
  tryCatch({
    # Check if all features in data are discretized (factors)
    all_discretized <- all(sapply(data, is.factor))
    if (!all_discretized) {
      stop("Error in 'synthetic_data_generation': All features in 'data' must be discretized (factors).")
    }
    
    if(mode == 'bias-aware'){
      # Call the bias-aware synthetic data generators
      if(bias_mode == 'expert_knowledge'){
        synthetic_data <- bias_aware_synthetic_data_generation(data, sensitive_attributes, size, bias_mode, expected_proportions, calculate_additional_data_method)
      }else{
        synthetic_data <- bias_aware_synthetic_data_generation(data, sensitive_attributes, size, bias_mode, expected_proportions = NULL, calculate_additional_data_method)
      }
      
    }else{
      # Call the bias-unaware synthetic data generators
      synthetic_data <- bias_unaware_synthetic_data_generation(data, size)
    }
    
  }, error = function(e) {
    stop("Error during processing in 'synthetic_data_generation': ", conditionMessage(e))
  })
  
  return(synthetic_data)
}

# Bias-Unaware Synthetic Data Generation
#
# This function generates synthetic data without considering any bias, using a Bayesian network learned from the original dataset.
#
# param data A dataframe containing the original dataset.
# param size Integer, the size of the synthetic dataset to generate (default is the size of the original dataset).
#
# return A dataframe containing the generated synthetic data.
#
# examples
# bias_unaware_synthetic_data_generation(df, size = 1000)
#
bias_unaware_synthetic_data_generation <- function(data, size = nrow(data)){
  
  # Input Validation
  tryCatch({
    if (!is.data.frame(data)) {
      stop("Error in 'bias_unaware_synthetic_data_generation': 'data' must be a dataframe.")
    }
    if (!is.numeric(size) || size <= 0) {
      stop("Error in 'bias_unaware_synthetic_data_generation': 'size' must be a positive integer.")
    }
  }, error = function(e) {
    stop("Input validation failed in 'bias_unaware_synthetic_data_generation': ", conditionMessage(e))
  })
  
  # Main Logic for Bias-Unaware Synthetic Data Generation
  synthetic_data <- tryCatch({
    
    bn <-structural.em(data, maximize = "hc", maximize.args = list(), fit = "mle",
                       fit.args = list(), impute = "parents", impute.args = list(), return.all = FALSE,
                       start = NULL, max.iter = 1000, debug = FALSE)
    
    bn_fitted <- bn.fit(bn, data)
    
    
    synthetic_data <- rbn(bn_fitted, n = size )
    
  }, error = function(e) {
    stop("Error during processing in 'bias_unaware_synthetic_data_generation': ", conditionMessage(e))
  })
  
  return(synthetic_data)
}

# Bias-Aware Synthetic Data Generation
#
# This function generates synthetic data while considering coverage bias for the specified sensitive attributes. 
# It ensures that the generated data maintains proper representation for intersectional demographic groups.
#
# param data A dataframe containing the original dataset.
# param sensitive_attributes A character vector of sensitive attributes.
# param size Integer, the size of the synthetic dataset to generate (default is the size of the original dataset).
# param bias_mode A character string specifying the bias mode (currently only supports 'coverage' and 'expert_knowledge').
# param calculate_additional_data_method Options are:
#   \itemize{
#     \item \code{"proportional_adjustment"}: Adjust counts proportionally. (when bias_mode == 'coverage')
#     \item \code{"preserve_excess"}: Retain over-covered categories and fill shortfalls. (when bias_mode == 'coverage')
#     \item \code{"match_to_maximum"}: Equalize all categories to the maximum count. (when bias_mode == 'coverage')
#     \item \code{"expert_knowledge"}: ----------- (when bias_mode == 'expert_knowledge')
# return A dataframe containing the generated synthetic data.
#
# examples
# bias_aware_synthetic_data_generation(df, sensitive_attributes = c("gender"), size = 1000, bias_mode = "coverage")
#
bias_aware_synthetic_data_generation <- function(data, sensitive_attributes, size = nrow(data), bias_mode, expected_proportions = NULL, calculate_additional_data_method){
  
  # Input Validation
  tryCatch({
    if (!is.data.frame(data)) {
      stop("Error in 'bias_aware_synthetic_data_generation': 'data' must be a dataframe.")
    }
    if (!is.character(sensitive_attributes) || length(sensitive_attributes) == 0) {
      stop("Error in 'bias_aware_synthetic_data_generation': 'sensitive_attributes' must be a non-empty character vector.")
    }
    if (!is.numeric(size) || size <= 0) {
      stop("Error in 'bias_aware_synthetic_data_generation': 'size' must be a positive integer.")
    }
    if (!(bias_mode %in% c('coverage','expert_knowledge'))) {
      stop("Error in 'bias_aware_synthetic_data_generation': Unsupported 'bias_mode'. Only 'coverage' and 'expert_knowledge' are supported.")
    }
  }, error = function(e) {
    stop("Input validation failed in 'bias_aware_synthetic_data_generation': ", conditionMessage(e))
  })
  
  
  # Main Logic for Bias-Aware Synthetic Data Generation
  tryCatch({
    
    # Bayesian Network learning
    bn <- structural.em(data, 
                        maximize = "hc", 
                        maximize.args = list(), 
                        fit = "mle",
                        fit.args = list(), 
                        impute = "parents", 
                        impute.args = list(), 
                        return.all = FALSE,
                        start = NULL, 
                        max.iter = 1000, 
                        debug = FALSE)
    
    bn_fitted <- bn.fit(bn, data)
    
    
    # if(bias_mode == 'coverage'){
    #   # Calculate data needed based on coverage
    #   # method = c("proportional_adjustment", "preserve_excess", "match_to_maximum")
    #   coverage_intersection <- calculate_coverage_metrics_intersection(data, sensitive_attributes)
    #   additional_data_needed <- calculate_additional_data_needed(coverage_intersection, desired_coverage_ratio = 1, expected_proportions = NULL, desired_data_size = size, bias_mode = bias_mode, method = calculate_additional_data_method)
    # }
    # 
    
    # Calculate data needed 
    statistics <- calculate_coverage_metrics_intersection(data, sensitive_attributes)
    additional_data_needed <- calculate_additional_data_needed(statistics, desired_coverage_ratio = 1, expected_proportions = expected_proportions, desired_data_size = size, bias_mode = bias_mode, method = calculate_additional_data_method)
    
    
    
    # Create Evidence
    if(bias_mode == 'coverage' && calculate_additional_data_method == 'proportional_adjustment'){
      evidence <- create_evidence_df(data, sensitive_attributes, desired_data_size = size)
    }
    if(calculate_additional_data_method == 'match_to_maximum' || calculate_additional_data_method == 'preserve_excess' || (bias_mode == 'expert_knowledge' && calculate_additional_data_method == 'proportional_adjustment')){
      evidence <- create_evidence_df_based_on_additional_data_needed(data, additional_data_needed)
    }
    print(evidence)
    
    # Format Evidence
    evidence_formatted <- format_evidence(additional_data_needed, evidence, sensitive_attributes)
    
    # Conditional Sampling
    synthetic_data <- conditional_sampling(bayesian_network_name = 'bn_fitted', 
                                           bayesian_network = bn, 
                                           fitted_bayesian_network = bn_fitted,
                                           evidence_formatted = evidence_formatted)
    
    
    # # check data size
    # if(nrow(synthetic_data) > size){
    #   synthetic_data <- dplyr::sample_n(synthetic_data,size)
    # }
    
    # check data size
    if(calculate_additional_data_method %in% c('proportional_adjustment','expert_knowledge') & nrow(synthetic_data) >= size){
      if(nrow(synthetic_data) > size){
        boosted_data <- dplyr::sample_n(synthetic_data,size)
      }else{
        boosted_data <- synthetic_data
      }
      
      
    } else {
      # We have to add synthetic data to original data
      boosted_data <- rbind(data,synthetic_data)
      
    }
    
    
    return(boosted_data)
    
    
  }, error = function(e) {
    stop("Error during processing in 'bias_aware_synthetic_data_generation': ", conditionMessage(e))
  })
  
  
}

# Create Evidence DataFrame for Synthetic Data Generation
#
# This function creates a dataframe of all possible combinations of attribute values for the specified sensitive attributes. 
# It also calculates the expected number of samples for each combination.
#
# param data A dataframe containing the original dataset.
# param attributes A character vector of sensitive attributes.
#
# return A dataframe containing the combinations of attribute values and the expected number of samples for each combination.
#
# examples
# create_evidence_df(df, attributes = c("gender", "age"))
#
create_evidence_df <- function(data, attributes, desired_data_size = nrow(data)) {
  
  cat('create_evidence_df')
  
  # Input Validation
  tryCatch({
    if (!is.data.frame(data)) {
      stop("Error in 'create_evidence_df': 'data' must be a dataframe.")
    }
    if (!is.character(attributes) || length(attributes) == 0) {
      stop("Error in 'create_evidence_df': 'attributes' must be a non-empty character vector.")
    }
    if (any(!attributes %in% colnames(data))) {
      stop("Error in 'create_evidence_df': Some attributes not found in the dataframe.")
    }
  }, error = function(e) {
    stop("Input validation failed in 'create_evidence_df': ", conditionMessage(e))
  })
  
  
  # Main Logic for Creating Evidence DataFrame
  tryCatch({
    # Generate all possible combinations of attribute values
    data_combinations <- data %>%
      select(all_of(attributes)) %>%
      distinct() %>%
      group_by(across(all_of(attributes))) %>%
      summarise(count = n(), .groups = 'drop')
    
    # Calculate the total number of combinations
    total_combinations <- nrow(data_combinations)
    
    # Calculate the expected number of samples per combination
    expected_samples <- ceiling((1 / total_combinations) * desired_data_size)
    
    # Add the expected samples and the combination name to each row
    data_combinations <- data_combinations %>%
      mutate(expected_n_samples = expected_samples) %>%
      ungroup()
    
    return(as.data.frame(data_combinations))
    
  }, error = function(e) {
    stop("Error during processing in 'create_evidence_df': ", conditionMessage(e))
  })
}



# Create Evidence DataFrame Based on Additional Data Needed
#
# This function creates a dataframe of all possible combinations of attribute values for the specified sensitive attributes
# and integrates the additional data needed to achieve the desired target proportions.
#
# param data A dataframe containing the original dataset.
# param attributes A character vector specifying the sensitive attributes for which combinations need to be evaluated.
# param additional_data_needed A dataframe containing the results of additional data needed calculations. 
#        It must include the same attributes and a column `additional_data_needed` with the data to be added.
#
# return A dataframe containing:
# \itemize{
#   \item \code{attributes}: Columns representing the attribute combinations.
#   \item \code{count}: The current count of records for each combination.
#   \item \code{additional_data_needed}: The additional data needed for each combination (defaulted to 0 if missing).
#   \item \code{total_required_samples}: The total required samples (sum of count and additional data needed).
# }
#
create_evidence_df_based_on_additional_data_needed <- function(data, additional_data_needed) {
  
  # Input Validation
  tryCatch({
    if (!is.data.frame(data)) {
      stop("Error in 'create_evidence_df_based_on_additional_data_needed': 'data' must be a dataframe.")
    }
    if (!is.list(additional_data_needed) || length(additional_data_needed) == 0) {
      stop("Error in 'create_evidence_df_based_on_additional_data_needed': 'additional_data_needed' must be a non-empty list.")
    }
  }, error = function(e) {
    stop("Input validation failed in 'create_evidence_df_based_on_additional_data_needed': ", conditionMessage(e))
  })
  
  # Initialize an empty evidence dataframe
  evidence_df <- data.frame()
  
  # Iterate over each element in the list additional_data_needed
  for (combination_name in names(additional_data_needed)) {
    
    # Extract the dataframe for the current combination
    current_df <- additional_data_needed[[combination_name]]
    
    # Split the combination_name (e.g., "name1&name2") to identify variable names
    variables <- unlist(strsplit(combination_name, "&"))
    
    # Check if all variables exist in the dataframe
    if (!all(variables %in% colnames(current_df))) {
      stop(paste("Error: Some variables", paste(variables, collapse = ", "), 
                 "are not present in the 'additional_data_needed' dataframe."))
    }
    
    # Calculate the expected samples based on 'additional_data_needed'
    evidence <- current_df[, c(variables, "count")]
    evidence$expected_n_samples <- current_df$additional_data_needed
    
    # Combine evidence for all combinations
    evidence_df <- rbind(evidence_df, evidence)
  }
  
  # Return the final evidence dataframe
  return(evidence_df)
}

# Format Evidence for Conditional Sampling
#
# This function formats the evidence (combinations of attribute values) for use in conditional sampling by the Bayesian network.
#
# param additional_data_needed A list of additional data needed for each group.
# param evidence A dataframe containing the combinations of attribute values.
# param sensitive_attributes A character vector of sensitive attributes.
#
# return A list of formatted evidence for conditional sampling, each with the evidence string and number of samples needed.
#
# examples
# format_evidence(additional_data_needed, evidence, sensitive_attributes = c("gender", "age"))
#
format_evidence <- function(additional_data_needed,evidence,sensitive_attributes){
  
  cat('format_evidence')
  
  # Input Validation
  tryCatch({
    if (!is.list(additional_data_needed)) {
      stop("Error in 'format_evidence': 'additional_data_needed' must be a list.")
    }
    if (!is.data.frame(evidence)) {
      stop("Error in 'format_evidence': 'evidence' must be a dataframe.")
    }
    if (!is.character(sensitive_attributes) || length(sensitive_attributes) == 0) {
      stop("Error in 'format_evidence': 'sensitive_attributes' must be a non-empty character vector.")
    }
    if (any(!sensitive_attributes %in% colnames(evidence))) {
      stop("Error in 'format_evidence': Some sensitive attributes are not found in the 'evidence' dataframe.")
    }
  }, error = function(e) {
    stop("Input validation failed in 'format_evidence': ", conditionMessage(e))
  })
  
  # Main Logic for Formatting Evidence
  tryCatch({
    evidence_formatted <- list()
    index <- 1
    
    for (comb_name in names(additional_data_needed)) {
      data_comb <- additional_data_needed[[comb_name]]
      
      for (row_index in 1:nrow(evidence)) {
        current_row <- evidence[row_index, ]
        
        if (current_row$expected_n_samples > 0) {
          
          # Build evidence string 
          evidence_items <- sapply(sensitive_attributes, function(col) {
            if (!is.na(current_row[[col]])) {  # Ensure column is not NA before adding to evidence
              paste(col, "==", shQuote(current_row[[col]], type="sh"), sep="")
            }
          })
          evidence_items <- evidence_items[!sapply(evidence_items, is.null)]  # Remove NULL entries
          evidence_string <- paste(evidence_items, collapse = " & ")
          
          # Store evidence and number of samples needed
          num_samples <- ceiling(current_row$expected_n_samples)
          evidence_formatted[[index]] <- list(evidence = evidence_string, num_samples = num_samples)
          index <- index + 1
        }
      }
    }  
    return(evidence_formatted)
    
  }, error = function(e) {
    stop("Error during processing in 'format_evidence': ", conditionMessage(e))
  })
  
}

# Conditional Sampling from Bayesian Network
#
# This function performs conditional sampling from a Bayesian network based on the provided evidence. It generates synthetic data
# while adhering to the specified demographic distributions.
#
# param bayesian_network_name A character string specifying the name of the fitted Bayesian network.
# param bayesian_network A Bayesian network object learned from the original dataset.
# param fitted_bayesian_network A fitted Bayesian network object used for sampling.
# param evidence_formatted A list of formatted evidence for conditional sampling.
#
# return A dataframe containing the generated synthetic data.
#
# examples
# conditional_sampling("bn_fitted", bn, bn_fitted, evidence_formatted)
#
conditional_sampling <- function(bayesian_network_name = 'bn_fitted', bayesian_network, fitted_bayesian_network, evidence_formatted){
  
  cat('conditional_sampling')
  # Input Validation
  tryCatch({
    if (!is.character(bayesian_network_name) || length(bayesian_network_name) != 1) {
      stop("Error in 'conditional_sampling': 'bayesian_network_name' must be a single character string.")
    }
    if (!inherits(bayesian_network, "bn")) {
      stop("Error in 'conditional_sampling': 'bayesian_network' must be a valid Bayesian network object.")
    }
    if (!inherits(fitted_bayesian_network, "bn.fit")) {
      stop("Error in 'conditional_sampling': 'fitted_bayesian_network' must be a valid fitted Bayesian network object.")
    }
    if (!is.list(evidence_formatted) || length(evidence_formatted) == 0) {
      stop("Error in 'conditional_sampling': 'evidence_formatted' must be a non-empty list.")
    }
  }, error = function(e) {
    stop("Input validation failed in 'conditional_sampling': ", conditionMessage(e))
  })
  
  # Prepare Bayesian network environment
  # This <<- is the trick used to then make bn_fitted accessible by the evaluation of the query below
  bn_fitted <<- fitted_bayesian_network
  bayesian_network_name <<- bayesian_network_name
  
  # Set up parallel backend
  require(parallel)
  num_cores <- detectCores() / 2  # Use half of the cores to leave enough memory
  cl <- makeCluster(num_cores)
  on.exit(stopCluster(cl))  # Ensure cluster is stopped when done or on error
  
  synthetic_data <- data.frame()
  node_names <- names(bayesian_network$nodes)
  nodes_string <<- paste("c('", paste(node_names, collapse = "', '"), "')", sep = "")
  max_retries <- 5  # Reduce the number of retries
  multiplication_factor <- 10 # Can be reduced
  
  # Export the fitted Bayesian network to the cluster
  clusterExport(cl, varlist = c(bayesian_network_name))
  
  # Main Logic for Conditional Sampling
  tryCatch({
    
    for (evidence in evidence_formatted) {
      
      print(evidence)
      
      # Reset variables for this evidence
      accumulated_samples <- list()
      total_samples_collected <- 0
      
      
      
      if (evidence$num_samples > 0) {
        
        evidence_string <<- paste0("(", as.character(evidence$evidence), ")")
        batch_size <<- 100  # Start with a smaller batch size
        num_batches <- ceiling(evidence$num_samples / batch_size)
        
        for (attempt in 1:max_retries) {
          print(paste0("attempt:", attempt))
          
          # Export necessary variables to the cluster
          clusterExport(cl, varlist = c("nodes_string", "evidence_string", "batch_size", "bayesian_network_name"))
          
          # Perform parallel sampling in batches
          query_results <- parLapply(cl, 1:num_batches, function(i) {
            eval(parse(text = paste("bnlearn::cpdist(", bayesian_network_name,
                                    ", nodes = ", nodes_string,
                                    ", evidence = ", evidence_string,
                                    ", n = ", batch_size,
                                    ")")))
          })
          
          # Combine the results from parallel execution
          if (length(query_results) > 0) {
            query_results <- do.call(rbind, query_results)
            accumulated_samples[[attempt]] <- query_results
            total_samples_collected <- total_samples_collected + nrow(query_results)
          }
          
          # Adjust batch size dynamically based on results
          if (total_samples_collected >= evidence$num_samples) {
            break
          } else {
            batch_size <- batch_size*multiplication_factor
          }
        }
        
        # Combine all collected samples
        if (total_samples_collected > 0) {
          combined_samples <- do.call(rbind, accumulated_samples)
          
          if (nrow(combined_samples) >= evidence$num_samples) {
            synthetic_data <- rbind(synthetic_data, combined_samples[seq_len(evidence$num_samples), ])
            print(paste0("Nrow added:", nrow(combined_samples[seq_len(evidence$num_samples), ])))
          } else {
            # Add only the remaining required samples to match evidence$num_samples
            remaining_samples <- evidence$num_samples - nrow(synthetic_data)
            synthetic_data <- rbind(synthetic_data, combined_samples[seq_len(remaining_samples), ])
            warning(sprintf("Only %d samples added due to insufficient combined samples", nrow(combined_samples[seq_len(remaining_samples), ])))
          }
          
          
        } else {
          warning(sprintf("Failed to generate any samples after %d attempts for evidence %s", max_retries, evidence_string))
        }
        
        
        
        
        # Cleanup
        gc()  # Force garbage collection
        rm(accumulated_samples)  # Remove large variables
      }
      
    }
  }, error = function(e) {
    stop("Error during conditional sampling in 'conditional_sampling': ", conditionMessage(e))
  })
  return(synthetic_data)
}



#  Function to run ITERATIVE SYNTHETIC REPLICATION (Loop on Learning + Sampling)
iterative_syngen <- function(data, n_reps, size= nrow(data), whitelist, blacklist, resample = FALSE, max_iter= 1000, seed = seed ) {
  
  
  replicate_results <- vector("list", n_reps)
  
  for (i in seq_len(n_reps)) {
    cat(sprintf("\n▶ Replicate %d / %d\n", i, n_reps))
    
    # ----------------------------------------------------------
    # Step 1. (optional) Bootstrap resampling
    # ----------------------------------------------------------
    if (isTRUE(resample)) {
      data_iter <- data[sample(seq_len(nrow(data)), replace = TRUE), ]
    } else {
      data_iter <- data
    }
    
    # ----------------------------------------------------------
    # Step 2. Learn Bayesian Network via Structural EM
    # ----------------------------------------------------------
    bn_model <- tryCatch({
      structural.em(
        x = data_iter,
        maximize = "hc",
        maximize.args = list(
          whitelist = whitelist,
          blacklist = blacklist,
          score = "bic"
        ),
        fit = "mle",
        impute = "parents",
        return.all = TRUE,
        max.iter = max_iter
      )
    }, error = function(e) {
      cat("⚠️ Error fitting BN in replicate", i, ":", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(bn_model)) next
    # Visualise the updated network (if you're using your helper functions)
    # bn_hc_mle_igraph <- convert_bn_to_igraph(bn_model)
    #visualise_network <- visualise_bn_with_visnetwork(bn_hc_mle_igraph, TRUE)
    #print(visualise_network)
    
    # ----------------------------------------------------------
    # Step 3. Generate synthetic dataset (unconditional sampling)
    # ----------------------------------------------------------
    syn_data <- tryCatch({
      rbn(bn_model$fitted, n = size)
    }, error = function(e) {
      cat("⚠️ Sampling failed in replicate", i, ":", e$message, "\n")
      return(NULL)
    })
    
    replicate_results[[i]] <- syn_data
  }
  
  attr(replicate_results, "scenario") <- "Synthetic_Replication"
  return(replicate_results)
}

# Function to run boosting validation experiment
boosting_validation_experiment <- function(data, group_value, drop_frac = 0.5, group_var = "ethncty", max_iter = 5, whitelist, blacklist) {
  
  # 1. Remove fraction of patients in the chosen group
  set.seed(123)
  group_idx <- which(data[[group_var]] == group_value)
  drop_idx  <- sample(group_idx, size = round(length(group_idx) * drop_frac))
  data_reduced <- data[-drop_idx, ]
  
  cat("Original", group_var, "=", group_value, ":", length(group_idx), 
      "| After removal:", sum(data_reduced[[group_var]] == group_value), "\n")
  
  # 2. Re-learn BN on reduced data
  bn_reduced <- structural.em(
    x = data_reduced,
    maximize = "hc",
    maximize.args = list(
      whitelist = whitelist,
      blacklist = blacklist,
      score = "bic"
    ),
    fit = "mle",
    impute = "parents",
    return.all = TRUE,
    max.iter = 1000
  )
  
  # 3. Conditional sample back missing patients using cpdist() with retry loop
  needed <- length(drop_idx)
  collected <- data.frame()
  attempt <- 1
  
  evidenceString <- paste0("(", group_var, " == '", group_value, "')")
  
  while (nrow(collected) < needed && attempt <= max_iter) {
    cpdist_call <- paste0(
      "cpdist(bn_reduced$fitted, nodes = names(data_reduced), evidence = ", 
      evidenceString, ", n = ", needed*10, ")"
    )
    tmp <- try(eval(parse(text = cpdist_call)), silent = TRUE)
    
    if (inherits(tmp, "try-error") || nrow(tmp) == 0) {
      cat("Iteration", attempt, ": no rows generated. Retrying...\n")
    } else {
      collected <- rbind(collected, tmp)
      cat("Iteration", attempt, ": generated", nrow(tmp), 
          "rows | total =", nrow(collected), "/", needed, "\n")
    }
    attempt <- attempt + 1
  }
  
  # 4. Adjust if too few or too many rows
  if (nrow(collected) < needed) {
    warning(sprintf("Only %d/%d samples generated after %d iterations", 
                    nrow(collected), needed, max_iter))
  }
  if (nrow(collected) > needed) {
    cat(sprintf("Oversampled: %d/%d rows. Downsampling...\n", nrow(collected), needed))
    set.seed(123)
    collected <- collected[sample(1:nrow(collected), needed), ]
  }
  
  # 5. Combine reduced + boosted
  data_boosted <- rbind(data_reduced, collected)
  
  cat("Boosted dataset size:", nrow(data_boosted), "\n")
  cat("Recovered", group_var, "=", group_value, ":", 
      sum(data_boosted[[group_var]] == group_value), "\n\n")
  
  return(list(
    reduced = data_reduced,
    boosted = data_boosted,
    bn_reduced = bn_reduced
  ))
}

# Function to iterate boosting_validation_experiment and return list of boosted datasets
boosting_replicates <- function(data, group_value, drop_frac, group_var, max_iter = 5, n_reps, whitelist, blacklist) {
  replicate_results <- vector("list", n_reps)
  
  for (i in seq_len(n_reps)) {
    replicate_results[[i]] <- boosting_validation_experiment(
      data = data,
      group_value = group_value,
      drop_frac = drop_frac,
      group_var = group_var,
      max_iter = max_iter, 
      whitelist, 
      blacklist
    )$boosted
  }
  
  return(replicate_results)
}
