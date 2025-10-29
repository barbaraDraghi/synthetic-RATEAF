# ------------------------------------------------------------------
# Define the functions used for Bias Investigation 
# ------------------------------------------------------------------

# Shannon Diversity Index
#
# This function calculates the Shannon Diversity Index for a given probability distribution.
#
# param p A numeric vector representing probabilities.
#
# return A numeric value representing the Shannon Diversity Index.
#
# examples
# p <- c(0.2, 0.3, 0.5)
# shannon_index(p)
#
shannon_index <- function(p) {
  
  
  # Input Validation
  tryCatch({
    if (!is.numeric(p)) {
      stop("Error in 'shannon_index': 'p' must be a numeric vector.")
    }
    if (any(p < 0) || any(p > 1)) {
      stop("Error in 'shannon_index': 'p' values must be between 0 and 1 (inclusive).")
    }
    if (sum(p) != 1) {
      stop("Error in 'shannon_index': 'p' must sum to 1.")
    }
    
  }, error = function(e) {
    stop("Input validation failed in 'shannon_index': ", conditionMessage(e))
  })
  
  # Shannon Diversity Index Calculation
  tryCatch({
    p <- p[p > 0] # Remove zero values to avoid log(0)
    result <- -sum(p * log(p))
    
    return(result)
    
  }, error = function(e) {
    stop("Error during Shannon Index calculation in 'shannon_index': ", conditionMessage(e))
  })
}


# Density Score Calculation
#
# This function calculates the density score based on the deviation from ideal density. 
# The score is transformed into a percentage, with 100% being ideal.
#
# param density A numeric vector representing density values.
#
# return A numeric value representing the density score as a percentage.
#
# examples
# density <- c(0.1, 0.2, 0.3, 0.4)
# density_score(density)
#
density_score <- function(density) {
  
  # Input Validation
  tryCatch({
    if (!is.numeric(density)) {
      stop("Error in 'density_score': 'density' must be a numeric vector.")
    }
  }, error = function(e) {
    stop("Input validation failed in 'density_score': ", conditionMessage(e))
  })
  
  
  # Density Score Calculation
  tryCatch({
    k <- length(density)
    if (k == 0) {
      stop("Error in 'density_score': 'density' vector is empty.")
    }
    
    ideal_density <- 1 / k
    deviations <- abs(density - ideal_density)
    relative_deviation <- deviations / ideal_density
    mean_deviation <- mean(relative_deviation)
    score <- (1 / (1 + mean_deviation)) * 100  # Transform to percentage
    
    return(max(0, score))  # Ensure non-negative
    
  }, error = function(e) {
    stop("Error during density score calculation in 'density_score': ", conditionMessage(e))
  })
}


# Calculate Coverage Metrics for Intersectional Demographics
#
# This function calculates coverage metrics for intersections of all demographic features in the dataset.
#
# param data A dataframe containing the demographic data.
# param demographics A character vector of demographic attributes.
#
# return A list containing the coverage metrics for the intersection of all demographic features.
#
# examples
# calculate_coverage_metrics_intersection(df, c("gender", "age"))
#
calculate_coverage_metrics_intersection <- function(data, demographics) {
  
  
  # Input Validation
  tryCatch({
    if (!is.data.frame(data)) {
      stop("Error in 'calculate_coverage_metrics_intersection': 'data' must be a dataframe.")
    }
    if (!is.character(demographics)) {
      stop("Error in 'calculate_coverage_metrics_intersection': 'demographics' must be a character vector.")
    }
    if (any(!demographics %in% colnames(data))) {
      stop("Error in 'calculate_coverage_metrics_intersection': Some 'demographics' columns are not present in the dataframe.")
    }
    if (nrow(data) == 0) {
      stop("Error in 'calculate_coverage_metrics_intersection': 'data' contains no rows.")
    }
  }, error = function(e) {
    stop("Input validation failed in 'calculate_coverage_metrics_intersection': ", conditionMessage(e))
  })
  
  
  # Coverage Metrics Calculation
  tryCatch({
    coverage_results <- list()
    
    # Calculate for the intersection of all features
    comb_name <- paste(demographics, collapse = "&")
    coverage <- data %>%
      group_by(across(all_of(demographics))) %>%
      summarize(count = n(), .groups = 'drop') %>%
      mutate(proportion_dataset = count / sum(count))
    
    # Validate if coverage calculation is valid
    if (nrow(coverage) == 0) {
      stop(paste("Error in 'calculate_coverage_metrics_intersection': Coverage calculation returned an empty result for the combination:", comb_name))
    }
    
    # Calculate target population proportions (1 / number of unique combinations)
    unique_combinations <- data %>%
      select(all_of(demographics)) %>%
      distinct() %>%
      nrow()
    
    proportion_target <- 1 / unique_combinations
    
    # Add target population proportions and calculate coverage ratio
    coverage <- coverage %>%
      mutate(proportion_target = proportion_target) %>%
      mutate(coverage_ratio = proportion_dataset / proportion_target) %>%
      mutate(across(where(is.numeric), round, 2))
    
    # Store the coverage results for the intersection of demographics
    coverage_results[[comb_name]] <- as.data.frame(coverage)
    
    return(coverage_results)
    
  }, error = function(e) {
    stop("Error during calculation in 'calculate_coverage_metrics_intersection': ", conditionMessage(e))
  })
}


# Calculate Additional Data Needed for Coverage
#
# This function calculates the additional data needed to achieve a desired coverage ratio or total dataset size
# using different strategies.
#
# param coverage_results A list of data frames with columns `count` and `proportion_target`.
# param desired_coverage_ratio Target coverage ratio (default is 1).
# param expected_proportions Data Frame containing, for each subgroup, the expected proportion
# param desired_data_size Optional. Target total dataset size.
# bias_mode A character string specifying the bias mode (currently only supports 'coverage' and 'expert_knowledge').
# param method Strategy to compute additional data. Options are:
#   \itemize{
#     \item \code{"proportional_adjustment"}: Adjust counts proportionally.
#     \item \code{"preserve_excess"}: Retain over-covered categories and fill shortfalls.
#     \item \code{"match_to_maximum"}: Equalize all categories to the maximum count.
#   }
#
# return A list of data frames with additional data calculations.
##
calculate_additional_data_needed <- function(statistics, 
                                             desired_coverage_ratio = 1, 
                                             expected_proportions = NULL,
                                             desired_data_size = NULL, 
                                             bias_mode,
                                             method = c("proportional_adjustment", "preserve_excess", "match_to_maximum")) {
  
  # Input Validation
  method <- match.arg(method)  # Ensures method is one of the three valid options
  tryCatch({
    if (!is.list(statistics)) {
      stop("Error in 'calculate_additional_data_needed': 'statistics' must be a list.")
    }
    if (!is.numeric(desired_coverage_ratio) || length(desired_coverage_ratio) != 1 || desired_coverage_ratio <= 0) {
      stop("Error in 'calculate_additional_data_needed': 'desired_coverage_ratio' must be a positive numeric value.")
    }
    if (!is.null(desired_data_size) && (!is.numeric(desired_data_size) || length(desired_data_size) != 1 || desired_data_size <= 0)) {
      stop("Error in 'calculate_additional_data_needed': 'desired_data_size' must be a positive numeric value.")
    }
    if (length(statistics) == 0) {
      stop("Error in 'calculate_additional_data_needed': 'statistics' list is empty.")
    }
  }, error = function(e) {
    stop("Input validation failed in 'calculate_additional_data_needed': ", conditionMessage(e))
  })
  
  
  # Main Logic for Calculating Additional Data Needed
  tryCatch({
    additional_data_needed <- list()
    
    # Compute the current total size of the dataset
    current_total_size <- sum(unlist(lapply(statistics, function(df) sum(df$count))))
    
    # Determine the scaling factor based on desired data size, if provided
    if (!is.null(desired_data_size)) {
      scaling_factor <- desired_data_size / current_total_size
    } else {
      scaling_factor <- 1  # No scaling needed if desired_data_size is not provided
    }
    
    
    #########################
    # Coverage-based approach
    #########################
    if(bias_mode == 'coverage')
    {
      for (combination in names(statistics)) {
        
        coverage <- statistics[[combination]]
        
        if (method == "proportional_adjustment") {
          # Logic 1: Proportional Scaling
          coverage <- coverage %>%
            mutate(
              scaled_count = count * scaling_factor,
              additional_data_needed = pmax(0, ceiling((proportion_target * sum(scaled_count) * desired_coverage_ratio) - scaled_count)),
              additional_data_needed_proportion = additional_data_needed / sum(scaled_count)
            )
          
        } else if (method == "preserve_excess") {
          # Logic 2: Retain Excess, Fill Shortfall
          coverage <- coverage %>%
            mutate(
              scaled_count = count,  # Retain current counts
              additional_data_needed = pmax(0, ceiling((proportion_target * sum(scaled_count) * desired_coverage_ratio) - scaled_count)),
              additional_data_needed_proportion = additional_data_needed / sum(scaled_count)
            )
          
        } else if (method == "match_to_maximum") {
          # Logic 3: Equalize to Maximum Count
          
          # Compute the maximum count across all combinations
          max_count <- max(unlist(lapply(statistics, function(df) df$count)))
          
          coverage <- coverage %>%
            mutate(
              additional_data_needed = pmax(0, max_count - count),  # Difference to reach max count
              additional_data_needed_proportion = additional_data_needed / max_count
            )
        }
        
        # Round numeric values for clarity
        coverage <- coverage %>%
          mutate(across(where(is.numeric), \(x) round(x, 2)))
        #mutate(across(where(is.numeric), round, 2))
        
        
        additional_data_needed[[combination]] <- as.data.frame(coverage)
      }
    }
    
    ###########################
    # knowledge-based approach
    ###########################
    if(bias_mode == 'expert_knowledge')
    {
      
      for (combination in names(statistics)) {
        
        statistics_extract <- statistics[[combination]]
        
        if(method == "proportional_adjustment"){
          
          # Extract all columns up to 'proportion_dataset'
          proportions <- statistics_extract[, 1:which(names(statistics_extract) == "proportion_dataset")]
          
          # Merge the two tables on the first column
          merged_table <- merge(proportions, expected_proportions, by = combination)
          
          
          # Calculate the total count after scaling
          total_count <- sum(merged_table$count) * scaling_factor
          
          # Compute target_scaled_count and data_needed independently
          merged_table <- merged_table %>%
            mutate(
              scaled_count = count * scaling_factor,
              additional_data_needed = round(target_proportions * total_count)
            )
          
          
          # Store the results
          additional_data_needed[[combination]] <- as.data.frame(merged_table)
        }
        
        
        
        if (method == "preserve_excess") {
          # Extract all columns up to 'proportion_dataset'
          proportions <- statistics_extract[, 1:which(names(statistics_extract) == "proportion_dataset")]
          
          # Merge the two tables on the first column
          merged_table <- merge(proportions, expected_proportions, by = combination)
          
          # Identify the group with the maximum count
          max_count <- max(merged_table$count)
          max_group <- merged_table[which.max(merged_table$count), ]
          
          # Calculate the new total count based on the preserved maximum count
          new_total_count <- max_count / max_group$target_proportions
          
          # Compute additional data needed for all groups
          merged_table <- merged_table %>%
            mutate(
              scaled_count = count,
              additional_data_needed = ifelse(
                count == max_count,
                0, # No additional data needed for the max group
                round((target_proportions * new_total_count) - count, 0)
              )
            )
          
          # Store the results
          additional_data_needed[[combination]] <- as.data.frame(merged_table)
        }
      }
      
      
    }
    
    
    
    return(additional_data_needed)
    
  }, error = function(e) {
    stop("Error during calculation in 'calculate_additional_data_needed': ", conditionMessage(e))
  })
}

