# =====================================================================
# functions.R — RATE-AF Synthetic Data Project
# Core utility and preprocessing functions for data preparation,
# discretisation, post-processing and visualisation.
# =====================================================================

# ------------------------------------------------
# Utility: Package management and script loading
# ------------------------------------------------

# Ensure necessary libraries are installed and loaded.
check_and_install_packages <- function(packages) {
  installed_packages <- rownames(installed.packages())
  for (pkg in packages) {
    if (!pkg %in% installed_packages) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# Safely sources external scripts with error handling and messages.
safe_source <- function(script_path) {
  tryCatch({
    source(script_path)
    cat("Successfully executed:", script_path, "\n")
  }, error = function(e) {
    cat("Error in executing script:", script_path, "\nError Message:", e$message, "\n")
  })
}


# ------------------------------------------------
# RATE-AF specific — Data loading and preparation
# ------------------------------------------------

# Loads and prepares RATE-AF baseline data from .dta or dataframe input.
rateAF_datapreparation_step1 <- function(data_path = NULL, df = NULL, recode_ethnicity = TRUE) {
  
  
  # Data Loading
  if (!is.null(data_path)) {
    # Case 1: Load data from .dta file
    if (!grepl("\\.dta$", data_path)) {
      stop("Error: The provided file is not a .dta file.")
    }
    data <- read_dta(data_path)
    
  } else if (!is.null(df)) {
    # Case 2: Work on an existing dataframe
    data <- df
  } else {
    stop("Error: Please provide either a data path or a dataframe.")
  }
  
  # Ensure data is in dataframe format
  data <- as.data.frame(data)  
  
  
  # Set correct data type
  # Pre-processed based on data specification
  {
    
    
    
    variable_names <- c("pid", 
                        "alloc", "age", "sex", 
                        "ethncty", 
                        "sodium", "urea", "creat", "egfr", "hb","ntprobnp", "ntprobnp_ln", 
                        "anticoag", "antiplat", "antihypertensive", "acei", "arb",
                        "acearb", 
                        "thzdlpdrtcs", "ccbs", "alphablckrs", "aldstrnantgnsts", "otherhtn", "inhlrsarwydis",
                        "yearsaf", 
                        "mdfdehrascr", 
                        "hf", 
                        "hfdiag", "nyha", "mihist", "pci", "cabg", "strokehist",
                        "tiahist", "smkngstatus", "mrthn8alchl", "hvasthma", "hvecopdemphysma", "hvlvrdisease","typ1dm", "typ2dm", "insln", "hsrcvdtrtmntthyrd", "hshdmjrblds", "ppm", 
                        "ecghr", 
                        "estmtdejctnfrctn",
                        "bpsysav", "bpdiaav", 
                        "signshf", 
                        "ht", "wt", "waist", 
                        "mmse", 
                        "death", "death_date","ecghr12m", "bpsysav12m", "bpdiaav12m", "ntprobnp12m", "ntprobnp12m_ln", "nyhaimprove", 
                        "nyha12m",
                        "mehra2class12m", "sae_yn", "sae_n_perpatient", "ae_yn", "gp_cv_yn", "gp_cv", "gp_noncv_yn","gp_noncv", "cvadj")
    
    # List of variable names and types
    variable_types <- c("double", 
                        "float", "float", "float", 
                        "byte",
                        "int", "double", "int", "int", "int","int", "float", 
                        "byte", "byte", "byte", "byte", "byte", 
                        "float", 
                        "byte", "byte", "byte", "byte", "byte", "byte", 
                        "float", 
                        "byte", 
                        "float", 
                        "byte", "byte", "byte", "byte", "byte", "byte", 
                        "byte", "byte", "byte", "byte", "byte", "byte", 
                        "byte", "byte", "byte", "byte", "byte", "byte", 
                        "int", "byte", "float", "float", 
                        "byte",
                        "int", "int", "int", 
                        "byte", 
                        "float", "float", "int", "float", "float", "long", "float","float", 
                        "byte", 
                        "float", "float", "float", "float", "float", "float", "float", "float","float")
    
    # Convert each variable to the desired data type
    for (i in seq_along(variable_types)) {
      var_name <- names(data)[i]
      var_type <- variable_types[i]
      
      # Check if the variable exists in the dataframe
      if (var_name %in% colnames(data)) {
        # Convert the variable to the desired data type
        if (var_type == "double") {
          data[[var_name]] <- as.double(data[[var_name]])
        } else if (var_type == "float") {
          data[[var_name]] <- as.double(data[[var_name]])
        #} else if (var_type == "byte") {
          #data[[var_name]] <- as.factor(data[[var_name]])
        } else if (var_type == "int") {
          data[[var_name]] <- as.integer(data[[var_name]])
        } else if (var_type == "long") {
          data[[var_name]] <- as.integer(data[[var_name]])
        }
      }
    }
    
    # Display the structure of the updated dataframe
    data <- as.data.frame(data)
    
    
    
  }
  
  #List of variable names to convert to factors
  variables_to_factorise <- c('alloc', 'sex', 'ethncty', 'anticoag', 'antiplat', 'antihypertensive', 'acei', 'arb', 'thzdlpdrtcs',
                              'ccbs', 'alphablckrs', 'aldstrnantgnsts', 'otherhtn', 'inhlrsarwydis', 'hf', 'hfdiag', 'nyha',
                              'mihist', 'pci', 'cabg', 'strokehist', 'tiahist', 'mrthn8alchl', 'hvasthma', 'hvecopdemphysma',
                              'hvlvrdisease', 'typ1dm', 'typ2dm', 'insln', 'hsrcvdtrtmntthyrd', 'hshdmjrblds', 'ppm', 'signshf',
                              'nyha12m', 'mehra2class12m', 'sae_yn', 'ae_yn', 'gp_cv_yn', 'gp_noncv_yn', 'cvadj','acearb', 'death', 'mdfdehrascr','smkngstatus','estmtdejctnfrctn')
  
  
  # Convert specified variables to factors
  data <-convert_to_factors(data,variables_to_factorise)
  
  # No of years since AF diagnosis: shouldn't it be positive?
  data$yearsaf <- abs(data$yearsaf)
  data$yearsaf <- as.integer(data$yearsaf)
  data$age <- as.integer(data$age)

  # Ensure 0 and 1 (instead of 1 and 2)
  data <- reshape_binary_variables(data)

  # Remove metadata (inherited from Stata)
  data <- remove_metadata_from_numeric_variables(data)
  
  # Remove pid
  data$pid <- NULL
  
  # Remove death date 
  data$death_date <- NULL
  
  # Re-code ethnicity
  if(recode_ethnicity){
    data <- update_ethnicity_coding(data)
  }

  
  return(data)
}

# Converts binary variables with levels (1,2) to (0,1).
reshape_binary_variables <- function(df) {
  for (col_name in names(df)) {
    # Extract the column
    col <- df[[col_name]]
    
    # Handle factors
    if (is.factor(col)) {
      levels_col <- as.numeric(levels(col))  # Extract numeric levels
      
      # Case 1: Levels are exactly 1 and 2
      if (all(levels_col == c(1, 2))) {
        # Map 1 -> 0 and 2 -> 1
        df[[col_name]] <- factor(as.numeric(col) - 1, levels = c(0, 1))
        
        # Case 2: Only one level (e.g., 2)
      } else if (length(levels_col) == 1 && levels_col == 2) {
        # Recode everything to 1
        df[[col_name]] <- factor(rep(1, length(col)), levels = c(0, 1))
      }
    } else if (is.numeric(col)) {
      # Handle numeric columns
      
      unique_vals <- unique(col)
      
      # Case 1: Numeric values are 1 and 2
      if (all(unique_vals %in% c(1, 2, NA))) {
        # Convert 1 -> 0 and 2 -> 1
        df[[col_name]] <- factor(col - 1, levels = c(0, 1))
        
        # Case 2: Only one unique value (e.g., 2)
      } else if (length(unique_vals) == 1 && unique_vals == 2) {
        # Recode everything to 1
        df[[col_name]] <- factor(rep(1, length(col)), levels = c(0, 1))
      }
    }
  }
  return(df)
}

# Converts selected variables to factors, preserving Stata labels if available.
convert_to_factors <- function(data, variables) {
  
  for (var in variables) {
    if (var %in% colnames(data)) {
      
      # Check if the variable is a byte (numeric with labels)
      if (!is.null(attributes(data[[var]])$labels)) {
        
        # Extract numeric levels from labels
        var_labels <- attributes(data[[var]])$labels
        levels_numeric <- as.numeric(var_labels)
        # Convert to factor while preserving levels
        data[[var]] <- factor(data[[var]], levels = levels_numeric)
      } else{
        # Default to factor for numeric/character variables without labels
        data[[var]] <- as.factor(data[[var]])
      } 
    } else {
      warning(paste("Variable", var, "not found in data."))
    }
  }
  return(data)
}

# Removes derived or constant variables to avoid redundancy.
remove_variables <- function(data){
  
  
  # --------------------------------------------- 
  # Remove derived/computed variables
  # ---------------------------------------------
  # mehra2class12m: I think this variable can be kept because we don't have the measure at 12m but only the baseline in the dataset (so we won't be able to compute it!)
  # nyhaimprove: I think this variable is computed using nyha and nyha12m, so we should omit nyhaimprove from the dataset and compute it afterward
  # ntprobnp_ln: ln(ntprobnp), can be computed afterward
  # ntprobnp12m_ln: ln(ntprobnp12m), can be computed afterward
  derived_variables <- c('ntprobnp_ln','ntprobnp12m_ln', 'nyhaimprove')
  data[,derived_variables] <- NULL
  
  
  
  
  
  
  # --------------------------------------------- 
  # Remove Constant Columns
  # ---------------------------------------------
  # Rule
  # - If the column has NA values, we can impute random values to introduce variation
  # - If the column doesn't have NA values, but it is a constant column, we can try to omit it (?)
  
  for (col in colnames(data)) {
    if (length(unique(data[[col]])) == 1) {
      # If the column doesn't have NA values and it's a constant column, we can try to omit it
      cat("Column", col, "is a constant column. Omitting it from the data.\n")
      data[[col]] <- NULL
    }
  }
  
  return(data)
  
  
}

# Ensures minimal variation by imputing missing values in near-constant columns.
deal_with_variation_lack <- function(data){
  
  # Rule
  # - If the column has NA values, we can impute random values to introduce variation
  
  for (col in colnames(data)) {
    if(length(unique(na.omit(data[[col]])))<2){
      if (anyNA(data[[col]])) {
        # Impute random values to introduce variation
        data[[col]] <- as.character(data[[col]])
        #data[[col]][is.na(data[[col]])] <- as.factor(sample(c(0, 1), size = length(data[[col]][is.na(data[[col]])]), replace = TRUE))
        data[[col]][is.na(data[[col]])] <- 0
        # Convert character column back to factor
        if (is.character(data[[col]])) {
          data[[col]] <- as.factor(data[[col]])
        }
      }
    }
    
  }
  
  return(data)
  
}

# Removes Stata metadata attributes (labels, formats) from numeric columns.
remove_metadata_from_numeric_variables <- function(data, attributes_to_remove = c("label", "format.stata")) {
  for (col_name in names(data)) {
    if (is.numeric(data[[col_name]])) {
      for (attr in attributes_to_remove) {
        if (!is.null(attr(data[[col_name]], attr))) {
          cat("Removing attribute", attr, "from", col_name, "\n")
          attr(data[[col_name]], attr) <- NULL
        }
      }
    }
  }
  return(data)
}

# Performs additional cleaning and type conversion for Bayesian analysis.
rateAF_datapreparation_step2 <- function(data){
  
  # Deal with variational lack
  data <- deal_with_variation_lack(data)
  
  # Remove variables
  data <- remove_variables(data)
  
  # Convert Integer to Numeric (for later mle-cg)
  integer_columns <- colnames(select_if(data, is.integer))
  data[, integer_columns] <- lapply(data[, integer_columns], function(x) as.numeric(as.character(x)))
  
 
  return(data)
}

# Adds a combined cardiovascular endpoint (CVEP) indicator from related columns.
add_cvep_column <- function(df) {
  # Check if the required columns are present
  if (!all(c("gp_cv_yn", "cvadj") %in% colnames(df))) {
    stop("The dataframe must contain 'gp_cv_yn' and 'cvadj' columns.")
  }
  
  # Ensure the input columns are factors with levels 0 and 1
  df$gp_cv_yn <- factor(df$gp_cv_yn, levels = c(0, 1))
  df$cvadj <- factor(df$cvadj, levels = c(0, 1))
  
  # Add the cvep column, explicitly handling NAs
  df$cvep <- factor(
    ifelse(
      is.na(df$gp_cv_yn) & is.na(df$cvadj), 0, # If both are NA, set to 0
      ifelse(
        (!is.na(df$gp_cv_yn) & as.numeric(as.character(df$gp_cv_yn)) == 1) |
          (!is.na(df$cvadj) & as.numeric(as.character(df$cvadj)) == 1), 1, 0
      )
    ),
    levels = c(0, 1)
  )
  
  return(df)
}

# Re-codes ethnicity variable into biologically meaningful categories.
update_ethnicity_coding <- function(data) {
  # Ensure ethncty column exists and is a factor
  if (!"ethncty" %in% colnames(data)) stop("Column 'ethncty' not found in the dataframe.")
  if (!is.factor(data$ethncty)) stop("'ethncty' column must be a factor.")
  
  # Convert factor to numeric temporarily
  ethncty_numeric <- as.numeric(as.character(data$ethncty))
  
  # Define biologically similar groupings with numeric labels
  group_mapping <- list(
    white   = list(codes = c(31, 32, 33, 34), label = 1),
    asian   = list(codes = c(39, 40, 41, 42, 43), label = 2),
    black   = list(codes = c(44, 45, 46), label = 3),
    mixed   = list(codes = c(35, 36, 37, 38), label = 4),
    other   = list(codes = c(47, 48, 98), label = 5),
    unknown = list(codes = c(99), label = 6)
  )
  
  # Flatten the group mapping for matching
  all_codes <- unlist(lapply(group_mapping, function(x) x$codes))
  all_labels <- rep(names(group_mapping), times = sapply(group_mapping, function(x) length(x$codes)))
  all_numeric_labels <- rep(1:length(group_mapping), times = sapply(group_mapping, function(x) length(x$codes)))
  
  # Map ethnicity_group and ethnicity_code
  ethnicity_group <- ifelse(!is.na(match(ethncty_numeric, all_codes)), 
                            all_labels[match(ethncty_numeric, all_codes)], 
                            "unmapped")
  
  ethnicity_code <- ifelse(!is.na(match(ethncty_numeric, all_codes)), 
                           all_numeric_labels[match(ethncty_numeric, all_codes)], 
                           7)
  
  # Add the new columns back to the dataframe as factors
  # By directly converting ethnicity_group and ethnicity_code to factors without explicitly assigning levels, 
  # only the levels present in the data will appear. This avoids adding unused levels with 0 counts in tables.
  #data$ethnicity_group <- factor(ethnicity_group)  
  #data$ethnicity_code <- factor(ethnicity_code)   
  data$ethncty <- factor(ethnicity_code)   
  # If you later decide you want all levels explicitly (including unused ones), 
  # you can reintroduce the levels argument in factor():
  # data$ethnicity_group <- factor(ethnicity_group, levels = c("white", "asian", "black", "mixed", "other", "unknown", "unmapped"))
  # data$ethnicity_code <- factor(ethnicity_code, levels = c(1, 2, 3, 4, 5, 6, 7))
  
  return(data)
}

# ------------------------------------------------
# RATE-AF specifics - PostProcessing
# ------------------------------------------------

# Recomputes derived clinical variables (e.g., NT-proBNP logs, NYHA improvement).
compute_derived_variables <- function(data) {
  
  # 1. ntprobnp_ln: Log of ntprobnp
  if ("ntprobnp" %in% colnames(data)) {
    data$ntprobnp_ln <- log1p(pmax(data$ntprobnp, 0))#log(data$ntprobnp)
  } else {
    warning("ntprobnp column not found in the dataset. ntprobnp_ln not computed.")
  }
  
  
  # 2. ntprobnp12m_ln: Log of ntprobnp12m
  if ("ntprobnp12m" %in% colnames(data)) {
    data$ntprobnp12m_ln <- log1p(pmax(data$ntprobnp12m, 0))#log(data$ntprobnp12m)
  } else {
    warning("ntprobnp12m column not found in the dataset. ntprobnp12m_ln not computed.")
  }
  
  # 3. nyhaimprove: Compute based on nyha and nyha12m
  if ("nyha" %in% colnames(data) & "nyha12m" %in% colnames(data)) {
    
    # Convert nyha and nyha12m to numeric locally without modifying original columns
    nyha_numeric <- as.numeric(as.character(data$nyha))
    nyha12m_numeric <- as.numeric(as.character(data$nyha12m))
    
    # Compute nyhaimprove
    nyhaimprove_numeric <- ifelse(
      is.na(nyha_numeric) | is.na(nyha12m_numeric),
      NA, # Assign NA if either value is missing
      ifelse(nyha12m_numeric <= nyha_numeric - 1, 1, 0) # 1 for improvement, 0 otherwise
    )
    
    # Convert nyhaimprove to a factor
    data$nyhaimprove <- factor(
      nyhaimprove_numeric,
      levels = c(0, 1)
    )
  } else {
    warning("nyha or nyha12m column not found in the dataset. nyhaimprove not computed.")
    data$nyhaimprove <- factor(NA, levels = c(0, 1), labels = c("No", "Yes"))
  }
  
  return(data)
}

# Restores specific constant columns removed during data generation.
add_constant_columns <- function(original_data,synthetic_data){

  # Ensure same number of rows
  if (nrow(original_data) == nrow(synthetic_data)) {
    
    synthetic_data$hf <- original_data$hf
    synthetic_data$typ1dm <- original_data$typ1dm
    
    
  } else {
    # If different number of rows, repeat rows of `original_data` to match rows of `synthetic_data`
    original_data = original_data[,c('hf','typ1dm')]
    original_repeated <- original_data[rep(1:nrow(original_data), length.out = nrow(synthetic_data)), ]
    
    # Add the repeated columns from `original_data` to `synthetic_data`
    synthetic_data <- cbind(synthetic_data, original_repeated)
  }
  
  
  
  return(synthetic_data)
 
}


# ---------------------------------------------------
# Data preprocessing for Bayesian Networks Modelling
# ---------------------------------------------------

# Discretises numeric variables adaptively for Bayesian network learning.
discretise_RATEAF_adaptive <- function(df, 
                                default_quantiles = 10,
                                max_bins = 20,
                                min_bin_freq = 0.03,
                                log_transform_vars = c("ntprobnp", "creatinine"),
                                keep_original = FALSE) {
  
  numeric_columns <- colnames(dplyr::select_if(df, is.numeric))
  pattern <- paste(log_transform_vars, collapse = "|")
  
  for (col in numeric_columns) {
    vals <- df[[col]]
    unique_vals <- unique(vals)
    
    # Binary
    if (all(unique_vals %in% c(0, 1)) && length(unique_vals) == 2) {
      df[[col]] <- factor(vals, levels = c(0,1), labels = c("0","1"))
      next
    }
    
    # Low-cardinality numeric
    if (length(unique_vals) <= 6) {
      df[[col]] <- factor(vals, levels = sort(unique_vals))
      next
    }
    
    # Adaptive bin count
    n_q <- ifelse(stringr::str_detect(tolower(col), pattern),
                  max_bins, default_quantiles)
    
    probs <- seq(0, 1, length.out = n_q + 1)
    qtiles <- quantile(vals, probs = probs, na.rm = TRUE)
    breaks <- unique(qtiles)
    
    new_col <- paste0(col, "_discretised")
    df[[new_col]] <- cut(vals, breaks = breaks, include.lowest = TRUE, right = FALSE)
    
    # Merge small bins
    freq <- table(df[[new_col]]) / nrow(df)
    small <- names(freq[freq < min_bin_freq])
    if (length(small) > 0 && length(freq) > 3) {
      df[[new_col]] <- as.character(df[[new_col]])
      df[[new_col]][df[[new_col]] %in% small] <- "merged_tail"
      df[[new_col]] <- factor(df[[new_col]])
    }
    
    if (!keep_original) df[[col]] <- NULL
    
    message(sprintf("Discretised %s into %d bins", col, nlevels(df[[new_col]])))
  }
  
  # Convert residual numeric columns
  df[] <- lapply(df, function(x) if (is.numeric(x)) factor(x) else x)
  
  return(df)
}




# Replaces NA values with a placeholder level ('3') for selected columns.
convert_na_to_null <- function(data, columns) {
  data[columns] <- lapply(data[columns], function(x) {
    if (is.factor(x)) {
      levels(x) <- c(levels(x), '3')
      x[is.na(x)] <- '3'
    } else if (is.character(x)) {
      x[is.na(x)] <- '3'
    } else {
      x <- as.character(x)
      x[is.na(x)] <- '3'
      x <- as.factor(x)
    }
    return(x)
  })
  return(data)
}



# ------------------------------------------------
# Post-processing of synthetic data
# ------------------------------------------------

# Reconstructs continuous variables from discretised bins in synthetic data.
post_processing <- function(original_data, synthetic_data) {
  
  # Identify columns ending with '_discretised'
  discretised_cols <- grep("_discretised", names(synthetic_data), value = TRUE)
  
  for (col in discretised_cols) {
    # Extract the base variable name
    base_name <- sub("_discretised", "", col)
    
    # Ensure the base variable exists in the original data
    if (!base_name %in% names(original_data)) {
      stop(paste("Column", base_name, "not found in original data."))
    }
    
    # Get original column
    original_column <- original_data[[base_name]]
    
    # convert discretised labels to numeric bin indices
    disc_col <- synthetic_data[[col]]
    
    if (is.factor(disc_col)) {
      # if factor levels are numeric-like ("1","2","3"...) keep them as numeric
      if (all(grepl("^[0-9]+$", levels(disc_col)))) {
        disc_values <- as.numeric(as.character(disc_col))
      } else {
        # otherwise use the internal integer codes (1..nlevels)
        disc_values <- as.numeric(disc_col)
      }
    } else {
      disc_values <- suppressWarnings(as.numeric(disc_col))
    }
    
    n_bins <- length(unique(na.omit(disc_values)))
    
    # Compute quantiles of the original column
    probs <- seq(0, 1, length.out = n_bins + 1)
    quantile_intervals <- quantile(original_column, probs = probs, na.rm = TRUE)
    quantile_breaks <- unique(quantile_intervals)
    quantiles_values <- as.numeric(quantile_breaks)
    
    # Check for invalid discretised values
    if (any(disc_values < 1 | disc_values > n_bins, na.rm = TRUE)) {
      warning(paste("Invalid discretised values detected in column", col, "- replacing with 1"))
      disc_values[disc_values < 1 | disc_values > n_bins] <- 1
    }
    
    # Map bins to numeric ranges
    min_vals <- quantiles_values[disc_values]
    max_vals <- quantiles_values[disc_values + 1]
    
    # Sample continuous values within each bin
    sampled_values <- rep(NA, length(disc_values))
    valid_idx <- which(!is.na(min_vals) & !is.na(max_vals))
    sampled_values[valid_idx] <- runif(
      n = length(valid_idx),
      min = min_vals[valid_idx],
      max = max_vals[valid_idx]
    )
    
    # Round if original column is integer-like
    if (is.integer(original_column) || all(original_column == floor(original_column), na.rm = TRUE)) {
      sampled_values <- round(sampled_values)
    }
    
    # Assign reconstructed continuous values
    synthetic_data[[base_name]] <- sampled_values
    
    # Remove discretised column
    synthetic_data[[col]] <- NULL
  }
  
  # Match original data type
  synthetic_data <- match_gt_data_types(original_data, synthetic_data)
  
  # Ensure numeric precision consistency
  numeric_cols <- sapply(original_data, is.numeric)
  common_numeric_cols <- intersect(names(original_data)[numeric_cols], names(synthetic_data)) 
  for (col_name in common_numeric_cols) {
    if (any(original_data[[col_name]] %% 1 != 0, na.rm = TRUE)) {
      precision <- max(sapply(original_data[[col_name]][!is.na(original_data[[col_name]])], function(x) {
        decimal_part <- sub("^[^.]*\\.", "", as.character(x))
        if (decimal_part == as.character(x)) 0 else nchar(decimal_part)
      }), na.rm = TRUE)
      synthetic_data[[col_name]] <- round(synthetic_data[[col_name]], digits = precision)
    }
  }
  
  # Add constant columns that were removed
  synthetic_data <- add_constant_columns(original_data, synthetic_data)
  
  # Compute derived variables
  synthetic_data <- compute_derived_variables(synthetic_data)
  
  # Align factor levels
  synthetic_data <- align_factors(original_data, synthetic_data)
  
  # Add any additional columns that have been removed
  unique_columns <- setdiff(colnames(original_data), colnames(synthetic_data))
  unique_original <- original_data[, unique_columns, drop = FALSE]
  
  if (nrow(original_data) == nrow(synthetic_data)) {
    synthetic_data <- cbind(synthetic_data, unique_original)
  } else {
    original_repeated <- original_data[rep(1:nrow(original_data), length.out = nrow(synthetic_data)), ]
    synthetic_data <- cbind(synthetic_data, original_repeated[, unique_columns, drop = FALSE])
  }
  
  # Final type alignment
  synthetic_data <- match_gt_data_types(original_data, synthetic_data)
  
  return(synthetic_data)
}

# Aligns factor levels between original and synthetic datasets.
align_factors <- function(original_data, synthetic_data) {
  # Iterate over the columns of the original dataframe
  for (col_name in names(original_data)) {
    # Check if the column is a factor and has more than 2 levels
    if (is.factor(original_data[[col_name]]) && length(levels(original_data[[col_name]])) > 2) {
      
      # Get the levels and unique values from the original data
      original_levels <- levels(original_data[[col_name]])
      original_unique_values <- unique(original_data[[col_name]])
      
      # Check if the column exists in the synthetic data and is a factor
      if (col_name %in% names(synthetic_data) && is.factor(synthetic_data[[col_name]])) {
        
        # Get the levels and unique values from the synthetic data
        synthetic_levels <- levels(synthetic_data[[col_name]])
        synthetic_unique_values <- unique(synthetic_data[[col_name]])
        
        # Check if the synthetic unique values match the original unique values
        if (!all(synthetic_unique_values %in% original_unique_values)) {
          cat("Mismatch in unique values for column:", col_name, "\n")
          cat("Original unique values: ", paste(original_unique_values, collapse = ", "), "\n")
          cat("Synthetic unique values: ", paste(synthetic_unique_values, collapse = ", "), "\n")
          
          # Filter out synthetic values not in original unique values
          synthetic_data[[col_name]] <- factor(
            synthetic_data[[col_name]],
            levels = original_levels
          )
          
          # Rescale the synthetic column to match original data's unique values
          synthetic_data[[col_name]][!synthetic_data[[col_name]] %in% original_unique_values] <- NA
        }
      }
    }
  }
  
  return(synthetic_data)
}


# Matches data types (numeric, integer, factor) between ground-truth and synthetic datasets.
match_gt_data_types <- function(gt_data, syn_data) {
  # Iterate over each column in gt_data
  for (col_name in colnames(gt_data)) {
    
    # Get the column type from gt_data
    target_type <- class(gt_data[[col_name]]) 
    
    # Convert the corresponding column in syn_data to the same type as gt_data
    if (col_name %in% colnames(syn_data)) {
      if (target_type == "numeric") {
        syn_data[[col_name]] <- as.numeric(as.character(syn_data[[col_name]]))
      } else if (target_type == "integer") {
        syn_data[[col_name]] <- as.integer(as.character(syn_data[[col_name]]))
      } else if (target_type == "factor") {
        # Convert to factor with the same levels as gt_data
        syn_data[[col_name]] <- factor(syn_data[[col_name]], levels = levels(gt_data[[col_name]]))
      } else {
        # Default conversion for other types
        syn_data[[col_name]] <- as(syn_data[[col_name]], target_type)
      }
    }
  }
  
  return(syn_data)
}


# ------------------------------------------------
# Bayesian Network visualisation
# ------------------------------------------------

# Converts a bnlearn Bayesian network model into an igraph object.
convert_bn_to_igraph <- function(bn) {
  bn_igraph <- bnlearn::as.igraph(bn$dag)
  return(bn_igraph)
}

# Visualises a Bayesian network using visNetwork.
visualise_bn_with_visnetwork <- function(bn_igraph, node_id_selection = TRUE) {
  
  
  
  nodes <- data.frame(id = V(bn_igraph)$name, label = V(bn_igraph)$name)
  edges <- data.frame(
    from = get.edgelist(bn_igraph)[,1],
    to = get.edgelist(bn_igraph)[,2],
    #total_occurrences = E(bn_igraph)$weight,
    width = 2
  )
  
  # Preparing edges with additional attributes for visualisation
  vis_network <- visNetwork(nodes, edges) %>%
    visEdges(arrows = 'to', color = list(color = '#848484', highlight = '#848484'), label = edges$title, labelHighlightBold = TRUE) %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = node_id_selection) %>%
    visLayout(randomSeed = 123, improvedLayout = TRUE) %>%
    visNodes(shape = 'dot', scaling = list(label = list(enabled = TRUE))) %>%
    visInteraction(navigationButtons = TRUE, keyboard = TRUE) %>%
    visPhysics(solver = 'forceAtlas2Based') %>%
    visEdges(smooth = FALSE)
  
  return(vis_network)
}
