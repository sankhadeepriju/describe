#' @title Categorical Analysis by Grouping
#'
#' @description Reports Chi-sq/Fisher's Exact test to compare distribution of subjects across a Grouping Variable.
#'
#' @param df A data set
#' @param categorical_col (character; column name) Categorical or Grouping column to compare
#' @param factor_columns (character vector of column names) Other categorical columns to create summary of the variables chosen with their distribution compared with respect to the previous grouping column
#'
#' @return A summary table
#'
#' @export analyze_category_by_category
#' @importFrom dplyr "%>%"
#' @name analyze_category_by_category



library(dplyr)
library(tidyverse)

analyze_category_by_category <- function(df, categorical_col, factor_columns) {
  factor_cols <- sort(setdiff(factor_columns, c("Name", "PatientID", categorical_col)))

  # Formula Generator
  generate_formula <- function(response_var, predictor_vars) {
    if (length(predictor_vars) == 0) {
      stop("Predictor variables cannot be empty.")
    }

    formula_str <- paste(response_var, "~", paste(predictor_vars, collapse = " + "))
    formula <- as.formula(formula_str)

    return(formula)
  }

  # Define the name of the response variable
  response_var <- categorical_col

  # Create the categorical/factor col formula
  formula <- generate_formula(response_var, factor_cols)
  res <- univariateTable(formula, data = df)

  summary_res <- summary(res)

  ## fetch p value by fisher/chi sq accordingly
  get_p_value <- function(df, cat1, cat2) {
    # Create the contingency table
    contingency_table <- table(df[[cat1]], df[[cat2]])

    # Calculate sample size and expected frequencies, suppressing warnings
    sample_size <- sum(contingency_table)
    expected_frequencies <- suppressWarnings(chisq.test(contingency_table)$expected)

    # Check conditions for Fisher's exact test
    if (sample_size < 30 || any(expected_frequencies < 5)) {
      # Use Fisher's exact test
      p_value <- paste(suppressWarnings(fisher.test(contingency_table)$p.value), "#")
    } else {
      # Use Chi-square test
      p_value <- paste(suppressWarnings(chisq.test(contingency_table, correct = FALSE)$p.value), "*")
    }

    return(p_value)
  }

  original_vec <- summary_res$`p-value`
  new_values <- sapply(X = factor_cols, FUN = function(x) { get_p_value(df = df,
                                                                              cat1 = response_var,
                                                                              cat2 = x)})
  # Define a function to process the vector
  format_numeric_string <- function(x) {
    num_part <- as.numeric(gsub("[^0-9.-]", "", x))
    symbol_part <- gsub("[0-9.-]", "", x)

    rounded_num <- ifelse(abs(num_part) < 1e-04, "< 1e-04", format(round(num_part, 4), nsmall = 4))
    return(paste0(rounded_num, symbol_part))
  }

  formatted_p_vals <- sapply(new_values, format_numeric_string)

  numeric_indexes <- which(original_vec != "")
  original_vec[numeric_indexes] <- formatted_p_vals

  summary_res$`p-value` <- original_vec

  colnames(summary_res)[3:4] <- paste0(colnames(summary_res)[3:4], "<br>", "N (%)")

  return(summary_res)
}

#data = readxl::read_excel("R/descriptives_data.xlsx")
#names(data)
#analyze_category_by_category(df = data, categorical_col = "Arm_of_randomization", factor_columns = c("Gender", "Comorbidities"))
