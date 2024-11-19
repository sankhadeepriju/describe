#' @title Numerical Analysis by Grouping
#'
#' @description Reports Group Comparison like t-test/Mann WHitney/ANOVA/Kruskal-Wallis to compare mean difference across a Grouping Variable.
#'
#' @param df A data set
#' @param categorical_col (character; column name) Categorical or Grouping column to compare
#' @param numeric_cols (character vector of column names) Other numerical columns to create summary table of the variable means copared with respect to the previous grouping column
#'
#' @return A summary table
#'
#' @export analyze_numeric_by_category
#' @name analyze_numeric_by_category



# Define the function
analyze_numeric_by_category <- function(df, categorical_col, numeric_cols) {
  # Check if the categorical column exists
  if (!categorical_col %in% names(df)) {
    stop("Categorical column not found in dataframe.")
  }

  unique_vals <- unique(df[[categorical_col]])

  results <- list()

  for (col in numeric_cols) {
    # Calculate summary statistics
    cat_data <- df %>% group_by(across(all_of(categorical_col))) %>% summarise(
      mean = round(mean(!!sym(col), na.rm = TRUE), digits = 2),
      sd = round(sd(!!sym(col), na.rm = TRUE), digits = 2),
      median = round(median(!!sym(col), na.rm = TRUE), digits = 2),
      Q1 = round(quantile(!!sym(col), probs = 0.25, na.rm = TRUE), digits = 2),
      Q3 = round(quantile(!!sym(col), probs = 0.75, na.rm = TRUE), digits = 2),
      count = n(),
      .groups = 'drop'
    )

    overall_data <- df %>% summarise(
      mean = round(mean(!!sym(col), na.rm = TRUE), digits = 2),
      sd = round(sd(!!sym(col), na.rm = TRUE), digits = 2),
      median = round(median(!!sym(col), na.rm = TRUE), digits = 2),
      Q1 = round(quantile(!!sym(col), probs = 0.25, na.rm = TRUE), digits = 2),
      Q3 = round(quantile(!!sym(col), probs = 0.75, na.rm = TRUE), digits = 2),
      count = n(),
    )

    # Test for normality for each category
    normality_pvalues <- list()
    for (val in unique_vals) {
      cat_vals <- df %>% filter(!!sym(categorical_col) == val) %>% pull(!!sym(col))
      n_val <- length(na.omit(cat_vals))

      if (n_val > 3 && n_val <= 5000) {
        normality_test <- shapiro.test(cat_vals)
        normality_pvalues[[val]] <- normality_test$p.value
      } else {
        normality_pvalues[[val]] <- NA
      }
    }

    # Determine the number of levels in the categorical variable
    num_levels <- length(unique(df[[categorical_col]]))

    # Check if all levels are normally distributed
    if (all(sapply(normality_pvalues, function(p) !is.na(p) && p >= 0.05))) {
      if (num_levels == 2) {
        # Perform two-sample t-test if only 2 levels and normal
        t_test_result <- t.test(df[df[[categorical_col]] == unique_vals[1], ][[col]], df[df[[categorical_col]] == unique_vals[2], ][[col]])
        p_value <- t_test_result$p.value
        test_type <- "Mean (SD)"
      } else {
        # Perform ANOVA if all levels are normally distributed and 3 or more levels
        anova_result <- aov(as.formula(paste(col, "~", categorical_col)), data = df)
        p_value <- summary(anova_result)[[1]][["Pr(>F)"]][1]
        test_type <- "Mean (SD)"
      }
    } else {
      if (num_levels == 2) {
        # Perform Wilcoxon test if only 2 levels and not normal
        wilcox_result <- wilcox.test(df[df[[categorical_col]] == unique_vals[1], ][[col]], df[df[[categorical_col]] == unique_vals[2], ][[col]])
        p_value <- wilcox_result$p.value
        test_type <- "Median (IQR)"
      } else {
        # Perform Kruskal-Wallis test if any level is not normal and 3 or more levels
        kruskal_result <- kruskal.test(as.formula(paste(col, "~", categorical_col)), data = df)
        p_value <- kruskal_result$p.value
        test_type <- "Median (IQR)"
      }
    }

    # Format summary statistics for all categories
    summary <- data.frame(
      Variable = col,
      Levels = test_type
    )

    for (val in unique_vals) {
      cat_stats <- cat_data %>% filter(across(all_of(categorical_col)) == val)
      cat_header <- paste0(val, " (n = ", cat_stats$count, ")")

      summary[[cat_header]] <- if (test_type == "ANOVA (mean (sd))") {
        paste(cat_stats$mean, "(", cat_stats$sd, ")")
      } else {
        paste(cat_stats$median, "(", cat_stats$Q1, ", ", cat_stats$Q3, ")")
      }
    }

    overall_header <- paste0("Overall (n = ", overall_data$count, ")")

    summary[[overall_header]] <- if (test_type == "ANOVA (mean (sd))") {
      paste(overall_data$mean, "(", overall_data$sd, ")")
    } else {
      paste(overall_data$median, "(", overall_data$Q1, ", ", overall_data$Q3, ")")
    }

    summary$p_value <- round(p_value, 3)

    # Check for missing values and add an extra row if needed
    missing_values <- df %>%
      group_by(across(all_of(categorical_col))) %>%
      summarise(missing = sum(is.na(!!sym(col))), .groups = 'drop')

    if (any(missing_values$missing > 0)) {
      missing_summary <- data.frame(
        Variable = col,
        Levels = "Missing Values"
      )

      for (val in unique_vals) {
        missing_summary[[paste0(val, " (n = ", cat_data %>% filter(across(all_of(categorical_col)) == val) %>% pull(count), ")")]] <-
          missing_values %>% filter(across(all_of(categorical_col)) == val) %>% pull(missing)
      }

      missing_summary[[overall_header]] <- sum(missing_values$missing)
      missing_summary$p_value <- NA

      results[[col]] <- rbind(summary, missing_summary)
    } else {
      results[[col]] <- summary
    }
  }

  # Combine all results into one dataframe
  final_results <- bind_rows(results)

  original_vec <- final_results$p_value
  new_values <- as.numeric(original_vec)
  numeric_values <- suppressWarnings(as.numeric(new_values))
  new_values[!is.na(numeric_values) & numeric_values < 0.0001] <- "<1e-04"

  final_results$p_value <- new_values

  final_results <- final_results %>%
    mutate(Variable = ifelse(duplicated(Variable), "", Variable),
           p_value = ifelse(is.na(p_value), "", p_value)) %>%
    rename(`p-value` = p_value,
           `Description` = Levels)

  return(final_results)
}

#dff <- readxl::read_excel("R/descriptives_data.xlsx")
#num_cols = names(dff)[sapply(dff, is.numeric)][-c(6,7,8,20,21)]    ## 6,7,8,20,21
#kk <- analyze_numeric_by_category(df = dff, categorical_col = "Arm_of_randomization", numeric_cols = num_cols)
