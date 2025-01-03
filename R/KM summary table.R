#' @title Kaplan Meier Survival Analysis Table (Standard Format)
#'
#' @description Reports a summary table from the KM analysis, with p-value (log rank test), median follow-up (using reverse KM) & median survival (with 95\% CI).
#'
#' @param data A data set
#' @param time_var (column name of time duration) A numeric column containing the duration of survival of all subjects.
#' @param event_var (column name of status/event variable) A binary numeric column containing the indicators whether the concerned event has occurred or not.
#' @param stratify_var (column name of stratification variable) A binary column containing the names of groups in which the study subjects belongs to.
#' @param duration_unit Pass on the choice of duration unit. Choose from "days", "weeks", "months", and "years". Default: "days".
#'
#' @return A standard format dataframe from the Kaplan Meier Survival Analysis. Reports the Median Follow-up of the study, the Variable, its Levels, total subjects in each level, no. of events in each level, Mediabn Survival (95\% CI), and p-value from Log-Rank test.
#'
#' @export km_summary
#' @name km_summary



library(survival)
library(dplyr)

# Kaplan-Meier Summary Table Function
km_summary <- function(data, time_var, event_var, stratify_var, duration_unit = "months") {

  # Strip labels, keep raw values only
  data[] <- lapply(data, as.vector)

  # Duration Conversion Factor
  conversion_factor <- switch(duration_unit,
                              "days" = 1,
                              "weeks" = 7,
                              "months" = 30.4375,
                              "years" = 365.25)

  # Convert durations to specified units
  data$duration <- data[[time_var]] / conversion_factor

  # Median Follow-up Calculation
  rev_km <- survfit(Surv(duration, 1 - data[[event_var]]) ~ 1, data = data)
  median_followup <- round(quantile(rev_km, probs = 0.5)$quantile, 2)

  # Check if stratify_var has only one level
  strat_levels <- levels(factor(data[[stratify_var]]))
  single_group <- length(strat_levels) == 1

  # KM Fit
  km_formula <- as.formula(paste("Surv(duration, ", event_var, ") ~ ", stratify_var))
  fit_km <- survfit(km_formula, data = data)

  # Log-Rank Test
  log_rank_p <- if (!single_group) {
    surv_diff <- survdiff(km_formula, data = data)
    round(surv_diff$pvalue, 4)
  } else {
    NA
  }

  # Prepare Summary Table
  summary_table <- data.frame(
    Stratification = character(),
    Level = character(),
    Total = integer(),
    Events = integer(),
    Median_CI = character(),  # Combined Median and 95% CI
    p_value = character(),    # Allow blank p-values
    stringsAsFactors = FALSE
  )

  for (i in seq_along(strat_levels)) {
    level_data <- subset(data, data[[stratify_var]] == strat_levels[i])

    # Total and Events
    total_obs <- ifelse(is.na(nrow(level_data)), 0, nrow(level_data))
    events <- ifelse(is.na(sum(level_data[[event_var]])), 0, sum(level_data[[event_var]]))

    # Median Survival and CI
    if (single_group) {
      median_surv <- summary(fit_km)$table["median"]
      lower_ci <- summary(fit_km)$table["0.95LCL"]
      upper_ci <- summary(fit_km)$table["0.95UCL"]
    } else {
      median_surv <- summary(fit_km)$table[i, "median"]
      lower_ci <- summary(fit_km)$table[i, "0.95LCL"]
      upper_ci <- summary(fit_km)$table[i, "0.95UCL"]
    }

    # Handle 'NR' for median survival
    median_surv <- ifelse(is.na(median_surv), "NR", round(median_surv, 2))
    lower_ci <- ifelse(is.na(lower_ci), "NR", round(lower_ci, 2))
    upper_ci <- ifelse(is.na(upper_ci), "NR", round(upper_ci, 2))

    # Combine Median and CI
    median_ci <- unname(ifelse(median_surv == "NR", "NR", paste0(median_surv, " [", lower_ci, ", ", upper_ci, "]")))

    # Append Row
    summary_table <- rbind(summary_table, data.frame(
      Stratification = if (single_group) stratify_var else stratify_var,
      Level = if (single_group) "" else strat_levels[i],
      Total = total_obs,
      Events = events,
      Median_CI = median_ci,  # Combined Median and CI
      p_value = if (single_group) "" else if (i == 1) log_rank_p else ""
    ))
  }

  # Add Median Follow-up as the Top Row (no total, events, or p-value)
  summary_table <- rbind(
    data.frame(
      Stratification = "",
      Level = "Median Follow-up",
      Total = paste0(median_followup, " ", duration_unit),
      Events = "",
      Median_CI = "",
      p_value = ""
    ),
    summary_table
  )

  summary_table <- summary_table %>%
    rename(Variable = Stratification,
           `Median Survival (95% CI)` = Median_CI,
           `p-value` = p_value) %>%
    mutate(Variable = ifelse(duplicated(Variable), "", Variable))

  return(summary_table)
}
