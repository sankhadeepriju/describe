#' @title Multivariate Cox Proportional Hazard Regression Summary
#'
#' @description Creates a Multivariate Cox Proportional Hazard Regression Model summary table;   Preferred: categorical columns are formatted as factors with sub levels.
#'
#' @param data A data set
#' @param start_date_var The start date of survival. Eg: Date of Registration/Randomization/Treatment start, etc.
#' @param event_date_var The date of occurrence of the Event of Interest Eg: Date of Death/Relapse/Progression/Recurrence, etc.
#' @param end_date_var The end date of survival. Eg: Date of Death/Last Follow Up (LFU), etc.
#' @param event_var binary variable: The variable that tells us if the Event of Interest has occurred or not. Usually, '1' for Event Occurred and '0' for No Event Occurrence.
#' @param event_marker single input: Choose the correct option among the 2 sub-categories of event_var. The event_var is likely to have binary options like "Dead" or "Alive", '1' or '0', "Progressed" or "Not Progressed", etc.
#' @param predictor_vars (character vector; covariate columns' names)
#' @param predictor_refs (named list of reference levels of each categorical covariate)
#'
#' @return A Logistic Model summary with ODDS RATIO (95\% CI), p-value & Model Performance, along with OR Forest Plot. Returns the summary by building the model with the outcome variable & all the predictor variables.
#'
#' @export mv_cox_ph_model_summary
#' @name mv_cox_ph_model_summary

mv_cox_ph_model_summary <- function(data, start_date_var, event_date_var, end_date_var, event_var, event_marker, predictor_vars, predictor_refs) {

  # Check if event_var is numeric and consists of 0 and 1 only
  if (is.numeric(data[[event_var]]) && all(data[[event_var]] %in% c(0, 1))) {
    data$new_event_column <- data[[event_var]]
    message("The event variable is already numeric and contains only 0 and 1.")
  } else {
    # If it doesn't meet the criteria, create a new column
    data$new_event_column <- ifelse(data[[event_var]] == event_marker, 1, 0)
    message("New event column created with 1 for subjects matching the Event marker.")
  }

  # Identify categorical predictors and convert them to factors
  cat_columns <- predictor_vars[sapply(data[, predictor_vars], function(col) {
    is.character(col) || haven::is.labelled(col)
  })]
  data[cat_columns] <- lapply(data[cat_columns], function(col) {
    if (is.numeric(col) && haven::is.labelled(col)) {
      # If the column is numeric but labelled, directly convert to factor
      factor(as.character(col))
    } else {
      # For character columns or others, replace "NA" and empty strings, then convert to factor
      col[col == "NA"] <- NA
      col[col == ""] <- NA
      factor(as.character(col))
    }
  })

  data <- data %>%
    mutate(across(where(is.factor), droplevels))

  # Convert predictor variables with selected reference levels
  for (predictor in cat_columns) {
    if (is.factor(data[[predictor]])) {
      data[[predictor]] <- relevel(data[[predictor]], ref = predictor_refs[[predictor]])
    }
  }

  # Check for factors with fewer than 2 levels and remove them
  valid_predictors <- predictor_vars[sapply(data[, predictor_vars], function(x) {
    if (is.factor(x)) {
      nlevels(x) > 1  # Keep only factors with more than one level
    } else {
      TRUE  # Keep non-factor variables
    }
  })]

  # Create a new time duration column
  data$survival_time <- ifelse(data$new_event_column == 1,
                               as.numeric(difftime(data[[event_date_var]], data[[start_date_var]], units = "days"))/30.4375,  # From start_date to event_date for event = 1
                               as.numeric(difftime(data[[end_date_var]], data[[start_date_var]], units = "days"))/30.4375)  # From start_date to end_date for event = 0
  survival_obj <- Surv(time = data$survival_time, event = data$new_event_column)

  # Create the Cox Proportional Hazards model formula
  formula <- as.formula(paste("survival_obj ~", paste(valid_predictors, collapse = " + ")))

  # Fit the Cox Proportional Hazards model
  model <- coxph(formula, data = data)

  # Extract model results using broom
  model_results <- tidy(model) %>%
    mutate(
      exp.coef = exp(estimate),
      conf.low = exp(confint(model)[, 1]),
      conf.high = exp(confint(model)[, 2]),
      HR_CI = paste0(round(exp.coef, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")"),
      p.value = round(p.value, 4)
    )

  # Construct final table similar to the logistic regression function
  final_table <- data.frame(Predictor = character(),
                            Category = character(),
                            Count = numeric(),
                            Event = numeric(),
                            HR_CI = character(),
                            p.value = character(),
                            stringsAsFactors = FALSE)

  for (predictor in valid_predictors) {
    if (is.factor(data[[predictor]])) {
      reference <- levels(data[[predictor]])[1]
      category_counts <- data %>%
        group_by(.data[[predictor]]) %>%
        summarise(
          n = n(),
          event = sum(.data[[event_var]] == 1)
        ) %>%
        filter(!is.na(.data[[predictor]]))

      final_table <- rbind(final_table, data.frame(
        Predictor = predictor,
        Category = reference,
        Count = category_counts[category_counts[[predictor]] == reference, ]$n,
        Event = category_counts[category_counts[[predictor]] == reference, ]$event,
        HR_CI = paste0("Reference: ", reference),
        p.value = ""
      ))

      for (category in levels(data[[predictor]])[-1]) {
        model_row <- model_results %>% filter(grepl(paste0(predictor, category), term))
        if (nrow(model_row) > 0) {
          final_table <- rbind(final_table, data.frame(
            Predictor = predictor,
            Category = category,
            Count = category_counts$n[category_counts[[predictor]] == category],
            Event = category_counts$event[category_counts[[predictor]] == category],
            HR_CI = model_row$HR_CI,
            p.value = model_row$p.value
          ))
        }
      }
    } else {
      model_row <- model_results %>% filter(term == predictor)
      final_table <- rbind(final_table, data.frame(
        Predictor = predictor,
        Category = "Continuous",
        Count = nrow(data),
        Event = sum(data[[event_var]] == 1),
        HR_CI = model_row$HR_CI,
        p.value = model_row$p.value
      ))
    }
  }

  final_table <- final_table %>%
    rename(Variable = Predictor, Levels = Category, N = Count, `Event` = Event, `HR (95% CI)` = HR_CI, `p value` = p.value) %>%
    mutate(Variable = ifelse(duplicated(Variable), "", Variable))

  # Add a header row for the event variable and its reference level
  header_row <- setNames(data.frame(
    Variable = paste0("Event Variable: ", event_var),
    Levels = paste0("Reference: ", "0"),
    N = NA,
    Event = NA,
    `HR (95% CI)` = NA,
    `p value` = NA,
    stringsAsFactors = FALSE
  ), colnames(final_table))

  # Prepend the header row to the final table
  final_table <- rbind(header_row, final_table)

  # Calculate model performance statistics
  model_summary <- summary(model)

  performance_stats <- data.frame(
    Statistic = c("Log-Likelihood", "AIC", "BIC", "R-squared"),
    Value = c(
      round(model_summary$loglik[2], 2),       # Log-likelihood
      round(AIC(model), 2),                    # AIC value
      round(BIC(model), 2),                    # BIC value
      round(1 - (model_summary$loglik[2] / model_summary$loglik[1]), 4) # Pseudo R-squared
    ),
    stringsAsFactors = FALSE
  )

  #### HR FOREST PLOT ####
  vec <- final_table$`HR (95% CI)`

  HR <- rep(NaN, length(vec))
  HR_lower <- rep(NaN, length(vec))
  HR_upper <- rep(NaN, length(vec))

  # Extract numeric values from HR CI
  numeric_indices <- grepl("\\d", vec)

  HR[numeric_indices] <- as.numeric(gsub("^(\\d+\\.\\d+).*", "\\1", vec[numeric_indices]))
  HR_lower[numeric_indices] <- as.numeric(gsub(".*\\((\\d*\\.?\\d+),.*", "\\1", vec[numeric_indices]))
  HR_upper[numeric_indices] <- as.numeric(gsub(".*\\(\\d*\\.?\\d+,\\s*(\\d*\\.?\\d+).*", "\\1", vec[numeric_indices]))

  mod_base_data <- tibble::tibble(mean  = HR,
                                  lower = HR_lower,
                                  upper = HR_upper,
                                  Variable = final_table$Variable,
                                  Levels = final_table$Levels,
                                  p_values = final_table$`p value`,
                                  HR = final_table$`HR (95% CI)`)

  adj_size = 0
  hr_fp <- mod_base_data |>
    forestplot(labeltext = c(Variable, Levels, HR, p_values),
               clip = c(0.1, 3),
               boxsize = 0.18,
               graphwidth = unit(9, "cm"),
               vertices = TRUE,
               xlab = "Hazard Ratio",
               xlog = TRUE,
               xticks = c(0.1, 0.2, 0.5, 1, 1.5, 2, 3),
               txt_gp = fpTxtGp(label = gpar(cex = 1.1 - adj_size),   # Data label size
                                ticks = gpar(cex = 1.1 - adj_size),   # Tick size
                                xlab = gpar(cex = 1.3 - adj_size)),
               title = "Forest Plot") |>
    fp_set_style(box = "royalblue",
                 line = "darkblue",
                 summary = "royalblue") |>
    fp_add_header(Variable = c("", "Variable"),
                  Levels = c("", "Levels"),
                  HR = c("HR", "(95% CI)"),
                  p_values = c("", "p-value")) |>
    fp_decorate_graph(box = gpar(lty = 2, col = "lightgray"),
                      graph.pos = 3) |>
    fp_add_lines() |>
    fp_set_zebra_style("#EFEFEF")

  # Return both the results table and the model performance table
  return(list(Model_Results = final_table, Model_Performance = performance_stats,
              model = model, 'HR FP' = hr_fp))
}
