#' @title Multivariate Logistic Regression Summary
#'
#' @description Creates a Multivariate Logistic Regression Model summary table;   Preferred: categorical columns are formatted as factors with sub levels.
#'
#' @param data A data set
#' @param response_var (character; binary column's name)
#' @param predictor_vars (character vector; covariate columns' names)
#' @param response_ref (reference level of the binary outcome)
#' @param predictor_refs (named list of references of each categorical covariate)
#'
#' @return A Logistic Model summary with ODDS RATIO (95% CI), p-value & Model Performance. Returns the summary by building the model with outcome variable & all the predictor variables.
#'
#' @export mv_logistic_regression_summary
#' @name mv_logistic_regression_summary




library(broom)
library(dplyr)

mv_logistic_regression_summary <- function(data, response_var, predictor_vars, response_ref, predictor_refs) {

  # Convert response variable to binary with selected reference
  data[[response_var]] <- relevel(as.factor(data[[response_var]]), ref = response_ref)
  event_group <- levels(data[[response_var]])[2] # Identify the non-reference group

  # Identify categorical predictors and convert them to factors
  cat_columns <- predictor_vars[sapply(data[, predictor_vars], function(col) {
    is.character(col) || haven::is.labelled(col)
  })]
  data[cat_columns] <- lapply(data[cat_columns], factor)
  data <- data %>%
    mutate(across(where(is.factor), droplevels))

  # Convert predictor variables with selected reference levels
  for (predictor in predictor_vars) {
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

  # Create formula and fit logistic regression model
  formula <- as.formula(paste(response_var, "~", paste(valid_predictors, collapse = " + ")))
  model <- glm(formula, data = data, family = binomial())

  # Extract model results
  model_results <- tidy(model) %>%
    mutate(
      odds_ratio = exp(estimate),
      conf.low = exp(confint.default(model)[, 1]),
      conf.high = exp(confint.default(model)[, 2]),
      OR_CI = paste0(round(odds_ratio, 4), " (", round(conf.low, 4), ", ", round(conf.high, 4), ")"),
      p.value = round(p.value, 4)
    )

  # Construct final table as in the original function
  final_table <- data.frame(Predictor = character(),
                            Category = character(),
                            Count = numeric(),
                            Event = numeric(),
                            OR_CI = character(),
                            p.value = character(),
                            stringsAsFactors = FALSE)

  for (predictor in predictor_vars) {
    if (is.factor(data[[predictor]])) {
      reference <- levels(data[[predictor]])[1]
      category_counts <- data %>%
        group_by(.data[[predictor]]) %>%
        summarise(
          n = n(),
          event = sum(.data[[response_var]] == event_group)
        )

      final_table <- rbind(final_table, data.frame(
        Predictor = predictor,
        Category = reference,
        Count = category_counts[category_counts[[predictor]] == reference, ]$n,
        Event = category_counts[category_counts[[predictor]] == reference, ]$event,
        OR_CI = paste0("Reference: ", reference),
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
            OR_CI = model_row$OR_CI,
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
        Event = sum(data[[response_var]] == event_group),
        OR_CI = model_row$OR_CI,
        p.value = model_row$p.value
      ))
    }
  }

  final_table <- final_table %>%
    rename(Variable = Predictor, Levels = Category, N = Count, `Event` = Event, `OR (95% CI)` = OR_CI, `p value` = p.value) %>%
    mutate(Variable = ifelse(duplicated(Variable), "", Variable))

  # Add a header row for the response variable and its reference level
  header_row <- setNames(data.frame(
    Variable = paste0("Response Variable: ", response_var),
    Levels = paste0("Reference: ", response_ref),
    N = NA,
    Event = NA,
    `OR (95% CI)` = NA,
    `p value` = NA,
    stringsAsFactors = FALSE
  ), colnames(final_table))

  # Prepend the header row to the final table
  final_table <- rbind(header_row, final_table)

  # Calculate model performance statistics
  model_summary <- summary(model)

  performance_stats <- data.frame(
    Statistic = c("AIC", "BIC", "Deviance", "R-squared"),
    Value = c(
      round(AIC(model), 2),       # AIC value
      round(BIC(model), 2),       # BIC value
      round(model_summary$deviance, 2), # Deviance value
      round(1 - (model_summary$deviance / model_summary$null.deviance), 4) # R-squared
    ),
    stringsAsFactors = FALSE
  )

  # Return both the results table and the model performance table
  return(list(Model_Results = final_table, Model_Performance = performance_stats, model = model))
}
