#' @title Logistic Regression Summary
#'
#' @description Creates a Logistic Regression Model; Preferred: categorical columns are formatter as factors with sub levels.
#'
#' @param data A data set
#' @param response_var (character; binary column's name)
#' @param predictor_vars (character vector; covariate columns' names)
#' @param response_ref (reference level of the binary outcome)
#' @param predictor_refs (named list of references of each categorical covariate)
#'
#' @return A Logistic Model summary with ODDS RATIO (95% CI), p-value & Model Performance
#'
#' @export logistic_regression_summary
#' @name logistic_regression_summary




library(broom)
library(dplyr)

logistic_regression_summary <- function(data, response_var, predictor_vars, response_ref, predictor_refs) {

  # Convert response variable to binary with selected reference
  data[[response_var]] <- relevel(as.factor(data[[response_var]]), ref = response_ref)

  # Identify categorical predictors and convert them to factors
  cat_columns <- predictor_vars[sapply(data[, predictor_vars], is.character)]
  data[cat_columns] <- lapply(data[cat_columns], factor)
  data <- data %>%
    mutate(across(where(is.factor), droplevels))

  # Convert predictor variables with selected reference levels
  for (predictor in predictor_vars) {
    if (is.factor(data[[predictor]])) {
      data[[predictor]] <- relevel(data[[predictor]], ref = predictor_refs[[predictor]])
    }
  }

  #Check for factors with fewer than 2 levels and remove them
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
                            OR_CI = character(),
                            p.value = character(),
                            stringsAsFactors = FALSE)

  for (predictor in predictor_vars) {
    if (is.factor(data[[predictor]])) {
      reference <- levels(data[[predictor]])[1]
      category_counts <- data %>%
        group_by(.data[[predictor]]) %>%
        summarise(n = n())

      final_table <- rbind(final_table, data.frame(
        Predictor = predictor,
        Category = reference,
        Count = category_counts[category_counts[[predictor]] == reference,]$n,
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
        OR_CI = model_row$OR_CI,
        p.value = model_row$p.value
      ))
    }
  }

  final_table <- final_table %>%
    rename(Variable = Predictor, Levels = Category, N = Count, `OR (95% CI)` = OR_CI, `p value` = p.value) %>%
    mutate(Variable = ifelse(duplicated(Variable), "", Variable))



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
