#' @title Univariate Logistic Regression Summary
#'
#' @description Creates multiple Single-Predictor Logistic Regression Model summary table;   Preferred: categorical columns are formatted as factors with sub levels.
#'
#' @param data A data set
#' @param response_var (character; binary column's name)
#' @param predictor_vars (character vector; covariate columns' names)
#' @param response_ref (reference level of the binary outcome)
#' @param predictor_refs (named list of references of each categorical covariate)
#'
#' @return A Logistic Model summary with ODDS RATIO (95\% CI), p-value & Model Performance, along with OR Forest Plot. Returns the summary by building the model with outcome variable & exactly 1 predictor variable, at a time. For example, if 4 predictors are selected, a total of 4 models would be built, with each predictor entering once in a simple logistic model, Y ~ X. The summary henceforth presented would be of the 4 models built.
#'
#' @export univariate_logistic_regression_summary
#' @name univariate_logistic_regression_summary




library(broom)
library(dplyr)

univariate_logistic_regression_summary <- function(data, response_var, predictor_vars, response_ref, predictor_refs) {
  # Convert response variable to binary with selected reference
  data[[response_var]] <- relevel(as.factor(data[[response_var]]), ref = response_ref)
  event_group <- levels(data[[response_var]])[2]

  # Identify categorical predictors and convert them to factors
  cat_columns <- predictor_vars[sapply(data[, predictor_vars], function(col) {
    is.character(col) || haven::is.labelled(col)
    })]
  data[cat_columns] <- lapply(data[cat_columns], function(col) {
    # Replace "NA" (as a string) with actual NA
    col[col == "NA"] <- NA
    col[col == ""] <- NA
    # Convert to factor
    factor(col)
  })

  data <- data %>%
    mutate(across(where(is.factor), droplevels))

  # Prepare a list to store model results
  all_results <- list()

  # Loop through each predictor variable to build individual models
  for (predictor in predictor_vars) {
    # Handle factor predictors: set reference levels if applicable
    if (is.factor(data[[predictor]])) {
      data[[predictor]] <- relevel(data[[predictor]], ref = predictor_refs[[predictor]])
    }

    # Check if the predictor has at least 2 levels (or is valid for regression)
    if (is.factor(data[[predictor]]) && nlevels(data[[predictor]]) < 2) next

    # Create formula and fit logistic regression model
    formula <- as.formula(paste(response_var, "~", predictor))
    model <- glm(formula, data = data, family = binomial())

    # Extract model results
    model_results <- tidy(model) %>%
      mutate(
        odds_ratio = exp(estimate),
        conf.low = exp(confint.default(model)[, 1]),
        conf.high = exp(confint.default(model)[, 2]),
        OR_CI = paste0(round(odds_ratio, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")"),
        p.value = round(p.value, 4)
      )

    # Create a results table for this predictor
    predictor_table <- data.frame(Predictor = character(),
                                  Category = character(),
                                  Count = numeric(),
                                  Event = numeric(),
                                  OR_CI = character(),
                                  p.value = character(),
                                  stringsAsFactors = FALSE)

    if (is.factor(data[[predictor]])) {
      # Handle categorical predictors
      reference <- levels(data[[predictor]])[1]
      category_counts <- data %>%
        group_by(.data[[predictor]]) %>%
        summarise(n = n(), event = sum(.data[[response_var]] == event_group)) %>%
        filter(!is.na(.data[[predictor]]))

      predictor_table <- rbind(predictor_table, data.frame(
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
          predictor_table <- rbind(predictor_table, data.frame(
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
      # Handle continuous predictors
      model_row <- model_results %>% filter(term == predictor)
      predictor_table <- rbind(predictor_table, data.frame(
        Predictor = predictor,
        Category = "Continuous",
        Count = nrow(data),
        Event = sum(data[[response_var]] == event_group),
        OR_CI = model_row$OR_CI,
        p.value = model_row$p.value
      ))
    }

    all_results[[predictor]] <- predictor_table
  }

  # Combine all individual results into a single final table
  final_table <- do.call(rbind, all_results) %>%
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


  #### OR FOREST PLOT ####
  vec <- final_table$`OR (95% CI)`

  OR <- rep(NaN, length(vec))
  OR_lower <- rep(NaN, length(vec))
  OR_upper <- rep(NaN, length(vec))

  # Use regular expressions to extract the numbers only in entries containing numeric values
  numeric_indices <- grepl("\\d", vec)

  OR[numeric_indices] <- as.numeric(gsub("^(\\d+\\.\\d+).*", "\\1", vec[numeric_indices]))
  OR_lower[numeric_indices] <- as.numeric(gsub(".*\\((\\d*\\.?\\d+),.*", "\\1", vec[numeric_indices]))
  OR_upper[numeric_indices] <- as.numeric(gsub(".*\\(\\d*\\.?\\d+,\\s*(\\d*\\.?\\d+).*", "\\1", vec[numeric_indices]))

  mod_base_data <- tibble::tibble(mean  = OR,
                                  lower = OR_lower,
                                  upper = OR_upper,
                                  Variable = final_table$Variable,
                                  Levels = final_table$Levels,
                                  p_values = final_table$`p value`,
                                  OR = final_table$`OR (95% CI)`)

  adj_size = 0
  or_fp <- mod_base_data |>
    forestplot(labeltext = c(Variable, Levels, OR, p_values),
               clip = c(0.1, 3),
               boxsize = 0.18,
               graphwidth = unit(9, "cm"),
               vertices = TRUE,
               xlab = "Odds Ratio",
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
                  OR = c("OR", "(95% CI)"),
                  p_values = c("", "p-value")) |>
    fp_decorate_graph(box = gpar(lty = 2, col = "lightgray"),
                      graph.pos = 3) |>
    fp_add_lines() |>
    fp_set_zebra_style("#EFEFEF") #|> prGridPlotTitle(title = "Forest Plot", gp = gpar(fontsize = 15, fontface = "bold"))

  # Return the results
  return(list('Final Table' = final_table,
              'OR FP' = or_fp))
}
