#' @title Scatter/Box/Bar Plot
#'
#' @description Plots using data, name of plot, & column names
#'
#' @param df A data set
#' @param type (character; choose from "scatter" or "boxplot" or "barplot")
#' @param var1 (character; column name)
#' @param var2 (character; column name)
#'
#' @return PLOTS
#'
#' @export draw
#' @name Plot



library(ggplot2)
library(dplyr)
library(scales)

draw <- function(df, type, var1, var2 = NULL, x_breaks = 5, y_breaks = 5) {
  if (!type %in% c("scatter", "boxplot", "barplot")) {
    stop("Invalid plot type. Choose from 'scatter', 'boxplot', 'barplot'.")
  }

  if (!var1 %in% colnames(df) || (!is.null(var2) && !var2 %in% colnames(df))) {
    stop("One or both variables are not present in the dataframe.")
  }

  if (type == "scatter") {
    if (is.null(var2)) stop("Scatter plot requires two variables (var1 and var2).")
    ggplot(df, aes_string(x = var1, y = var2)) +
      geom_point(size = 1.8) +
      scale_x_continuous(breaks = round(seq(min(df[[var1]], na.rm = TRUE), max(df[[var1]], na.rm = TRUE), length.out = x_breaks), 2)) +
      scale_y_continuous(breaks = round(seq(min(df[[var2]], na.rm = TRUE), max(df[[var2]], na.rm = TRUE), length.out = y_breaks), 2)) +
      theme_minimal()

  } else if (type == "boxplot") {
    counts <- df %>%
      group_by(!!sym(var1)) %>%
      summarize(count = n())

    labels <- paste0(counts %>% pull(!!sym(var1)), "\n(n = ", counts %>% pull(count), ")")

    ggplot(df, aes_string(x = var1, y = var2)) +
      geom_boxplot(fill = "#bbe7fe", color = "black") +
      scale_x_discrete(labels = labels) +
      theme_minimal()

  } else if (type == "barplot") {
    count_data <- df %>%
      count(!!sym(var1)) %>%
      arrange(desc(n))

    ggplot(count_data, aes_string(x = var1, y = "n", fill = "n")) +
      geom_bar(stat = "identity") +
      scale_fill_gradientn(colors = c("lightblue", "darkblue"),
                           values = scales::rescale(c(0, max(count_data$n))),
                           name = "Count") +
      theme_minimal()
  }
}
