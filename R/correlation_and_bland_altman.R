#' @title Correlation & Bland Altman Plot
#'
#' @description Calculates Correlation for 2 variables (eg: 2 methods of treatment) & produces a Bland Altman Plot
#'
#' @param data A data set
#' @param v1 (character; numeric column's name) To calculate correlation & conduct appropriate test for the same.
#' @param v2 (character; numeric column's name) To calculate correlation & conduct appropriate test for the same.
#' @param var1 (numeric vector) To plot the Bland Altman plot by calculating the means and differences of the 2 variables chosen.
#' @param var2 (numeric vector) To plot the Bland Altman plot by calculating the means and differences of the 2 variables chosen.
#'
#' @return Correlation Summary with 95\% CI (use function correlation_test) & Bland Altman Plot (use function blandaltmanplot)
#'
#' @export correlation_test
#' @export blandaltmanplot
#' @name correlation_test_and_blandaltmanplot



library(DescTools)
library(ggplot2); library(BlandAltmanLeh); library(blandr); library(dplyr); library(readxl)
library(plotly)

#### Correlation ####
correlation_test <- function(v1, v2, data) {

  var1 = data[[v1]]
  var2 = data[[v2]]

  # Check if the variables are numeric
  if (!is.numeric(var1) || !is.numeric(var2)) {
    stop("Both variables must be numeric.")
  }

  # Helper function to calculate correlation and CIs
  calculate_correlation <- function(var1, var2, method) {
    if (method == "pearson") {
      # Use cor.test for Pearson correlation
      cor_test <- cor.test(var1, var2, method = "pearson")
      cor_coef <- cor_test$estimate
      ci_lower <- cor_test$conf.int[1]
      ci_upper <- cor_test$conf.int[2]
      p_value <- round(cor_test$p.value, 4)
      p_value <- ifelse(p_value < 1e-04, "<1e-04", p_value)
    } else if (method == "spearman") {
      # Use SpearmanRho from DescTools for Spearman's rho
      cor_test <- SpearmanRho(var1, var2, conf.level = 0.95)
      cor_coef <- cor_test[1]
      ci_lower <- cor_test[2]
      ci_upper <- cor_test[3]
      # Calculate p-value using cor.test for Spearman's rho
      cor_test_full <- cor.test(var1, var2, method = "spearman")
      p_value <- round(cor_test_full$p.value, 4)
      p_value <- ifelse(p_value < 1e-04, "<1e-04", p_value)
    } else if (method == "kendall") {
      # Use KendallTau from DescTools for Kendall's tau
      cor_test <- KendallTauB(var1, var2, conf.level = 0.95)
      cor_coef <- cor_test[1]
      ci_lower <- cor_test[2]
      ci_upper <- cor_test[3]
      # Calculate p-value using cor.test for Kendall's tau
      cor_test_full <- cor.test(var1, var2, method = "kendall")
      p_value <- round(cor_test_full$p.value, 4)
      p_value <- ifelse(p_value < 1e-04, "<1e-04", p_value)
    }
    return(list(cor = cor_coef, p_value = p_value,
                ci_lower = ci_lower, ci_upper = ci_upper))
  }

  # Check normality using Shapiro-Wilk test
  normal1 <- shapiro.test(var1)$p.value > 0.05
  normal2 <- shapiro.test(var2)$p.value > 0.05

  # Count number of ties in each variable
  ties_var1 <- length(var1) - length(unique(var1))
  ties_var2 <- length(var2) - length(unique(var2))

  # Choose correlation method based on normality and ties
  if (normal1 & normal2) {
    # Both variables are normal, use Pearson
    method <- "pearson"
  } else if (ties_var1 > 0.2 * length(var1) | ties_var2 > 0.2 * length(var2)) {
    # If there are a lot of ties (more than 20% of the data), use Kendall
    method <- "kendall"
  } else {
    # Otherwise, use Spearman
    method <- "spearman"
  }

  # Perform the correlation test
  cor_result <- calculate_correlation(var1, var2, method)

  # Add results to the table
  result_table <- data.frame(Variable1 = v1,
                             Variable2 = v2,
                             "Correlation Coefficient" = cor_result$cor,
                             "CI Lower" = cor_result$ci_lower,
                             "CI Upper" = cor_result$ci_upper,
                             "p value" = cor_result$p_value,
                             "Test Used" = method,
                             check.names = FALSE)
  # Return the result table
  return(result_table)
}


#### Bland Altman Plot ####
blandaltmanplot = function(var1, var2){

  # Ensure both variables are numeric
  if (!is.numeric(var1) || !is.numeric(var2)) {
    return(NULL)
  }
  # Compute differences between the methods
  differences <- var1 - var2
  means <- (var1 + var2) / 2
  mean_data = data.frame(means, differences)

  # Check normality of differences using Shapiro-Wilk test
  shapiro_test <- shapiro.test(differences)
  normality_result <- shapiro_test$p.value

  # If normal, use blandr.draw, else calculate median and percentiles
  if (normality_result > 0.05) {
    # Normal case: Use blandr.draw for a crisp plot
    bland_altman = blandr.draw(var1, var2, plotTitle = "Bland-Altman Plot",
                               ciDisplay = FALSE, ciShading = FALSE)

    # Calculate the mean difference and limits of agreement
    mean_diff <- mean(differences)
    upper_loa <- mean_diff + 1.96 * sd(differences)
    lower_loa <- mean_diff - 1.96 * sd(differences)

    # Add annotations to the plot
    bland_altman <- bland_altman +
      theme(
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),  # Title font size
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)
      ) +
      annotate("text", x = max(means), y = mean_diff,
               label = paste("Mean =", round(mean_diff, 2)),
               vjust = -0.8, hjust = 0.5, color = "blue", size = 5) +
      annotate("text", x = max(means), y = upper_loa,
               label = paste("Upper LoA =", round(upper_loa, 2)),
               vjust = -0.35, hjust = 0.7, color = "blue", size = 5) +
      annotate("text", x = max(means), y = lower_loa,
               label = paste("Lower LoA =", round(lower_loa, 2)),
               vjust = 1.2, hjust = 0.7, color = "blue", size = 5)

    # Plot the final annotated Bland-Altman plot
    print(bland_altman)

  } else {

    # Non-normal case: Use median and percentiles
    median_diff <- median(differences)
    lower_limit <- quantile(differences, 0.025)
    upper_limit <- quantile(differences, 0.975)

    # Plot Bland-Altman using median and percentiles
    bland_altman <- ggplot(mean_data,aes(x = means, y = differences)) +
      geom_point() +
      geom_hline(yintercept = median_diff, color = "black", linetype = "dashed", size = 0.55) +
      geom_hline(yintercept = lower_limit, color = "black", linetype = "dashed", size = 0.55) +
      geom_hline(yintercept = upper_limit, color = "black", linetype = "dashed", size = 0.55) +
      xlab("Means") +
      ylab("Differences") +
      ggtitle("Bland - Altman Plot") +
      theme(
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),  # Title font size
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)
        ) +
      annotate("text", x = max(means), y = mean_diff,
               label = paste("Median =", round(median_diff, 2)),
               vjust = -0.8, hjust = 0.5, color = "blue", size = 5) +
      annotate("text", x = max(means), y = upper_loa,
               label = paste("97.5th LoA =", round(upper_limit, 2)),
               vjust = -0.35, hjust = 0.7, color = "blue", size = 5) +
      annotate("text", x = max(means), y = lower_loa,
               label = paste("2.5th LoA =", round(lower_loa, 2)),
               vjust = 1.2, hjust = 0.7, color = "blue", size = 5)

    bland_altman
  }
}
