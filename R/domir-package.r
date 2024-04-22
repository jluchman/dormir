#' @name domir-package
#' @title Tools to Support Relative Importance Analysis
#' @author Joseph Luchman \email{jluchman@gmail.com}
#' @keywords internal
#' "_PACKAGE"
#'
#' @description
#' Methods to apply dominance analysis-based relative importance analysis for
#' predictive modeling functions.
#'
#' @details
#' This package supports the determination of importance for inputs
#' (i.e., independent variables, predictors, features, parameter estimates;
#' called 'names' in the package) using dominance analysis
#' (Azen & Budescu, 2004; Budescu, 1993).
#'
#' Dominance analysis resolves the indeterminancy of ascribing
#' the value returned by a predictive modeling function to inputs/names when
#' it is not possible to do so analytically. The most common use case for the
#' application of dominance analysis is in comparing inputs/names in terms of
#' their contribution to a predictive model's fit statistic or metric.
#'
#' Dominance analysis is a common, and generally well accepted, method for
#' determining the relative importance of inputs/names that is, in part,
#' a conceptual extension of the well-known Shapley value
#' decomposition (e.g., Grömping, 2007; Lipovetsky & Conklin, 2001).
#'
#' @references
#' \itemize{
#' \item Azen, R., & Budescu, D. V. (2003). The dominance analysis approach
#' for comparing predictors in multiple regression. Psychological Methods,
#' 8(2), 129-148. doi:10.1037/1082-989X.8.2.129
#' \item Budescu, D. V. (1993). Dominance analysis: A new approach to the
#' problem of relative importance of predictors in multiple regression.
#' Psychological Bulletin, 114(3), 542-551. doi:10.1037/0033-2909.114.3.542
#' \item Grömping, U. (2007). Estimators of relative importance in linear
#' regression based on variance decomposition. The American Statistician,
#' 61(2), 139-147. doi:10.1198/000313007X188252
#' \item Lipovetsky, S, & and Conklin, M. (2001). Analysis of regression in
#' game theory approach. Applied Stochastic Models in Business and Industry,
#' 17(4), 319-330. doi:10.1002/asmb.446
#'}

NULL
