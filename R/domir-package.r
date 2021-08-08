#' The \code{domir} package is intended to provide a set of flexible wrapper and helper functions for conducting relative importance analysis with a focus on dominance analysis. That is, the intention of this package is to provide tools that allow relative importance analysis across a wide variety of practical data analytic situations.
#' 
#' @details
#' Relative importance analysis is a methodology focused on comparing independent variables (IVs)/features/predictors as well as parameter estimates to one another in terms of how they predict some dependent variable/response/outcome in the context of a predictive model. 
#' 
#' The intention of this package is to focus on what I will call "model evaluation" or post hoc examination of IV predictive utility in the context of a vetted, selected predictive model.  That is, the methods to apply here will assume that the user has previously applied model selection methods and that the IVs of the predictive model for relative importance analysis have non-trivial effects in improving model fit. The methods in this package are not intended for model selection - though I will acknowledge that many importance methods are (at least implicitly) focused on model selection/identifying which IVs have a trivial effect and removing them.
#' 
#' The only method implemented at current in \code{domir} is dominance analysis (DA) method \code{domin}. \code{domin} is a flexible wrapper function that can be used with many modeling functions. \code{domin} is an extension of the Stata command by the same name (see Luchman, 2021) to the R environment.
#' 
#' @section Dominance Analysis:
#' 
#' As a relative importance method, DA determines the relative importance of independent variables or parameter estimates (PEs) in an prediction model based on contribution to an overall model fit statistic/metric (see Budescu, 1993; Groemping, 2007 for a discussions). DA is an ensemble method in which importance determinations about IV/PEs are made by producing weighted averages of fit statistic results across many models, though implementations of the method usually require the ensemble contain fit statistics for each possible combination of models where the IV/PEs are either included or excluded and parameter estimates/prediction rules for these IV/PEs are estimated from the data for each model.
#' 
#' The implementation of DA in the \code{domin} function computes dominance statistics by moving through three processes. 
#' 
#' First, all IV/PEs (inlcuding sets of IV/PEs) for the DA are identified and all possible combinations of these entries are collected. The all possible combinations ensemble with \eqn{p} IV/PEs in focal model submitted to \code{domin} results in \eqn{2^p} combinations of IV/PEs. That is, each combination of \eqn{p} IV/PEs alternating between included versus excluded (i.e., the base of 2 to the \eqn{p} exponent number of IV/PEs).  This process can take 
#' 
#' Second, all the models implied by the IV/PE combinations collected in the first step are estimated, a fit statistic extractor function is applied, and the fit statistics for each model are recorded.  
#' 
#' Third, for each model in the DA, all marginal/incremental fit statistic values are computed.  This requires identifying the submodel that contains all the IV/PEs in the current model with the exception of the focal IV/PE and computing and recording the difference between these two models.  From these marginal fit statistics, conditional, complete, and general dominance statistics are computed.
#' 
#' @section Dominance Statistics:
#' 
#' There are three sets of statistics produced by DA.
#' 
#' \strong{General Dominance Statistics} are an additive decomposition of the fit statistic associated with the full, selected model that are ascribed to each IV/PE.  These represent the average of the conditional dominance statistics described below.  An IV \emph{generally dominates} another IV when its general dominance statistic is bigger.
#' 
#' \strong{Conditional Dominance Statistics} are the average marginal contribution to the fit statistic for an IV/PE for all models including a specific number of IV/PEs.  The number of conditional dominance statistics is equal to \eqn{p}. An IV \emph{conditionally dominates} another IV when its conditional dominance statistics are bigger by number of IV/PEs in the model (i.e., by column in the conditional dominance matrix).
#' 
#' \strong{Complete Dominance Designations} are logicals indicating whether, for each comparable submodel examining the marginal contribution of two IV/PEs, one IV/PE is \emph{always} larger than the other.  The rows of the matrix indicate an IV/PE dominating the IV/PE in the columns.  \code{TRUE} values indicate such a dominance relationship.  \code{FALSE} indicates the opposite, that the IV/PE in the row is dominated by the IV/PE in the column and \code{NA} indicates no such relationship can be determined.
#' 
#' Complete dominance is the strongest form of dominance, followed by conditional, followed by general.  
#'
#' @name domir-package
#' @aliases domir
#' @docType package
#' @title Tools to Support Relative Importance Analysis
#' @author Joseph Luchman \email{jluchman_at_gmail_com}
#' @references
#' \itemize{
#' \item Azen, R., & Budescu, D. V. (2003). The dominance analysis approach for comparing predictors in multiple regression. Psychological Methods, 8(2), 129-148. doi:10.1037/1082-989X.8.2.129
#' \item Budescu, D. V. (1993). Dominance analysis: A new approach to the problem of relative importance of predictors in multiple regression. Psychological Bulletin, 114(3), 542-551. doi:10.1037/0033-2909.114.3.542
#' \item Groemping, U. (2007). Estimators of relative importance in linear regression based on variance decomposition. The American Statistician, 61(2), 139-147. doi:10.1198/000313007X188252
#' \item Luchman, J. N., Lei, X., & Kaplan, S. A. (2020). Relative Importance Analysis With Multivariate Models: Shifting the Focus from Independent Variables to Parameter Estimates. Journal of Applied Structural Equation Modeling, 4(2), 1-20. doi:10.47263/JASEM.4(2)02
#'  \item Luchman, J. N. (2021). Determining relative importance in Stata using dominance analysis: domin and domme. Stata Journal 21(2), 510-538. doi:10.1177/1536867X211025837
#'}
NULL