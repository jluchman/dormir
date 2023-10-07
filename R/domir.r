#' @title Dominance analysis methods
#' @name domir
#' @description
#' Parses input object to obtain list of names, determines all required
#' combinations of subsets of the name list, submits name list subsets to
#' a function as the input type, and computes dominance decomposition
#' statistics based on the returned values from the function.
#'
#' @param .obj A `formula` or `formula_list`.
#'
#' Parsed to produce list of names. Combinations of subsets the name list are
#' [`sapply`]-ed to `.fct`.
#'
#' The name list subsets submitted to `.fct` are formatted to be of
#' the same [`class`] as `.obj` and are submitted to
#' `.fct` as the first, unnamed argument.
#'
#' @param .fct A [`function`] or string function name.
#'
#' Applied to all subsets of elements as received from `.obj`.
#' Must return a length 1/scalar, numeric, atomic vector.
#'
#' @param .set A `list`.
#'
#' Must be comprised of elements of the same class as `.obj`.
#' Elements of the list can be named.
#'
#' @param .wst Not yet used.
#'
#' @param .all A `formula` or `formula_list`.
#'
#' Must be the same class as `.obj`.
#'
#' @param .adj Logical.
#'
#' If `TRUE` then a model including only an intercept is submitted to `.fct`
#' and the value returned is subtracted from the values returned from all
#' subsets in the dominance analysis.
#'
#' @param .cdl Logical.
#'
#' If `FALSE` then conditional dominance matrix is not computed and
#' method to produce general dominance statistics changes.
#'
#' @param .cpt Logical.
#'
#' If `FALSE` then complete dominance matrix is not computed.
#'
#' @param .rev Logical.
#'
#' If `TRUE` then standardized vector, ranks, and complete dominance
#' designations are reversed in their interpretation.
#'
#' @param ... Passes arguments to other methods during method dispatch;
#' passes arguments to the function in `.fct` during function execution.
#'
#' @return Returns an object of [`class`] "domir" composed of:
#' \describe{
#'  \item{`General_Dominance`}{Vector of general dominance values.}
#'  \item{`Standardized`}{Vector of general dominance values normalized
#'  to sum to 1.}
#'  \item{`Ranks`}{Vector of ranks applied to the general dominance values.}
#'  \item{`Conditional_Dominance`}{Matrix of conditional dominance values.
#'  Each row represents an element in `.obj`;
#'  each column represents a number of elements from `.obj` in a subset.}
#'  \item{`Complete_Dominance`}{Logical matrix of complete dominance
#'  designations.
#'  The `.obj` elements represented in each row indicates dominance status;
#'  the `.obj` elements represented in each column indicates
#'  dominated-by status.}
#'  \item{`Value`}{Value returned by `.fct` with all elements (i.e.,
#'  from `.obj`, `.all`, and `.adj`.}
#'  \item{`Value_All`}{Value of `.fct` associated with elements included
#'  in `.all`;
#'  when elements are in `.adj`, will be adjusted for `Value_Adjust`.}
#'  \item{`Value_Adjust`}{Value of `.fct` associated with elements in `.adj`.}
#'  \item{`Call`}{The matched call.}
#' }
#'
#' @details
#' ## Element Parsing
#'
#' `.obj`s is parsed into a name list that is used to determine
#' the required number of combinations of subsets of the name list
#' included the dominance analysis.  How the name list is obtained
#' depends on `.obj`'s class.
#'
#' ### `formula`
#'
#' The `formula` creates a name list using all terms in the formula.
#' The terms are obtained using [`terms.formula`]. All processing
#' that is normally applied to the right hand side of a formula is
#' implemented (see [`formula`]).
#'
#' A response/left hand side is not required but, if present, is
#' included in all `formula`s passed to `.fct`.
#'
#' ### `formula_list`
#'
#' The [`formula_list`] creates a name list out of response-term pairs.
#' The terms are obtained using `terms.formula` applied to each individual
#' formula in the list.
#'
#' ### Additional Details
#'
#' By default, names obtained from `.obj` are all considered separate
#' 'value-generating names' with the same priority.
#' Each value-generating name will be a separate element when
#' computing combination subsets and will be compared to all other
#' value-generating names.
#'
#' `formula`s and `formula_list` elements are assumed to have an intercept
#' except if explicitly removed with a `- 1` in the `formula`(s) in `.obj`.
#' If removed, the intercept will be removed in all `formula`(s) in each
#' `sapply`-ed subset to `.fct`.
#'
#' If [`offset`]s are included, they are passed, like intercepts, while
#' `sapply`-ing subsets to `.fct`. Currently, only the `formula_list` method
#' allows `offsets`.
#'
#' ## Changing Element Parsing
#'
#' All methods' default behavior that considers all value-generating names
#' to be of equal priority can be overriden using `.set` and `.all` arguments.
#'
#' Names in `.set` and `.all` must also be present in `.obj`.
#'
#' ### `.set`
#'
#' `.set` binds together value-generating names such that
#' they are of equal priority and are never separated when submitted to
#' `.fct`.
#' Thus, the elements in `.obj` bound together contribute jointly to the
#' returned value and are considered, effectively, a single
#' value-generating name.
#'
#' If list elements in `.set` are named, this name will be used in all
#' returned results as the name of the set of value-generating names bound
#' together.
#'
#' `.set` thus considers the value-generating names an 'inseparable set' in the
#' dominance analysis and are always included or excluded together.
#'
#' ### `.all`
#'
#' `.all` gives immediate priority to value-generating names.
#' The value-generating names in `.all` are bound together, are
#' ascribed their full amount of the returned value from `.fct`, and
#' are not adjusted for contribution of other value-generating names.
#'
#' The value of `.fct` ascribed to the value-generating names bound
#' together in `.all` is returned separately from, and not directly
#' compared to, the other value-generating names.
#'
#' The `formula` method for `.all` does not allow a left hand side.
#'
#' `.all` includes the value-generating names in 'all subsets' submitted to
#' the dominance analysis which effectively removes the value associated with
#' this set of names.
#'
#' ### `.adj`
#'
#' `.adj` indicates that an intercept-only model should be supplied to `.fct`.
#' This intercept-only subset is given most immediate priority and the
#' value of `.fct` ascribed to it is removed from all other
#' value-generating names and groups including those in `.all`.
#'
#' The `formula` method will submit an intercept-only formula to `.fct`.
#' The `formula_list` method creates a separate, intercept-only subset for each
#' of the `formula`s in the list.
#' Both the `formula` and `formula_list` methods will respect the user's
#' removal of an intercept. The `formula_list` method will also respect the
#' user's inclusion of an `offset` and will include them in the submission to
#' `.fct`.
#'
#' `.adj` then 'adjusts' the returned value for a non-0 value-returning
#' null model when no value generating names are included.
#'
#' ### Additional Details
#'
#' All methods submit combinations of subsets of names as an
#' object of the same class as `.obj`.
#' A `formula` in `.obj` will submit all combinations of subsets of names
#' as `formula`s to `.fct`.
#' A `formula_list` in `.obj` will submit all combinations of subsets of names
#' as `formula_list`s to `.fct`.
#' In the case that `.fct` requires a different `class` (i.e.,
#' a vector of names, a [`Formula::Formula`] see [`fmllst2Fml`]) the
#' subsets of names will have to be processed in `.fct` to
#' obtain the correct `class`.
#'
#' The all subsets of names will be submitted to `.fct` as the first, unnamed
#' argument.
#'
#' ## `.fct` as Analysis Pipeline
#'
#' The function `sapply`-ed and to which the combinations of subsets of
#' names will be applied.
#'
#' `.fct` is expected to be a complete analysis pipeline that receives a
#' subset of names of the same `class` as `.obj`, uses the names in the
#' `class` as submitted to generate a returned value of the appropriate
#' type to dominance analyze. Typically, this returned value is a
#' fit statistic extracted from a predictive model.
#'
#' At current, only atomic (i.e., non-`list`), numeric scalars (i.e.,
#' vectors of length 1) are allowed as returned values.
#'
#' The `.fct` argument is strict about names submitted and returned value
#' requirements for functions used and applies a series of checks to
#' ensure the submitted names and returned value adhere to these requirements.
#' The checks include whether the `.obj` can be submitted to `.fct` without
#' producing an error and whether the
#' returned value from `.fct` is a length 1, atomic, numeric vector.
#' In most circumstances, the user will have to make their own named or
#' anonymous function to supply as `.fct` to satisfy the checks.
#'
#' # Notes
#'
#' ## `formula` method
#'
#' Prior to version 1.1.0, the `formula` method allowed a `formula`
#' to be submitted to `.adj`.
#' Submitting an intercept-only `formula` as opposed to a
#' logical has been depreciated and submitting a `formula` with more than an
#' intercept is defunct.
#'
#' The `formula` and `formula_list` methods can be used to pass responses,
#' intercepts, and in some cases, `offset`s to all combinations of subsets
#' of names.
#' If the user seeks to include other model components integral to
#' estimation
#' (i.e., a random effect term in [`lme4::glmer()`]) include them as
#' [`update`][update.formula] to the submitted `formula` or `formula_list`
#' imbedded in `.fct`.
#'
#' Second-order or higher terms (i.e., interactions like`~ a*b`) are parsed
#' by default but not used differently from first-order terms for producing
#' subsets. The values ascribed to such terms may not be valid unless
#' the user ensures that second-order and
#' higher terms are used appropriately in `.fct`.
#'
#' @export
#' @examples
#' ## Linear model returning r-square
#' lm_r2 <-
#'   function(fml, data) {
#'     lm_res <- lm(fml, data = data)
#'     summary(lm_res)[["r.squared"]]
#'  }
#'
#' domir(mpg ~ am + vs + cyl, lm_r2, data = mtcars)
#'
#'
#' ## Linear model including set
#' domir(
#'   mpg ~ am + vs + cyl + carb + gear + disp + wt,
#'   lm_r2,
#'   .set = list(~ carb + gear, ~ disp + wt),
#'   data = mtcars
#' )
#'
#'
#' ## Multivariate regression with multivariate r-square and
#' ## all subsets variable
#' mlm_rxy <-
#'   function(fml, data, dvnames) {
#'     mlm_res <- lm(fml, data = data)
#'     mlm_pred <- predict(mlm_res)
#'     cancor(mlm_pred, data[dvnames])$cor[[1]]^2
#'   }
#'
#' domir(
#'   cbind(wt, mpg) ~ vs + cyl + am + carb,
#'   mlm_rxy,
#'   .all = ~ carb,
#'   data = mtcars,
#'   dvnames = c("wt", "mpg")
#' )
#'
#'
#' ## Named sets
#' domir(
#'   mpg ~ am + gear + cyl + vs + qsec + drat,
#'   lm_r2,
#'   data = mtcars,
#'   .set =
#'     list( trns = ~ am + gear,
#'           eng = ~ cyl + vs, misc = ~ qsec + drat
#'     )
#' )
#'
#'
#' ## Linear model returning AIC
#' lm_aic <-
#'   function(fml, data) {
#'     lm_res <- lm(fml, data = data)
#'     AIC(lm_res)
#'  }
#'
#' domir(
#'   mpg ~ am + carb + cyl,
#'   lm_aic,
#'   .adj = TRUE,
#'   .rev = TRUE,
#'   data = mtcars
#'  )
#'
#'
#' ## 'systemfit' with 'formula_list' method returning AIC
#' if (requireNamespace("systemfit", quietly = TRUE)) {
#'   domir(
#'     formula_list(mpg ~ am + cyl + carb, qsec ~ wt + cyl + carb),
#'     function(fml) {
#'       res <- systemfit::systemfit(fml, data = mtcars)
#'       AIC(res)
#'     },
#'     .adj = TRUE, .rev = TRUE
#'   )
#' }
#'
domir <- function(.obj, ...) {
  UseMethod("domir")
}
#' @rdname domir
#' @exportS3Method
domir.formula <- function(
    .obj, .fct,
    .set = NULL, .wst = NULL, .all = NULL, .adj = FALSE,
    .cdl = TRUE, .cpt = TRUE, .rev = FALSE,
    ...) {
  # placeholder argument for within set evaluated names
  if (!is.null(.wst)) {
    .NotYetUsed(".wst")
  }
  # process and check 'formula' ----
  # use 'formula_parse()' to obtain formula inputs
  fml_parsed <- formula_parse(.obj)
  # no intercept-only check
  if (length(fml_parsed$rhs_names) == 0)
    stop("'.obj' cannot be intercept-only.", call. = FALSE)
  # logical vector to assist removing names from sub-setting processes
  rmv_frm_subst <- rep(FALSE, times = length(fml_parsed$rhs_names))
  # confirm .fct works as applied to .obj - returns value for full model
  full_model <- function_checker(.obj, .fct, ...)
  # estimate '.adj' value ----
  adj_model <- est_adj_model_fml(fml_parsed, .fct, .adj, ...)
  # process '.all' ----
  # ~~ put 'meta_domir' here - theb estimate 'all_model' afterward ~~ devise way to change 'rmv_frm_subst' and 'select_lgl from parsed formula ~~ ----
# TODO: disallow offsets and intercept removal - extend to 'formula_list' ----
  if (!is.null(.all)) {
    # '.all' must be a 'formula' check
    if (!inherits(.all, "formula")) {
      stop("'.all' must be a 'formula'.", call. = FALSE)
    }
    # parse '.all' using 'formula_parse()'
    all_parsed <- formula_parse(.all)
    # no intercept-only check
    if (length(all_parsed$rhs_names) == 0)
      stop("'.all' cannot be intercept-only.", call. = FALSE)
    #  no offsets in '.all' check
    if (!is.null(all_parsed$offset))
      stop("Offsets not allowed in '.all'.", call. = FALSE)
    # no removed intercepts in '.all' check
    if (!all_parsed$intercept_lgl)
      stop("Removing intercepts not allowed in '.all'.", call. = FALSE)
    # names in '.all' are in '.obj' check
    valid_all_lgl <- all_parsed$rhs_names %in% fml_parsed$rhs_names
    if (!all(valid_all_lgl)) {
      bad_all_names <- 
        paste(all_parsed$rhs_names[!valid_all_lgl], collapse = " ")
      stop("Name(s) ", bad_all_names,
           " in '.all' not found on the right hand side of '.obj'.",
           call. = FALSE)
    }
    # select '.all' names in .$select_lgl and remove as candidate subset
    pos_all_name <- which(fml_parsed$rhs_names %in% all_parsed$rhs_names)
    fml_parsed$select_lgl[pos_all_name] <- TRUE
    rmv_frm_subst[pos_all_name] <- TRUE
  }
  # # Process '.set' ----
  # if (!is.null(.set)) {
  # 
  #   if (!is.list(.set)) {
  #     stop("'.set' must be a 'list'.", call. = FALSE)
  #   }
  # 
  #   # apply checks to all elements of '.set' argument
  #   .mapply(fml_check, list(.set,
  #           paste0("Position ", 1:(length(.set)), " of '.set'."),
  #           rep(TRUE, times = length(.set))),
  #           NULL)
  # 
  #   # name vectors from each element of '.set'
  #   Set_names <-
  #     lapply(.set,
  #            function(x) attr(stats::terms(x), "term.labels"))
  # 
  #   # recursively check names in '.set'-s and remove from rhs
  #   for (set in 1:length(Set_names)) {
  # 
  #     rhs_names <-
  #       rhs_name_remover(Set_names[[set]], rhs_names,
  #                        paste0("Position ", set, " of '.set'."))
  # 
  #   }
  # 
  #   # apply labels to '.set'
  #   if (!is.null(names(.set)))
  #     Set_labels <- names(.set)
  # 
  #   else Set_labels <- paste0("set", 1:length(.set))
  # 
  #   missing_Set_labels <- which(Set_labels == "")
  # 
  #   if (length(missing_Set_labels) > 0)
  #     Set_labels[missing_Set_labels] <- paste0("set", missing_Set_labels)
  # 
  #   if (any(Set_labels %in% rhs_names)) {
  # 
  #     repeat_names <- Set_labels[which(Set_labels %in% rhs_names)]
  # 
  #     stop("Set element names ",
  #          paste(repeat_names, collapse = " "),
  #          " are also the names of elements in '.obj'.\n",
  #          "Please rename these '.set' elements.", call. = FALSE)
  #   }
  # 
  # }
  # 
  # else Set_names <- NULL
  # 
  # # Too few subsets error
  # if ((length(rhs_names) + length(Set_names)) < 2)
  #   stop("At least two subsets are needed for a dominance analysis.",
  #        call. = FALSE)
  # 
  # # Define meta-function to coordinate .fct calls ----
  # meta_domir <-
  #   function(Selector_lgl,
  #            RHS, LHS,
  #            .fct, .all, .adj,
  #            intercept, args_2_fct) {
  # 
  #     # logical vector to select subset of names
  #     selected_names <-
  #       unlist(RHS[Selector_lgl])
  # 
  #     # reconstruct formula to submit to '.fct'
  #     formula_to_use <-
  #       stats::reformulate(
  #         c(selected_names, .all, .adj),
  #         response = LHS, intercept = intercept)
  # 
  #     # submit formula and arguments to '.fct'
  #     returned_scalar <-
  #       do.call(.fct,
  #               append(formula_to_use, args_2_fct) )
  # 
  #     return(returned_scalar)
  # 
  #   }
  # 
  # # Check '.fct' inputs ----
  # 
  # # does '.fct' work?
  # test_model <-
  #   tryCatch(
  #     do.call(eval(.fct), append(.obj, list(...))),
  #     error = function(err)
  #       stop("'.fct' produced an error when applied to '.obj'.\n",
  #            "Also, check arguments passed to '.fct'",
  #            call. = FALSE)
  #     )
  # 
  # # is '.fct's returned value an atomic numeric scalar?
  # if (!is.numeric(test_model) || !is.vector(test_model) ||
  #     !is.atomic(test_model) || length(test_model) != 1)
  #   stop("result of '.fct' is not an atomic numeric, scalar/",
  #        "vector of length of 1 value.", call. = FALSE)
  # 
  # # is '.fct's returned value regarded as a list?
  # if (is.list(test_model))
  #   stop("result of '.fct' is a list.  It must be flattened/unlisted.",
  #        call. = FALSE)
  # 
  # # Define arguments to `dominance_scalar` ----
  # args_list <-
  #   list(RHS = append(rhs_names, Set_names),
  #        LHS = lhs_names,
  #        .fct = .fct,
  #        .all = All_names,
  #        .adj = Adj_names,
  #        intercept = intercept_lgl,
  #        args_2_fct = list(...))
  # 
  # cons_args <-
  #   list(RHS = rhs_names,
  #        LHS = lhs_names,
  #        .fct = .fct,
  #        .all = NULL,
  #        .adj = Adj_names,
  #        intercept = intercept_lgl,
  #        args_2_fct = list(...))
  # 
  # # Call `dominance_scalar` ----
  # return_list <-
  #   dominance_scalar(meta_domir, args_list, cons_args, test_model,
  #                .cdl, .cpt, .rev)
  # 
  # # Finalize returned values and attributes ----
  # 
  # if (is.null(.set))
  #   IV_Labels <- rhs_names
  # 
  # else IV_Labels <- c(rhs_names, Set_labels)
  # 
  # names(return_list$General_Dominance) <-
  #   IV_Labels
  # 
  # names(return_list$General_Dominance_Ranks) <-
  #   IV_Labels
  # 
  # if (.cdl)
  #   dimnames(return_list$Conditional_Dominance) <-
  #   list(names(return_list$General_Dominance),
  #        paste0("subset_size_", 1:length(return_list$General_Dominance)))
  # 
  # if (.cpt)
  #   dimnames(return_list$Complete_Dominance) <-
  #   list(paste0("Dmnates_", names(return_list$General_Dominance)),
  #        paste0("Dmnated_", names(return_list$General_Dominance)))
  # 
  # if (.rev == FALSE)
  #   Standardized <-
  #   return_list$General_Dominance /
  #   (
  #     return_list$Value -
  #       ifelse(length(return_list$Adj_result) > 0, return_list$Adj_result, 0)
  #   )
  # 
  # else Standardized <-
  #   -return_list$General_Dominance /
  #   -(
  #     return_list$Value -
  #       ifelse(length(return_list$Adj_result) > 0, return_list$Adj_result, 0)
  #   )
  # 
  # return_list <- list(
  #   "General_Dominance" = return_list$General_Dominance,
  #   "Standardized" = Standardized,
  #   "Ranks" = return_list$General_Dominance_Ranks,
  #   "Complete_Dominance" = return_list$Complete_Dominance,
  #   "Conditional_Dominance" = return_list$Conditional_Dominance,
  #   "Value" = return_list$Value,
  #   "Value_All" =
  #     return_list$All_result -
  #     ifelse(is.null(return_list$Adj_result), 0, return_list$Adj_result),
  #   "Value_Adjust" = return_list$Adj_result,
  #   "Call" = match.call()
  # )
  # 
  # # apply class 'domir'
  # class(return_list) <- c("domir")
  # 
  # return(return_list)
  list(
    over = fml_parsed, adj = adj_model
  )
}

#' @rdname domir
#' @exportS3Method 
domir.formula_list <- function(
    .obj, .fct, 
    .set = NULL, .wst = NULL, .all = NULL, .adj = FALSE,
    .cdl = TRUE, .cpt = TRUE, .rev = FALSE,
    ...) {
  
  # TODO: Unify domir.formula_list and domir.formula interfaces ----

  if (!is.null(.wst)) {
    
    .NotYetUsed(".wst")
    
  }
  
  # Process 'formula_list' ----
  # Obtain RHS, LHS, intercepts, and offsets in list
  list_parsed <- lapply(.obj, formula_parse)
  # no empty formulas check
  rhs_counts <- 
    sapply(list_parsed, 
         function(elem) {
           length(elem$rhs_names)
         })
  
  if (any(rhs_counts == 0)) {
    stop(paste("Each formula in '.obj' must have one or more terms.",
                "Formulas", paste(which(rhs_counts == 0), collapse = " "), "have no terms."),
               call. = FALSE)
  }
  
  # Record locations the .$select_lgl elements for all sub-lists
  element_count <- 
    sum( 
      sapply(
        1:length(list_parsed),
        function(elem) {
          length( list_parsed[[elem]][["select_lgl"]] )
        } ) )
  selector_locations <- 
    vector(mode = "list", length = element_count)
  pos = 1
  for ( elem in 1:length(list_parsed) ) {
    for ( loc in 1:length( list_parsed[[elem]][["select_lgl"]] ) ) {
      selector_locations[[pos]] <- c(elem, 5, loc)
      pos <- pos + 1
    } }
  
  # LHS ~ RHS pair list
  LHS_RHS_pairlist <- # use this when printing - also to check against sets
    unlist( lapply(
      list_parsed, 
      function(elem) {
        paste0(elem$lhs_names, "~", elem$rhs_names)
      } ) )
  
  remove_loc <- rep(FALSE, times = length(LHS_RHS_pairlist))
  
  # Apply '.adj' ----
  if (.adj) {
    
    # Generate '.adj' 'formula_list'
#TODO: reformulate here must include offset ----
    adj_fml_lst <- 
      lapply(
        list_parsed, 
        function(elem) {
          stats::reformulate("1", response = elem$lhs_names, 
                      intercept = elem$intercept_lgl)
        }
      ) 
    class(adj_fml_lst) <- c("formula_list", "list")
  }
  
  # Process '.all' ----
  if (!is.null(.all)) {
    
    # '.all' must be 'formula_list'
    if (!inherits(.all, "formula_list")) {
      stop("'.all' must be a 'formula_list'.", call. = FALSE)
    }
    
    # Obtain RHS, LHS, intercept, and offsets from '.all'
    all_parsed <- 
      lapply(.all, formula_parse)
    
    # no empty formulas check
    rhs_counts_all <- 
      sapply(all_parsed, 
             function(elem) {
               length(elem$rhs_names)
             })
    
    if (any(rhs_counts_all == 0)) {
      stop(paste("Each formula in '.all' must have one or more terms.",
                 "Formulas", paste(which(rhs_counts_all == 0), collapse = " "), "have no terms."),
           call. = FALSE)
    }
    
    # Check validity of LHS-RHS pairs in '.all'
    all_pairs <-
      lapply(
        all_parsed,
        function(elem) {
          paste0(elem$lhs_names, "~", elem$rhs_names)
        } )
    is_valid_pair <- 
      unlist(all_pairs) %in% LHS_RHS_pairlist
    if ( !all(is_valid_pair) ) 
      stop("Pair ", paste(unlist(all_pairs)[!is_valid_pair], collapse = " "), 
           " in '.all' not found among those in '.obj'.", call. = FALSE)

    # Activate '.all' pairs and remove as candidate subset
    for ( elem in unlist(all_pairs) ) {
      list_parsed[[ selector_locations[[ which(LHS_RHS_pairlist %in% elem) ]] ]] <- 
        TRUE
    }
    remove_loc[ which(LHS_RHS_pairlist %in% unlist(all_pairs) )] <- TRUE

  }
  
  # Process '.set' ----
  if (!is.null(.set)) {
    
    # '.set' must be list
    if (!is.list(.set)) {
      stop("'.set' must be a 'list'.", call. = FALSE)
    }
    
    # '.set' must be comprised of 'formula_list's
    set_fmllst <- 
      sapply(.set, 
             function(elem) {
               inherits(elem, "formula_list")
             } )
    
    if (!all(set_fmllst)) 
      stop("'.set' element ", paste( which(!set_fmllst), collapse = " " ), 
           " not of class 'formula_list'.", call. = FALSE)
    
    # Obtain RHS, LHS, intercept, and offsets from '.set'
    sets_parsed <- 
      lapply(.set, 
             function(set) {
               lapply(set, formula_parse)
             } )
    
    # no empty formulas check
    rhs_counts_sets <- 
      lapply(sets_parsed, 
             function(set) {
               sapply(set, 
                      function(elem) {
                        length(elem$rhs_names)
                      })
             })
    
    if (any(unlist(rhs_counts_sets) == 0)) {
      stop(paste("Each formula in '.set' must have one or more terms.",
                 "The following:", 
                 paste(sapply(seq_len(length(rhs_counts_sets)), 
                              function(elem) {
                                no_term <- 
                                  which(rhs_counts_sets[[elem]] == 0)
                                if (length(no_term) > 0) 
                                  paste("list", elem, "formulas", 
                                        paste0(paste(no_term, collapse = " "), ","))
                                else NULL
                              }),
                       collapse = " "), " have no terms."),
           call. = FALSE)
    }
    
    # Check validity of LHS-RHS pairs in '.set'
    set_pairs <-
      lapply(
        sets_parsed,
        function(sublist) {
          unlist( lapply(
            sublist, 
            function(elem) {
              paste0(elem$lhs_names, "~", elem$rhs_names)
            } ) ) } )
    is_valid_pair <- 
      unlist(set_pairs) %in% LHS_RHS_pairlist
    if ( !all(is_valid_pair) ) 
      stop("Pair ", paste(unlist(set_pairs)[!is_valid_pair], collapse = " "), 
                 " in '.set' not found among those in '.obj'.", call. = FALSE)

    # Group together locations of the .$select_lgl elements for '.set's
    selector_locations_sets <- 
      vector(mode = "list", length = length(.set))
    pos = 1
    for ( elem in 1:length(sets_parsed) ) {
      selector_locations_sets[[pos]] <- 
        lapply( set_pairs[[elem]] , 
                function(loc) {
                  selector_locations[[ which(LHS_RHS_pairlist %in% loc) ]]
                } )
      remove_loc[which(LHS_RHS_pairlist %in% set_pairs[[elem]])] <- TRUE
      pos <- pos + 1
    }
    
  }
  else selector_locations_sets <- NULL
  
  # Update selector locations with .set
  selector_locations <- 
    append(selector_locations[!remove_loc], selector_locations_sets)
  
  # apply labels to '.set'
  if (!is.null(names(.set))) 
    Set_labels <- names(.set)
  
  else Set_labels <- paste0("set", 1:length(.set))
  
  missing_Set_labels <- which(Set_labels == "")
  
  if (length(missing_Set_labels) > 0)
    Set_labels[missing_Set_labels] <- paste0("set", missing_Set_labels)
  
  if (any(Set_labels %in% LHS_RHS_pairlist)) {
    
    repeat_names <- Set_labels[which(Set_labels %in% LHS_RHS_pairlist)]
    
    stop("Set element names ",
         paste(repeat_names, collapse = " "), 
         " are also the names of elements in '.obj'.\n",
         "Please rename these '.set' elements.", call. = FALSE)
  }
  
  # Define meta-function to coordinate .fct calls ----
  meta_domir <- 
    function(Selector_lgl, 
             list_parsed, .fct, selector_locations,
             args_2_fct, ...) {

      # distribute logical vector to elements of formula_list
      for (elem in selector_locations[Selector_lgl]) {
        if ( is.list(elem) ) {
          for (pos in 1:length(elem)) {
            list_parsed[[ elem[[pos]] ]] <- TRUE
          } }
        else list_parsed[[elem]] <- TRUE
      }
      
      # reconstruct the selected formula_list
      fml_lst <- 
        lapply(
          list_parsed, 
          function(elem) {
            if ( all(!elem$select_lgl) )
              stats::reformulate(c("1", elem$offset), response = elem$lhs_names, 
                          intercept = elem$intercept_lgl)
            else
              stats::reformulate(
                c(elem$RHS_name[elem$select_lgl], elem$offset),
                response = elem$lhs_names, intercept = elem$intercept_lgl) } )
      
      class(fml_lst) <- c("formula_list", "list")
      
      # submit formula_list to '.fct'
      returned_scalar <-
        do.call(.fct,
                append(list(fml_lst), args_2_fct) )
      
      return(returned_scalar)
      
    }
  
  # Confirm .fct works as applied to .obj - fits full-model
  test_model <- function_checker(.obj, .fct, ...)
  
  # Estimate .adj model if applicable
  if (.adj) {
    adj_model <- function_checker(adj_fml_lst, .fct, ...)
  }
  else adj_model <- NULL
  
  # Estimate .all model if applicable
  if (!is.null(.all))
    all_model <- 
    meta_domir(rep(FALSE, times = length(selector_locations)),
               list_parsed, .fct, selector_locations, 
               args_2_fct = list(...))
  else all_model <- NULL 
  
  # Define arguments to `dominance_scalar` ----
  args_list <- 
    list(RHS = selector_locations,
         list_parsed = list_parsed,
         .fct = .fct,
         selector_locations = selector_locations,
         .all = all_model, .adj = adj_model,
         args_2_fct = list(...))
  
  return_list <-
    dominance_scalar(
      meta_domir,
      args_list,
      NULL, test_model,
      .cdl, .cpt, .rev)
  
  # Finalize returned values and attributes ----
  
  if (is.null(.set)) 
    IV_Labels <- LHS_RHS_pairlist[!remove_loc]
  
  else IV_Labels <- c(LHS_RHS_pairlist[!remove_loc], Set_labels)
  
  names(return_list$General_Dominance) <- 
    IV_Labels
  
  names(return_list$General_Dominance_Ranks) <- 
    IV_Labels
  
  if (.cdl)
    dimnames(return_list$Conditional_Dominance) <- 
    list(names(return_list$General_Dominance), 
         paste0("subset_size_", 1:length(return_list$General_Dominance)))
  
  if (.cpt)
    dimnames(return_list$Complete_Dominance) <- 
    list(paste0("Dmnates_", names(return_list$General_Dominance)),  
         paste0("Dmnated_", names(return_list$General_Dominance)))
  
  if (!.rev) 
    Standardized <- 
    return_list$General_Dominance / 
    (
      return_list$Value - 
        ifelse(length(return_list$Adj_result) > 0, return_list$Adj_result, 0)
    ) 
  
  else Standardized <- 
    -return_list$General_Dominance / 
    -(
      return_list$Value - 
        ifelse(length(return_list$Adj_result) > 0, return_list$Adj_result, 0)
    )
  
  return_list <- list(
    "General_Dominance" = return_list$General_Dominance,
    "Standardized" = Standardized,
    "Ranks" = return_list$General_Dominance_Ranks,
    "Complete_Dominance" = return_list$Complete_Dominance,
    "Conditional_Dominance" = return_list$Conditional_Dominance,
    "Value" = return_list$Value,
    "Value_All" = 
      return_list$All_result - 
      ifelse(is.null(return_list$Adj_result), 0, return_list$Adj_result),
    "Value_Adjust" = return_list$Adj_result,
    "Call" = match.call()
  )
  
  # apply class 'domir'
  class(return_list) <- c("domir") 
  
  return_list
    
}
# `domir` helper functions ----
#' @title Internal formula parsing function
#' @description Internal formula parsing function to facilitate re-construction
#' of a formula using `reformulate()`
#'
#' Not intended to be called by the user.
#'
#' @keywords internal
#'
#' @name formula_parse
#'
#' @rdname formula_parse
#'
#' @export
formula_parse <- function(.obj) {
  # obtain 'right hand side'/RHS name vector from `.obj`
  rhs_names <- attr(stats::terms(.obj), "term.labels")
  # intercept logical - confirms inclusion of intercept
  intercept_lgl <- as.logical(attr(stats::terms(.obj), "intercept"))
# TODO: offset wrong when no IVs given use of rownames ----
  # obtain offsets
  if (!is.null(attr(stats::terms(.obj), "offset"))) {
    offset_locs <- attr(stats::terms(.obj), "offset")
    offset <- rownames(attr(stats::terms(.obj), "factors"))[offset_locs]
  } else {
    offset <- NULL
  }
  # obtain 'left hand side'/LHS (if included)
  if (attr(stats::terms(.obj), "response") == 1) {
    lhs_names <- attr(stats::terms(.obj), "variables")[[2]]
  } else {
    lhs_names <- NULL
  }
  # set logical selection vector; defaults to not selected
  select_lgl <- rep(FALSE, times = length(rhs_names))
  list(rhs_names = rhs_names,
       lhs_names = lhs_names,
       intercept_lgl = intercept_lgl,
       offset = offset,
       select_lgl = select_lgl)
}

function_checker <- function(.obj, .fct, ...) { # make function checker a generic for 'formula' and 'list'? -- or make formula a list then simplify it before submitting?
  obj_submit <- 
    switch(class(.obj)[[1]],
      "formula_list" = list(.obj),
      "formula" = .obj
    )
  # does '.fct' work?
  test_model <- 
    tryCatch(
      do.call(eval(.fct), append(obj_submit, list(...))), # need to use list(.obj) for list method and just .obj for fml?
      error = function(err) 
        stop("'.fct' produced an error when applied to '.obj'.\n", 
             "Have you checked the arguments passed to '.fct'?", 
             call. = FALSE)
    )
  
  # is '.fct's returned value an atomic numeric scalar?
  if (!is.numeric(test_model) || !is.vector(test_model) || 
      !is.atomic(test_model) || length(test_model) != 1) 
    stop("Result of '.fct' is not an atomic, numeric, scalar object ",
         "(vector with a 'length()' value of 1).", call. = FALSE)
  
  return(test_model)
  
}

est_adj_model_fml <- function(fml_prs, .fct, .adj, ...) {
  if (.adj) {
    adj_fml <-
      stats::reformulate(c("1", fml_prs$offset),
                         response = fml_prs$lhs_names,
                         intercept = fml_prs$intercept_lgl)
    adj_model <- function_checker(adj_fml, .fct, ...)
  } else {
    adj_model <- NULL
  }
  adj_model
}
#' @title Print method for `domir`
#' @description Reports formatted results from `domir` class object.
#' @param x an object of class "domir".
#' @param ... further arguments passed to [`print.default`].
#' @return The submitted "domir" object, invisibly.
#' @details The print method for class `domir` objects reports out the 
#' following results:
#' \itemize{
#'  \item{Value when all elements are included in `obj`.}
#'  \item{Value for the elements included in `.all`, if any.}  
#'  \item{Value for the elements included in `.adj`, if any.}
#'  \item{Matrix describing general dominance values, standardized 
#'  general dominance values, and the ranking of the general 
#'  dominance values.}
#'  \item{Matrix describing the conditional dominance values, if computed}
#'  \item{Matrix describing the complete dominance designations, if evaluated}
#'  \item{If following [`summary.domir`], matrix describing the strongest 
#'  dominance designations between all elements.}}
#'  
#'  The `domir` print method alters dimension names for readability and they 
#'  do not display as stored in the `domir` object.
#'  
#' @exportS3Method 

# `domir` printing methods ----
print.domir <- function(x, ...) {
  
  cat("Overall Value:     ", x[["Value"]], "\n")
  
  if (length(x[["Value_All"]]) > 0) 
    cat("All Subset Value:  ", x[["Value_All"]], "\n")
  
  if (length(x[["Value_Adjust"]]) > 0) 
    cat("Adjustment Value:  ", 
        x[["Value_Adjust"]], "\n")
  
  cat("\n")
  
  cat("General Dominance Values:\n")
  
  Display_Std <- 
    t(rbind(x[["General_Dominance"]], x[["Standardized"]], x[["Ranks"]]))
  
  dimnames(Display_Std) <- 
    list(names(x[["Ranks"]]), c("General Dominance", "Standardized", "Ranks"))
  
  print(Display_Std, ...)
  
  cat("\n")
  
  if (length(x[["Conditional_Dominance"]] > 0)) {
    
    cat("Conditional Dominance Values:\n")
    
    colnames(x[["Conditional_Dominance"]]) <- 
      paste("Subset Size:", 1:ncol(x[["Conditional_Dominance"]]))
    
    print(x[["Conditional_Dominance"]], ...)
    
    cat("\n")
    
  }
  
  if (length(x[["Complete_Dominance"]] > 0)) {
    
    cat("Complete Dominance Designations:\n")
    
    colnames(x[["Complete_Dominance"]]) <- 
      gsub("^Dmnated_", "Dmnated?", colnames(x[["Complete_Dominance"]]))
    
    rownames(x[["Complete_Dominance"]]) <- 
      gsub("^Dmnates_", "Dmnates?", rownames(x[["Complete_Dominance"]]))
    
    print(x[["Complete_Dominance"]])
    
    cat("\n")
    
  }
  
  if (length(x[["Strongest_Dominance"]] > 0)) {
    
    cat("Strongest Dominance Designations:")
    
    print(x[["Strongest_Dominance"]])
    
    cat("\n")
    
  }
  
  invisible(x)
  
}


#' @title Summary method for `domir`
#' @description Reports dominance designation results from the `domir` 
#' class object.
#' @param object an object of class "domir".
#' @param ... further arguments passed to or from other methods. 
#' Not used currently.
#' @return The submitted "domir" object with an additional 
#' `Strongest_Dominance` element added.
#' \describe{
#'  \item{\code{Strongest_Dominance}}{Matrix comparing the element in the first 
#'  row to the element in the third row.  The second row denotes the strongest 
#'  designation between the two elements.}
#' }
#' 
#' @details The summary method for class `domir` objects is used for obtaining 
#' the strongest dominance designations (i.e., general, conditional, or 
#' complete) among all pairs of dominance analyzed elements.
#' 
#' @exportS3Method

summary.domir <- function(object, ...) {
  
  if (length(object[["Strongest_Dominance"]]) == 0) {
    
    reverse <- as.list(object$Call)$.rev
    if (is.null(reverse)) reverse <- FALSE
    
    reverse_cdl <- 
      ifelse(reverse, 
             rep(-1, times = length(object$General_Dominance)), 
             rep(1, times = length(object$General_Dominance)))
    
    reverse_gnl <- ifelse(reverse, -1, 1)
    
    pairs <- utils::combn(names(object$General_Dominance), 2)
    
    pairs <- rbind(pairs[1,], rep("", times = ncol(pairs)), pairs[2,])
    
    location <- 0
    
    for (IV1 in 1:(length(object$General_Dominance) - 1)) {
      
      for (IV2 in (IV1+1):length(object$General_Dominance)) {
        
        location <- location + 1
        
        if (length(object[["Complete_Dominance"]] > 0)) {
          
          if (!is.na(object$Complete_Dominance[IV1, IV2])) {
            pairs[2, location] <- 
              ifelse(object$Complete_Dominance[IV1, IV2], 
                     "completely dominates",
                     "is completely dominated by")
            
            next
            
          }
          
        }
        
        if (length(object[["Conditional_Dominance"]] > 0)) {
          
          if (all(object$Conditional_Dominance[IV1,]*reverse_cdl > 
                  object$Conditional_Dominance[IV2,]*reverse_cdl)) { 
            
            pairs[2, location] <- "conditionally dominates"
            
            next
            
          }
          
          else if (all(object$Conditional_Dominance[IV1,]*reverse_cdl < 
                       object$Conditional_Dominance[IV2,]*reverse_cdl)) {
            
            pairs[2, location] <- "is conditionally dominated by"
            
            next
            
          }
          
        }
        
        pairs[2, location] <- 
          ifelse(object$General_Dominance[[IV1]]*reverse_gnl > 
                   object$General_Dominance[[IV2]]*reverse_gnl, 
                 "generally dominates", 
                 ifelse(object$General_Dominance[[IV1]]*reverse_gnl < 
                          object$General_Dominance[[IV2]]*reverse_gnl,
                        "is generally dominated by", 
                        "has no dominance designation with"))
        
      }
      
    }
    
    rownames(pairs) <- rep("", times = 3)
    
    colnames(pairs) <- rep("", times = ncol(pairs))
    
    res <- append(object, list(Strongest_Dominance = pairs))
    
    class(res) <- c("domir")
    
    return(res)
    
  }
  
  else return(object)
  
}
