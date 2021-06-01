Tools to Support Relative Importance Analysis
================

# Overview

The intention of the `{domir}` package is to provide tools that allow
relative importance analysis across a wide variety of practical data
analytic situations. With `{domir}`, if you have a statistical/machine
learning model and an extractor function to obtain a fit metric, you can
conduct a relative importance analysis.

More specifically, `{domir}` is intended to implement a set of flexible
wrapper and helper functions for conducting relative importance analysis
and the current implementation of the package has a focus on dominance
analysis with the `domin` function as a relative importance analysis
method. Currently, `domin` is the only importance function that is
implemented.

For readers looking to familiarize themselves more with the dominance
analysis methodology, a more extensive conceptual discussion of
dominance analysis (which focuses on the Stata version of `domin`) as a
method is available
[here](https://github.com/jluchman/domin/blob/master/README.md).

# Installation

To install the most recent stable version of `domir` from CRAN use:

`install.packages("domir")`

To install the working development version of `{domir}` using the
`devtools` package use:

`devtools::install_github("https://github.com/jluchman/domir")`

# What `{domir}` Does

Before discussing details of the `{domir}` package, I provide some
examples of what `{domir}` can do.

The focus of this section is on outlining how `domir::domin` extends
existing packages and on the structure of the function.

## Comparison with Existing Relative Importance Packages

Fundamentally, `domir::domin` is an extension of the “lmg” type for the
`calc.relimpo` function in the `{relaimpo}` package as well as the
`dominanceAnalysis` function in the `{dominanceanalysis}` package.

`domir::domin` can replicate the results produced by the above packages
but, as will be seen, requires a “deconstructed” the model to be
submitted to it. This difference in structure does make `domin` more
complex but also allows the function a great deal more flexibility in
terms of the kinds of models and fit metrics that can be dominance
analyzed.

Before discussing some of the elements that make `domin` flexible,
consider the following example that shows how `domin` is similar to
existing packages. All three of the dominance analysis results to come
are based on the following linear model:

`lm(mpg ~ am + vs + cyl, data = mtcars)`

The variance explained *R*<sup>2</sup> will be the focal fit metric.

### `{domir}`’s `domin`

``` r
domin(mpg ~ am + vs + cyl, 
      lm, 
      list(summary, "r.squared"), 
      data = mtcars)
```

    ## Overall Fit Statistic:      0.7619773 
    ## 
    ## General Dominance Statistics:
    ##     General Dominance Standardized Ranks
    ## am          0.1774892    0.2329324     3
    ## vs          0.2027032    0.2660226     2
    ## cyl         0.3817849    0.5010450     1
    ## 
    ## Conditional Dominance Statistics:
    ##        IVs: 1    IVs: 2      IVs: 3
    ## am  0.3597989 0.1389842 0.033684441
    ## vs  0.4409477 0.1641982 0.002963748
    ## cyl 0.7261800 0.3432799 0.075894823
    ## 
    ## Complete Dominance Designations:
    ##             Dmnated?am Dmnated?vs Dmnated?cyl
    ## Dmnates?am          NA         NA       FALSE
    ## Dmnates?vs          NA         NA       FALSE
    ## Dmnates?cyl       TRUE       TRUE          NA

In `domin`, the `lm` model is submitted in pieces. Specifically, the key
inputs were the formula (`mpg ~ am + vs + cyl`) and the model to be
called using the formula (`lm`). In this way, `domin` is a `Map`- or
`apply`-like function as it receives an object on which to operate
(i.e., the formula) and a function to which to apply to it.

In being like `apply`, `domin` is agnostic to the fit metric to use for
the model called and it must be supplied with a list outlining extractor
function information (`list(summary, "r.squared")`; described further in
the [Details](#Details) section).

Like `apply`, other arguments (`data = mtcars`) can also be passed to
each call of `lm`.

The focus of `domin`’s `print`-ed results focuses on the numerical
results from “General Dominance Statistics” and “Conditional Dominance
Statistics” and, a logical matrix of “Complete Dominance Designations”.

### `{relaimpo}`’s `calc.relimp` with “lmg”

``` r
relaimpo::calc.relimp(mpg ~ am + vs + cyl, 
                      data = mtcars, 
                      type = "lmg")
```

    ## Response variable: mpg 
    ## Total response variance: 36.3241 
    ## Analysis based on 32 observations 
    ## 
    ## 3 Regressors: 
    ## am vs cyl 
    ## Proportion of variance explained by model: 76.2%
    ## Metrics are not normalized (rela=FALSE). 
    ## 
    ## Relative importance metrics: 
    ## 
    ##           lmg
    ## am  0.1774892
    ## vs  0.2027032
    ## cyl 0.3817849
    ## 
    ## Average coefficients for different model sizes: 
    ## 
    ##            1X       2Xs       3Xs
    ## am   7.244939  4.316851  3.026480
    ## vs   7.940476  2.995142  1.294614
    ## cyl -2.875790 -2.795816 -2.137632

`{relaimpo}`’s `calc.relimp` accepts only `lm` models and the variance
explained *R*<sup>2</sup> as a fit metric. As a result, the function
does not ask for model or fit metric type.

The function’s printed results provide the “lmg” relative importance
statistics (i.e., General Dominance Statistics) and, in addition, report
the average `lm` coefficients across all models.

### `{dominanceanalysis}`’s `dominanceAnalysis`

``` r
dominanceanalysis::dominanceAnalysis(lm(mpg ~ am + vs + cyl, 
                                        data = mtcars))
```

    ## 
    ## Dominance analysis
    ## Predictors: am, vs, cyl 
    ## Fit-indices: r2 
    ## 
    ## * Fit index:  r2 
    ##     complete conditional general
    ## am                              
    ## vs                            am
    ## cyl    am,vs       am,vs   am,vs
    ## 
    ## Average contribution:
    ##   cyl    vs    am 
    ## 0.382 0.203 0.177

`{dominanceanalysis}`’s `dominanceAnalysis` implements dominance
analysis for specific models, of which `lm` is a supported model.
`dominanceAnalysis` accepts a fitted `lm` object as input and uses the
explained variance *R*<sup>2</sup> as the fit metric.

`dominanceAnalysis`’s printed output is focused on qualitative dominance
designations but also reports the, magnitude sorted, average
contribution (i.e., General Dominance Statistic) values.

## How `{domir}` Extends on Previous Packages

The intention of `{domir}` is to extend relative importance to new data
analytic situations the user might encounter where a dominance analysis
could be valuable.

The sections below outline some pertinent examples of specific models
that the `domin` function can accommodate.

### Linear Model Revisited

`domin` is fit metric agnostic and, as such, one component of its
flexibility is in allowing the user to apply any applicable fit metric
for a model for the purposes of relative importance analysis.

In this example, the explained variance *R*<sup>2</sup> is swapped with
an alternative, but nonetheless applicable, fit metric: the McFadden
pseudo-*R*<sup>2</sup> as implemented by the `{pscl}` package.

Note the use of the pipes to `capture.output` and `invisible`. These are
not not strictly necessary but if not used will print far more output
than is needed as `pscl::pR2` is a rather verbose function and will
print a message for each model fitted.

``` r
(mcf_da_lm <- 
   domin(mpg ~ am + vs + cyl, 
         lm, 
         list(pscl::pR2, "McFadden"), 
         data = mtcars)) |> 
  capture.output() |> 
  invisible()

mcf_da_lm
```

    ## Overall Fit Statistic:      0.2243283 
    ## 
    ## General Dominance Statistics:
    ##     General Dominance Standardized Ranks
    ## am         0.04848726    0.2161442     3
    ## vs         0.04970277    0.2215627     2
    ## cyl        0.12613826    0.5622931     1
    ## 
    ## Conditional Dominance Statistics:
    ##         IVs: 1     IVs: 2      IVs: 3
    ## am  0.06969842 0.05507782 0.020685547
    ## vs  0.09088103 0.05629333 0.001933959
    ## cyl 0.20243215 0.13272881 0.043253806
    ## 
    ## Complete Dominance Designations:
    ##             Dmnated?am Dmnated?vs Dmnated?cyl
    ## Dmnates?am          NA         NA       FALSE
    ## Dmnates?vs          NA         NA       FALSE
    ## Dmnates?cyl       TRUE       TRUE          NA

Note that this fit metric produces effectively the same answers, in
terms of qualitative importance inferences about the terms, as that from
the explained variance *R*<sup>2</sup>.

### Ordered Logistic Regression

`domin` acts like an `apply` function for models and does not have built
in methods. This is another component of its flexibility as it can
accommodate functions that, to this point, have not been supported in
relative importance analysis. One pertinent example is the `polr`
function from the `{MASS}` package also using `pscl::pR2` as fit metric.

``` r
mtcars2 <- data.frame(mtcars, carb2 = as.factor(mtcars$carb))

(da_polr <- 
    domin(carb2 ~ am + vs + mpg, 
          MASS::polr, 
          list(pscl::pR2, "McFadden"), 
          data = mtcars2)) |> 
  capture.output() |> 
  invisible()

da_polr
```

    ## Overall Fit Statistic:      0.2647682 
    ## 
    ## General Dominance Statistics:
    ##     General Dominance Standardized Ranks
    ## am         0.04221668    0.1594477     3
    ## vs         0.09264306    0.3499026     2
    ## mpg        0.12990844    0.4906497     1
    ## 
    ## Conditional Dominance Statistics:
    ##          IVs: 1     IVs: 2     IVs: 3
    ## am  0.001505741 0.05272927 0.07241503
    ## vs  0.161029601 0.10315565 0.01374394
    ## mpg 0.151278401 0.14042103 0.09802589
    ## 
    ## Complete Dominance Designations:
    ##             Dmnated?am Dmnated?vs Dmnated?mpg
    ## Dmnates?am          NA         NA       FALSE
    ## Dmnates?vs          NA         NA          NA
    ## Dmnates?mpg       TRUE         NA          NA

### Decision Trees

`domin` can also accept models that do not produce model coefficients
like `rpart::rpart`.

``` r
domin(mpg ~ am + vs + cyl, 
      rpart::rpart, 
      list(\(model) 
            list(R2 = 1-model$cptable[nrow(model$cptable), 3]), 
            "R2"),
      data = mtcars)
```

    ## Overall Fit Statistic:      0.7324601 
    ## 
    ## General Dominance Statistics:
    ##     General Dominance Standardized Ranks
    ## am          0.1199330    0.1637400     3
    ## vs          0.1605074    0.2191346     2
    ## cyl         0.4520197    0.6171254     1
    ## 
    ## Conditional Dominance Statistics:
    ##        IVs: 1     IVs: 2    IVs: 3
    ## am  0.3597989 0.00000000 0.0000000
    ## vs  0.4409477 0.04057437 0.0000000
    ## cyl 0.7324601 0.33208674 0.2915124
    ## 
    ## Complete Dominance Designations:
    ##             Dmnated?am Dmnated?vs Dmnated?cyl
    ## Dmnates?am          NA         NA       FALSE
    ## Dmnates?vs          NA         NA       FALSE
    ## Dmnates?cyl       TRUE       TRUE          NA

Note that an anonymous function can be used as a valid submission to the
`fitstat` argument. In this case, the anonymous function transforms and
extracts the proportion of error from the `rpart` object. If the model
object returns its own fit statistic, it can be extracted using an
anonymous function.

### Multinomial Logistic (softmax) Regression with Extra Features

`domin`, similar to other packages, can combine multiple terms into a
single set as well as use one or more terms as covariate(s) in all model
subsets.

This example outlines another model, `multinom` from the `{nnet}`
package,  
another function that has not been accommodated in relative importance
packages, that uses sets and all/covariate terms.

In addition, `complete = FALSE` which saves a little computation time
and suppresses reporting complete dominance designations.

``` r
(da_mnl <- 
  domin(carb2 ~ mpg, 
      nnet::multinom, 
      list(pscl::pR2, "McFadden"), 
      sets = list(c("am", "vs"), c("cyl", "disp")),
      all = c("gear"),
      complete = FALSE,
      data = mtcars2)) |> 
  capture.output() |> 
  invisible()

da_mnl
```

    ## Overall Fit Statistic:      0.9282015 
    ## All Subsets Fit Statistic:  0.1393919 
    ## 
    ## General Dominance Statistics:
    ##      General Dominance Standardized Ranks
    ## mpg          0.2958544    0.3187394     2
    ## set1         0.1770852    0.1907832     3
    ## set2         0.3158700    0.3403033     1
    ## 
    ## Conditional Dominance Statistics:
    ##         IVs: 1    IVs: 2    IVs: 3
    ## mpg  0.4452671 0.2553281 0.1869679
    ## set1 0.2886101 0.1365589 0.1060867
    ## set2 0.5769312 0.2753437 0.0953351
    ## 
    ## Components of sets:
    ## set1 : am vs 
    ## set2 : cyl disp 
    ## 
    ## All subsets variables: gear

``` r
da_mnl$Subset_Details$Full_Model
```

    ## [1] "carb2 ~ mpg + am + vs + cyl + disp"

The `domin` automatically combines the entries in the `formula_overall`,
`sets`, and `all` arguments. The full model formula can be obtained from
the `domin` object in the `.$Subset_Details$Full_Model` element.

### Zero-Inflated Poisson with Wrapper Function

Although `domin` can work directly with modeling functions that accept
standard formula, more complex formulas such as those used by models
such as `zeroinfl` models from the package `{pscl}` can also be
accommodated using wrapper functions.

The below wrapper function`zinfl_wrap` uses the entries in the formula
to create a symmetric count and zero-inflation formulas that will be
submitted to `zeroinfl` model.

In an effort to illustrate what each model submitted to `zeroinfl` looks
like, the model formula for all 7 models is printed before each run.

``` r
zinfl_wrap <- function(model, ...) {
  zip_terms <- model |> terms() |> attr("term.labels") |> paste(collapse = " + ")
  zip_formula_rhs <- zip_terms |> rep(times = 2) |> paste(collapse = " | ")
  zip_formula_lhs <- (model |> all.vars())[[1]]
  zip_formula <- c(zip_formula_lhs, zip_formula_rhs) |> paste(collapse = " ~ ") |> as.formula()
  print(deparse(zip_formula))
  pscl::zeroinfl(zip_formula, ...)
}

domin(art ~ fem + mar + kid5, 
      zinfl_wrap,
      list(\(model) {capture.output(result <- pscl::pR2(model)); result}, "McFadden"), 
      data=pscl::bioChemists)
```

    ## [1] "art ~ fem | fem"
    ## [1] "art ~ mar | mar"
    ## [1] "art ~ kid5 | kid5"
    ## [1] "art ~ fem + mar | fem + mar"
    ## [1] "art ~ fem + kid5 | fem + kid5"
    ## [1] "art ~ mar + kid5 | mar + kid5"
    ## [1] "art ~ fem + mar + kid5 | fem + mar + kid5"

    ## Overall Fit Statistic:      0.009101817 
    ## 
    ## General Dominance Statistics:
    ##      General Dominance Standardized Ranks
    ## fem       0.0059812901   0.65715343     1
    ## mar       0.0008482014   0.09319035     3
    ## kid5      0.0022723252   0.24965622     2
    ## 
    ## Conditional Dominance Statistics:
    ##            IVs: 1      IVs: 2      IVs: 3
    ## fem  0.0054489923 0.006059012 0.006435866
    ## mar  0.0005852711 0.000925923 0.001033410
    ## kid5 0.0008854100 0.002350047 0.003581519
    ## 
    ## Complete Dominance Designations:
    ##              Dmnated?fem Dmnated?mar Dmnated?kid5
    ## Dmnates?fem           NA        TRUE         TRUE
    ## Dmnates?mar        FALSE          NA        FALSE
    ## Dmnates?kid5       FALSE        TRUE           NA

Further discussion of how to generate wrapper commands is outlined below
in the [Details](#Details) section.

------------------------------------------------------------------------

# Details

Having provided some examples of what `domin` can do, this section moves
on to outline details of how `domin` works by way of the structure of
the function.

`domin` estimates models for all possible subsets of the terms submitted
to it by repeatedly calling different models and collecting their
results. `domin` takes inspiration from the `apply` family of functions
and works in a similar way - invoking repeated function calls from three
‘building block’ arguments:

1.  a formula
2.  a modeling function
3.  list of instructions to call an extractor function that obtains a
    model fit metric

You can think of `domin` as repeatedly invoking the following process:

`modeling_function(formula) |> fit_metric()`

Hence, the modeling function is called using the formula input (i.e., a
subset of the formuala elements) as its argument and the results of the
modeling function are are ‘piped’ to the fit metric function that is
used for dominance statistic computation.

In the sections below, each of the three arguments to `domin` is
discussed in greater detail.

## 1) Formula

The first, and most important, argument for `domin` is the formula
input. Understanding how the formula input is constructed and submitted
to each call of the modeling function is important for the effective use
of `domin`.

The formula components are the most important pieces of `domin` as they
directly define the terms used to dominance analyze the model and, thus,
the scope of all subsets of models.

The formula input is derived from three arguments in `domin`.

1.  the `formula_overall` argument
2.  the `sets` argument
3.  the `all` argument

***section in progress***

### `formula_overall`

The `formula_overall` argument must take the form of `response ~ terms`
as that is a standard format for many modeling functions such as `lm`
and `glm`.

The entries on the right hand side of the formula are parsed using the
`terms` utilities in the `stats` library. The desired behavior is that
variable names separated by `+` are divided up and used as terms in
computing all combinations.

It is also important to note that all special formula processing is
applied to the formula including the use of `I()`, `:`, `*`, and
`offset()`.  
The actual list of terms on which `domin` computes all combinations is
obtained from `attr(,"term.labels")`. If the user wants to test to see
what `domin` will do with the formula submitted use:
`formula(.) |> terms() |> attr("term.labels")` to see which labels it
will produce (note also that `domin` does not have a method to
accommodate second or higher order terms like `{relaimpo}` and will
issue a warning when second or higher order terms are used).

### `sets`

After processing, the formula object is combined with entries from the
`sets` argument. The entries in the `sets` argument are minimally
processed and are combined

and is not ‘data frame aware’. That is, shorthand such as `~ .` will not
work to select variables in a data frame even if a `data` argument is
supplied to the `domin` function.

``` r{formula_example}
formula(mpg ~ .) |> terms(data=mtcars) |> formula()
```

…note all and sets go here as well

many modeling functions that accept a formula that follow the standard
`response ~ terms` format.

`domin` works closely with the `do.call` function and the structure of
the `domin` function is similar as a result. Users are encouraged to
read the documentation for the `do.call` function to understand how
`domin` is implemented and for bug/error checking.

wrapper function that can be used with The format used by domin can be
extended to other functions focused on independent
variable/predictor/feature-based RI can be accommodated with a
additional wrapper functions based on the formula it creates and submits
to modeling functions. Some examples are offered below.

## Modeling Function

## Fit Statistic Extractor Function

# Consideraions

This section outlines a few key considerations

`domin` does *not* check to ensure that the sample underlying the
modeling is consistent across modeling runs on which the dominance
statistics are computed. If the modeling function uses a `na.action`
that omits missing responses, and different variables have different
missing observations, the sample included for each modeling run will
vary. I recommend filtering the data to include *only* the sample that
does not have any missing data on the variables included.

In my view, `domin` is a decision tool for model explanation (i.e.,
understanding a fit, selected model) and is not as useful for model
selection or inference. Hence, if the user
