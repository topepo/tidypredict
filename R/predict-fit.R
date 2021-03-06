#' Returns a Tidy Eval formula to calculate fitted values
#'
#' It parses a model or uses an already parsed model to return a
#' Tidy Eval formula that can then be used inside a dplyr command.
#'
#' @param model An R model or a list with a parsed model. 
#'
#' @examples
#'
#' model <- lm(mpg ~ wt + cyl * disp, offset = am, data = mtcars)
#' tidypredict_fit(model)
#'
#' @export
tidypredict_fit <- function(model) {
  UseMethod("tidypredict_fit")
}

#' @export
tidypredict_fit.list <- function(model) {
  mt <- model$general$model
  fit <- NULL
  if(mt == "lm" | mt == "glm" | mt == "earth") 
    fit <- build_fit_formula(model)
  if(mt == "randomForest" | mt == "ranger") 
    fit <- build_fit_formula_rf(model)
  if(is.null(fit)) stop("Model type not supported")
  fit
}