#' Verify the presence and class of a data.frame variable
#' 
#' @param object \code{data.frame} or a \code{\link[IRanges]{DataFrame}}
#' @param description description of class being tested to print in error
#' @param var.name column or slot name to check for
#' @param var.class expected class of column or slot

.validColClass <- function(object, description, var.name, var.class) {
  if (!var.name %in% colnames(object)) 
    stop(description, " must contain ", var.name, call. = FALSE)
  
  if (var.class != class(object[[var.name]]))
    stop(var.name, " in ", description, " should be ", var.class, call. = FALSE)
}
