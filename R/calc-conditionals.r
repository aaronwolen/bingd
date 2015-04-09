#' Calculate conditional probabilities for annotated features
#' 
#' @param object \code{AnnotatedGWAS} object
#' @param risk.thresh the p-value threshold applied to identify ``risky'' markers
#' @param adjust optional parameter to adjust baseline and conditional
#' probabilities for linkage disequilibrium
#' @inheritParams calc.priors
#' 
#' @exportMethod calc.conditionals
#' 
#' @return data.frame containing probabilities of observing each combination of
#' annotated features.

setGeneric("calc.conditionals", 
  function(object, risk.thresh = NULL, adjust = NULL, verbose = FALSE) {
  standardGeneric("calc.conditionals")
})

setMethod("calc.conditionals", "AnnotatedGWAS", 
  function(object, risk.thresh = NULL, adjust = NULL, verbose = FALSE) {
  
  if (!is.consolidated(object)) 
    stop("Features must be consolidated with consolidate().", call. = FALSE)
  
  if (is.null(risk.thresh)) risk.thresh <- calc.threshold(pvalue(object), verbose)
  
  labels <- names(fcols(object))
  vars <- DataFrame(conditional = pvalue(object) < risk.thresh, fcols(object))
  
  # Baseline probabilities
  base <- transform(plyr::count(vars[labels]), 
                    prob = freq /sum(freq))
  
  # Probability given conditional variable
  risk <- transform(plyr::count(subset(vars, conditional)[labels]),
                    prob = freq / sum(freq))
  
  # Impute combined probabilities for unobserved combinations 
  if (nrow(base) > nrow(risk)) risk <- impute.conditionals(base, risk, labels)
  
  probs <- merge(base, risk, by = labels, suffixes = c(".base", ".risk"))
  
  # Calculate enrichment
  if (is.null(adjust)) {
    probs$prob <- probs$prob.risk - probs$prob.base 
  } else {
    probs$prob <- probs$prob.risk * adjust - (adjust - 1) * probs$prob.base  
  }
  
  # Label variable combinations
  out <- data.frame(label = label.groups(probs[labels]), probs)
  out <- structure(out, class = c("gwas.conditionals", class(out)))
  
  if (verbose) report(out, "Conditional probabilities")
  return(out)
})


# Impute combined conditional probabilities
impute.conditionals <- function(base, risk, labels) {
  
  # Identify variable combinations absent in risk
  base.ids <- apply(base[labels], 1, paste, collapse = ".")
  risk.ids <- apply(risk[labels], 1, paste, collapse = ".")
  null.combos <- which(!base.ids %in% risk.ids)
  
  warning("Probabilties were estimated for the variable combinations listed below because none were present in the provided data", call. = FALSE)
  print(base[null.combos, labels])
  
  imputed <- as.numeric()
  for (i in seq_along(null.combos)) {
    new.row <- as.list(base[null.combos[i], labels])
    sub.exp <- paste(paste(names(new.row), new.row, sep = "=="), collapse = "|")
    var.probs <- subset(risk, eval(parse(text = sub.exp)))$prob
    new.prob <- c(do.call("prod", as.list(var.probs)) * length(var.probs))
    risk <- rbind(risk, c(new.row, freq = 0, prob = new.prob))
  }
  return(risk)
}
