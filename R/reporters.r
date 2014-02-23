# Report intermediate statistics

report <- function(x, header, digits = 2) UseMethod("report")

report.character <- function(x, header, digits = 2) {
  print.header(header)
  message(x, "\n")
}

report.gwas.priors <- function(x, header, digits = 2) {
  print.header(header)  
  
  s.names <- paste0(format(names(x), justify = "right"), ":")
  s.out <- lapply(x, format, digits = digits)

  message(Map(paste, s.names, s.out, "\n"))
}


report.gwas.conditionals <- function(x, header, digits = 2) {
  print.header(header)
  
  out <- data.frame(x[, c("freq.base", "prob.base", "freq.risk", "prob.risk")])
  
  f.prob <- function(x) format(x, digits = digits)
  f.freq <- function(x) format(x, big.mark = ",")
  
  out <- transform(out,
                   freq.base = f.freq(freq.base),
                   prob.base = f.prob(prob.base),
                   freq.risk = f.freq(freq.risk),
                   prob.risk = f.prob(prob.risk))
  
  out <- mapply(function(o, n) format(c(n, o)),
                out, colnames(out), SIMPLIFY = FALSE)
  
  out <- do.call("paste", c(as.list(out), sep = " | "))
  
  categories <- gsub("\\.", " ", x$label)
  categories <- format(c("Category", categories), justify = "right")
  out <- paste(categories, out, sep = " | ")
  
  message(paste(out, collapse = "\n"))
}


print.header <- function(x) message("\n", x, "\n", rep("-", nchar(x)))