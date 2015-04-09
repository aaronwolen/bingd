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
  
  out.cols <- c("n.base", "prob.base", "n.risk", "prob.risk", "prob")
  out <- data.frame(x[, out.cols])
  
  f.prob <- function(x) format(x, digits = digits)
  f.freq <- function(x) format(x, big.mark = ",")
  
  out <- transform(out,
                   n.base = f.freq(n.base),
                   prob.base = f.prob(prob.base),
                   n.risk = f.freq(n.risk),
                   prob.risk = f.prob(prob.risk),
                        prob = f.prob(prob))
  
  out <- mapply(function(o, n) format(c(n, o)),
                out, colnames(out), SIMPLIFY = FALSE)
  
  out <- do.call("paste", c(as.list(out), sep = " | "))
  
  categories <- gsub("\\.", " ", x$label)
  categories <- format(c("Category", categories), justify = "right")
  out <- paste(categories, out, sep = " | ")
  
  message(paste(out, collapse = "\n"))
}


print.header <- function(x) message("\n", x, "\n", rep("-", nchar(x)))