# Generate labels based on column names of logical variables in a data.frame
# 
# example:
# label.groups(data.frame(a = c(T,F, T), b = c(F, T, T)))

label.groups <- function(df, sep = ".and.", empty = "neither") {
  
  if (!all(sapply(df, is.logical))) 
    stop("All columns in df must be of type logical.")
  
  # Construct groups
  groups <- do.call("expand.grid", rep(list(c(TRUE, FALSE)), ncol(df)))
  colnames(groups) <- colnames(df)

  # Group labels
  labels <- sapply(names(groups), function(x) ifelse(groups[[x]], x, NA))
  labels <- replace(labels, is.na(labels), "")
  labels <- do.call("paste", as.data.frame(labels, stringsAsFactors = FALSE))
  labels <- gsub("\\s{2,}", " ", labels)
  labels <- gsub("^\\s+|\\s+$", "", labels)
  labels[labels == ""] <- empty
  labels <- gsub(" ", sep, labels, fixed = TRUE)
  
  # Apply group labels
  row.ids   <- do.call("paste0", as.list(df))
  group.ids <- labels
  names(group.ids) <- do.call("paste0", as.list(groups))
  
  group.labels <- structure(group.ids[row.ids], names = NULL)
  
  return(group.labels)
}
