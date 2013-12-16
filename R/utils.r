
#' Fix a GRanges object's chromosome information
#' 
#' Edit a GRanges object's chromosome information by adding chromosome lengths 
#' from UCSC and properly ordering the seqnames factor levels. Adding lengths
#' requires an internet connection. If one is not available then
#' \code{add.length} should be set to \code{FALSE}.
#' 
#' @param object object of class \code{GRanges}
#' @param add.length Add missing chromosome length information
#' @param db Optional bioconductor organism package (e.g., org.Hs.eg.db) to
#'   provide chromosome lengths so that an internet connection isn't necessary
#'   
#' @importFrom GenomicFeatures getChromInfoFromUCSC
#' @export

fix_chrs <- function(object, genome = "hg19", add.length = TRUE, db = NA) {

  stopifnot(class(object) == "GRanges")
  
  # Add chr prefix if missing
  seqlevels(object) <- paste0("chr", sub("^chr", "", seqlevels(object)))
  
  # Add chr lengths if a org.db file is provided
  if (!missing(db)) add.length <- TRUE
  
  if (add.length) {
    if (is.na(db)) {
      # Obtain chr lengths from UCSC  
      genome <- match.arg(genome, choices = ucscGenomes()[, "db"])
      chrs <- suppressMessages(getChromInfoFromUCSC(genome))
    } else {
      require(db, character.only = TRUE, quietly = TRUE)
      chrs <- eval(parse(text = gsub("\\.db", "CHRLENGTHS", db)))
      chrs <- data.frame(chrom = paste0("chr", names(chrs)), length = chrs)
    }
  } else {
    chrs <- data.frame(chrom = seqlevels(object), length = seqlengths(object))    
  }
  
  chrs <- subset(chrs, chrom %in% seqlevels(object))
  chrs <- chrs[match(order_chrs(chrs$chrom), chrs$chrom),]
  rownames(chrs) <- 1:nrow(chrs)
  
  seqlevels(object) <- as.character(chrs$chrom)
  seqlengths(object) <- chrs$length
  genome(object) <- genome
  
  return(sort(object))
}


#' Properly order a vector of chromosomes
#' 
#' @param chrs unordered character vetor of chromosome labels

order_chrs <- function(chrs) {

  stopifnot(is.null(dim(chrs)))
  chrs <- as.character(chrs)
  
  prefix <- ifelse(any(grepl("chr", chrs)), TRUE, FALSE)
  chrs <- sub("chr", "", chrs)
  
  nums <- grepl("[0-9]", chrs)
  alphas <- grepl("[A-Z]", chrs, ignore.case = TRUE)
  
  if (sum(c(nums, alphas)) != length(chrs))  {
    stop("This is weird. One or more of your chromosomes could not be\n",
         "classified as a numbered or lettered chromosome\n.", call. = FALSE)
  }
  
  nums <- sort(as.numeric(chrs[nums]))
  alphas <- toupper(chrs[alphas])
  alphas <- alphas[order(match(alphas, c("X", "Y", "M")))]
  
  out <- c(nums, alphas)
  if (prefix) out <- paste0("chr", out)
  
  return(out)
}
