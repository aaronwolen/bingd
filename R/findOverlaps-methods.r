#' Find features overlapping GWAS markers
#' 
#' Determine which GWAS markers overlap each \code{FeatureList} feature.
#'  
#' @return 
#' \code{findOverlaps} returns either a \code{\link[IRanges]{HitsList-class}}
#' object when \code{select="all"} (the default), or a 
#' \code{\link[IRanges]{CompressedIntegerList}} when \code{select} is not
#' \code{"all"}.
#' 
#' @param query \code{GWAS} object
#' @param subject \code{FeatureList} object
#' @inheritParams IRanges::findOverlaps
#' @inheritParams GenomicRanges::findOverlaps
#' @inheritParams featureOverlaps
#' 
#' @rdname findOverlaps

setMethod("findOverlaps", c(query = "GWAS", subject = "FeatureList"),
  function(query, subject, maxgap = 0L, minoverlap = 1L,
           type = c("any", "start", "end", "within"),
           select = c("all", "first", "last", "arbitrary"), 
           ignore.strand = FALSE) {
    
    type <- match.arg(type)
    select <- match.arg(select)

    f.index <- stack(LocalPath(subject))
    f.paths <- structure(f.index$values, names = rownames(f.index))
  
    chrs <- GenomeInfoDb::seqlevels(query)
    
    result <- mclapply(f.paths, function(p) 
                       findOverlaps(query = query, 
                                    subject = load.feature(p, chrs, TRUE),
                                    maxgap = maxgap, minoverlap = minoverlap,
                                    type = type, select = select, 
                                    ignore.strand = ignore.strand),
                       mc.cores = foreach::getDoParWorkers())
    
    if (select == "all") {
      result <- as(SimpleList(result), "HitsList")  
    } else {
      result <- IntegerList(result)
    }
    
    return(result)
})
