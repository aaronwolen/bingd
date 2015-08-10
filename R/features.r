# Load a specific feature from an AnnotationHub object
#
# object: AnnotationHub object
# which: numeric index or feature label
# seqlevels: limit the seqlevels for the returned object
# GNCList: should feature be converted into a GNCList
load_feature <- function(object, which, seqlevels = NULL, GNCList = FALSE) {
  
  if (!is(object, "AnnotationHub")) stop("object must be an AnnotatoinHub object.")
  if (missing(which)) stop("Must provide index or name of AnnotationHub record.")
  
  title <- AnnotationHub::mcols(object[which])$title
  if (grepl("DNase.hotspot.fdr0.01", title, fixed = TRUE)) {

    if (grepl("01.peaks.bed", title, fixed = TRUE)) {
      extra.cols <- c(name="character", 
                      signalValue = "numeric", 
                      pValue = "numeric")
    } else if (grepl("01.broad.bed", title, fixed = TRUE)) {
      extra.cols <- c(name="character", 
                      pValue = "numeric")
    } else {
      stop("Encountered unanticipated RoadMap DNase BED file format")
    }
    
    er <- AnnotationHub::cache(object[which])
    obj <- rtracklayer::import(er, format="bed", extraCols = extra.cols)
  } else {
    obj <- object[[which]]  
  }
  
  
  if (!is(obj, "GRanges")) stop("Features must be GRanges objects.")
  
  if (!is.null(seqlevels)) obj <- GenomeInfoDb::keepSeqlevels(obj, seqlevels)
  if (GNCList) obj <- GenomicRanges::GNCList(obj)
  obj
}
