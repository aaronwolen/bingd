# Load a specific feature from an AnnotationHub object
#
# object: AnnotationHub object
# which: numeric index or feature label
# seqlevels: limit the seqlevels for the returned object
# GNCList: should feature be converted into a GNCList
load_feature <- function(object, which, seqlevels = NULL, GNCList = FALSE) {
  
  if (!is(object, "AnnotationHub")) stop("object must be an AnnotatoinHub object.")
  if (missing(which)) stop("Must provide index or name of AnnotationHub record.")
  
  obj <- object[[which]]
  if (!is(obj, "GRanges")) stop("Features must be GRanges objects.")
  
  if (!is.null(seqlevels)) obj <- GenomeInfoDb::keepSeqlevels(obj, seqlevels)
  if (GNCList) obj <- GenomicRanges::GNCList(obj)
  obj
}
