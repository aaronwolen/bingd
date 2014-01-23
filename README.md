bingd
=====

Bayesian integration of genomic data.

Features
========

bingd assumes all features are peak ranges, represented as `GRanges` objects and saved as `RData` files. 

Feature metadata
----------------

[AnnotationHub] provides detailed metadata for each feature in its repository. bingd's `hub.metadata(online = TRUE)` function is used to retrieve the latest metadata when an internet connection is available. An extra column, `LocalPath`, is appended to the metadata, which provides the full path to a cached copy of each feature; uncached features are denoted by an `NA` in this column. If `online = FALSE` then `hub.metadata()` derived metadata from locally cached AnnotationHub features and provides only `Title` and `LocalPath`. 

Searching for features
----------------------

Genomic feature files can be searched for locally using ______ or in the AnnotationHub repository using `hub.search()`, which can be conducted online or offline.

Both functions return a list of `DataFrame`s that, at minimum, include `Title` and `LocalPath` columns.

<!-- links -->
[AnnotationHub]: http://master.bioconductor.org/packages/release/bioc/html/AnnotationHub.html
