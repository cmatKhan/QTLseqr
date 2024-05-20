# comparisons ----

#' @rdname comparisons
#' @family BSAExperiment-methods
#' @usage NULL
#' @export
setGeneric("comparisons<-", function(object, ..., value)
  standardGeneric("comparisons<-")
)

#' @rdname comparisons
#' @usage NULL
#' @export
setGeneric("comparisons", function(object, comparisons)
  standardGeneric("comparisons")
)

# totalDepth ----

#' @rdname totalDepth
#' @family BSAExperiment-methods
#' @usage NULL
#' @export
setGeneric("totalDepth<-", function(object, ..., value)
  standardGeneric("totalDepth<-")
)

#' @rdname totalDepth
#' @usage NULL
#' @export
setGeneric("totalDepth", function(object, sample_set = NULL)
  standardGeneric("totalDepth")
)

# altDepth ----

#' @rdname altDepth
#' @family BSAExperiment-methods
#' @usage NULL
#' @export
setGeneric("altDepth<-", function(object, ..., value)
  standardGeneric("altDepth<-")
)

#' @rdname altDepth
#' @usage NULL
#' @export
setGeneric("altDepth", function(object, sample_set = NULL)
  standardGeneric("altDepth")
)

# refDepth ----

#' @rdname refDepth
#' @family BSAExperiment-methods
#' @usage NULL
#' @export
setGeneric("refDepth", function(object, sample_set = NULL)
  standardGeneric("refDepth")
)

# variantFilter ----

#' @rdname variantFilter
#' @family BSAExperiment-methods
#' @usage NULL
#' @export
setGeneric("variantFilter<-", function(object, ..., value)
  standardGeneric("variantFilter<-")
)

#' @rdname variantFilter
#' @usage NULL
#' @export
setGeneric("variantFilter", function(object, value)
  standardGeneric("variantFilter")
)

# sampleFilter ----

#' @rdname sampleFilter
#' @family BSAExperiment-methods
#' @usage NULL
#' @export
setGeneric("sampleFilter<-", function(object, ..., value)
  standardGeneric("sampleFilter<-")
)

#' @rdname sampleFilter
#' @usage NULL
#' @export
setGeneric("sampleFilter", function(object, value)
  standardGeneric("sampleFilter")
)

# matrixFilter ----

#' @rdname matrixFilter
#' @family BSAExperiment-methods
#' @usage NULL
#' @export
setGeneric("matrixFilter<-", function(object, ..., value)
  standardGeneric("matrixFilter<-")
)

#' @rdname matrixFilter
#' @usage NULL
#' @export
setGeneric("matrixFilter", function(object, value)
  standardGeneric("matrixFilter")
)

# applyFilter ----

#' @rdname applyFilter
#' @family BSAExperiment-methods
#' @usage NULL
#' @export
setGeneric("applyFilter", function(object,
                                   variant_filter = TRUE,
                                   sample_filter = TRUE,
                                   matrix_filter = TRUE)
  standardGeneric("applyFilter")
)

# populationDepths ----

#' @rdname populationDepths
#' @family BSAExperiment-methods
#' @usage NULL
#' @export
setGeneric(name = "populationDepths",
           signature = "object",
           def = function(object,
                          population_1_sample,
                          population_2_sample
                          ) {
             standardGeneric("populationDepths")
           })

# snpIndex ----

#' @rdname snpIndex
#' @family BSAExperiment-methods
#' @usage NULL
#' @export
setGeneric(name = "snpIndex",
           signature = "object",
           def = function(object, sample_set = NULL) {
             standardGeneric("snpIndex")
           })

# deltaSnpIndex ----

#' @rdname deltaSnpIndex
#' @family BSAExperiment-methods
#' @usage NULL
#' @export
setGeneric(name = "deltaSnpIndex",
           signature = "object",
           def = function(object,
                          population_1_sample,
                          population_2_sample
                          ) {
             standardGeneric("deltaSnpIndex")
           })
