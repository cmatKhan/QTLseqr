#' @rdname BSAExperiment-methods
#' @export
setGeneric("comparisons", function(object, comparisons)
  standardGeneric("comparisons")
)

#' @rdname BSAExperiment-methods
#' @export
setGeneric("comparisons<-", function(object, ..., value)
  standardGeneric("comparisons<-")
)

#' @rdname BSAExperiment-methods
#' @export
setGeneric(name = "createComparisonsFrame",
           signature = "object",
           def = function(object,
                          grouping_variable,
                          population_variable,
                          population_base_condition,
                          var1_name,
                          var2_name,
                          base_cond_in_each_group) {
             standardGeneric("createComparisonsFrame")
           })

#' @rdname BSAExperiment-methods
#' @export
setGeneric(name = "populationDepths",
           signature = "object",
           def = function(object,
                          population_1_sample,
                          population_2_sample
                          ) {
             standardGeneric("populationDepths")
           })


#' @rdname BSAExperiment-methods
#' @export
setGeneric(name = "deltaAltFrequency",
           signature = "object",
           def = function(object,
                          population_1_sample,
                          population_2_sample
                          ) {
             standardGeneric("deltaAltFrequency")
           })

#' @rdname BSAResults-method
