#----Coordinate object----
#' The Coordinate Class
#'
#' The Coordinate object stores Coordinate data of cornplot.
#'
#' @slot coordinate A data frame contains the x, y, and spot name of geoseq.
#' @slot text The location name of geoseq.
#' @slot text_height The y axis of text on cornplot.
#' @slot tick_start_of_layers The start site of y tick.
#' @slot x The range of x in cornplot.
#' @slot y The range of y in cornplot.
#' @slot title_height The height of title in cornplot.
#' @slot legend_position The legend position in cornplot.
#' @slot shift The shift of y tick.
#' @slot size The dot size of spot.
#' @slot text_size The text size of cornplot.
#'
#' @name Coordinate-class
#' @rdname Coordinate-class
#' @exportClass Coordinate
#'
Coordinate = setClass(
  Class = 'Coordinate',
  slots = c(
    coordinate = 'data.frame',
    text = 'character',
    text_height = 'numeric',
    tick_start_of_layers = 'numeric',
    x = 'numeric',
    y = 'numeric',
    title_height = 'numeric',
    legend_position = 'numeric',
    shift = 'numeric',
    size = 'numeric',
    text_size = 'numeric'
  )
)

#' Create a Coordinate object
#'
#' Create a Coordinate object contain the basic cornplot information.
#'
#' @param coordinate A data frame contains the x, y, and spot name of geoseq.
#' @param text The location name of geoseq.
#' @param text_height The y axis of text on cornplot.
#' @param tick_start_of_layers The start site of y tick.
#' @param x The range of x in cornplot.
#' @param y The range of y in cornplot.
#' @param title_height The height of title in cornplot.
#' @param legend_position The legend position in cornplot.
#' @param shift The shift of y tick.
#' @param size The dot size of spot.
#' @param text_size The text size of cornplot.
#'
#' @return A Coordinate object
#'
#' @export
#'
CreateCoordinateObject = function(
    coordinate,
    text,
    text_height,
    tick_start_of_layers,
    x,
    y,
    title_height,
    legend_position,
    shift,
    size,
    text_size
){

  object = new(
    Class = 'Coordinate',
    coordinate = coordinate,
    text = text,
    text_height = text_height,
    tick_start_of_layers = tick_start_of_layers,
    x = x,
    y = y,
    title_height = title_height,
    legend_position = legend_position,
    shift = shift,
    size = size,
    text_size = text_size
  )

  return(object)
}
#----Coordinate object----

#----GEOseq object----
#' The GEOseq Class
#'
#' The GEOseq object stores GEOseq data.
#'
#' @slot exp_mat A data frame contains the log normalized GEO-seq data.
#'
#' @name GEOseq-class
#' @rdname GEOseq-class
#' @exportClass GEOseq
#'
GEOseq = setClass(
  Class = 'GEOseq',
  slots = c(
    exp_mat = 'ANY'
  )
)

#' Create a GEOseq object
#'
#' Create a GEOseq object contain the expression of GEO-seq data.
#'
#' @param exp_mat A data frame contains the log normalized GEO-seq data.
#'
#' @return A GEOseq object
#'
#' @export
#'
CreateGEOseqObject = function(
    exp_mat
){

  object = new(
    Class = 'GEOseq',
    exp_mat = exp_mat
  )

  return(object)
}
#----GEOseq object----

#----snRNASpatial object----
#' The snRNASpatial Class
#'
#' The snRNASpatial object stores Tangram mapped snRNASpatial data.
#'
#' @slot exp_mat A data frame contains the Tangram mapped snRNASpatial data.
#'
#' @name snRNASpatial-class
#' @rdname snRNASpatial-class
#' @exportClass snRNASpatial
#'
snRNASpatial = setClass(
  Class = 'snRNASpatial',
  slots = c(
    exp_mat = 'ANY'
  )
)

#' Create a snRNASpatial object
#'
#' Create a snRNASpatial object contain Tangram mapped snRNASpatial data.
#'
#' @param exp_mat A data frame contains Tangram mapped snRNASpatial data.
#'
#' @return A snRNASpatial object
#'
#' @export
#'
CreatesnRNASpatialObject = function(
    exp_mat
){

  object = new(
    Class = 'snRNASpatial',
    exp_mat = exp_mat
  )

  return(object)
}
#----snRNASpatial object----

#----snATACSpatial object----
#' The snATACSpatial Class
#'
#' The snATACSpatial object stores Tangram mapped snATACSpatial data.
#'
#' @slot exp_mat A data frame contains the Tangram mapped snATACSpatial data.
#' @slot exp_mat2 Extra snATACSpatial data.
#'
#' @name snATACSpatial-class
#' @rdname snATACSpatial-class
#' @exportClass snATACSpatial
#'
snATACSpatial = setClass(
  Class = 'snATACSpatial',
  slots = c(
    exp_mat = 'ANY',
    exp_mat2 = 'ANY'
  )
)

#' Create a snATACSpatial object
#'
#' Create a snATACSpatial object contain Tangram mapped snATACSpatial data.
#'
#' @param exp_mat A data frame contains Tangram mapped snATACSpatial data.
#' @param exp_mat2 Extra snATACSpatial data.
#'
#' @return A snATACSpatial object
#'
#' @export
#'
CreatesnATACSpatialObject = function(
    exp_mat,
    exp_mat2
){

  object = new(
    Class = 'snATACSpatial',
    exp_mat = exp_mat,
    exp_mat2 = exp_mat2
  )

  return(object)
}
#----snATACSpatial object----

#----CelltypeSpatial object----
#' The CelltypeSpatial Class
#'
#' The CelltypeSpatial object stores Tangram mapped CelltypeSpatial data.
#'
#' @slot exp_mat A data frame contains the Tangram mapped CelltypeSpatial data.
#'
#' @name CelltypeSpatial-class
#' @rdname CelltypeSpatial-class
#' @exportClass CelltypeSpatial
#'
CelltypeSpatial = setClass(
  Class = 'CelltypeSpatial',
  slots = c(
    exp_mat = 'ANY'
  )
)

#' Create a CelltypeSpatial object
#'
#' Create a CelltypeSpatial object contain Tangram mapped CelltypeSpatial data.
#'
#' @param exp_mat A data frame contains Tangram mapped CelltypeSpatial data.
#'
#' @return A CelltypeSpatial object
#'
#' @importFrom methods new
#'
#' @export
#'
CreateCelltypeSpatialObject = function(
    exp_mat
){

  object = new(
    Class = 'CelltypeSpatial',
    exp_mat = exp_mat
  )

  return(object)
}
#----CelltypeSpatial object----

#----Spatial object----
#' The spatial Class
#'
#' The spatial object is a wrapper of multi-omic gastrulation data.
#'
#' @slot Coordinate Coordinate object.
#' @slot GEOseq GEOseq object.
#' @slot snRNASpatial snRNA spatial object.
#' @slot snATACSpatial snATAC spatial object.
#' @slot CelltypeSpatial Celltype spatial object.
#' @slot stage The stage of embryo data (E6.5-E7.5).
#'
#' @name spatial-class
#' @rdname spatial-class
#' @exportClass spatial
#'
spatial = setClass(
  Class = 'spatial',
  slots = c(
    Coordinate = 'Coordinate',
    GEOseq = 'GEOseq',
    snRNASpatial = 'snRNASpatial',
    snATACSpatial = 'snATACSpatial',
    CelltypeSpatial = 'CelltypeSpatial',
    stage = 'character'
  )
)

#' Create a spatial object
#'
#' Create a spatial object from Coordinate_object, GEOseq_object, snRNASpatial_object, snATACSpatial_object and CelltypeSpatial_object.
#'
#' @param Coordinate_object Coordinate object.
#' @param GEOseq_object GEOseq object.
#' @param snRNASpatial_object snRNA spatial object.
#' @param snATACSpatial_object snATAC spatial object.
#' @param CelltypeSpatial_object celltype spatial object
#' @param stage The stage of embryo data (E6.5-E7.5).
#'
#' @return A spatial object
#'
#' @export
#'
#'
#'
CreateSpatialObject = function(Coordinate_object,
                               GEOseq_object,
                               snRNASpatial_object,
                               snATACSpatial_object,
                               CelltypeSpatial_object,
                               stage){

  object = new(
    Class = 'spatial',
    Coordinate = Coordinate_object,
    GEOseq = GEOseq_object,
    snRNASpatial = snRNASpatial_object,
    snATACSpatial = snATACSpatial_object,
    CelltypeSpatial = CelltypeSpatial_object,
    stage = stage

  )

  return(object)

}
#----Spatial object----

#----show method----
setMethod(f = 'show',
          signature = 'Coordinate',
          definition = function(object){
            cat('A Coordinate object contains the coordinate of cornplot\n')
          })
setMethod(f = 'show',
          signature = 'GEOseq',
          definition = function(object){
            cat('A GEOseq object contains a', nrow(object@exp_mat), '*', ncol(object@exp_mat), 'gene-spot matrix\n')
          })
setMethod(f = 'show',
          signature = 'snRNASpatial',
          definition = function(object){
            cat('A snRNASpatial object contains a', nrow(object@exp_mat), '*', ncol(object@exp_mat), 'gene-spot matrix\n')
          })
setMethod(f = 'show',
          signature = 'snATACSpatial',
          definition = function(object){
            cat('snATACSpatial slot contains a', nrow(object@exp_mat), '*', ncol(object@exp_mat), 'Signac peak-spot matrix, and a', nrow(object@exp_mat2), '*', ncol(object@exp_mat2), 'Scenicplus peak-spot matrix\n')
          })
setMethod(f = 'show',
          signature = 'CelltypeSpatial',
          definition = function(object){
            cat('CelltypeSpatial slot contains a', nrow(object@exp_mat), '*', ncol(object@exp_mat), 'Celltype-spot matrix\n')
          })
setMethod(f = 'show',
          signature = 'spatial',
          definition = function(object){
            cat('A spatial object used for cornplot\n',
                'GEOseq slot contains a', nrow(object@GEOseq@exp_mat), '*', ncol(object@GEOseq@exp_mat), 'gene-spot matrix\n',
                'snRNASpatial slot contains a', nrow(object@snRNASpatial@exp_mat), '*', ncol(object@snRNASpatial@exp_mat), 'gene-spot matrix\n',
                'snATACSpatial slot contains a', nrow(object@snATACSpatial@exp_mat), '*', ncol(object@snATACSpatial@exp_mat), 'Signac peak-spot matrix, and a', nrow(object@snATACSpatial@exp_mat2), '*', ncol(object@snATACSpatial@exp_mat2), 'Scenicplus peak-spot matrix\n',
                'CelltypeSpatial slot contains a', nrow(object@CelltypeSpatial@exp_mat), '*', ncol(object@CelltypeSpatial@exp_mat), 'celltype-spot matrix\n')
          })

#----show method----



