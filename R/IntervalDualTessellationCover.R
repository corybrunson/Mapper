#' @title Interval Dual Tessellation Cover
#'
#' @docType class
#' @author Jason Cory Brunson
#' @description \code{\link{R6}} class representing an interval dual
#'   tessellation cover.
#' @family cover

#' @details
#'
#' A (regular) interval tessellation \eqn{T} of \eqn{Z} (the real line, or of a
#' segment of it) consists of intervals of fixed length, shared endpoints, and
#' empty overlap. In order to completely partition \eqn{Z}, each interval is
#' considered closed at one end and open at the other. (**NB:** This has not
#' been carefully implemented.) The dual \eqn{T'} of \eqn{T} is obtained by
#' defining a new interval cover with endpoints at the centers of the intervals
#' of \eqn{T}. The union \eqn{C} of \eqn{T} and \eqn{T'} is an _interval dual
#' tessellation cover_.
#'
#' An interval dual tessellation cover \eqn{C} is a special case of a
#' \code{\link{FixedIntervalCover}}, for a 1-dimensional lens and with
#' \code{percent_overlap = 50}. \eqn{C} is _gomic_ (a **g**eneric, **o**pen,
#' **m**inimal, **i**nterval **c**over) in the sense of Carriere and Oudot
#' (2016).
#'
#' \eqn{C} has constant "depth", in the sense that every point in \eqn{Z} is a
#' member of exactly two sets in \eqn{C}. As discussed by Carriere and Oudot,
#' because the mapper construction collapses sets in the pullback cover to
#' vertices of the nerve complex, only 1-dimensional features of an original
#' point cloud \eqn{X} can be recovered from a one-dimensional lens
#' \eqn{f:X->Z}. A cover of depth 2 is sufficient to recover these features, so
#' no tessellation covers of greater depth are available.
#'
#' \code{IntervalDualTessellationCover} takes the extrema of \eqn{f(X)} to be
#' the centers of the first set in \eqn{T} and of the last set in \eqn{T} (if
#' \code{number_intervals} is even) or \eqn{T'} (if \code{number_intervals} is
#' odd). Thus \code{number_intervals} must be at least 2.
#' 

#' @field number_intervals An integer. The number of bins used to cover \eqn{Z};
#'   must be at least 2.
#' @export
IntervalDualTessellationCover <- R6::R6Class(
  classname = "IntervalDualTessellationCover",
  inherit = CoverRef,
  private = list(.number_intervals = NULL)
)

#' @export
IntervalDualTessellationCover$set(
  which = "public", name = "initialize",
  value = function(...) {
    super$initialize(typename = "Interval Dual Tessellation")
    params <- list(...)
    if ("number_intervals" %in% names(params)) {
      self$number_intervals <- params[["number_intervals"]]
    }
  }
)

## Active binding to set `number_intervals`
## number_intervals ----
IntervalDualTessellationCover$set(
  which = "active", name = "number_intervals",
  value = function(value) {
    if (missing(value)) { private$.number_intervals } else {
      if (length(value) != 1) {
        stop("An interval tessellation covers a 1-dimensional space.")
      }
      if (value < 2) {
        stop("An interval dual tessellation needs at least 2 sets.")
      }
      private$.number_intervals <- value
      self
    }
  }
)

## Validate parameter settings
## validate ----
IntervalDualTessellationCover$set(
  which = "public", name = "validate",
  value = function(filter) {
    stopifnot(! is.null(private$.number_intervals))
    stopifnot(length(self$number_intervals) == 1)
    stopifnot(self$number_intervals >= 2)
    f_dim <- ncol(filter())
    if (f_dim != 1) {
      stop("An interval tessellation requires a 1-dimensional lens.")
    }
  }
)

## format ----
IntervalDualTessellationCover$set(
  which = "public", name = "format",
  value = function(...) {
    titlecase <- function(x) {
      s <- strsplit(x, " ")[[1]]
      paste(
        toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " "
      )
    }
    sprintf(
      "%s Cover: (number intervals = [%s])",
      titlecase(private$.typename),
      self$number_intervals
    )
  }
)

## This function is specific to the interval-type covers
## interval_bounds ----
IntervalDualTessellationCover$set(
  which = "public", name = "interval_bounds",
  value = function(filter, index = NULL) {
    self$validate(filter)
    
    ## Get filter range and length
    fv <- filter()
    f_ran <- range(fv)
    f_len <- diff(f_ran)
    
    ## Calculate hyper-parameters
    i_len <- 2 / (self$number_intervals - 1) * f_len
    i_rad <- i_len/2
    n1 <- ceiling(self$number_intervals / 2)
    n2 <- floor(self$number_intervals / 2)
    
    ## If no index is given, construct the entire cover
    if (missing(index) || is.null(index)) {
      ## primary tessellation centers
      a_centers <- f_ran[1] + 0:(n1 - 1) * i_len
      ## dual tessellation centers
      b_centers <- f_ran[1] + i_rad + 0:(n2 - 1) * i_len
      ## bounds
      ls_centers <- c(a_centers, b_centers)
      ls_bounds <- cbind(ls_centers - i_rad, ls_centers + i_rad)
    } else {
      stopifnot(index %in% self$index_set)
      s_dual <- match(substr(index, 2L, 2L), c("A", "B")) - 1
      s_order <- as.integer(sub("^\\([AB]: ([0-9]+)\\)$", "\\1", index))
      s_center <- f_ran[1] + s_dual*i_rad + (s_order - 1)
      ls_bounds <- s_center + c(-i_rad, i_rad)
    }
    
    ## Return bounds
    return(ls_bounds)
  }
)

## Setup dual index sets
## construct_index_set ----
IntervalDualTessellationCover$set(
  which = "public", name = "construct_index_set",
  value = function(...) {
    
    ## Calculate hyper-parameters
    n1 <- ceiling(self$number_intervals / 2)
    n2 <- floor(self$number_intervals / 2)
    
    ## Primary, then dual, indices
    a_ind <- paste0("(A: ", 1:n1, ")", sep = "")
    b_ind <- paste0("(B: ", 1:n2, ")", sep = "")
    
    self$index_set <- c(a_ind, b_ind)
  }
)

## Given the current set of parameter values, construct the level sets whose
## union covers the filter space
## construct_cover ----
IntervalDualTessellationCover$set(
  which = "public", name = "construct_cover",
  value = function(filter, index = NULL) {
    stopifnot(is.function(filter))
    self$validate(filter)
    
    ## Get filter range and length
    fv <- filter()
    f_ran <- range(fv)
    f_len <- diff(f_ran)
    
    ## If the index set hasn't been made yet, construct it.
    self$construct_index_set()
    
    ## If no index specified, return the level sets either by construction
    if (missing(index) || is.null(index)) {
      stopifnot(! index %in% self$index_set)
      set_bnds <- self$interval_bounds(filter)
      self$level_sets <- constructIsoAlignedLevelSets(fv, as.matrix(set_bnds))
      return(invisible(self)) ## return invisibly 
    } else {
      if (! is.null(self$level_sets) && index %in% names(self$level_sets)) {
        return(self$level_sets[[index]])
      } else {
        p_idx <- which(index == self$index_set)
        set_bnds <- self$interval_bounds(filter, index)
        level_set <- constructIsoAlignedLevelSets(fv, set_bnds)
        return(level_set)
      }
    }
    
  }
)

## Constructs a 'neighborhood', which is an (n x k+1) subset of pullback ids
## representing the set of n unique (k+1)-fold intersections are required to
## construct the nerve.
## neighborhood ----
IntervalDualTessellationCover$set(
  which = "public", name = "neighborhood",
  value = function(filter, k) {
    stopifnot(!is.null(private$.index_set))
    
    ## Take advantage of bipartite overlap structure; could even allow buffers
    ## at endpoints and simply not calculate within-layer intersections
    
    if (k == 1L) {
      
      ## Calculate hyper-parameters
      n1 <- ceiling(self$number_intervals / 2)
      n2 <- floor(self$number_intervals / 2)
      odd <- as.integer(n2 < n1)
      
      ## Each set only (potentially) overlaps with two sets in the dual layer
      return(cbind(
        rep(self$index_set[seq(n1)], c(1, rep(2, n1 - 2), 2 - odd)),
        rep(self$index_set[seq(n1 + 1, n1 + n2)], c(rep(2, n2 - 1), 1 + odd))
      ))
      
    } else if (k > 1L) {
      
      ## The cover has two partitioning layers, so no three sets overlap
      return(matrix("", nrow = 0, ncol = k + 1L))
      
    } else {
      ## Default to the generic neighborhood calculator
      return(super$neighborhood(filter, k))
    }
    
  }
)
