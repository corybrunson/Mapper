#' @title Interval Dual Tessellation Cover
#'
#' @docType class
#' @author Jason Cory Brunson
#' @description [R6] class representing an interval dual tessellation cover.
#' @family cover

#' @details
#'
#' The regular interval tessellation \eqn{T} of \eqn{Z} (the real line, or a
#' segment of it) consists of intervals of fixed length, shared endpoints, and
#' empty overlap. In order to completely partition \eqn{Z}, each interval is
#' considered closed at one end and open at the other. (**This remains to be
#' carefully implemented.**) The dual tessellation \eqn{T'} of \eqn{T} is
#' obtained by taking the endpoints of the intervals at the centers of the
#' intervals of \eqn{T}. The union \eqn{C} of \eqn{T} and \eqn{T'} is an
#' _interval dual tessellation cover_.
#'
#' An interval dual tessellation cover \eqn{C} is a special case of a
#' `[FixedIntervalCover]`, for a 1-dimensional lens and with `percent_overlap =
#' 50`. \eqn{C} is _gomic_ (a **g**eneric, **o**pen, **m**inimal, **i**nterval
#' **c**over) in the sense of Carriere and Oudot (2016). \eqn{C} has constant
#' "depth", in the sense that every point in \eqn{Z} is a member of exactly two
#' sets in \eqn{C}. As discussed by Carriere and Oudot, because the mapper
#' construction collapses sets in the pullback cover to vertices of the nerve
#' complex, only 1-dimensional features of a point cloud \eqn{X} can be
#' recovered from a 1-dimensional lens \eqn{f:X->Z}. A cover of depth 2 is
#' sufficient to recover these features, so no tessellation covers of greater
#' depth are available. Every non-empty intersection of cover sets has the same
#' length.
#'
#' `IntervalDualTessellationCover` takes the center of \eqn{f(X)} to be the
#' center of the overlap between a pair of intervals in \eqn{T} and in \eqn{T'}.
#' Thus, the tessellations have negation symmetry and the same number of sets.
#' When `width` is sufficiently large (twice the range of the data), the total
#' number of intervals is minimized at 2.
#' 

#' @field width A positive number. The width of each interval in the cover, as a
#'   fraction of the range of the lensed data.
#' @export
IntervalDualTessellationCover <- R6::R6Class(
  classname = "IntervalDualTessellationCover",
  inherit = CoverRef,
  private = list(.width = NULL)
)

#' @export
IntervalDualTessellationCover$set(
  which = "public", name = "initialize",
  value = function(...) {
    super$initialize(typename = "Interval Dual Tessellation")
    params <- list(...)
    if ("width" %in% names(params)) {
      self$width <- params[["width"]]
    }
  }
)

## Active binding to set `width`
## width ----
IntervalDualTessellationCover$set(
  which = "active", name = "width",
  value = function(value) {
    if (missing(value)) { private$.width } else {
      if (length(value) != 1) {
        stop("A tessellation cover takes only one `width` value.")
      }
      if (value <= 0) {
        stop("The width of a tessellation must be positive.")
      }
      private$.width <- value
      self
    }
  }
)

## Validate parameter settings
## validate ----
IntervalDualTessellationCover$set(
  which = "public", name = "validate",
  value = function(filter) {
    stopifnot(! is.null(private$.width))
    stopifnot(length(self$width) == 1)
    f_dim <- ncol(filter())
    if (f_dim != 1) {
      stop("A linear tessellation requires a 1-dimensional lens.")
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
      "%s Cover: (width = %s)",
      titlecase(private$.typename), self$width
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
    f_ran <- range(filter())
    f_len <- diff(f_ran)
    f_cent <- f_ran[1] + f_len/2
    
    ## Calculate hyper-parameters
    i_len <- self$width*f_len
    i_num <- c(ceiling(.5*f_len/i_len - .75), ceiling(.5*f_len/i_len - .25))
    
    ## If no index is given, construct the entire cover
    if (missing(index) || is.null(index)) {
      
      ## primary tessellation endpoints
      a_pts <- f_cent + (seq(-i_num[1] - 1, i_num[2]) + .25)*i_len
      ## secondary tessellation endpoints
      b_pts <- f_cent + (seq(-i_num[2], i_num[1] + 1) - .25)*i_len
      ## dual tessellation bounds
      ls_bounds <- cbind(
        c(a_pts[-length(a_pts)], b_pts[-length(b_pts)]),
        c(a_pts[-1], b_pts[-1])
      )
      
    } else {
      stopifnot(index %in% self$index_set)
      
      s_dual <- as.logical(match(substr(index, 2L, 2L), c("A", "B")) - 1L)
      s_order <- as.integer(sub("^\\([AB]: ([0-9\\-]+)\\)$", "\\1", index))
      ls_bounds <- f_cent + if (s_dual) {
        (s_order + c(0, 1) - .25)*i_len
      } else {
        (s_order + c(-1, 0) + .25)*i_len
      }
      
    }
    
    ## Return bounds
    return(ls_bounds)
  }
)

## Setup dual index sets
## construct_index_set ----
IntervalDualTessellationCover$set(
  which = "public", name = "construct_index_set",
  value = function(filter) {
    if (is.function(filter)) {
      self$validate(filter)
      fv <- filter()
    } else {
      fv <- filter
    }
    
    ## Get filter length
    f_len <- diff(range(fv))
    
    ## Calculate hyper-parameters
    i_len <- self$width*f_len
    i_num <- c(ceiling(.5*f_len/i_len - .75), ceiling(.5*f_len/i_len - .25))
    
    ## Primary and dual indices have same counts
    a_ind <- paste0("(A: ", seq(-i_num[1], i_num[2]), ")", sep = "")
    b_ind <- paste0("(B: ", seq(-i_num[2], i_num[1]), ")", sep = "")
    
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
    
    ## Get filter
    fv <- filter()
    
    ## If the index set hasn't been made yet, construct it
    self$construct_index_set(fv)
    
    ## If no index specified, return the level sets either by construction
    if (missing(index) || is.null(index)) {
      #stopifnot(! index %in% self$index_set)
      
      set_bnds <- self$interval_bounds(filter)
      
      self$level_sets <- constructIsoAlignedLevelSets(fv, as.matrix(set_bnds))
      return(invisible(self))
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
    
    if (k == 1L) {
      
      ## Get filter length
      f_len <- diff(range(filter()))
      
      ## Calculate hyper-parameters
      i_len <- self$width*f_len
      i_num <- c(ceiling(.5*f_len/i_len - .75), ceiling(.5*f_len/i_len - .25))
      i_tot <- 1 + sum(i_num)
      ab <- as.integer(i_num[1] >= i_num[2])
      
      ## Each set only overlaps with at most two sets, all in the dual layer
      return(cbind(
        rep(self$index_set[seq(i_tot)], c(2-ab, rep(2, i_tot-2), 1+ab)),
        rep(self$index_set[i_tot + seq(i_tot)], c(1+ab, rep(2, i_tot-2), 2-ab))
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
