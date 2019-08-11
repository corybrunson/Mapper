#' @title Square Dual Tesselation Cover
#'
#' @docType class
#' @author Jason Cory Brunson
#' @description [R6] class representing a square dual tessellation cover.
#' @family cover

#' @details
#'
#' A (regular) tessellation \eqn{T} of \eqn{Z} (the plane, or )
#'
#' A (regular) square tessellation \eqn{T} of \eqn{Z} (the real line, or a
#' path-connected region of it) consists of squares of fixed side length
#' arranged as tiles, with shared vertices and edges and empty overlap. In order
#' to completely partition \eqn{Z}, each square is considered closed along half
#' of its perimeter and open along the other half. (**NB:** This has not been
#' carefully implemented.) The dual \eqn{T'} of \eqn{T} is obtained by defining
#' a new square cover with vertices at the centroids of the squares of \eqn{T}
#' and edges perpendicular to those of \eqn{T}. The union \eqn{C} of \eqn{T} and
#' \eqn{T'} is a _square dual tessellation cover_.
#'
#' A square dual tessellation cover \eqn{C} is the union of two
#' `[FixedIntervalCover]`s, for a 2-dimensional lens and with the same uniform
#' `number_intervals` and  `percent_overlap = 0`, each offset from the other by
#' half of the side length along each dimension.
#'
#' \eqn{C} has constant "depth", in the sense that every point in \eqn{Z} is a
#' member of exactly two sets in \eqn{C}. A cover of depth 2 is only sufficient
#' to recover 1-dimensional features of a point cloud.
#'
#' `SquareDualTessellationCover` takes the centroid of \eqn{f(X)} to be the
#' centroid of the overlap between a pair of squares in \eqn{T} and in \eqn{T'}.
#' Thus, the tessellations have rotation symmetry and the same number of sets
#' (before dropping empty sets). When `width` is sufficiently large (twice the
#' range of the data), the total number of squares is minimized at 2.
#' 

#' @field width A positive number. The side length of each square in the cover,
#'   as a fraction of the larger coordinate range of the lensed data.
#' @export
SquareDualTessellationCover <- R6::R6Class(
  classname = "SquareDualTessellationCover",
  inherit = CoverRef,
  private = list(.width = NULL)
)

#' @export
SquareDualTessellationCover$set(
  which = "public", name = "initialize",
  value = function(...) {
    super$initialize(typename = "Square Dual Tessellation")
    params <- list(...)
    if ("width" %in% names(params)) {
      self$width <- params[["width"]]
    }
  }
)

SquareDualTessellationCover$set(
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

SquareDualTessellationCover$set(
  which = "public", name = "validate",
  value = function(filter) {
    stopifnot(! is.null(private$.width))
    stopifnot(length(self$width) == 1)
    f_dim <- ncol(filter())
    if (f_dim != 2) {
      stop("A square tessellation requires a 2-dimensional lens.")
    }
  }
)

SquareDualTessellationCover$set(
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

SquareDualTessellationCover$set(
  which = "public", name = "interval_bounds",
  value = function(filter, index = NULL) {
    self$validate(filter)
    
    ## Get filter range and length
    fv <- filter()
    f_ran <- apply(fv, 2, range)
    f_len <- diff(f_ran)
    f_cent <- f_ran[1, ] + f_len/2
    
    ## Calculate hyper-parameters
    i_len <- self$width*max(f_len)
    i_num <- rbind(
      ceiling(.5*f_len/i_len - .75),
      ceiling(.5*f_len/i_len - .25)
    )
    i_tot <- apply(i_num, 2, sum) + 1
    
    ## If no index is given, construct the entire cover
    if (missing(index) || is.null(index)) {
      
      ## primary tessellation coordinate bounds
      a_xbnd <- f_cent[1] + (seq(-i_num[1, 1] - 1, i_num[2, 1]) + .25)*i_len
      a_ybnd <- f_cent[2] + (seq(-i_num[1, 2] - 1, i_num[2, 2]) + .25)*i_len
      ## primary tessellation vertices (xmin, ymin, xmax, ymax):
      ## left-to-right, bottom-to-top
      a_bounds <- cbind(
        rep(a_xbnd[-length(a_xbnd)], each = length(a_ybnd) - 1),
        rep(a_ybnd[-length(a_ybnd)], times = length(a_xbnd) - 1),
        rep(a_xbnd[-1], each = length(a_ybnd) - 1),
        rep(a_ybnd[-1], times = length(a_xbnd) - 1)
      )
      
      ## secondary tessellation vertices (rotate & shift primary vertices)
      b_bounds <- t(2*rep(f_cent, times = 2) -
                      t(a_bounds))[seq(prod(i_tot), 1), c(3, 4, 1, 2)]
      
      ## dual tessellation bounds
      ls_bounds <- rbind(a_bounds, b_bounds)
    } else {
      stopifnot(index %in% self$index_set)
      s_dual <- as.logical(match(substr(index, 2L, 2L), c("A", "B")) - 1L)
      s_order <- as.integer(strsplit(
        sub("^\\([AB]: ([0-9,\\-]+)\\)$", "\\1", index),
        split = ","
      )[[1]])
      ls_bounds <- rep(f_cent, times = 2) + if (s_dual) {
        (s_order + c(0, 0, 1, 1) - .25)*i_len
      } else {
        (s_order + c(-1, -1, 0, 0) + .25)*i_len
      }
    }
    
    ## Return bounds
    return(ls_bounds)
  }
)

SquareDualTessellationCover$set(
  which = "public", name = "construct_index_set",
  value = function(fv) {
    if (is.function(fv)) {
      self$validate(fv)
      fv <- fv()
    }
    
    ## Get filter length
    f_len <- diff(apply(fv, 2, range))
    
    ## Calculate hyper-parameters
    i_len <- self$width*max(f_len)
    i_num <- rbind(
      ceiling(.5*f_len/i_len - .75),
      ceiling(.5*f_len/i_len - .25)
    )
    
    ## Primary, then dual, indices
    i_rc <- as.vector(t(outer(
      seq(-i_num[1, 1], i_num[2, 1]),
      seq(-i_num[1, 2], i_num[2, 2]),
      paste, sep = ","
    )))
    a_ind <- paste0("(A: ", i_rc, ")", sep = "")
    b_ind <- paste0("(B: ", i_rc, ")", sep = "")
    
    self$index_set <- c(a_ind, b_ind)
  }
)

## -+- Modularize this function to reduce duplication? -+-
SquareDualTessellationCover$set(
  which = "public", name = "construct_cover",
  value = function(filter, index = NULL) {
    stopifnot(is.function(filter))
    self$validate(filter)
    
    ## Get filter
    fv <- filter()
    
    ## If the index set hasn't been made yet, construct it.
    self$construct_index_set(fv)
    
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

SquareDualTessellationCover$set(
  which = "public", name = "neighborhood",
  value = function(filter, k) {
    stopifnot(!is.null(private$.index_set))
    
    if (k == 1L) {
      
      ## Get filter length
      f_len <- diff(apply(filter(), 2, range))
      
      ## Calculate hyper-parameters
      i_len <- self$width*max(f_len)
      i_num <- rbind(
        ceiling(.5*f_len/i_len - .75),
        ceiling(.5*f_len/i_len - .25)
      )
      i_tot <- apply(i_num, 2, sum) + 1
      
      ## Each set only overlaps with at most four sets, all in the dual layer:
      ## (A: x,y) intersects only (B: x,y), (B: x,y-1), (B: x-1,y), (B: x-1,y-1)
      ## Take advantage of naming convention by generating pairs of IDs directly
      a_coord <- rbind(
        x = rep(seq(-i_num[1, 1], i_num[2, 1]), each = i_tot[2]),
        y = rep(seq(-i_num[1, 2], i_num[2, 2]), times = i_tot[1])
      )[, rep(seq(prod(i_tot)), each = 4), drop = FALSE]
      b_coord <- a_coord - c(0, 0, 0, 1, 1, 0, 1, 1)
      w_coord <- which(b_coord["x", ] >= -i_num[1, "x"] &
                         b_coord["y", ] >= -i_num[1, "y"])
      ab_ind <- matrix(paste("(", c("A", "B"), ": ", rbind(
        apply(a_coord[, w_coord, drop = FALSE], 2, paste, collapse = ","),
        apply(b_coord[, w_coord, drop = FALSE], 2, paste, collapse = ",")
      ), ")", sep = ""), ncol = 2, byrow = TRUE)
      
      return(ab_ind)
    } else if (k > 1L) {
      
      ## The cover has two partitioning layers, so no three sets overlap
      return(matrix("", nrow = 0, ncol = k + 1L))
      
    } else {
      ## Default to the generic neighborhood calculator
      return(super$neighborhood(filter, k))
    }
    
  }
)
