#' @title Triangle Dual Tesselation Cover
#'
#' @docType class
#' @author Jason Cory Brunson
#' @description [R6] class representing a triangle--hexagon dual tessellation
#'   cover.
#' @family cover

#' @details
#'
#' The regular triangle tessellation \eqn{T} of \eqn{Z} (the plane, or a
#' path-connected region of it; Schlafli symbol \eqn{{3,6}}) consists of
#' congruent equiangular triangles arranged as tiles, with shared vertices and
#' edges and empty overlap. In order to completely partition \eqn{Z}, each
#' triangle is considered closed along part of its perimeter and open along the
#' remainder. The dual \eqn{T'} of \eqn{T} (the regular hexagonal tessellation;
#' Schlafli symbol \eqn{6,3}) is obtained by taking vertices at the centroids of
#' the tiles of \eqn{T} and edges perpendicular to those of \eqn{T}. The union
#' \eqn{C} of \eqn{T} and \eqn{T'} is a _triangle--hexagon dual tessellation
#' cover_.
#'
#' A triangle--hexagon dual tessellation cover \eqn{C} has constant "depth", in
#' the sense that every point in \eqn{Z} is a member of exactly two sets in
#' \eqn{C}. A cover of depth 2 is only sufficient to recover 1-dimensional
#' features of a point cloud. Every non-empty intersection of cover sets is
#' congruent up to boundary.
#'
#' `TriangleDualTessellationCover` takes the centroid of \eqn{f(X)} to be the
#' centroid of the overlap between a pair of tiles in \eqn{T} and in \eqn{T'}.
#' When `width` is sufficiently large, the total number of tiles is minimized at
#' 6 (though some may be empty).
#' 

#' @field width A positive number. The side length of each triangle tile and the
#'   side-to-side length of each hexagon tile, as a fraction of the larger
#'   coordinate range of the lensed data.
#' @export
TriangleDualTessellationCover <- R6::R6Class(
  classname = "TriangleDualTessellationCover",
  inherit = CoverRef,
  private = list(.width = NULL)
)

#' @export
TriangleDualTessellationCover$set(
  which = "public", name = "initialize",
  value = function(...) {
    super$initialize(typename = "Triangle-Hexagon Dual Tessellation")
    params <- list(...)
    if ("width" %in% names(params)) {
      self$width <- params[["width"]]
    }
  }
)

TriangleDualTessellationCover$set(
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

TriangleDualTessellationCover$set(
  which = "public", name = "validate",
  value = function(filter) {
    stopifnot(! is.null(private$.width))
    stopifnot(length(self$width) == 1)
    f_dim <- ncol(filter())
    if (f_dim != 2) {
      stop("A planar tessellation requires a 2-dimensional lens.")
    }
  }
)

TriangleDualTessellationCover$set(
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

## Calculate x and y positions of triangle tiles (for set indices and bounds)
triangle_xy <- function(f_ran, tri_num) {
  stopifnot(all.equal(dim(f_ran), c(2L, 2L)))
  stopifnot(all.equal(dim(tri_num), c(2L, 2L)))
  
  tri_tot <- apply(tri_num, 2, sum) + 1
  cbind(
    x = rep(seq(-tri_num[1, 1], tri_num[2, 1]), each = tri_tot[2]),
    y = rep(seq(-tri_num[1, 2], tri_num[2, 2]), times = tri_tot[1])
  )
}

## Convert triangle x and y positions to ternary positions (as matrices)
triangle_xy2abc <- function(xy) {
  if (is.vector(xy)) xy <- t(xy)
  stopifnot(ncol(xy) == 2L)
  
  ## Tile orientation
  inv <- apply(xy, 1, sum) %% 2L
  ## Matrix of ternary positions
  cbind(
    # in direction of 1/6*pi
    a = as.integer((xy[, 2] - xy[, 1] + inv) / 2),
    # in direction of 1/2*pi
    b = xy[, 2],
    # in direction of 5/6*pi
    c = as.integer((xy[, 2] + xy[, 1] + inv) / 2)
  )
}

## Convert triangle ternary positions to x and y positions (as matrices)
triangle_abc2xy <- function(abc) {
  if (is.vector(abc)) abc <- t(abc)
  stopifnot(ncol(abc) == 3L)
  
  ## Matrix of x,y positions
  cbind(
    # x coordinate
    x = abc[, 1] - abc[, 3],
    # y coordinate
    y = abc[, 2]
  )
}

## Calculate bounds for the primary (triangle) tessellation (for plotting)
TriangleDualTessellationCover$set(
  which = "public", name = "primary_bounds",
  value = function(filter, index = NULL) {
    self$validate(filter)
    
    ## Get filter lengths
    fv <- filter()
    f_ran <- apply(fv, 2, range)
    f_len <- diff(f_ran)
    f_cent <- f_ran[1, ] + f_len/2
    
    ## Calculate hyper-parameters
    y_bar <- (1/9 - 1/16) / sqrt(3) / (1/3 - 1/4) # origin overlap centroid
    s_len <- f_len[1]*self$width # triangle side length
    d_hts <- f_len[2] / s_len / sqrt(3) # number of triangle heights in data/2
    tri_num <- ceiling(cbind(
      # numbers of triangles on each side of origin triangle
      x = 1/self$width,
      y = d_hts - c(sqrt(3)/2 - y_bar, y_bar)*2/sqrt(3)
    ))
    
    ## Calculate ternary positions from Cartesian positions
    tri_abc <- triangle_xy2abc(triangle_xy(f_ran, tri_num))
    
    ## Tile orientation
    inv <- apply(tri_abc, 1, sum) %% 2
    
    ## Calculate unit normal and origin projection for each side of each tile
    tri_bnd <- cbind(
      u1 = -sqrt(3)/2 * (-1)^inv,
      v1 = -1/2 * (-1)^inv,
      x1 = f_cent[1] + (sqrt(3)/4*y_bar + floor(tri_abc[, 1]/2) * 3/2)*s_len,
      y1 = f_cent[2] + (1/4*y_bar + floor(tri_abc[, 1]/2) * sqrt(3)/2)*s_len,
      u2 = 0,
      v2 = 1 * (-1)^inv,
      x2 = f_cent[1],
      y2 = f_cent[2] +
        (-sqrt(3)/2 + y_bar + ceiling(tri_abc[, 2]/2) * sqrt(3))*s_len,
      u3 = sqrt(3)/2 * (-1)^inv,
      v3 = -1/2 * (-1)^inv,
      x3 = f_cent[1] - (sqrt(3)/4*y_bar + floor(tri_abc[, 3]/2) * 3/2)*s_len,
      y3 = f_cent[2] + (1/4*y_bar + floor(tri_abc[, 3]/2) * sqrt(3)/2)*s_len
    )
    
    return(tri_bnd)
  }
)

## Calculate x and y positions of hexagon tiles (for set indices and bounds)
hexagon_xy <- function(f_ran, hex_num) {
  stopifnot(all.equal(dim(f_ran), c(2L, 2L)))
  stopifnot(all.equal(dim(hex_num), c(2L, 2L)))
  
  hex_tot <- apply(hex_num, 2, sum) + 1
  hex_xy <- cbind(
    x = rep(seq(-hex_num[1, 1], hex_num[2, 1]), each = hex_tot[2]),
    y = rep(seq(-hex_num[1, 2], hex_num[2, 2]), times = hex_tot[1])
  )
  hex_xy[apply(hex_xy %% 2, 1, sum) != 1, ]
}

## Convert hexagon ternary positions to x and y positions (as matrices)
hexagon_abc2xy <- function(abc) {
  if (is.vector(abc)) abc <- t(abc)
  stopifnot(ncol(abc) == 3L)
  
  ## Matrix of x,y positions
  y <- floor((-abc[, 1] + abc[, 3]) / 3)
  y_mod2 <- y %% 2
  cbind(
    # x coordinate
    x = 2 * floor((abc[, 2] - y_mod2) / 2) + y_mod2,
    # y coordinate
    y = y
  )
}

TriangleDualTessellationCover$set(
  which = "public", name = "construct_index_set",
  value = function(filter) {
    if (is.function(filter)) {
      self$validate(filter)
      fv <- filter()
    } else {
      fv <- filter
    }
    
    ## Get filter lengths
    f_ran <- apply(fv, 2, range)
    f_len <- diff(f_ran)
    
    ## Calculate hyper-parameters
    y_bar <- (1/9 - 1/16) / sqrt(3) / (1/3 - 1/4) # origin overlap centroid
    s_len <- f_len[1]*self$width # triangle side length
    d_hts <- f_len[2] / s_len / sqrt(3) # number of triangle heights in data/2
    tri_num <- ceiling(cbind(
      # numbers of triangles on each side of origin triangle
      x = 1/self$width,
      y = d_hts - c(sqrt(3)/2 - y_bar, y_bar)*2/sqrt(3)
    ))
    hex_num <- ceiling(cbind(
      # numbers of hexagons on each side of origin hexagon
      x = 1/self$width,
      y = d_hts - c(1/4 - y_bar, 1/4 + y_bar)*2/sqrt(3)
    ))
    
    ## Primary (triangle) indices
    tri_xy <- triangle_xy(f_ran, tri_num)
    ## Dual (hexagon) indices
    hex_xy <- hexagon_xy(f_ran, hex_num)
    
    ## Index sets
    self$index_set <- c(
      paste0("(A: ", apply(tri_xy, 1, paste, collapse = ","), ")"),
      paste0("(B: ", apply(hex_xy, 1, paste, collapse = ","), ")")
    )
  }
)

## Given the current set of parameter values,
## construct the level sets whose union covers the filter space,
## such that every point lies in exactly one set of each tessellation
TriangleDualTessellationCover$set(
  which = "public", name = "construct_cover",
  value = function(filter, index = NULL) {
    stopifnot(is.function(filter))
    self$validate(filter)
    
    ## Calculate filter and summary statistics
    fv <- filter()
    f_ran <- apply(fv, 2, range)
    f_len <- diff(f_ran)
    f_cent <- f_ran[1, ] + f_len/2
    
    ## Calculate additional hyper-parameters
    y_bar <- (1/9 - 1/16) / sqrt(3) / (1/3 - 1/4) # origin overlap centroid
    s_len <-
      # triangle side length
      f_len[1]*self$width +
      # ensure all lensed points are covered
      sqrt(.Machine$double.eps)
    
    ## If the index set hasn't been made yet, construct it
    self$construct_index_set(fv)
    
    ## If no index specified, return the level sets either by construction
    if (missing(index) || is.null(index)) {
      #stopifnot(! index %in% self$index_set)
      
      ## Center lensed point cloud
      f_vec <- sweep(fv, 2, f_cent, "-")
      
      ## Triangle memberships
      f_tri_abc <- floor(cbind(
        a = (f_vec %*% c(sqrt(3), 1)/2) / s_len + sqrt(3)/2 - 1/2*y_bar,
        b = f_vec[, 2] / s_len + sqrt(3)/2 - y_bar,
        c = (f_vec %*% c(-sqrt(3), 1)/2) / s_len + sqrt(3)/2 - 1/2*y_bar
      ) * 2/sqrt(3))
      f_tri_xy <- triangle_abc2xy(f_tri_abc)
      f_tri_set <-
        paste0("(A: ", apply(f_tri_xy, 1, paste, collapse = ","), ")")
      stopifnot(all(f_tri_set %in% self$index_set))
      
      ## Hexagon memberships
      ## (zeroth sub-triangle column is to the left of the centroid)
      f_hex_abc <- floor(cbind(
        a = (f_vec %*% c(1, -sqrt(3))/2) / s_len + sqrt(3)/2*y_bar,
        b = f_vec[, 1] / s_len + 1/2,
        c = (f_vec %*% c(1, sqrt(3))/2) / s_len + sqrt(3)/2*y_bar
      ) * 2)
      f_hex_xy <- hexagon_abc2xy(f_hex_abc)
      f_hex_set <-
        paste0("(B: ", apply(f_hex_xy, 1, paste, collapse = ","), ")")
      stopifnot(all(f_hex_set %in% self$index_set))
      
      ## Store membership vectors as level sets
      lvl_set <- c(
        lapply(
          self$index_set[grep("\\(A\\:", self$index_set)],
          function(x) which(f_tri_set %in% x)
        ),
        lapply(
          self$index_set[grep("\\(B\\:", self$index_set)],
          function(x) which(f_hex_set %in% x)
        )
      )
      names(lvl_set) <- self$index_set
      self$level_sets <- lvl_set
      return(invisible(self))
    } else {
      
      if (! is.null(self$level_sets) && index %in% names(self$level_sets)) {
        
        return(self$level_sets[[index]])
      } else {
        stopifnot(index %in% self$index_set)
        
        p_idx <- which(index == self$index_set)
        stop("Cannot yet compute the level set for a single tile; ",
             "use `.$construct_cover()` to compute all level sets.")
        #lvl_set <- ...
        
        #return(lvl_set)
      }
    }
    
  }
)

TriangleDualTessellationCover$set(
  which = "public", name = "neighborhood",
  value = function(filter, k) {
    stopifnot(!is.null(private$.index_set))
    
    if (k == 1L) {
      
      ## Get filter length
      f_ran <- apply(filter(), 2, range)
      f_len <- diff(f_ran)
      
      ## Calculate hyper-parameters
      y_bar <- (1/9 - 1/16) / sqrt(3) / (1/3 - 1/4) # origin overlap centroid
      s_len <- f_len[1]*self$width # triangle side length
      d_hts <- f_len[2] / s_len / sqrt(3) # number of triangle heights in data/2
      tri_num <- ceiling(cbind(
        # numbers of triangles on each side of origin triangle
        x = 1/self$width,
        y = d_hts - c(sqrt(3)/2 - y_bar, y_bar)*2/sqrt(3)
      ))
      hex_num <- ceiling(cbind(
        # numbers of hexagons on each side of origin hexagon
        x = 1/self$width,
        y = d_hts - c(1/4 - y_bar, 1/4 + y_bar)*2/sqrt(3)
      ))
      
      ## Primary (triangle) indices
      tri_xy <- triangle_xy(f_ran, tri_num)
      tri_ran <- apply(tri_xy, 2, range)
      ## Dual (hexagon) indices
      hex_xy <- hexagon_xy(f_ran, hex_num)
      
      ## Each hexagon intersects an x,y grid of (at most) 6 triangles:
      ## (B: x,y) intersects (A: {x-1, x, x+1}*{y, y+1})
      b_xy <- hex_xy[rep(seq(nrow(hex_xy)), each = 6), , drop = FALSE]
      a_xy <- cbind(
        x = b_xy[, 1] + rep(seq(-1, 1), each = 2),
        y = b_xy[, 2] + rep(seq(0, 1), times = 3)
      )
      w_xy <- which(a_xy[, 1] >= tri_ran[1, 1] & a_xy[, 1] <= tri_ran[2, 1] &
                      a_xy[, 2] >= tri_ran[1, 2] & a_xy[, 2] <= tri_ran[2, 2])
      ab_xy <- matrix(paste("(", c("A", "B"), ": ", rbind(
        apply(a_xy[w_xy, , drop = FALSE], 1, paste, collapse = ","),
        apply(b_xy[w_xy, , drop = FALSE], 1, paste, collapse = ",")
      ), ")", sep = ""), ncol = 2, byrow = TRUE)
      
      return(ab_xy)
    } else if (k > 1L) {
      
      ## The cover has two partitioning layers, so no three sets overlap
      return(matrix("", nrow = 0, ncol = k + 1L))
      
    } else {
      ## Default to the generic neighborhood calculator
      return(super$neighborhood(filter, k))
    }
    
  }
)
