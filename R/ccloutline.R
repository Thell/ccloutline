#' @title Connected component labelling and eulerian cycle outlines.
#'
#' @description Given a 2 dimensional matrix generate the cartesian
#' coordinates to plot a grid graph path outlining the outer and inner
#' component edges.
#'
#' @name ccloutline
#' @docType package
#' @author Thell Fowler \email{thell@tbfowler.name}
#' @keywords package ccl border outline
#' @useDynLib ccloutline
#' @importFrom Rcpp sourceCpp
NULL

# groups dispersed randomly, seed for repeatability.
.random_source_gen <- function(n_row, n_col, group.count, seed = 1, prob=NULL) {
  print(match.call())
  set.seed(seed)
  dat <- sample.int( group.count, n_col * n_row,
                     group.count < (n_col * n_row), prob )
  array( dat, dim = c(n_row,n_col) )
}

.checkerboard_source_gen <- function(n_row, n_col) {
  print(match.call())
  sapply( 1:n_col, function(i) {
    sapply( 1:n_row, function(j) {
      ((i + j) %% 2) + 1
    })
  })
}

#' Connected component extraction.
#'
#' A matrix of binary values representing components.
#'
#' @format A matrix 9 rows and 17 variables:
#' @source \url{https://en.wikipedia.org/wiki/Connected-component_labeling#Graphical_example_of_two-pass_algorithm}
"wiki_sample_mat"
