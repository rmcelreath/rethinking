#' Turn "rasterized" kernel into "rasterized" density function
#'
#' Turn "rasterized" kernel into "rasterized" density function
#'
#' @param x,y numerical vectors of equal length giving the (x, y) coordinates
#'        of a kernel function \code{k} with \code{k(x) = y}.
#' @return a rescaled version of \code{y} that turns the kernel into a pdf.
#' @examples
#' x <- seq(0, 1, by = 0.05)
#' y1 <- (x - 0.5)^2
#' y2 <- densify(x, y1)
#' if (require(lattice)) {
#'   xyplot(y1 + y2 ~ x, type = "l")
#' }

densify <- function(x, y) {
  if ( diff(range(diff(x))) > 0.1 * min(abs(diff(x))) )
    stop("`x` must be (approximately) equally spaced.")

  width <- mean(diff(x))
  y / sum(y) / width
}