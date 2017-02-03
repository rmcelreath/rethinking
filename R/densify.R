

densify <- function(x, y) {
  if ( diff(range(diff(x))) > 0.1 * min(abs(diff(x))) )
    stop("`x` must be (approximately) equally spaced.")

  width <- mean(diff(x))
  y / sum(y) / width
}