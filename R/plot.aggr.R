#' @title Plot Amostra Agregada
#'
#' @export

plot.aggr <- function(y, main='Sample', xlab='x', ylab='y', ...) {
  x <- (1:nrow(y$fun))/nrow(y$fun)
  plot(x=x, y=x, main=main, xlab=xlab, ylab=ylab, type='n',
       ylim=range(y$fun), ...)
  for (i in 1:ncol(y$fun)) lines(x, y$fun[,i], col=i)
}
