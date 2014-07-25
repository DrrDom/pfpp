#' Load Rdata container into the local variable.
#'
#' @param file.name name of RData file to load.
#' @return object contained in RData file.
#' @seealso \code{\link{load}} which this function wraps.
#' @export
#' @examples
#' a <- local.load("model.RData")
local.load <- function(file.name) local(get(load(file.name)))