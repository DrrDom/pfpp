#' Load Rdata container into a local variable.
#'
#' @param file.name name of RData file to load.
#' @return object contained in RData file.
#' @seealso \code{\link{load}} which this function wraps.
#' @export
#' @examples
#' a <- 1:10
#' save(a, file="a.RData")
#' b <- local.load("a.RData")
local.load <- function(file.name) local(get(load(file.name)))



#' @title Statistic for binary classification task.
#' @description Returns TP, TN, FP, FN, sensitivity, specificity, balanced accuracy and accuracy for binary classification task.
#' @param pred vector of predicted labels (0, 1).
#' @param obs vector of observed labels (0, 1).
#' @return names vector contained TP, TN, FP, FN, sensitivity, specificity, balanced accuracy and accuracy.
#' @export
#' @examples
#' set.seed(42)
#' obs <- sample(0:1,10,T)
#' pred <- sample(0:1,10,T)
#' getBinaryStat(pred, obs)
getBinaryStat <- function (pred, obs) {
  if (length(obs) != length(pred)) {
    return(NULL)
  }
  tbl <- table(obs, pred)
  lst <- list()
  lst$TN <- tbl[1,1]
  lst$TP <- tbl[2,2]
  lst$FN <- tbl[2,1]
  lst$FP <- tbl[1,2]
  lst$sensitivity <- lst$TP / (lst$TP + lst$FN)
  lst$specificity <- lst$TN / (lst$TN + lst$FP)
  lst$balanced.acc <- (lst$sensitivity + lst$specificity) / 2
  lst$acc <- mean(obs == pred)
  return(unlist(lst))
}



#' @title Determination coefficient for external test set or cross-validation.
#' @description Returns determination coefficient for external test set or cross-validation.
#' @param ts.pred predicted values for external test set or cross-validation folds.
#' @param ts.obs observed values for external test set or cross-validation folds.
#' @param ws.obs.mean mean value for training set. If determination coefficient for cross-validation folds is calculated, then this parameter can be left as is.
#' @return value of determination coefficient calculated according to formula 1-PRESS/SS.
#' @export
#' @examples
#' set.seed(42)
#' cv.obs <- runif(50, 0, 10)
#' cv.pred <- cv.obs + rnorm(50)
#' R2test(cv.pred, cv.obs)
R2test <- function(ts.pred, ts.obs, ws.obs.mean = mean(ts.obs)) {
  return(1 - sum((ts.pred - ts.obs)**2)/sum((ts.obs - ws.obs.mean)**2))
}



#' @title Propability of difference between two correlation coefficients.
#' @description Retuns propability of significance of difference between two correlation coefficients.
#' @param x.cor first correlation coefficient.
#' @param nx number of objects used in calculation of first correlation coefficient.
#' @param y.cor second correlation coefficient.
#' @param ny number of objects used in calculation of second correlation coefficient.
#' @return propability of significance of differance between two correlation coefficients.
#' @export
#' @examples
#' set.seed(42)
#' cv.obs <- runif(50, 0, 10)
#' cv.pred <- cv.obs + rnorm(50)
#' R2test(cv.pred, cv.obs)
com.cor <- function(x.cor, nx, y.cor, ny){
  z1=0.5*log((1+x.cor)/(1-x.cor))
  z2=0.5*log((1+y.cor)/(1-y.cor))
  w=sqrt(1/(nx-1)+2/(nx-1)^2+1/(ny-1)+2/(ny-1)^2)
  2*pnorm(-abs((z1-z2)/w),0,1)
}


#' @title Non-matching comparison
#' @description Returns a logical vector indicating if there is a non-match or match for its left operand. Inverse action relatively to \%in\% function
#' @usage x \%nin\% y
#' @param x vector of values to be matched
#' @param y vector of values to be matched against.
#' @return logical vector of TRUE (non-match) and FALSE (match).
#' @seealso \code{\link{\%in\%}}.
#' @export
#' @examples
#' set.seed(42)
#' a <- sample(0:1, 20, T)
#' b <- sample(0:1, 20, T)
#' a %nin% b
#' !(a %in% b)
`%nin%` <- Negate(`%in%`)
