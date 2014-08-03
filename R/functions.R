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
#' a <- c(1,2,15)
#' b <- 1:10
#' a %nin% b
#' !(a %in% b)
`%nin%` <- Negate(`%in%`)



#' @title Read dat-file of MDA1 program
#' @description Read dat-file contained descriptors of set of chemical compounds.
#' @param fname.dat name of dat-file created by MDA1 program.
#' @return data.frame of descriptors of set of chemical compounds.
#' @export
#' @examples
#' df <- read.mda.dat("all_simplexes.dat")
read.mda.dat <- function(fname.dat) {
  # read cds file
  fname.cds <- sub(".dat$", ".cds", fname.dat)
  lines <- scan(fname.cds, sep = "\n", quote = "", what = character())
  lines <- strsplit(lines, " ")
  nCompounds <- as.numeric(lines[[1]][2])
  nDescriptors <- as.numeric(lines[[1]][1])
  # read dat file
  df <- readBin(fname.dat, "numeric", n = nCompounds * nDescriptors, size = 4)
  df <- as.data.frame(matrix(df, nCompounds, nDescriptors, byrow = TRUE))
  # set names for compounds and descriptors
  colnames(df) <- sapply(lines[2:(nDescriptors+1)], "[", 3)
  rownames(df) <- sapply(lines[(nDescriptors+2):length(lines)], "[", 2)
  return(df)
}



#' @title Cross-validation result of optimal caret model
#' @description Returns result of cross-validation step for the optimal model obtained by caret.
#' @param caret.model object of classs train returned by caret train function.
#' @return vector of predicted values of all folds during cross-validation.
#' @export
#' @examples
#' cv.pred <- getCV(pls.model)
getCV <- function(caret.model) {
  df <- caret.model$pred
  best <- caret.model$bestTune
  ids <- apply(df[ ,names(best), drop=FALSE], 1, function(r) all(r == best[1,]) )
  df <- df[ids, ]
  df <- df[order(df$rowIndex), c("pred")]
  return(df)
}



#' @title Named list
#' @description Creates a named list from given objects.
#' @usage named.list(...)
#' @param ... comma separated list of objects.
#' @return list of given objects with corresponding names.
#' @seealso \code{\link{list}}.
#' @export
#' @examples
#' a <- 1:10
#' b <- LETTERS[1:5]
#' named.list(a, b)
named.list <- function(...) {
  names <- as.list(substitute(list(...)))[-1]
  setNames(list(...), names)
}



#' @title Adjust dataset
#' @description Add new columns to data.frame, fill them with specified value and trim columns which are absent in col.names.
#' @param df data.frame which should be modified.
#' @param col.names columns names which should be present in output data.frame.
#' @param value value which should fill missing columns.
#' @return data.frame containing only specified columns in that order. Previously missing columns will be fill with the specified value.
#' @details this function is useful in adjust features of one dataset (e.g. test set) to another one (e.g. training set). Because some models need exactly the same features in the same order in test set as in training set to make correct prediction.
#' @export
#' @examples
#' df <- data.frame(A=1:10, B=11:20, C=21:30)
#' adjust.dataset(df, c("A","F","B"))
adjust.dataset <- function(df, col.names, value = 0) {
  df[,setdiff(col.names, colnames(df))] <- value
  df[,col.names]
}



#' @title Expand two-level factor
#' @description Create data.frame containing columns corresponding to first part of two-level factor and values equal to the second part of the latter. .
#' @param v vector of two-level factor variable. Each element represent as at least two-letter string.
#' @param split character, which will be used to split \code{v} on two parts.
#' @param var.name name of variable, which will be used for generation names of newly created variables.
#' @param comb.sep character, which will be used as a separator in newly created variable names between \code{var.name} and first part of two-level factor \code{v}.
#' @return data.frame.
#' @export
#' @examples
#' vec <- c("A1","B4",NA,"A4","C2")
#' expand.two.level.factor(vec)
expand.two.level.factor <- function(v, split="", var.name="var", comb.sep=":") {
  df <- as.data.frame(do.call(rbind, strsplit(v, split)), stringsAsFactors = F)
  colnames(df) <- c(var.name, "value")
  f <- as.formula(paste("~", colnames(df)[1], "-1", sep=""))
  m <- model.matrix(f, df)
  p <- paste("^(", var.name, ")(.*)$", sep="")
  s <- paste("\\1", comb.sep, "\\2", sep="")
  colnames(m) <- sub(p, s, colnames(m))
  m <- t(m)
  m[m == 1] <- na.omit(df[,2])
  m <- t(m)
  m <- as.data.frame(m, stringsAsFactors = F)
  m <- m[rownames(df),]
  rownames(m) <- rownames(df)
  m[is.na(m)] <- 0
  m
}
