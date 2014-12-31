#' Load RData container into a local variable.
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
#' com.cor(0.90, 100, 0.82, 150)
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
#' @description Create data.frame containing columns corresponding to first part of two-level factor and values equal to the second part of the latter.
#' @param vec vector of two-level factor or character variable. Each element is represented as at least two-letter string.
#' @param split character, which will be used to split \code{vec} on two parts.
#' @param var.name name of variable, which will be used for generation names of newly created variables.
#' @param comb.sep character, which will be used as a separator in newly created variable names between \code{var.name} and first part of two-level factor \code{vec}.
#' @param fill default value to represent missing (\code{NA}) values in input vector.
#' @return data.frame.
#' @export
#' @examples
#' vec <- c("A1","B4",NA,"A4","C2")
#' expand.two.level.factor(vec)
expand.two.level.factor <- function(vec, split="", var.name="var", 
                                    comb.sep=":", fill=NA) {
  if (!is.character(vec)) 
    vec <- as.character(vec)
  df <- as.data.frame(do.call(rbind, strsplit(vec, split)), stringsAsFactors = F)
  colnames(df) <- c(var.name, "value")
  f <- as.formula(paste("~", colnames(df)[1], "-1", sep=""))
  m <- model.matrix(f, df)
  p <- paste("^(", var.name, ")(.*)$", sep="")
  s <- paste("\\1", comb.sep, "\\2", sep="")
  colnames(m) <- sub(p, s, colnames(m))
  m <- t(m)
  m[m == 1] <- na.omit(df[,2])
  m <- t(m)
  mm <- matrix(fill, nrow=length(vec), ncol=ncol(m))
  rownames(mm) <- as.character(1:nrow(mm))
  mm[rownames(m),] <- m
  colnames(mm) <- colnames(m)
  rownames(mm) <- NULL
  as.data.frame(mm, stringsAsFactors = FALSE)
}



#' @title Expand one-level factor
#' @description Creates a design matrix like \code{model.matrix} and fill unused missing (\code{NA}) values with rows of 0.
#' @param vec factor or character vector..
#' @param var.name variable name, which will be used for generation names of newly created variables.
#' @param comb.sep character, which will be used as a separator in newly created variable names between \code{var.name} and \code{vec} levels.
#' @param fill default value to represent missing (\code{NA}) values in input vector.
#' @return data.frame.
#' @seealso \code{\link{model.matrix}}.
#' @export
#' @examples
#' vec <- c("A","B",NA,"A","C")
#' expand.one.level.factor(vec)
expand.one.level.factor <- function(vec, var.name="var", 
                                    comb.sep=":", fill=NA) {
  if (!is.character(vec)) 
    vec <- as.character(vec)
  m <- model.matrix(~ vec - 1)
  p <- paste("^(vec)(.*)$", sep="")
  s <- paste(var.name, comb.sep, "\\2", sep="")
  colnames(m) <- sub(p, s, colnames(m))
  # init zero matrix of target size
  mm <- matrix(fill, nrow=length(vec), ncol=ncol(m))
  rownames(mm) <- as.character(1:nrow(mm))
  mm[rownames(m),] <- m
  colnames(mm) <- colnames(m)
  rownames(mm) <- NULL
  as.data.frame(mm, stringsAsFactors = FALSE)
}



#' @title Read tab-delimited text file with sirms descriptors (output of sirms.py)
#' @description Load tab-delimited text file with sirms descriptors in \code{data.frame}.
#' @param fname file name of tab-delimited text file with sirms descriptors.
#' @param ... further arguments for \code{data.table::fread} function.
#' @return data.frame.
#' @export
sirms.read <- function(fname, ...) {
  if (!require(data.table)) stop("Error. Install data.table package.")
  x <- as.data.frame(data.table::fread(fname, sep="\t", header = TRUE, ...))
  rownames(x) <- x[,1]
  x <- x[,-1]
  return(x)
}



#' @title Remove constant variables
#' @description remove constant variables from a data.frame.
#' @param df data.frame from which constant variable will be removed.
#' @return data.frame with removed constant variables.
#' @export
#' @examples
#' df <- data.frame(a = rnorm(5), b = 1:5, c = rep(1, 5))
#' df <- remove.const.vars(df)
remove.const.vars <- function(df) {
  df <- df[, sapply(df, function(v) length(unique(v)) > 1)]
  return(df)
}



#' @title Creation of folds for group-out cross-validation based on Monte Carlo method.
#' @description Function takes input vector of variables which will be used to group 
#' input objects. It creates folds with approaximately equal number of objects in each fold. 
#' All object of the same group put in the single fold.
#' @param v input vector of group variables.
#' @param nfolds maximum number of folds. But if it is not possible to create 
#' the specified number of folds its number will be automatically decreased.
#' @param max_iter maximum number of Monte Carlo iterations.
#' @param start_opt_param strting parameter of Monter Carlo optimization.
#' @param error_limit if final error will be greater than this values results 
#' will be discarded.
#' @param seed seed for random number generator to reproduce results.
#' @return vector of folds numbers for each input object.
#' @export
#' @examples
#' df <- data.frame(A=runif(100), B=runif(100), C=sample(LETTERS[1:10], 100, T))
#' res <- create_folds_mc(df$C, 5)
create_folds_mc <- function(v, nfolds = "auto", max_iter = 1000, start_opt_param = 100, 
                            error_limit = NA, seed = 0) {
  
  opt_func <- function(tmp, opt) {
    return(var(aggregate(tmp, list(opt), sum)[,2]))
  }
  
  replace_value <- function(input_vector, value, pos) {
    input_vector[pos] <- value
    return(input_vector)
  }
  
  if (nfolds == "auto") {
    ngroups <- round(length(v) / max(table(v)))
    nfolds <- min(5, ngroups)
  }   
  
  tmp <- table(v)
  set.seed(seed)
  tmp <- sample(tmp)
  
  opt <- rep(1:nfolds, length(tmp) / (nfolds) + 1)
  length(opt) <- length(tmp)
  
  dec_param <- start_opt_param / max_iter
  niter <- 1
  cur_var <- opt_func(tmp, opt)
  
  while (niter < max_iter) {
    # avoid lost of fold number
    new_opt <- replace_value(opt, sample.int(nfolds, 1), sample.int(length(opt), 1))
    while (length(unique(new_opt)) < nfolds)
      new_opt <- replace_value(opt, sample.int(nfolds, 1), sample.int(length(opt), 1))
    new_var <- opt_func(tmp, new_opt)
    if (new_var < cur_var) {
      opt <- new_opt
      cur_var <- new_var
    } else {
      if (exp( -(new_var - cur_var) / start_opt_param ) > runif(1)) {
        opt <- new_opt
        cur_var <- new_var
      }
    }
    start_opt_param <- start_opt_param - dec_param
    niter <- niter + 1
  }
  
  names(opt) <- names(tmp)
  
  # if error is greater than specified limit then discard this split
  if (!is.na(error_limit) & sqrt(cur_var) > error_limit) {
    return(NULL)
  } else {
    return(opt[as.character(v)])
  }
  
}



#' @title Calculate Tanimoto similarity coefficient beetwen two sets of folds for cross-validation.
#' @description Calculate Tanimoto similarity coefficient beetwen two sets of folds for 
#' cross-validation considering objects grouping in each folds
#' @param g1,g2 input vectors with folds numbers for each objects.
#' @return mean Tanimoto between input vectors considering objects grouping in folds
#' @export
#' @examples
#' set.seed(42)
#' g1 <- sample(1:5, 100, T)
#' g2 <- sample(1:5, 100, T)
#' groupwise_tanimoto(g1, g2)
groupwise_tanimoto <- function(g1, g2) {
  
  tanimoto <- function(a, b) {
    return(length(intersect(a, b)) / length(union(a, b)))
  }
  
  g1 <- split(seq_along(g1), g1)
  g2 <- split(seq_along(g2), g2)
  result <- sapply(g1, function(i) {
    max(sapply(g2, function(j) {
      tanimoto(i, j)
    }))
  })
  mean(result)
}



#' @title Select dissimilar sets of folds.
#' @description Select most dissimilar folds sets from the input list based on groupwise Tanimoto measure.
#' @param folds_list list of sets of folds, each element is a vector of folds numbers.
#' @param max_sets maximum number of sets which will be remained
#' @param max_sim maximum similarity value, only sets which have less or equal similarity to all other sets will be remained.
#' @return list of sets of folds which are dissimilar to each other
#' @export
#' @examples
#' folds_list <- lapply(1:10, function(i) sample.int(5, 100, T))
#' selected_folds <- select_folds(folds_list, 5, 0.7)
select_folds <- function(folds_list, max_sets = 10, max_sim = 0.8, print_sim_matrix = FALSE) {
  
  res_group <- combn(seq_along(folds_list), 2, 
                     function(i) round(pfpp::groupwise_tanimoto(folds_list[[i[1]]], folds_list[[i[2]]]), 3))
  
  m <- matrix(0.0, length(folds_list), length(folds_list))
  m[lower.tri(m)] <- res_group
  colnames(m) <- seq_along(folds_list)
  rownames(m) <- seq_along(folds_list)
  
  while(nrow(m) > max_sets | max(m) > max_sim) {
    i <- which.max(apply(m, 1, max))
    m <- m[-i, -i]
  }
  
  if (print_sim_matrix)
    print(m)
  
  return(folds_list[as.integer(rownames(m))])
  
}



#' @title Create folds of test sets to use in caret.
#' @description Transform list of sets of folds to list which may be passed to caret train function.
#' @param folds_list list of sets of folds, each element is a vector of folds numbers.
#' @return list of training set object for each fold in format which suit the requirements of caret train function.
#' @export
#' @examples
#' folds_list <- lapply(1:10, function(i) sample.int(5, 100, T))
#' caret_folds_list <- create_caret_folds(folds_list)
create_caret_folds <- function(folds_list) {
  caret_folds <- lapply(folds_list, function(f) {
    lapply(unique(f), function(i) which(f != i))
  })
  caret_folds <- unlist(caret_folds, recursive = FALSE)
  names(caret_folds) <- unlist(lapply(seq_along(folds_list), function(i) {
    paste0("Fold", 1:length(unique(folds_list[[i]])), ".Rep", i)
  }))
  return(caret_folds)
}



#' @title Calculate RMSE.
#' @description Returns root mean squared deviation.
#' @param obs observed values.
#' @param pred predicted values.
#' @return value of RMSE.
#' @export
#' @examples
#' obs <- runif(100)
#' pred <- runif(100)
#' rmse(obs, pred)
rmse <- function(obs, pred) {
  sqrt(sum((obs - pred) ^ 2) / (length(obs) - 1))
}


