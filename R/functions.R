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
  lst$FN <- tbl[1,2]
  lst$FP <- tbl[2,1]
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
  
  if (grepl("\\.", df$Resample[1])) {
    df <- cbind(df, do.call(rbind, strsplit(df$Resample, "\\.")))
    colnames(df)[(ncol(df) - 1):ncol(df)] <- c("fold", "rep")
    # cast and remove row index column
    df <- reshape2::dcast(df, rowIndex ~ rep, value.var = "pred")[, -1]
  } else {
    df <- df[order(df$rowIndex), c("pred")]
  }
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
#' @param id.col column number or name which contains IDs and which will be converted 
#' into row.names. To omit this option set it to \code{NULL} or \code{NA}.
#' @param ... further arguments for \code{data.table::fread} function.
#' @return data.frame.
#' @export
sirms.read <- function(fname, id.col = 1, ...) {
  if (!require(data.table)) stop("Error. Install data.table package.")
  x <- as.data.frame(data.table::fread(fname, sep="\t", header = TRUE, ...))
  if (is.null(id.col) | is.na(id.col)) 
    return(x)
  if (is.numeric(id.col) & id.col > 0) {
    rownames(x) <- x[, id.col]
    x <- x[, -id.col]
  }
  if (is.character(id.col)) {
    id.col <- which(colnames(x) == id.col)
    rownames(x) <- x[, id.col]
    x <- x[, -id.col]
  }
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
    paste0("Fold", 
           formatC(1:length(unique(folds_list[[i]])), 
                   width = 1 + floor(log10(max(1:length(unique(folds_list[[i]]))))),
                   flag = "0"),
           ".Rep", i)
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



#' @title Combine datasets
#' @description Combine two datasets, add new columns absent in each of them and fill them with specified value.
#' @param df1 first data.frame to combine.
#' @param df2 second data.frame to combine.
#' @param value value which should fill missing columns.
#' @return combined data.frame.
#' @details this function is useful to rbind two data.frames with different columns.
#' @export
#' @examples
#' df1 <- data.frame(A=1:10, B=11:20, C=21:30)
#' df2 <- data.frame(E=1:10, C=11:20, B=21:30)
#' combine.datasets(df1, df2)
combine.datasets <- function(df1, df2, value = 0) {
  df1[,setdiff(colnames(df2), colnames(df1))] <- value
  df2[,setdiff(colnames(df1), colnames(df2))] <- value
  df2 <- df2[, colnames(df1)]
  rbind(df1, df2)
}



#' Copy data to clipboard suitable for Excel pasting.
#'
#' @param x vector, data.frame or matrix to copy.
#' @param row.names copy with row names.
#' @param col.names copy with column names.
#' @return nothing, just copy object to clipboard.
#' @details column names in copied object will be shifted one column left.
#' @export
#' @examples
#' a <- 1:10
#' names(a) <- LETTERS[1:10]
#' write.excel(a)
#' b <- data.frame(V1 = 1:10, V2 = rnorm(10), row.names = LETTERS[1:10])
#' write.excel(b)
write.excel <- function(x, row.names = TRUE, col.names = TRUE, ...) {
  write.table(x, "clipboard", sep = "\t", row.names = row.names, 
              col.names = col.names, quote = FALSE, ...)
}



#' Returns statistics of regression caret models with and without 
#' application of applicability domain.
#'
#' @param model model trained with caret package.
#' @return data.frame which conatains model's parameters, values of R square 
#' and RMSE with and without application of applicability domain.
#' @details Applicability domain is estimated by fragment control approach. 
#' If value of a variable is outside the range of that variable for training set 
#' then the corresponding case is outside applicability domain and is removed 
#' while computing model statistics
#' @export
#' @examples
#' reg <- get.DA.reg(model)
get.DA.reg <- function(model) {
  
  frag_control <- function(train_df, test_df) {
    train_names <- colnames(train_df)[!sapply(train_df, function(i) all(i == 0))]
    which(apply(test_df, 1, function(i) {
      n <- names(i)[i != 0]
      length(setdiff(n, train_names)) == 0
    }))
  }
  
  # add columns fold and rep for predicted data.frame
  df <- model$pred
  par_names <- colnames(model$bestTune)
  df <- cbind(df, do.call(rbind, strsplit(df$Resample, "\\.")))
  colnames(df)[(ncol(df) - 1) : ncol(df)] <- c("fold", "rep")
  
  # calc stat for each rep
  tmp <- split(df, df[,c(par_names, "rep")])
  res <- as.data.frame(t(sapply(tmp, function(tmp_df) {
    c("R2"=pfpp::R2test(tmp_df$pred, tmp_df$obs), "RMSE"=pfpp::rmse(tmp_df$obs, tmp_df$pred))
  })))
  
  # add column with parameters names to the result
  res <- cbind(res, do.call(rbind, strsplit(rownames(res), "\\.")))
  colnames(res)[(ncol(res) - length(par_names)):ncol(res)] <- c(par_names, "rep")
  
  # find ids in DA
  ids_da <- lapply(model$control$indexOut, function(case_ids) {
    ids <- frag_control(model$trainingData[-case_ids,], model$trainingData[case_ids,])
    case_ids[ids]
  })
  # hack, because $indexOut has names like (Resample01, ...) while $index (Fold1.Rep01, ...)
  names(ids_da) <- names(model$control$index)
  ids_da <- split(ids_da, list(sub("^.*(Rep[0-9]*)$", "\\1", names(ids_da))))
  ids_da <- lapply(ids_da, unlist)
  
  # calc stat with DA
  tmp <- split(df, df[,c(par_names, "rep")])
  res_da <- as.data.frame(t(sapply(tmp, function(tmp_df) {
    ids <- tmp_df$rowIndex %in% ids_da[[as.character(unique(tmp_df$rep))]]
    c("R2_DA" = pfpp::R2test(tmp_df$pred[ids], tmp_df$obs[ids], mean(tmp_df$obs)), 
      "RMSE_DA" = pfpp::rmse(tmp_df$obs[ids], tmp_df$pred[ids]),
      "DA_coverage" = sum(ids) / length(ids))
  })))
  
  # add column with parameters names to the result
  res_da <- cbind(res_da, do.call(rbind, strsplit(rownames(res_da), "\\.")))
  colnames(res_da)[(ncol(res_da) - length(par_names)):ncol(res_da)] <- c(par_names, "rep")
  
  output <- cbind(aggregate(res[,1:2], list(res[, par_names]), mean), 
                  aggregate(res[,1:2], list(res[, par_names]), sd),
                  aggregate(res_da[,1:2], list(res_da[, par_names]), mean),
                  aggregate(res_da[,1:2], list(res_da[, par_names]), sd),
                  aggregate(res_da[,3], list(res_da[, par_names]), mean))
  
  colnames(output)[1:length(par_names)] <- par_names
  output <- output[,!grepl("Group\\.", colnames(output))]
  colnames(output)[(length(par_names) + 1):ncol(output)] <- 
    c("R2", "RMSE", "R2_SD", "RMSE_SD", "R2_DA", "RMSE_DA", "R2_DA_SD", 
      "RMSE_DA_SD", "Coverage")
  
  output[,1:length(par_names)] <- sapply(output[,1:length(par_names)], function(i) as.numeric(as.character(i)))
  
  return(output)
}



#' Save text data as binary to save used encoding.
#'
#' @param df is `data.frame` or `matrix` to save.
#' @param sep separator of fields.
#' @param row.names write row names.
#' @param col.names write column names.
#' @return nothing.
#' @details This function was created as a workaround for Windows users to properly 
#' save text files with UTF-8 encoding.
#' @export
write.as.binary <- function(df, file.name, sep = "\t", row.names = TRUE, col.names = TRUE) {
  output <- character()
  if (col.names) {
    output[1] <- paste0(colnames(df), collapse = sep)
    if (row.names) {
      output[1] <- paste0("rownames", sep, output[1])
    }
  }
  tmp <- apply(df, 1, function(i) paste0(i, collapse = "\t"))
  if (row.names) {
    tmp <- sapply(1:nrow(df), function(i) paste0(rownames(df)[i], sep, tmp[i]))
  }
  output <- c(output, tmp)
  writeLines(output, file.name, useBytes = TRUE)
}



#' Read sirms descriptor saved in svm format.
#'
#' @param file_name path to filename to load/
#' @return data.frame.
#' @details there should be files with colnames and rownames in the same dir 
#' along with target file_name. They have to have the same name and extensions 
#' .colnames and .rownames.
#' @export
sirms.read.svm <- function(file_name) {
  # read main file
  x <- readLines(file_name)
  x <- strsplit(x, " |:")
  x <- sapply(seq_along(x), function(i) {
    tmp <- cbind(i, matrix(x[[i]], ncol = 2, byrow = TRUE))
    class(tmp) <- "numeric"
    tmp[, 2] <- tmp[, 2] + 1
    tmp
  })
  x <- do.call(rbind, x)
  # read colnames and rownames
  row_names <- as.character(read.table(paste0(tools::file_path_sans_ext(file_name), ".rownames"))[, 1])
  col_names <- as.character(read.table(paste0(tools::file_path_sans_ext(file_name), ".colnames"))[, 1])
  # create output data.frame
  df <- matrix(0, nrow = length(row_names), ncol = length(col_names))
  df[x[,c(1,2)]] <- x[, 3]
  df <- as.data.frame(df)
  colnames(df) <- col_names
  rownames(df) <- row_names
  return(df)
}



#' Read data from text file with fragments contributions.
#'
#' @param file_name name of input text file with fragments contributions.
#' @param addN logical whether or not to add number of occurencies of each fragments 
#' across the dataset.
#' @param addM logical whether or not to add number of molecules comprising 
#' each type of fragment.
#' @param sep separator between molecule and fragment name
#' @return data.frame.
#' @details Old version of input file had '#' as a separator.
#' @export
#' @examples
#' df <- spci.load_data(file_name)
spci.load_data <- function(file_name, addN = TRUE, addM=TRUE, sep="###") {
  
  # read data
  df <- as.data.frame(data.table::fread(file_name, sep = "\t", header = T, na.strings = "", autostart = 3))
  rnames <- colnames(df)[-1]
  colnames(df) <- NULL
  rownames(df) <- df[,1]
  df <- as.data.frame(t(df[,-1]))
  df <- cbind(do.call(rbind, strsplit(rnames, sep)), df)
  colnames(df)[1:2] <- c("MolID", "FragID")
  df[,1:2] <- sapply(df[,1:2], as.character)
  
  # add full names
  full_names <- as.data.frame(t(sapply(unique(df[,"FragID"]), function(frag_name) {
    if (addN) {
      N <- sum(df[,"FragID"] == frag_name)
    }
    if (addM) {
      M <- length(unique(df[df[,"FragID"] == frag_name, "MolID"]))
    }
    if (addN & addM) {
      return(c(paste0(frag_name, " (M = ", M, ", N = ", N, ")"), M, N))
    } else {
      if (addN) {
        return(c(paste0(frag_name, " (N = ", N, ")"), M, N))
      } else {
        if (addM) {
          return(c(paste0(frag_name, " (M = ", M, ")"), M, N))
        }
      }
    }
  })))
  full_names[, 1] <- as.character(full_names[, 1])
  full_names[, 2:ncol(full_names)] <- sapply(full_names[, 2:ncol(full_names)], function(i) as.integer(as.character(i)))
  colnames(full_names) <- c("full_names", "M", "N")
  df <- cbind(full_names[df[,"FragID"], ], df)
  df
}



#' Add consensus (average) predictions to the data.frame with fragments contributions.
#'
#' @param df input data.frame.
#' @param model_names vector of model names.
#' @return data.frame.
#' @details Model names have to correspond to column names of input data.frame. 
#' These column names are usually consist of two parts model name of atomic property name 
#' (svm_charge or rf_overall).
#' @export
#' @examples
#' df <- spci.add_consensus(df, c("svm", "rf", "gbm", "pls"))
spci.add_consensus <- function(df, model_names) {
  model_names <- setdiff(model_names, "consensus")
  regexp <- paste0("^", model_names, "_", collapse = "|")
  selected_colnames <- grep(regexp, colnames(df), value = TRUE)
  prop_names <- unique(sub(".*_(.*)", "\\1", selected_colnames))
  consensus <- sapply(prop_names, function(prop_name) {
    rowMeans(df[, grep(prop_name, selected_colnames, value = TRUE), drop = FALSE])
  })
  df[paste0("consensus_", colnames(consensus))] <- consensus
  df
}



#' Filter out fragments whose occurence is below specified limit and keep 
#' only specified models.
#'
#' @param df input data.frame.
#' @param minM number of the lowest possible number of molecules comprising each 
#' fragment to keep them.
#' @param minN number of the lowest possible number of occurencies of fragmnets 
#' across the whole dataset.
#' @param model_names vector with model names to use for filtering
#' @return data.frame.
#' @details Model names have to correspond to column names of input data.frame. 
#' These column names are usually consist of two parts model name of atomic property name 
#' (svm_charge or rf_overall).
#' First five columns should contain full_names, M, N, MolId and FragID and they 
#' will be always kept.
#' @export
#' @examples
#' df <- spci.filter_data(df, 10, 10, c("svm", "rf", "gbm", "pls", "consensus"))
spci.filter_data <- function(df, minM, minN, model_names) {
  regexp <- paste0("^", model_names, "_", collapse = "|")
  selected_model_columns <- grep(regexp, colnames(df), value = TRUE)
  selected_columns <- c(colnames(df)[1:5], selected_model_columns)
  df[df$N >= minN & df$M >= minM, selected_columns]
}



#' Reorder factor levels representing fragment names according to specified 
#' value and metric.
#'
#' @param df input data.frame.
#' @param order_by_var name of the column which will be used for reordering
#' @param FUN function used for summarization
#' @return data.frame.
#' @details all values of the specified column will be grouped by fragment names and 
#' summarized with specified function which output will used for ordering (descending).
#' @export
#' @examples
#' df <- spci.reorder_data(df, "consensus_overall")
spci.reorder_data <- function(df, order_by_var, FUN=median) {
  get_col <- function(x,y) {x[,if(is.name(substitute(y))) deparse(substitute(y)) else y, drop = FALSE][[1]]}
  if (order_by_var == "_")
    return(df)
  d <- df %>%
    group_by(full_names)
  sorted_levels <- eval(substitute(summarise(d, res = FUN(var)), list(var = as.name(order_by_var)))) %>%
    arrange(res)
  df$full_names <- factor(df$full_names, levels = get_col(sorted_levels, 1))
  return(df)
}



#' Melt data.frames with fragments for graphical representation .
#'
#' @param df input data.frame.
#' @return melted data.frame.
#' @details First five columns should contain full_names, M, N, MolId and FragID.
#' @export
#' @examples
#' df <- spci.melt_data(df, "consensus_overall")
spci.melt_data <- function(df) {
  # init contribution names
  contrib_names <- c("overall", "hydrophobic", "electrostatic", "hydrogen bonding", "dispersive")
  names(contrib_names) <- c("overall", "logp", "charge", "hb", "refractivity")
  mdf <- reshape2::melt(df, id.vars=1:5)
  mdf <- cbind(mdf[, 1:5], 
               do.call(rbind, strsplit(as.character(mdf$variable), "_")), 
               mdf$value)
  colnames(mdf)[6:8] <- c("Model", "Property", "Contribution")
  mdf[, c("Model", "Property")] <- sapply(mdf[, c("Model", "Property")], as.character)
  mdf[, "Property"] <- contrib_names[mdf[, "Property"]]
  return(mdf)
}
