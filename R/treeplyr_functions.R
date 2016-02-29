#' Function for making an object of class 'treedata'
#' 
#' This function generates an object of class 'treedata' that ensures that the ordering of tip labels and data remain intact. The object can be manipulated using \code{dplyr} functions.
#' 
#' @param tree An object of class 'phylo'
#' @param data A data frame or matrix
#' @param name_column The column of \code{data} that contains the names to be matched to the tree. By default, it is set to "detect" which finds the column with the most matches to the tree (including the rownames).
#' @return An object of class "\code{treedata}". The tree is pruned of tips not represented in the data, and the data is filtered for taxa not in the tree. The data is returned as a data frame tble that is compatible with \code{dplyr} functions. 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' @export
make.treedata <- function(tree, data, name_column="detect") {
  if(class(tree)!="phylo") stop("tree must be of class 'phylo'")
  if(is.vector(data)){
    data <- as.matrix(data)
    colnames(data) <- "trait"
  }
  if(is.null(colnames(data))){
    colnames(data) <- paste("trait", 1:ncol(data), sep="")
  }
  coln <- colnames(data)
  if(name_column=="detect"){
    if(is.null(rownames(data))){
      tmp.df <- data.frame(data)
      offset <- 0
    } else {
      tmp.df <- data.frame(rownames(data), data)
      offset <- 1
    }
    matches <- sapply(tmp.df, function(x) sum(x %in% tree$tip.label))
    if(all(matches==0)) stop("No matching names found between data and tree")
    name_column <- which(matches==max(matches))-offset
  }
  dat <- tbl_df(as.data.frame(lapply(1:ncol(data), function(x) type.convert(apply(data[,x, drop=FALSE], 1, as.character)))))
  #dat <- tbl_df(as.data.frame(lapply(1:ncol(data), function(x) type.convert(as.character(data[,x])))))
  #dat <- data
  colnames(dat) <- coln
  #dat <- apply(dat, 2, type.convert)
  if(name_column==0){
    clnm <- colnames(dat)
    dat <- dat[,clnm, drop=FALSE]
    dat.label <- as.character(rownames(data))
  } else {
    if(is.numeric(name_column)){
      clnm <- (1:ncol(data))[-name_column]
    } else {
      clnm <- colnames(dat)[-which(colnames(dat)==name_column)]
    }
    dat <- dat[, clnm, drop=FALSE]   
    dat.label <- as.character(as.data.frame(data)[[name_column]])
  }
  data_not_tree <- setdiff(dat.label, tree$tip.label)
  tree_not_data <- setdiff(tree$tip.label, dat.label)
  phy <- drop.tip(tree, tree_not_data)
  dat <- filter(dat, dat.label %in% phy$tip.label)
  dat.label <- dat.label[dat.label %in% phy$tip.label]
  if(any(duplicated(dat.label))){
    warning("Duplicated data in dataset, selecting first unique entry for each species")
    dat <- filter(dat, !duplicated(dat.label))
    dat.label <- dat.label[!duplicated(dat.label)]
  }
  o <- match(dat.label, phy$tip.label)
  dat <- arrange(dat, o)
  td <- list(phy=phy, dat=dat)
  class(td) <- c("treedata", "list")
  attributes(td)$tip.label <- phy$tip.label
  attributes(td$dat)$row.names <- phy$tip.label
  attributes(td)$dropped <- list(dropped_from_tree=data_not_tree,dropped_from_data=tree_not_data)
  rownames(td$dat) <- attributes(td)$tip.label
  return(td)
}

#' Function for mutating an object of class 'treedata'
#' 
#' This function can be used to add new variables to a treedata object; see \code{\link{mutate}}.
#' 
#' @param tdObject A "\code{treedata}" object
#' @param ... Arguments to mutate the treedata object
#' @return An object of class "\code{treedata}" with new data added. 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat[,-(3:5)], name_column=1)
#' tdmutate <- mutate(td, anolis$dat[,3])
#' @export
mutate.treedata <- function(tdObject, ...){
  if(is.null(list(substitute(...))[[1]])) stop("No expressions provided to add to the treedata object")
  dat <- mutate_impl.dplyr(tdObject$dat, named_dots.dplyr(...))
  tdObject$dat <- dat
  rownames(tdObject$dat) <- attributes(tdObject)$tip.label
  return(tdObject)
}

#' Function for selecting columns from an object of class 'treedata'
#' 
#' This function can be used to select a subset of variables (columns) from a treedata object; see \code{\link{select}}.
#' 
#' @param tdObject A "\code{treedata}" object
#' @param ... Additional arguments to select function
#' @return An object of class "\code{treedata}" with specified variables selected. 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' tdselect <- select(td, SVL)
#' @export
select.treedata <- function(tdObject, ...){
  if(is.null(list(substitute(...))[[1]]))  stop("No criteria provided for selection")
  vars <- select_vars(names(tdObject$dat), ..., env = parent.frame())
  dat <- select_impl.dplyr(tdObject$dat, vars)
  #dat <- lapply(dat, function(x){names(x) <- tdObject$dat[,1]; x})
  tdObject$dat <- dat
  rownames(tdObject$dat) <- attributes(tdObject)$tip.label
  return(tdObject)
}

#' Function for filtering rows from an object of class 'treedata'
#' 
#' This function can be used to select a subset of species (rows) from a treedata object; see \code{\link{filter}}.
#' 
#' @param tdObject A "\code{treedata}" object
#' @param ... Arguments used to filter the data.
#' @return An object of class "\code{treedata}" with specified species selected. 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' tdfilter <- filter(td, island=="Cuba")
#' @export
filter.treedata <- function(tdObject, ...){
  #if(is.null(list(substitute(...))[[1]]))  stop("No criteria provided for filtering")
  tdObject$dat <- mutate(tdObject$dat, tip.label=attributes(tdObject)$tip.label)
  tdObject$dat <- filter(tdObject$dat, ...)
  attributes(tdObject)$tip.label <- tdObject$dat$tip.label
  tdObject$dat <- select(tdObject$dat, 1:(ncol(tdObject$dat)-1))
  rownames(tdObject$dat) <- attributes(tdObject)$tip.label
  tdObject$phy <- drop.tip(tdObject$phy, tdObject$phy$tip.label[!(tdObject$phy$tip.label %in% attributes(tdObject)$tip.label)])
  return(tdObject)
}

#' @name summarise.treedata
#' @aliases summarize.treedata
#' @title Function for summarizing an object of class 'treedata'
#' 
#' This function can be used to summarize a treedata object. BUT IT DOES NOT WORK
#' 
#' @param tdObject A "\code{treedata}" object
#' @param ... additional arguments to summarize
#' @return A summary of the treedata object 
#' @examples 
#' data(anolis)
#' @export
summarize.treedata <- function(tdObject, ...){
  if(is.null(list(substitute(...))[[1]]))  stop("No expression provided to summarize data")
  phy <- tdObject$phy
  dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
  res <- summarise_(tdObject$dat, .dots = dots)
  return(res)
}

#' @export
summarise.treedata <- function(tdObject, ...){
  if(is.null(list(substitute(...))[[1]]))  stop("No expression provided to summarize data")
  phy <- tdObject$phy
  dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
  res <- summarise_(tdObject$dat, .dots = dots)
  return(res)
}

#' @export
summarise_.treedata <- function(tdObject, ..., .dots){
  #if(is.null(list(substitute(...))[[1]])) stop("No expression provided to summarize data")
  dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
  env <- new.env(parent = parent.frame(), size = 1L)
  env$phy <- tdObject$phy
  env$dat <- as.data.frame(tdObject$dat)
  env$tip.label <- tdObject$phy$tip.label
  rownames(env$dat) <- env$phy$tip.label
  for(i in 1:length(dots)){
    dots[[i]]$env <- env
  }
  res <- summarise_(tdObject$dat, .dots = dots)
  return(res)
}

summarise_.grouped_treedata <- function(tdObject, ..., .dots){
  #if(is.null(list(substitute(...))[[1]])) stop("No expression provided to summarize data")
  dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
  nind <- (1:nrow(tdObject$dat))
  group_levels <- attributes(tdObject$dat)$labels[[1]]
  group_by_var <- colnames(attributes(tdObject$dat)$labels)[1]
  phys <- lapply(attributes(tdObject$dat)$indices, function(x) drop.tip(tdObject$phy, nind[!(nind %in% (x+1))]))
  tds <- lapply(group_levels, function(x) filter_(tdObject, paste(group_by_var, "=='", x,"'", sep="")))
  envs <- lapply(1:length(group_levels), function(x){e <- new.env(parent=parent.frame(), size=1L);
                                                      e$phy <- phys[[x]];
                                                      e$dat <- tds[[x]]$dat;
                                                      e})
  OUT <- NULL
  for(j in 1:length(group_levels)){
    for(i in 1:length(dots)){
      dots[[i]]$env <- envs[[j]]
    }
    res <- summarise_(tds[[j]]$dat, .dots = dots)
    OUT <- rbind(OUT, res)
  }  
  return(OUT)
}

#' @export
group_by_.treedata <- function(tdObject, ..., .dots, add=FALSE){
  groups <- dplyr:::group_by_prepare(tdObject$dat, ..., .dots = .dots, add = add)
  dat <- grouped_df(groups$data, groups$groups)
  tdObject$dat <- dat
  class(tdObject) <- c("grouped_treedata", "treedata", "list")
  return(tdObject)
}

#' @rdname reorder.treedata
#' @export
reorder <- function(tdObject, ...){
  UseMethod("reorder")
}

#' Reorder a 'treedata' object
#' @description Reorders a 'treedata' object. Both the tips and the data are automatically reordered to match.
#' @param tdObject An object of class 'treedata'
#' @param order Method for reordering
#' @param index.only XXX
#' @param ... Additional arguments to reorder
#' @return An object of class 'treedata'
#' 
#' @examples
#' x <- 10
#' @export
reorder.treedata <- function(tdObject, order="postorder", index.only=FALSE, ...){
  dat.attr <- attributes(tdObject$dat)
  phy <- reorder(tdObject$phy, order, index.only)#, ...)
  index <- match(tdObject$phy$tip.label, phy$tip.label)
  tdObject$dat <- tdObject$dat[index,]
  attributes(tdObject$dat) <-dat.attr
  attributes(tdObject$dat)$row.names <- attributes(tdObject)$tip.label <- phy$tip.label
  tdObject$phy <- phy
  if(index.only){
    return(index)
  } else {
    return(tdObject)
  }
}

#' @rdname treeply.treedata
#' @export
treeply <- function(tdObject, ...){
  UseMethod("treeply")
}

#' Run a function on the phylogeny of a 'treedata' object
#' @description Applies a function to the phylogeny in a 'treedata' object. If the order of tips are changed, or if tips are dropped, then the data are automatically reordered to match the tree.
#' @param tdObject An object of class 'treedata'
#' @param FUN A function that operates on an object of class 'phylo'
#' @param ... Additional arguments

#' @return An object of class 'treedata'
#' 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' td_OU <- treeply(td, geiger::rescale.phylo, model="OU", 10)
#'   
#' par(mfrow=c(1,2))
#' plot(td$phy)
#' plot(td_OU$phy)
#' 
#' @export
treeply.treedata <- function(tdObject, FUN, ...){
  if(!class(tdObject)[1]=="treedata") stop("Object is not of class 'treedata'")
  FUN <- match.fun(FUN)
  out_FUN <- FUN(tdObject$phy, ...)
  if(class(out_FUN)[1] == "phylo"){
    tdObject$phy <- out_FUN
    tdObject$dat <- tdObject$dat[match(tdObject$phy$tip.label, attributes(tdObject)$tip.label),]
    attributes(tdObject)$tip.label <- tdObject$phy$tip.label
    return(tdObject)
  } else {
    warning("Function output was not of class 'phylo', returning output only")
    return(out_FUN)
  } 
}

#' Apply a function over all treedata object columns and return a list of results, analogously to the 
#' normal apply function
#' 
#' @param tdObject A treedata object
#' @param MARGIN the margin over which the data is applied (e.g. 1 = rows, 2 = columns)
#' @param FUN A function to apply over the data frame
#' @param ... Additional parameters passed on to FUN
#' 
#' @details Note that if the parameter "phy" is specified in the additional parameters (i.e. '...'), then 
#' it will be substituted with the treedata object $phy. 
#' 
#' @export
tdapply <- function(tdObject, MARGIN, FUN, ...){
  if(!class(tdObject)[1]=="treedata") stop("Object is not of class 'treedata'")
  FUN <- match.fun(FUN)
  phy <- tdObject$phy
  env <- new.env(parent=parent.frame(), size=1L)
  env$phy <- tdObject$phy
  res <- eval(substitute(apply(tdObject$dat, MARGIN, FUN, ...)), env)
  return(res)
}


#' @rdname treedply.treedata
#' @export
treedply <- function(tdObject, ...){
  UseMethod("treedply")
}

#' Run a function on a 'treedata' object
#' 
#' @param tdObject A treedata object
#' @param ... A function call.
#' 
#' @details This function allows arbitrary R functions that use trees and data to be run on a 'treedata' objects. 
#' 
#' @return Function output
#' 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat)
#' treedply(td, geiger::fitContinuous(phy, getVector(dat, SVL), model="BM", ncores=1))
#'  
#' treedply(td, phytools::phylosig(phy, getVector(dat, awesomeness), "lambda", test=TRUE))
#' treedply(td, phytools::phenogram(phy, getVector(dat, SVL), ftype="off", spread.labels=FALSE))
#' @export
treedply.treedata <- function(tdObject, ...){
  if(!is.call(substitute(...))){
    call <- list(...)[[1]]
  } else {
    call <- substitute(...)
  }
  env <- new.env(parent = parent.frame(), size = 1L)
  env$phy <- tdObject$phy
  env$dat <- tdObject$dat
  out <- eval(call, env)
  if(is.null(out)){
    invisible()
  } else {
    return(out)
  }
}

#' A function for returning a named vector from a data frame or matrix with row names
#'
#' @param dat A data frame or matrix with row names
#' @param ... The name of the column to select
#' 
#' @return A named vector
#' @export
getVector <- function(dat, ...){
  args <- as.character(substitute(list(...)))[-1L]
  arg_sub <- type.convert(args)
  if(is.numeric(arg_sub) | is.integer(arg_sub)) args <- arg_sub
  vecs <- lapply(args,function(x) dat[[x]])
  vecs <- lapply(vecs, function(x) setNames(x, attributes(dat)$row.names))
  if(length(vecs)==1){
    vecs = vecs[[1]]
  } else {names(vecs) <- args}
  return(vecs)
}

#' @export
getVector_ <- function(dat, ...){
  args <- eval(...)
  arg_sub <- type.convert(args)
  if(is.numeric(arg_sub) | is.integer(arg_sub)) args <- arg_sub
  vecs <- lapply(args,function(x) dat[[x]])
  vecs <- lapply(vecs, function(x) setNames(x, attributes(dat)$row.names))
  if(length(vecs)==1){
    vecs = vecs[[1]]
  } else {names(vecs) <- args}
  return(vecs)
}


#' @export
print.treedata <- function(x, ...){
  cat("$phy \n")
  print(x$phy)
  
  cat("\n$dat \n")
  print(x$dat)
}

#' @export
summary.treedata <- function(x, ...){
 
  cat('A treeplyr treedata object', "\n")
  cat(paste('The dataset contains ', ncol(x$dat), ' traits'), "\n")
  types <- setNames(suppressWarnings(detectAllCharacters(as.matrix(x$dat))), colnames(x$dat))
  cat("Continuous traits: ", names(types)[which(types=="continuous")], "\n")
  cat("Discrete traits: ", names(types)[which(types=="discrete")], "\n")
  cat(paste("The following traits have missing values:", paste(names(types)[apply(x$dat, 2, function(y) any(is.na(y)))], collapse=", "), "\n"))
  cat(paste("These taxa were dropped from the tree:", paste(attributes(x)$dropped$dropped_from_tree, collapse=", "), "\n"))
  cat(paste("These taxa were dropped from the data:", paste(attributes(x)$dropped$dropped_from_data, collapse=", "), "\n"))

  
  cat("$phy \n")
  print(x$phy)
  cat("\n$dat \n")
  print(x$dat)
}



#' Function for checking whether a treedata object contains only numeric columns and for forcing data columns into numeric format
#' 
#' This function can be used to check if a treedata object contains numeric columns and, if desired, drop all non-numeric columns.
#' 
#' @param tdObject A "\code{treedata}" object
#' @param return.numeric If TRUE, then a treedata object with all numeric columns will be returned; non-numeric columns will be removed.
#' @return If return.numeric, then an object of class "\code{treedata}" with only numeric columns. 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' tdnumeric <- forceNumeric(td)
#' @export
forceNumeric <- function(tdObject, return.numeric=TRUE) {
  valid <- which(sapply(tdObject$dat, is.numeric))
  if(length(valid) < ncol(tdObject$dat)){
    if(length(valid)==0){
      warning("Dataset does not contain any numeric data") 
      tdObject <- select(tdObject)
    } else {
      not.valid <- colnames(tdObject$dat)[which(!sapply(tdObject$dat, is.numeric))]
      warning(paste("Not all data continuous, dropping non-numeric data columns:", paste(not.valid, collapse=" ")))
      tdObject <- select(tdObject, valid)
    }
  }
  if(return.numeric){
    return(tdObject)
  } else {
    names(tdObject$dat)
  }
}

#' Function for checking whether a treedata object contains only factors and for forcing data columns into factor format
#' 
#' This function can be used to check if a treedata object contains factors and, if desired, convert all columns automatically to factors.
#' 
#' @param tdObject A "\code{treedata}" object
#' @param return.factor If TRUE, then a treedata object with all factors will be returned; columns will be forced into factors using \code{factor} and any with no repeated elements will be removed.
#' @return If return.factor, then an object of class "\code{treedata}" with all columns as factors. 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' tdforcefactor <- forceFactor(td)
#' @export
forceFactor <- function(tdObject, return.factor=TRUE) {
  classes <- sapply(tdObject$dat, class)
  valid <- which(classes=="factor")
  if(length(valid) < ncol(tdObject$dat)){
    #Which data are numeric
    are.numeric <- which(classes %in% c("numeric", "integer"))
    if(length(are.numeric) > 0){
      warning("Data contain numeric entries, which will be converted to factors")
      #convert them to factors
      tdObject$dat[, are.numeric] <- lapply(as.data.frame(tdObject$dat[,are.numeric]), factor)
      ##Check to see if converted data has any columns that appear continuous
      classes <- sapply(tdObject$dat, class)
      valid <- which(classes=="factor")
      too.many.levels <- which(sapply(tdObject$dat,function(x) length(unique(x)))==nrow(tdObject$dat))
      if(length(too.many.levels) > 0){
        warning(paste("Conversion failed for data columns", paste(colnames(tdObject$dat)[too.many.levels], collapse=" "), "as these data have no shared states. These data will be removed", sep=" "))
        valid <- valid[which(!(valid %in% too.many.levels))]
      }
    }

    if(length(valid)==0){
      warning("Data does not contain any discrete data") 
      tdObject <- select(tdObject)
    } else {
      tdObject$dat <- select(tdObject$dat, valid)
    }
  }
  if(return.factor){
    return(tdObject)
  } else {
    names(tdObject$dat)
  }
}

# Helper functions that allow string arguments for  dplyr's data modification functions like arrange, select etc. 
# Author: Sebastian Kranz

# Examples are below

#' Modified version of dplyr's filter that uses string arguments
#' @param .data dplyr object
#' @param ... Additional arguments
#' @export
s_filter = function(.data, ...) {
  eval.string.dplyr(.data,"filter", ...)
}

#' Modified version of dplyr's select that uses string arguments
#' @param .data dplyr object
#' @param ... Additional arguments
#' @export
s_select = function(.data, ...) {
  eval.string.dplyr(.data,"select", ...)
}

#' Modified version of dplyr's arrange that uses string arguments
#' @param .data dplyr object
#' @param ... Additional arguments
#' @export
s_arrange = function(.data, ...) {
  eval.string.dplyr(.data,"arrange", ...)
}

#' Modified version of dplyr's arrange that uses string arguments
#' @param .data dplyr object
#' @param ... Additional arguments
#' @export
s_mutate = function(.data, ...) {
  eval.string.dplyr(.data,"mutate", ...)
}

#' Modified version of dplyr's summarise that uses string arguments
#' @param .data dplyr object
#' @param ... Additional arguments
#' @export
s_summarise = function(.data, ...) {
  eval.string.dplyr(.data,"summarise", ...)
}

#' Modified version of dplyr's group_by that uses string arguments
#' @param .data dplyr object
#' @param ... Additional arguments
#' @export
s_group_by = function(.data, ...) {
  eval.string.dplyr(.data,"group_by", ...)
}

#' Modified version of dplyr's group_by that uses string arguments
#' @param .data dplyr object
#' @param ... Additional arguments
#' @export
s_treeply = function(.data, ...) {
  eval.string.dplyr(.data,"treeply", ...)
}

#' Modified version of dplyr's group_by that uses string arguments
#' @param .data dplyr object
#' @param ... Additional arguments
#' @export
s_treedply = function(.data, ...) {
  eval.string.dplyr(.data,"treedply", ...)
}

# Internal function used by s_filter, s_select etc.
eval.string.dplyr = function(.data, .fun.name, ...) {
  args = list(...)
  args = unlist(args)
  code = paste0(.fun.name,"(.data,", paste0(args, collapse=","), ")")
  df = eval(parse(text=code,srcfile=NULL))
  df  
}

#' @export
select_.treedata <- function(.data, ..., .dots){
  dots <- all_dots(.dots, ...)
  vars <- select_vars_(names(.data$dat), dots)
  dat <- .data$dat[, vars, drop = FALSE]
  row.names(dat) <- attributes(.data)$tip.label
  .data$dat <- dat
  return(.data)
}

#' @export
mutate_.treedata <- function(.data, ..., .dots){
  dots <- all_dots(.dots, ..., all_named=TRUE)
  dat <- mutate_impl.dplyr(.data$dat, dots)
  row.names(dat) <- attributes(.data)$tip.label
  .data$dat <- dat
  return(.data)
}

#' @export
mutate_.grouped_treedata <- function(.data, ..., .dots){
  dots <- all_dots(.dots, ..., all_named=TRUE)
  dat <- mutate_impl.dplyr(.data$dat, dots)
  row.names(dat) <- .data$phy$tip.label
  .data$dat <- dat
  return(.data)
}

#' @export
filter_.treedata <- function(.data, ..., .dots){
  dots <- all_dots(.dots, ..., all_named=TRUE)
  .data$dat <- mutate(.data$dat, tip.label=attributes(.data)$tip.label)
  dat <- filter_impl.dplyr(.data$dat, dots)
  .data$dat <- dat
  attributes(.data)$tip.label <- .data$dat$tip.label
  .data$dat <- select(.data$dat, 1:(ncol(.data$dat)-1))
  .data$phy <- drop.tip(.data$phy, .data$phy$tip.label[!(.data$phy$tip.label %in% attributes(.data)$tip.label)])
  attributes(.data$dat)$row.names <- .data$phy$tip.label
  return(.data)
}

#' @export
filter_.grouped_treedata <- function(.data, ..., .dots){
  dots <- all_dots(.dots, ..., all_named=TRUE)
  cl <- class(.data$dat)
  .data$dat$tip.label <- .data$phy$tip.label
  dat <- filter_impl.dplyr(.data$dat, dots)
  .data$dat <- dat
  attributes(.data)$tip.label <- .data$dat$tip.label
  .data$dat <- select(.data$dat, 1:(ncol(.data$dat)-1))
  class(.data$dat) <- cl
  .data$phy <- drop.tip(.data$phy, .data$phy$tip.label[!(.data$phy$tip.label %in% attributes(.data)$tip.label)])
  return(.data)
}

select_impl.dplyr <- function (df, vars) 
{
    .Call("dplyr_select_impl", PACKAGE = "dplyr", df, vars)
}

filter_impl.dplyr <- function (df, dots) 
{
    .Call("dplyr_filter_impl", PACKAGE = "dplyr", df, dots)
}

mutate_impl.dplyr <- function (df, dots) 
{
    .Call("dplyr_mutate_impl", PACKAGE = "dplyr", df, dots)
}

named_dots.dplyr <- function (...) 
{
    auto_name.dplyr(dots.dplyr(...))
}

dots.dplyr <- function (...) 
{
  eval(substitute(alist(...)))
}

auto_name.dplyr <- function (x) 
{
  names(x) <- auto_names.dplyr(x)
  x
}

auto_names.dplyr <- function (x) 
{
  nms <- names2.dplyr(x)
  missing <- nms == ""
  if (all(!missing)) 
    return(nms)
  deparse2 <- function(x) paste(deparse(x, 500L), collapse = "")
  defaults <- vapply(x[missing], deparse2, character(1), USE.NAMES = FALSE)
  nms[missing] <- defaults
  nms
}

names2.dplyr <- function (x) 
{
  names(x) %||% rep("", length(x))
}

#' Add regimes to a treedata object
#' @export
paint_clades <- function(tdObject, nclades=1, name="clades", interactive=TRUE, type="nodes", ids=NULL, plot=TRUE, ...){
  if(interactive){
    regimes <- bayou::identifyBranches(tdObject$phy, nclades)
    cat("branches ", paste(regimes$sb, collapse=", "), "selected\n")
    cols <- rainbow(length(regimes$sb) + 1)
    L <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    legend(min(L$xx), max(L$yy), legend=1:(nclades+1), lwd=3, col=cols)
  } else 
    {
    if(type == "nodes"){
      if(any(!ids %in% (length(tdObject$phy$tip.label)+1):max(tdObject$phy$edge))) stop("Please specify valid internal node ids")
      nclades = length(ids)
      regimes <- list(sb = match(ids, tdObject$phy$edge[,2]), loc = rep(0, nclades))

    } 
    if(type == "branches"){
      if(any(!ids %in% (1:nrow(tdObject$phy$edge)))) stop("Please specify valid branch ids")
      nclades = length(ids)
      regimes <- list(sb = ids, loc = rep(0, nclades))
    }      
      if(plot){
        sb <- regimes$sb
        n <- nclades
        loc <- regimes$loc
        dum <- numeric(length(tdObject$phy$tip.label))
        names(dum) <- tdObject$phy$tip.label
        cache <- bayou:::.prepare.ou.univariate(tdObject$phy, dum)
        pars <- list(sb = sb, loc = loc, t2 = 1:n + 1)
        tr <- bayou:::.toSimmap(bayou:::.pars2map(pars, cache), cache)
        cols <- rainbow(length(sb) + 1)
        names(cols) <- 1:(length(sb) + 1)
        phytools:::plotSimmap(tr, pts = FALSE, fsize = 0.5, colors = cols)
        L <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        legend(min(L$xx), max(L$yy), legend=1:(nclades+1), lwd=3, col=cols)
      }
      
    }
  regimes$t2 <- 2:(length(regimes$sb)+1)
  regimes$k <- nclades
  regimes$ntheta <- nclades+1
  tdObject <- mutate(tdObject, bayou:::.tipregime(regimes, tdObject$phy))
  colnames(tdObject$dat)[ncol(tdObject$dat)] <- name
  return(tdObject)
}
