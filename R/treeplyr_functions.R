#' Function for making an object of class \code{treedata}
#'
#' This function generates an object of class \code{treedata} that ensures that the ordering of tip labels
#' and data remain intact. The object can be manipulated using \code{dplyr} functions.
#'
#' @param tree An object of class 'phylo'
#' @param data A data frame or matrix
#' @param name_column An optional argument that specifies the column of \code{data} that contains the
#' names to be matched to the tree. By default, it is set to "detect" which finds the column with the
#' most matches to the tree (including
#' the rownames).
#' @param as.is Whether convert to factors. When FALSE (default), convert character vectors to factors. 
#' @return An object of class "\code{treedata}". The tree is pruned of tips not represented in the data,
#' and the data is filtered for taxa not in the tree. The data is returned as a data frame tble that is
#' compatible with \code{dplyr} functions.
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat)
#' @export
make.treedata <- function(tree, data, name_column="detect", as.is=FALSE) {
  if(!inherits(tree, "phylo")) stop("tree must be of class 'phylo'")
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
  } else{
    if(is.character(name_column)){
      name_column <- which(name_column==coln)[1]
    }
  }
  dat <- as_tibble(as.data.frame(lapply(1:ncol(data), function(x)
    type.convert(apply(data[,x, drop=FALSE], 1, as.character), as.is=as.is))))
  colnames(dat) <- coln
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
  ...my.order... <- match(dat.label, phy$tip.label)
  dat <- arrange(dat, ...my.order...)
  td <- list(phy=phy, dat=dat)
  class(td) <- c("treedata", "list")
  attributes(td)$tip.label <- phy$tip.label
  #attributes(td$dat)$row.names <- phy$tip.label
  attributes(td)$dropped <- list(dropped_from_tree=data_not_tree,dropped_from_data=tree_not_data)
  #Removing rownames for future versions of dplyr
  #rownames(td$dat) <- attributes(td)$tip.label
  return(td)
}

#' Function for mutating an object of class \code{treedata}
#'
#' This function can be used to add new variables to a treedata object; see \code{\link{mutate}}.
#'
#' @aliases mutate.treedata mutate.grouped_treedata
#' @param .data An object of class \code{treedata}
#' @param ... Arguments to mutate the treedata object
#' @return An object of class \code{treedata} with new data added.
#' @seealso \code{\link{mutate}}
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat)
#' tdmutate <- mutate(td, lnSVL = log(SVL), badassery = awesomeness + hostility)
#' @export
mutate.treedata <- function(.data, ...){
  dat <- mutate(.data$dat, ...)
  .data$dat <- dat
  return(.data)
}


#' @title mutate_.treedata
#' @description Now defunct and replaced with mutate.treedata and mutate_.grouped_treedata.
#'
#' @name mutate_-defunct
#' @seealso \code{\link{treeplyr-defunct}}
#' @keywords internal
#' @export
mutate_.treedata <- function(.data, ...) {
  .Defunct(msg = "'mutate_' has been removed from this package; you can use mutate instead")
}


#' Choose rows by their ordinal position in the tbl for an object of class \code{treedata}
#'
#' This function can be used to drop tips from tree and data; see \code{\link{slice}}.
#'
#' @param .data An object of class \code{treedata}
#' @param ... Integer row values
#' @return An object of class \code{treedata}.
#' @seealso \code{\link{slice}}
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat)
#' tdslice <- slice(td, 1:5)
#' tdslice
#' @export
slice.treedata <- function(.data, ...){
  .data$dat$labelTEMP0123 <- .data$phy$tip.label
  dat <- slice(.data$dat, ...)
  .data <- make.treedata(.data$phy, dat)
  return(.data)
}



#' @title slice_.treedata
#' @description Now defunct and replaced with slice.treedata
#'
#' @name slice_-defunct
#' @seealso \code{\link{treeplyr-defunct}}
#' @keywords internal
#' @export
slice_.treedata <- function(.data, ...) {
  .Defunct(msg = "'slice_' has been removed from this package; you can use slice instead")
}

#' Function for selecting columns from an object of class \code{treedata}
#'
#' This function can be used to select a subset of variables (columns) from a treedata object;
#' see \code{\link{select}}.
#'
#' @aliases select.grouped_treedata
#' @param .data An object of class \code{treedata}
#' @param ... Additional arguments to select columns
#' @return An object of class \code{treedata} with specified variables selected.
#' @seealso \code{\link{select}}
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat)
#' tdselect <- select(td, SVL, awesomeness)
#' @export
select.treedata <- function(.data, ...){
  .data$dat <- select(.data$dat, ...)
  return(.data)
}


#' @title select_.treedata
#' @description Now defunct and replaced with select.treedata
#'
#' @name select_-defunct
#' @seealso \code{\link{treeplyr-defunct}}
#' @keywords internal
#' @export
select_.treedata <- function(.data, ...) {
  .Defunct(msg = "'select_' has been removed from this package; you can use select instead")
}


#' Function for filtering rows from an object of class \code{treedata}
#'
#' This function can be used to select a subset of species (rows) from a treedata object;
#' see \code{\link{filter}}.
#'
#' @aliases filter.grouped_treedata
#' @param .data An object of class \code{treedata}
#' @param ... Additional arguments to filter by
#' @return An object of class \code{treedata} with the dataset filtered by the specified criteria.
#' @seealso \code{\link{filter}}
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' tdfilter <- filter(td, island=="Cuba", SVL > 3.5)
#' @export
filter.treedata <- function(.data, ...){
  tip.labels <- attributes(.data)$tip.label
  .data$dat <- mutate(.data$dat, tip.label=tip.labels)
  dat <- filter(.data$dat, ...)
  .data$dat <- dat
  attributes(.data)$tip.label <- .data$dat$tip.label
  nc <- ncol(.data$dat)
  .data$dat <- select(.data$dat, 1:(nc-1))
  .data$phy <- drop.tip(.data$phy, .data$phy$tip.label[!(.data$phy$tip.label %in% attributes(.data)$tip.label)])
  return(.data)
}


#' @title filter_.treedata
#' @description Now defunct and replaced with filter.treedata
#'
#' @name filter_-defunct
#' @seealso \code{\link{treeplyr-defunct}}
#' @keywords internal
#' @export
filter_.treedata <- function(.data, ...) {
  .Defunct(msg = "'filter_' has been removed from this package; you can use select instead")
}

#' Function for summarizing an object of class \code{treedata}
#'
#' This function can be used to summarize a treedata object.
#'
#' @aliases summarize.treedata summarize.grouped_treedata summarise.grouped_treedata
#' @param .data An object of class \code{treedata}
#' @param ... Additional expressions by which to summarize data in the \code{treedata} object
#' @details Summarizing \code{treedata} objects allows expressions using the objects \code{phy}. The \code{treedata}
#' object can also be grouped, with summary statistics being applied to the pruned groups and phylogenies.
#' @return An object of class \code{tbl_df} with the requested summary data.
#' @seealso \code{\link{summarize}}, \code{\link{group_by}}
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat)
#' summarize(td, ntips = length(phy$tip.label), meanSVL = mean(SVL), sdSVL = sd(SVL))
#' tdGrouped <- group_by(td, ecomorph)
#' summarize(tdGrouped, ntips = length(phy$tip.label),
#'    totalBL = sum(phy$edge.length), meanSVL = mean(SVL), sdSVL = sd(SVL))
#' @export
summarise.treedata <- function(.data, ...){
  env <- new.env(parent = parent.frame(), size = 1L)
  env$phy <- .data$phy
  env$dat <- as.data.frame(.data$dat)
  env$tip.label <- .data$phy$tip.label
  dots <- quos(...)
  for(i in 1:length(dots)){
    attributes(dots[[i]])$.Environment <- env
  }
  summarise(env$dat, !!!dots)
}



#' @rdname summarise.treedata
#' @export
summarise.grouped_treedata <- function(.data, ...){
  #if(is.null(list(substitute(...))[[1]])) stop("No expression provided to summarize data")
  dots <- quos(...)
  #lazyeval::all_dots(.dots, ..., all_named = TRUE)
  nind <- (1:nrow(.data$dat))
  group_levels <- group_data(.data$dat)[[1]]
  group_by_var <- group_vars(.data$dat)
  group_index <- group_indices(.data$dat)
  phys <- lapply(group_levels, function(x) drop.tip(.data$phy, which(!(group_index %in% as.numeric(x)))))
  dat <- as.data.frame(.data$dat)
  rownames(dat) <- attributes(.data)$tip.label
  dats <- lapply(phys, function(x) make.treedata(x, dat)$dat)
  envs <- lapply(group_levels, function(x){e <- new.env(parent=parent.frame(), size=1L);
                                                      e$phy <- phys[[x]];
                                                      e$dat <- dats[[x]];
                                                      e})
  OUT <- NULL
  for(i in seq_along(envs)){
    edots <- dots
    for(j in 1:length(edots)){
      attributes(edots[[j]])$.Environment <- envs[[i]]
    }
   OUT <- rbind(OUT, summarise(envs[[i]]$dat, !!! edots))
  }
  OUT <- data.frame(group_levels, OUT)
  colnames(OUT)[1] <- group_by_var
  return(OUT)
}

#' Function for grouping an object of class \code{treedata}
#'
#' This function can be used to group a treedata object by some factor.
#'
#' @aliases group_by.treedata
#' @param .data An object of class \code{treedata}
#' @param ... The name of the grouping factor.
#' @param .add By default, when .add = FALSE, group_by will override existing groups.
#' To instead add to the existing groups, use .add = TRUE
#' @param x An object of class \code{treedata}
#' @details Groups the data frame and phylogeny by one of the factors in the data table.
#' @return An object of class \code{grouped_treedata}.
#' @seealso \code{\link{summarize}}
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat)
#' tdGrouped <- group_by(td, ecomorph)
#' summarize(tdGrouped, ntips = length(phy$tip.label),
#'    totalBL = sum(phy$edge.length), meanSVL = mean(SVL), sdSVL = sd(SVL))
#' @export
group_by.treedata <- function(.data, ..., .add=FALSE){
  groups <- group_by_prepare(.data$dat, ..., .add = .add)
  dat <- grouped_df(groups$data, groups$group_names)
  .data$dat <- dat
  class(.data) <- c("grouped_treedata", "treedata", "list")
  return(.data)
}


#' @title group_by_.treedata
#' @description Now defunct and replaced with group_by.treedata.
#'
#' @name group_by_-defunct
#' @seealso \code{\link{treeplyr-defunct}}
#' @keywords internal
#' @export
group_by_.treedata <- function(.data, ..., add = FALSE) {
    .Defunct(msg = "'group_by_' has been removed from this package; you can use group_by instead")
}

#' @rdname reorder.treedata
#' @export
reorder <- function(tdObject, ...){
  UseMethod("reorder")
}

#' Reorder a \code{treedata} object
#'
#' Reorders a \code{treedata} object. Both the tips and the data are automatically reordered to match.
#'
#' @param tdObject An object of class \code{treedata}
#' @param order Method for reordering
#' @param index.only Whether a index is returned rather than the reordered treedata object
#' @param ... Additional arguments to reorder.phylo
#' @return An object of class \code{treedata}
#'
#' @seealso \code{\link{reorder.phylo}}
#'
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat)
#' td <- reorder(td, "postorder")
#' @export
reorder.treedata <- function(tdObject, order="postorder", index.only=FALSE, ...){
  dat.attr <- attributes(tdObject$dat)
  phy <- reorder(tdObject$phy, order, index.only)#, ...)
  index <- match(tdObject$phy$tip.label, phy$tip.label)
  tdObject$dat <- tdObject$dat[index,]
  attributes(tdObject$dat) <-dat.attr
  #attributes(tdObject$dat)$row.names
  attributes(tdObject)$tip.label <- phy$tip.label
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

#' Run a function on the phylogeny of a \code{treedata} object
#'
#' Applies a function to the phylogeny in a \code{treedata} object. If the order of
#' tips are changed, or if tips are dropped, then the data are automatically reordered
#' to match the tree.
#'
#' @param tdObject An object of class \code{treedata}
#' @param FUN A function that operates on an object of class 'phylo'
#' @param ... Additional arguments

#' @return An object of class \code{treedata}
#'
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat)
#' td2 <- treeply(td, drop.tip, 1:50)
#'
#' par(mfrow=c(1,2))
#' plot(td$phy)
#' plot(td2$phy)
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
#' @details Note that if the parameter \code{phy} is specified in the additional parameters (i.e. '...'),
#' then it will be substituted with the \code{treedata} object \code{$phy}.
#'
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat)
#' td %>% forceNumeric(.) %>% tdapply(., 2, phytools::phylosig, tree=phy)
#' @export
tdapply <- function(tdObject, MARGIN, FUN, ...){
  if(!class(tdObject)[1]=="treedata") stop("Object is not of class 'treedata'")
  FUN <- match.fun(FUN)
  phy <- tdObject$phy
  env <- new.env(parent=parent.frame(), size=1L)
  env$phy <- tdObject$phy
  dat <- as.matrix(tdObject$dat)
  rownames(dat) <- phy$tip.label
  #rownames(dat) <- tdObject$phy$tip.label
  res <- eval(substitute(apply(dat, MARGIN, FUN, ...)), env)
  return(res)
}


#' @rdname treedply.treedata
#' @export
treedply <- function(tdObject, ...){
  UseMethod("treedply")
}

#' Run a function on a \code{treedata} object
#'
#' @param tdObject A treedata object
#' @param ... A function call.
#'
#' @details This function allows arbitrary R functions that use trees and data to be run on
#' \code{treedata} objects.
#'
#' @return Function output
#'
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat)
#' treedply(td, geiger::fitContinuous(phy, getVector(td, SVL), model="BM", ncores=1))
#' treedply(td, phytools::phylosig(phy, getVector(td, awesomeness), "lambda", test=TRUE))
#' treedply(td, phytools::phenogram(phy, getVector(td, SVL), ftype="off", spread.labels=FALSE))
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
#' @param td A treedata object
#' @param ... The name of the column to select
#'
#' @return A named vector
#' @export
getVector <- function(td, ...){
  dat <- td$dat
  args <- as.character(substitute(list(...)))[-1L]
  arg_sub <- type.convert(args)
  if(is.numeric(arg_sub) | is.integer(arg_sub)) args <- arg_sub
  vecs <- lapply(args,function(x) dat[[x]])
  vecs <- lapply(vecs, function(x) setNames(x, attributes(td)$tip.label))
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
summary.treedata <- function(object, ...){

  cat('A treeplyr treedata object', "\n")
  cat(paste('The dataset contains ', ncol(object$dat), ' traits'), "\n")
  types <- setNames(suppressWarnings(detectAllCharacters(as.matrix(object$dat))), colnames(object$dat))
  cat("Continuous traits: ", names(types)[which(types=="continuous")], "\n")
  cat("Discrete traits: ", names(types)[which(types=="discrete")], "\n")
  cat(paste("The following traits have missing values:", paste(names(types)[apply(object$dat, 2, function(y) any(is.na(y)))], collapse=", "), "\n"))
  cat(paste("These taxa were dropped from the tree:", paste(attributes(object)$dropped$dropped_from_tree, collapse=", "), "\n"))
  cat(paste("These taxa were dropped from the data:", paste(attributes(object)$dropped$dropped_from_data, collapse=", "), "\n"))


  cat("$phy \n")
  print(object$phy)
  cat("\n$dat \n")
  print(object$dat)
}



#' Function for checking whether a treedata object contains only numeric columns and for forcing
#' data columns into numeric format
#'
#' This function can be used to check if a treedata object contains numeric columns and, if
#' desired, drop all non-numeric columns.
#'
#' @param tdObject A \code{treedata} object
#' @param return.numeric If TRUE, then a treedata object with all numeric columns will be returned;
#' non-numeric columns will be removed.
#' @return If return.numeric, then an object of class "\code{treedata}" with only numeric columns.
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat)
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

#' Function for checking whether a treedata object contains only factors and for forcing
#' data columns into factor format
#'
#' This function can be used to check if a treedata object contains factors and, if desired,
#' convert all columns automatically to factors.
#'
#' @param tdObject A \code{treedata} object
#' @param return.factor If TRUE, then a treedata object with all factors will be returned;
#' columns will be forced into factors using \code{factor} and any with no repeated elements
#' will be removed.
#' @return If return.factor, then an object of class "\code{treedata}" with all columns as
#' factors.
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat)
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

#' @rdname mutate.treedata
#' @export
mutate.grouped_treedata <- function(.data, ...){
  dat <- mutate(.data$dat, ...)
  .data$dat <- dat
  return(.data)
}


#' @title mutate_.grouped_treedata
#' @description Now defunct and replaced with mutate.groupedtreedata.
#'
#' @name mutate_-defunct
#' @seealso \code{\link{treeplyr-defunct}}
#' @keywords internal
#' @export
mutate_.grouped_treedata <- function(.data, ..., add = FALSE) {
  .Defunct(msg = "'mutate_' has been removed from this package; you can use group_by instead")
}

#' @rdname filter.treedata
#' @export
filter.grouped_treedata <- function(.data, ...){
  cl <- class(.data$dat)
  .data$dat$tip.label <- .data$phy$tip.label
  dat <- filter_(.data$dat, ...)
  .data$dat <- dat
  attributes(.data)$tip.label <- .data$dat$tip.label
  .data$dat <- select(.data$dat, 1:(ncol(.data$dat)-1))
  class(.data$dat) <- cl
  .data$phy <- drop.tip(.data$phy, .data$phy$tip.label[!(.data$phy$tip.label %in% attributes(.data)$tip.label)])
  return(.data)
}

#' @title filter_.grouped_treedata
#' @description Now defunct and replaced with filter.groupedtreedata.
#'
#' @name filter_-defunct
#' @seealso \code{\link{treeplyr-defunct}}
#' @keywords internal
#' @export
filter_.grouped_treedata <- function(.data, ...) {
  .Defunct(msg = "'filter_' has been removed from this package; you can use filter instead")
}

#' Add regimes to a treedata object
#'
#' This function paints clades on the phylogeny and adds a data column that specifies to which clade each species
#' belongs
#' @param tdObject A \code{treedata} object
#' @param nclades The number of clades that will be specified if used interactively
#' @param name The name of the resulting data column
#' @param interactive If \code{TRUE}, then a plot will appear that will allow the user to click on \code{nclades}
#' branches. The selections will then be coverted into the data table.
#' @param type Either "nodes" or "branches" specifying if the ids provided specify the branch id (assuming a
#' post-ordered tree) or the node number. Ignored if \code{interactive = TRUE}.
#' @param ids A vector of node numbers of branch numbers that specify clades. Ignored if \code{interactive=TRUE}.
#' @param plot If \code{TRUE} and \code{interactive = FALSE} then a simmap plot is produced.
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat)
#' td <- reorder(td, "postorder")
#' td.painted <- paint_clades(td, interactive=FALSE, type="nodes",
#'                                    ids=c(184, 160, 135, 122), plot=TRUE)
#' td.painted <- group_by(td.painted, clades)
#' summarise(td.painted,
#'              psig1 = phytools::phylosig(setNames(SVL, phy$tip.label), tree=phy),
#'                    meanSVL = mean(SVL))
#' @export
paint_clades <- function(tdObject, nclades=1, name="clades", interactive=TRUE, type="nodes", ids=NULL, plot=TRUE){
  tdObject <- reorder(tdObject, "postorder")
  if(interactive){
    regimes <- .identifyBranches(tdObject$phy, nclades)
    cat("branches ", paste(regimes$sb, collapse=", "), "selected\n")
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
        cache <- .prepare.branches(tdObject$phy)
        pars <- list(sb = sb, loc = loc, t2 = 1:n + 1)
        tr <- .toSimmap(.pars2map(pars, cache), cache)
        cols <- .plotRegimes(tr, lwd=1, pal=rainbow)
        L <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        legend(min(L$xx), max(L$yy), legend=names(cols), lwd=3, col=cols)
      }

    }
  regimes$t2 <- 2:(length(regimes$sb)+1)
  regimes$k <- nclades
  regimes$ntheta <- nclades+1
  tdObject <- mutate(tdObject, .tipregime(regimes, tdObject$phy))
  colnames(tdObject$dat)[ncol(tdObject$dat)] <- name
  return(tdObject)
}

#' @rdname group_by.treedata
#' @export
ungroup.grouped_treedata <- function(x, ...){
  x$dat <- ungroup(x$dat, ...)
  return(x)
}

#' @export
'[[.treedata' <- function(x, ..., exact=TRUE){
  y <- x$dat
  res <- '[[.data.frame'(y, ..., exact=exact)
  if(length(res) != nrow(y)){
    stop("Use '[' for selecting multiple columns")
  }
  return(setNames(res, x$phy$tip.label))
}

#' @export
'[.treedata' <- function(x, i, j, drop=FALSE, tip.label=FALSE){
  if(!tip.label){
    y <- as.data.frame(x$dat)
  } else {
    y <- as.data.frame(cbind(tip.label= x$phy$tip.label, x$dat))
    if(!missing(j)){
      if(is.numeric(j)){
        j <- c(1, j+1)
      } else {
        j <- c('tip.label', j)
      }
    }
  }
  if (drop)
    warning("drop ignored", call. = FALSE)
  nr <- nrow(y)
  if (nargs() <= 2L) {
    if (!missing(i)) {
      .check_names_df(y, i)
      result <- .subset(y, i)
    }
    else {
      result <- y
    }
    attr(result, "row.names") <- .set_row_names(nr)
    return(.as_data_frame.data.frame(result))
  }
  if (!missing(j)) {
    .check_names_df(y, j)
    y <- .subset(y, j)
  }
  if (!missing(i)) {
    if (length(y) == 0) {
      nr <- length(attr(y, "row.names")[i])
    }
    else {
      y <- lapply(y, `[`, i)
      nr <- length(y[[1]])
    }
  }
  attr(y, "row.names") <- .set_row_names(nr)
  res <- .as_data_frame.data.frame(y)
}

#' .check_names_df function from tibble package
#' @keywords internal
.check_names_df <- function(x,j){
  if (is.character(j) && any(wrong <- !j %in% names(x))) {
    names <- j[wrong]
    stop(sprintf("undefined columns: %s", paste(names, collapse = ", ")))
  }
}

#' as_data_frame.data.frame function from the tibble package
#' @keywords internal
.as_data_frame.data.frame <- function (x, ...)
{
  class(x) <- c("tbl_df", "tbl", "data.frame")
  x
}
