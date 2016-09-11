# do clean 1st
#' @export
corrected <- function(data, negative_controls) {
  
  if(is(data, "LumiBatch")){
    expr <- Biobase::exprs(data)
  } else if(is(data, "matrix")) { 
    expr <- data
  } else {
    stop("data must be either a LumiBatch object or a matrix")
  }
  
  no_variance = names(which(apply(expr, 2, sd) == 0))
  if(length(no_variance > 0)){
    stop("Found samples with no variance (e.g. gene expression values all zero). Please investigate samples: ",
         paste(unlist(no_variance), collapse=" "), ".")
  }
  
  negative = names(which(apply(expr, 2, min) < 0))
  if(length(negative > 0)){
    stop("Found samples with negative gene expression values. Please investigate samples: ",
         paste(unlist(negative), collapse=""), ".")
  }
  
  
  total <- rbind(expr, negative_controls[, colnames(expr)])
  probe_t <- c(rep("regular", nrow(expr)), rep("negative", nrow(negative_controls)))
  expr <- limma::nec(total, probe_t)
  expr <- expr[probe_t == "regular", ]
  
  if(is(data, "LumiBatch")){
    Biobase::exprs(data) <- expr
  } else {
    data <- expr
  }
  
  data
}


#' @export
normalized <- function(data) {
  tmp <- preprocessCore::normalize.quantiles(exprs(data), copy=TRUE)
  rownames(tmp) <- rownames(data)
  colnames(tmp) <- colnames(data)
  exprs(data) <- tmp

  data
}

# standards by gunther
# fval: present in at least (fval)% of samples
# also filters bad probes
#' @export
filtered <- function(data, pval=0.01, fval=0.01) {
  threshold <- round(fval*ncol(data))

  # counts number of samples with probe detectable at pval level
  present <- detectionCall(data, Th=pval, type="probe")
  remove <- present < threshold

  probeq <- illuminaHumanv4.db::illuminaHumanv4PROBEQUALITY
  probes <- lumi::nuID2IlluminaID(rownames(data), lib.mapping=NULL,
                                  species ="Human", idType='Probe')
  probe_quality <- unlist(AnnotationDbi::mget(as.character(probes), probeq, ifnotfound=NA))
  bad_quality <- (probe_quality == "Bad") | (probe_quality == "No match")

  data <- data[-which(remove | bad_quality), ]
  data
}

#' @export
probe_aggregated <- function(data) {
  annotation(data) <- "lumiHumanAll"
  data <- genefilter::nsFilter(data, var.filter=FALSE)$eset ## aggregate across annotated probes

  data
}

#' @export
gene_names <- function(data) {
  nu_ids <- rownames(data)
  mapping <- nuID2RefSeqID(nu_ids, lib.mapping='lumiHumanIDMapping', returnAllInfo =TRUE)

  mapping
}

#' @export
preprocessed <- function(data, negative_controls, pval=0.01, fval=0.01) {
  data <- corrected(data, negative_controls)
  data <- normalized(data)
  data <- filtered(data, pval=pval, fval=fval)
  data <- probe_aggregated(data)

  data
}
