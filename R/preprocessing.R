suggested_packages <- c(
  "lumi",
  "Biobase",
  "limma",
  "genefilter",
  "preprocessCore",
  "AnnotationDbi",
  "lumiHumanAll.db",
  "illuminaHumanv4.db"
)

req_pkg <- function(pkg_name) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop(pkg_name, " (biconductor) needed for this function to work. Please install it.",
         call. = FALSE)
  }
}

req_pkgs <- function() {
  plyr::a_ply(suggested_packages, 1, req_pkg)
}

# do clean 1st
#' @export
corrected <- function(data, negative_controls) {
  req_pkgs()

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
  req_pkgs()
  tmp <- preprocessCore::normalize.quantiles(exprs(data), copy=TRUE)
  rownames(tmp) <- rownames(data)
  colnames(tmp) <- colnames(data)
  exprs(data) <- tmp

  data
}

# standards by gunther
# fval: present in at least (fval*100)% of samples
# also filters bad probes
#' @export
filtered <- function(data, pval=0.01, fval=0.01, verbose=FALSE) {
  req_pkgs()
  threshold <- round(fval*ncol(data))

  # counts number of samples with probe detectable at pval level
  present <- lumi::detectionCall(data, Th=pval, type="probe")
  remove <- present < threshold

  probeq <- illuminaHumanv4.db::illuminaHumanv4PROBEQUALITY
  probes <- lumi::nuID2IlluminaID(rownames(data), lib.mapping=NULL,
                                  species ="Human", idType='Probe')
  probe_quality <- unlist(AnnotationDbi::mget(as.character(probes), probeq, ifnotfound=NA))
  bad_quality <- (probe_quality == "Bad") | (probe_quality == "No match")

  if(verbose){
    cat("Removing", sum(remove), "probes that were not present in at least",
        fval*100, "% of the samples.\n")
    cat("Removing", sum(bad_quality), "bad quality probes.\n")
  }

  data <- data[-which(remove | bad_quality), ]
  data
}

#' @export
probe_aggregated <- function(data, verbose=FALSE) {
  req_pkgs()
  probes_in = length(Biobase::featureNames(data))

  annotation(data) <- "lumiHumanAll"
  data <- genefilter::nsFilter(data, var.filter=FALSE)$eset ## aggregate across annotated probes

  probes_out = length(Biobase::featureNames(data))
  if(verbose){
    cat("Filtered out", probes_in-probes_out,
      "probes.\n")
  }

  data
}

#' @export
gene_names <- function(data) {
  req_pkgs()
  nu_ids <- rownames(data)
  mapping <- lumi::nuID2RefSeqID(nu_ids, lib.mapping='lumiHumanIDMapping', returnAllInfo =TRUE)

  mapping
}

#' Standard NOWAC preprocessing
#'
#' Performs the standard NOWAC preprocessing steps on your data. This happens
#' after outlier removal. Performs the following steps in the following order:
#' \itemize{
#'  \item Background correction. Uses negative control probes (probes that are
#'  designed to detect noise) to do a normal-exponential convolution correction
#'  \item Normalization. Quantile/shoehorn normalization.
#'  \item Probe filtering. Filters out porobes with bad annotation quality,
#'  probes below detection threshold, and probes that arent' present in a certain
#'  fraction of the population.
#'  \item Probe aggregation. Combine probes that detect the same gene to a single
#'  value.
#' }
#'
#' @export
preprocessed <- function(data, negative_controls, pval=0.01, fval=0.01, verbose=FALSE) {
  data <- corrected(data, negative_controls)
  data <- normalized(data)
  data <- filtered(data, pval=pval, fval=fval, verbose=verbose)
  data <- probe_aggregated(data, verbose=verbose)

  data
}
