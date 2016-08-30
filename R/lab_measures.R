#' Grounds for excluding a sample based on lab measures
#'
#' We have several quality control measurements from the lab; when we are
#' uncertain whether an array looks bad we can inspect these measures and
#' get a sense of whether something has gone wrong. We use these cutoffs to
#' say that an array has gone wrong:
#' \itemize{
#'  \item{\strong{RIN < 7:}}{ An array with RIN less than 7 is considered bad.}
#'  \item{\strong{260/280 RNA ratio < 2:}}{ An array with 260/280 less than 2 is considered bad.}
#'  \item{\strong{260/230 RNA ratio < 1.7:}}{ An array with 260/230 less than 1.7 is considered bad.}
#'  \item{\strong{Ng/ul RNA outside of (50, 500):}}{ Well-behaved arrays should have
#'  50 < Ng/ul RNA < 500.}
#' }
#' \strong{NB} that perfectly well-behaved arrays might have lab measures outside
#' of these limits, we never exclude arrays based only on lab measures.
#' @export
lab_thresholds <- c("Bad: RIN value < 7", "Bad: 260/280 RNA ratio < 2",
                    "Bad: 260/230 RNA ratio < 1.7", "Good: 50 < Ng/ul RNA < 500")

#' @export
lab_variables <- c("Ng/ul_RNA", "260/280_RNA", "260/230_RNA", "RIN",
               "Ng/ul_cRNA", "260/280_cRNA", "260/230_cRNA")

# TODO: this is too dangerous right now
# check_lab_thresholds <- function(samplename, info) {
#   if (!all(lab_variables %in% colnames(info))) stop("variables have different names in info object than expected, you must check thresholds by hand.")
#
#   inside <- c(info[1] > 50 & info[1] < 500, info[2] < 2, info[3] < 1.7,
#                info[4] < 7, info[5] > 50 & info[5] < 500, info[6] < 2, info[7] < 1.7)
#
#   ret <- ifelse(inside, "OK", "OUTSIDE")
#   names(ret) <- lab_variables
#
#   ret
# }

