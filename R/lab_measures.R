#' @export
lab_thresholds <- c("RIN value < 7", "260/280 ratio < 2",
                    "260/230 ratio < 1.7", "50 < RNA < 500")

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

