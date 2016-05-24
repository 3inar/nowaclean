#' Blood type probes
#'
#' Returns a list of probes that cover HLA genes. These are usually strongly
#' expressed and it may sometimes be desirable to remove them before
#' analysis. The function \code{blood_probes()} is a wrapper to this
#' function with a friendlier name.
#'
#' @return A character vector of nuIDs. If you for some reason want to
#' find the Illumina IDs behind these nuIDs you can use \code{names()} on the
#' returned vector.
#'
#' @source The Genomics Core Facility, NTNU \url{https://www.ntnu.edu/dmf/gcf}
#' @export
#' @examples
#' ## nuIDs:
#' hla_hist()
#'
#' ## Illumina IDs:
#' names(hla_hist())
hla_hist <- function() {
  # this is serialized from a character vector with dput()
  structure(c("BgSEayCHZU.6DuFt7U", "0cVOsh0SeXuuB.p1eU", "xhfxSSl1UUjleZ6MXA",
              "3xABh5KCFUHEtVJn1U", "KYiXiWoXopFIr0oicQ", "T3p2gk3kQf4pR0eRVE",
              "xVb14TRJJdQRCnUIXI", "HzqJNpx567t9VBpGis", "oCkSogqMFcRUhVKehM",
              "HgCiovEdKeQkuEuVKk", "9V4f795ApF5O7e727I", "3M5JGKp5CJV0V4nqJQ",
              "uUi6XlzRBRNX1VRCkM", "cdeqd0u5TIoSQeuPv0", "BpQH6JK7U1SV7svIko",
              "Qd31eNfoolYHj3one4", "NFtdNMC3eb3pThValQ", "9S4nlVA0L8uV01Pzt0",
              "TfXanqXzTU0Si0gKnU", "9XqXh7q3rqTTz6hTfQ", "cXS4CS030k.1JXMlFU",
              "BHeh9JQSj14nguSOEM", "xSMQcar7oifRLkmmLU", "KBVEl_OkmXbfQf_59c",
              "lsRoeDeRIkI54ui2qk", "lBLUXe7cU4VX10R4Xs", "Bzkii9KFxeipi5ei_c",
              "NUU4V31054Hk9fVQ0U", "QjIAoqJx3SnkJpBLlQ", "B157U7P3256g16V7VE",
              "TPeolujqponqv65L9Q", "BqlfKW_of4xUoIMSYs", "xXcolO813tVTItCXLk",
              "3NHP_7nWbnMn.qmleU", "W1ASS5ebVnUaLfUgRc", "37SdG6X8AFx9VlZKJ0",
              "xnL3utKaRQ4gKDkhXU", "KGhlVCGTOHRRnt3hTg"),
            .Names = c("160132",  "5080333", "6420253", "4260736", "5340382",
                       "5910220",  "1240070",  "4560047", "7040008", "2350066",
                       "3400438",  "7200398", "540563",  "4900731", "2470553",
                       "3450338",  "1190039", "1050360", "7160474",  "270168",
                       "2570564",  "1770504", "510079", "620544", "4200725",
                       "1030747",  "5220070", "2070088", "1980592", "1050746",
                       "6380022",   "3120576", "2900451", "2900204", "270136",
                       "1470673",  "2650370",  "5080692"))
}

#' @describeIn hla_hist
#' A more descriptive name for hla_hist
#' @export
blood_probes <- function() {
  hla_hist()
}
