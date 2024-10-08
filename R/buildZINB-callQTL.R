#' Build zinb model.
#'
#' @param counts a vector for gene expression.
#' @importFrom gamlss gamlssML
#' @return Four parameters
#' @export
#' @examples
#' data(testGene)
#' gene <- unlist(testGene[1, ])
#' result <-buildZINB(gene)
buildZINB <- function(counts) {
    options(show.error.messages = FALSE)
    zinb_gamlssML <- try(gamlssML(counts, family = "ZINBI"), silent = TRUE)
    options(show.error.messages = TRUE)

    if (inherits(zinb_gamlssML, "try-error")) {
        message("MLE of ZINB failed! Please choose another model")
        return(list(theta = NA, mu = NA, size = NA, prob = NA))
    } else {
        theta <- zinb_gamlssML$nu
        mu <- zinb_gamlssML$mu
        size <- 1/zinb_gamlssML$sigma
        prob <- size/(size + mu)
    }
    return(list(theta = theta, mu = mu, size = size, prob = prob))
}
