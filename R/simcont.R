#' Simulate data for a prediction model of a continuous outcome
#'
#' This function generates a dataframe with 6 patient covariates and a continuous outcome
#' simulated from a model that uses the covariates.
#' @param Npat Number of patients to simulate.
#' @return The function returns a dataframe with:
#'
#' x1, x2, x3, x4: patient covariates.
#'
#' t= treatment assignment (0 for control, 1 for active).
#'
#' y.control= the outcome if the patient takes the control treatment.
#'
#' y.active= the outcome if the patient takes the active treatment.
#'
#' benefit= the treatment benefit, i.e. y.active-y.control.
#'
#' y.observed= the observed outcome.
#' @importFrom MASS mvrnorm
#' @importFrom stats rbinom rnorm
#' @examples
#' dat1=simcont(100)$dat
#' head(dat1)
#' @export
simcont<-function (Npat=100)
{
  expit = function(x) {
    exp(x)/(1 + exp(x))
  }
  x1x2 <- mvrnorm(n = Npat, c(0, 0), Sigma = matrix(c(1, 0.3,
                                                      0.3, 1), 2, 2))
  simdat <- data.frame(x1 = x1x2[, 1])
  simdat$x2 <- x1x2[, 2]
  simdat$x3 <- rbinom(Npat, 1, prob = 0.2)
  simdat$x4 <- rnorm(Npat, 0, 1)

  pt <- 0.5
  simdat$t <- rbinom(Npat, 1, prob = pt)
  simdat$y.control <- with(simdat, 5 + 0.3 * x1 + 0.3 * x2 +
                             0.3 * x3 + 0.3 * x4
                           + rnorm(Npat,0, 1))

  simdat$benefit <- with(simdat, 0.5+ 0.3*x1-x3)
  simdat$y.active <- with(simdat, y.control + benefit)
  simdat$y.observed <- with(simdat, y.control + t * benefit)
  return(list(dat = simdat))}
