#' Simulate data for a binary outcome
#'
#' This function generates a dataframe with 6 patient covariates and a binary outcome simulated from
#' a model that uses the covariates.
#' @param Npat Number of patients to simulate.
#' @return The function returns a dataframe with:
#'
#' x1, x2, x3, x4: patient covariates.
#'
#' t= treatment assignment (0 for control, 1 for active).
#'
#' logit.control= the logit of the probability of an outcome in the control treatment.
#'
#' logit.active= the logit of the probability of an outcome in the active treatment.
#'
#' benefit= treatment benefit in log odds ratio.
#'
#' py=the probability of the outcome for each patient, under the treatment actually administered.
#'
#' logit.py= the logit of py.
#'
#' y.observed= the observed outcome
#' @importFrom MASS mvrnorm
#' @importFrom stats rbinom rnorm
#' @examples
#' dat1=simbinary(100)$dat
#' head(dat1)
#' @export
simbinary<-function (Npat = 100)
{
  x1x2 = mvrnorm(n = Npat, c(0, 0), Sigma = matrix(c(1, 0.2,
                                                     0.2, 1), 2, 2))
  simdat = data.frame(x1 = x1x2[, 1])
  simdat$x2 = x1x2[, 2]
  simdat$x3 = rbinom(Npat, 1, prob = 0.2)
  simdat$x4 = rbinom(Npat, 1, prob = 0.1)

  pt <- 0.5
  simdat$t <- rbinom(Npat, 1, prob = pt)
  simdat$logit.control <- with(simdat, -2 + 0.5 * x1 + 0.4 *
                                 x2 + 0.3 * x3 + 0.5 * x4  )
  simdat$benefit = with(simdat, -0.3+0.3*x1+0.2*x3+0.1*x4)
  simdat$logit.active = with(simdat, logit.control + benefit)
  simdat$logit.py = with(simdat, logit.control*(t==0) + logit.active*(t==1))
  simdat$py <- expit(simdat$logit.py)
  simdat$y.observed = rbinom(Npat, 1, prob = simdat$py)
  return(list(dat = simdat))


}



#' Expit
#'
#' Calculates the expit of a real number
#' @param x A real number
#' @return exp(x)/(1+exp(x))
#' @examples
#' expit(2.3)
#' @export
expit=function(x){ exp(x)/(1 + exp(x))}


#' Logit
#'
#' Calculates the logit of a real number between 0 and 1
#' @param x A real number between 0 and 1
#' @return log(x/(1-x))
#' @examples
#' logit(0.2)
#' @export
logit=function(x){l=log(x/(1-x)); return(l)}
