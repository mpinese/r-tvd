# tvd1d.R: 1-D TVD algorithms
# 
# Copyright 2014 Mark Pinese
#
# This file is distributed under the terms of the Eclipse Public 
# License v1.0, available at:
# https://www.eclipse.org/org/documents/epl-v10.html


#' Perform Total Variation Denoising on a 1-Dimensional Signal
#' 
#' When supplied a noisy sequential signal in vector y, performs
#' TVD with regularization parameter lambda, and returns a
#' denoised version of y.
#'
#' 1D TVD is a filtering technique for a sequential univariate signal that attempts
#' to find a vector x_tvd that approximates a noisy vector y, as:
#'   \deqn{x_{tvd} = argmin_{x}(E(x, y) + \lambda V(x))}{x_tvd = argmin_x(E(x, y) + \lambda*V(x))}
#' where E(x, y) is a loss function measuring the error in approximating
#' y with x, and V(x) is the total variation of x:
#'   \deqn{V(x) = sum(|x_{i+1} - x_{i}|)}
#'
#' TVD is particularly well-suited to recovering piecewise constant 
#' signals.  The degree of approximation is controlled by the parameter
#' lambda: for lambda = 0, x_tvd = y, and as lambda increases, x_tvd
#' contains increasingly fewer value transitions, until, for a high
#' enough value, it is nearly constant.
#'
#' Currently only implements Condat's fast squared-error loss TVD
#' algorithm (method "Condat"), which is restricted to vectors of
#' length 2^32 - 1 and shorter.
#' 
#' @param y a numeric vector of sequential noisy data values
#' @param lambda the total variation penalty coefficient
#' @param method a string indicating the algorithm to use for
#'   denoising.  Currently only supports method "Condat"
#'
#' @return a numeric vector of the same length as y, containing
#'   denoised data.
#'
#' @importFrom Rcpp evalCpp
#'
#' @export
#' @references Condat, L. (2013) A Direct Algorithm for 1-D Total Variation Denoising
#'   IEEE Signal Processing Letters 20(11): 1054-1057.  \url{doi:10.1109/LSP.2013.2278339}
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
#'
#' @examples
#' ## Generate a stepped signal
#' x = rep(c(1, 2, 3, 4, 2, 4, 3, 2, 1), each = 100)
#'
#' ## Create a noisy version of the signal
#' y = x + rnorm(length(x), sd = 0.5)
#'
#' ## Denoise the signal by Condat's method
#' x.denoised = tvd1d(y, lambda = 10, method = "Condat")
#'
#' ## Plot the original signal, the noisy signal, and the denoised signal
#' plot(y, col = "black", pch = 19, cex = 0.3)
#' lines(x, col = "blue", lwd = 3)
#' lines(x.denoised, col = "red", lwd = 3)
#' legend("topleft", legend = c("Original", "Noisy", "Denoised"), 
#'   col = c("blue", "black", "red"), lty = c("solid", "solid", "solid"), 
#'   lwd = c(2, 0, 1), pch = c(NA, 19, NA), pt.cex = c(NA, 0.3, NA), inset = 0.05)
tvd1d <- function(y, lambda, method = c("Condat"))
{
	method = match.arg(method)

	if (method == "Condat")
	{
		if (length(y) > 2^32-1)
		{
			stop("Method \"Condat\" is currently limited to length(y) < 2^32")
		}
		x = tvd_1d_condat_worker(y, lambda)
	}

	x
}


SquaredError <- function(x, y)
{
	sum((x - y)^2)		# TODO: Maybe look into inlining with Rcpp for speed and memory efficiency?
}


#' Identify the Optimal lambda Penalty for 1-D TVD
#' 
#' @export
#' @examples
#' ## Generate a stepped signal
#' x = rep(c(1, 2, 3, 4, 2, 4, 3, 2, 1), each = 100)
#'
#' ## Create a noisy version of the signal
#' y = x + rnorm(length(x), sd = 0.25)
#'
#' ## Determine the best value for lambda
#' cv.result = cv.tvd1d(y)
#'
#' ## Plot the effect of lambda on prediction error
#' plot(cv.result$losses[1,] ~ cv.result$lambdas, log = "x", ylim = range(c(cv.result$losses[1,] + cv.result$losses[2,], cv.result$losses[1,] - cv.result$losses[2,])))
#' arrows(cv.result$lambdas, cv.result$losses[1,] - cv.result$losses[2,], cv.result$lambdas, cv.result$losses[1,] + cv.result$losses[2,], length = 0.05, angle = 90, code = 3)
#' abline(v = cv.result$lambda.min, col = "red")
#' abline(v = cv.result$lambda.1se, col = "orange")
#'
#' ## Denoise the signal using the best lambdas
#' x.denoised.min = tvd1d(y, lambda = cv.result$lambda.min)
#' x.denoised.1se = tvd1d(y, lambda = cv.result$lambda.1se)
#'
#' ## Plot the original signal, the noisy signal, and the denoised signal
#' plot(y, col = "black", pch = 19, cex = 0.3)
#' lines(x, col = "blue", lwd = 3)
#' lines(x.denoised.min, col = "red", lwd = 3)
#' lines(x.denoised.1se, col = "orange", lwd = 3)
#' legend("topleft", legend = c("Original", "Noisy", "Denoised Min", "Denoised 1SE"), 
#'   col = c("blue", "black", "red", "orange"), lty = c("solid", "solid", "solid", "solid"), 
#'   lwd = c(2, 0, 1, 1), pch = c(NA, 19, NA, NA), pt.cex = c(NA, 0.3, NA, NA), inset = 0.05)
cv.tvd1d <- function(y, lambdas = NULL, method = c("Condat"), nfolds = 10, lossfunc = c("MatchMethod", "SquaredError"))
{
	method = match.arg(method)
	lossfunc = match.arg(lossfunc)

	if (lossfunc == "MatchMethod")
	{
		# Determine the loss function to use based on the TVD method selected
		if (method == "Condat")
		{
			lossfunc = "SquaredError"
		}
	}

	if (is.null(lambdas))
	{
		# Construct a reasonable series of lambdas
		# TODO: Make this more intelligent!
		lambdas = 10^seq(-2, 3, 0.1)
	}

	losses = sapply(lambdas, cv.tvd1d.singlelambda, y = y, method = method, nfolds = nfolds, lossfunc = lossfunc)
	mean_losses = losses[1,]
	se_losses = losses[2,]

	best_loss_min = min(mean_losses)
	best_lambda_min = lambdas[which.min(mean_losses)]
	index_1se = max(which(mean_losses - se_losses <= best_loss_min))
	best_loss_1se = mean_losses[index_1se]
	best_lambda_1se = lambdas[index_1se]

	return(list(lambdas = lambdas, losses = losses, lambda.min = best_lambda_min, lambda.1se = best_lambda_1se))
}


cv.tvd1d.singlelambda <- function(lambda, y, method, nfolds, lossfunc)
{
	folds = rep(1:nfolds, times = ceiling(length(y) / nfolds))[1:length(y)]

	cv.losses = sapply(1:nfolds, cv.tvd1d.singlelambda.singlefold, y = y, lambda = lambda, method = method, nfolds = nfolds, lossfunc = lossfunc)
	c(mean(cv.losses), sd(cv.losses) / sqrt(nfolds))
}


cv.tvd1d.singlelambda.singlefold <- function(foldid, y, lambda, method, nfolds, lossfunc)
{
	indices.test = seq.int(foldid, length(y), nfolds)
	indices.train = (1:length(y))[-indices.test]
	y.test = y[indices.test]
	y.train = y[indices.train]

	x.train = tvd1d(y.train, lambda, method)
	x.test.pred = approx(indices.train, x.train, indices.test, method = "linear", rule = 2)$y
	do.call(lossfunc, list(x = x.test.pred, y = y.test))
}
