# install.packages("timedeppar")
# install.packages("readxl")

# libraries
set.seed(42)
library(readxl)
library(timedeppar)

##############
# reading data
data   = read_excel("~/Scrivania/MOD B/esercizi/github/Project/data/syntheticData.xlsx", skip = 1)
data   = as.data.frame(data[c(1, 2, 4, 5, 7, 8)])
colnames(data) = c("t", "y", "t", "y", "t", "y")
data   = na.omit(data)
data1  = data[1:2]
data2  = data[3:4]
data3  = data[5:6]
n      = dim(data)[1]
n_interval = 100             # default: 50
n_iter     = 5000            # default: 10000
omega      = 2 * pi / 11                 # fixed (period T = 11)
count      = 0                           # fixed
max_count  = n_iter * (n_interval + 1)   # fixed
resol  = 3                   # default: 5
delta  = floor(max_count * resol / 100)
#
# time-dependent parameters
mu         = 0
tau        = 1
sigma_xi   = 0.03
#################


#######################################
# maker of vector of initial parameters
param_ini    = list(A        = 10,
					phi      = 2,
				    sigma_y  = 0.05,
				    xi       = data.frame(t  = 1:(n - 1),
										  xi = 1))

# turning off working with log-scale inside infer.timedeppar
param_logic  = c(A           = FALSE,
				 phi         = FALSE,
				 sigma_y     = FALSE,
				 xi          = FALSE)

# ranges initialization (both time-dependent and not)
param_ranges = list(A        = c(5,    15),
					phi      = c(0,    2 * pi),
					sigma_y  = c(0.01, 10),
					xi_mean  = c(-10,  10),
					xi_sd    = c(0.01, 10),
					xi_gamma = c(0.05, 20))

# vector of initial values of parameters associated to time-dependent xi
param_ou_ini = c(xi_mean     = mu,
				 xi_sd       = sigma_xi * sqrt(2 / tau),
				 xi_gamma    = 1 / tau)
#######################################


########################
# loglikelihood function
loglikeli <- function(param) {
	xis       = param$xi$xi
	A         = param$A
	phi       = param$phi
	sigma_y   = param$sigma_y
	#
	# WORKING WITH xis:
	#y_likeli  = y - (A * cos(omega * (sum_{j=i}^{2}(xi_j) + 1) + phi)) / sigma_y
	y_likeli  = (y - A * cos(omega * c(1, cumsum(xis) + 1) + phi)) / sigma_y
	y_likeli  = sum(dnorm(y_likeli, 0, 1, TRUE))
	#
	# feedback
	count <<- count + 1
	if (count %% delta == 0) {
		prog = count * resol / delta
		print(paste(prog, "%. loglikeli: ", y_likeli))
	}
	return(y_likeli)
}

logprior    <- function(param) {
	A         = param_ini$A
	phi       = param_ini$phi
	sigma_y   = param_ini$sigma_y
	#
	prior     = dnorm(param[["A"]],        A,        A / 2, TRUE) +   # A
				dnorm(param[["phi"]],      phi,      phi,   TRUE) +   # phi
				dnorm(param[["sigma_y"]],  sigma_y,  1,     TRUE)     # sigma_y
	#
	return(prior)
}

ou_logprior <- function(param) {
	xi_mean   = param_ou_ini[["xi_mean"]]
	xi_sd     = param_ou_ini[["xi_sd"]]
	xi_gamma  = param_ou_ini[["xi_gamma"]]
	#
	prior     = dnorm(param[["xi_mean"]],  xi_mean,  2,     TRUE) +   # sigma_xi
				dnorm(param[["xi_sd"]],    xi_sd,    1,     TRUE) +   # tau
				dnorm(param[["xi_gamma"]], xi_gamma, 1,     TRUE)     # mu
	#
	return(prior)
}

CALL <- function() {
	res <- infer.timedeppar(loglikeli           = loglikeli,
							loglikeli.keepstate = FALSE,
							param.ini           = param_ini,
							param.range         = param_ranges,
							param.log           = param_logic,
							param.logprior      = logprior,
							param.ou.ini        = param_ou_ini,
							param.ou.logprior   = ou_logprior,
							task                = "start",
							n.iter              = n_iter,
							verbose             = 0,
							control             = list(n.interval = n_interval))
	return(res)
}
#
y = data1$y
res1 = CALL()
#
#y = data2$y
#res2 = CALL()
#
#y = data3$y
#res3 = CALL()


PLOT <- function(res, d) {
	params = res[["param.maxpost"]]
	A      = params$A
	phi    = params$phi
	xis    = xis.extractor(param = params)
	#xis    = xis[1:40]
	ts     = c(1, cumsum(xis) + 1)
	#ts     = ts[1:40]
	ys     = A * cos(omega * ts + phi)
	#ys     = ys[1:40]
	#
	plot(ts, ys,
		 type = "l",
		 col  = "blue",
		 ylim = c(min(c(ys, d$y)),
				  max(c(ys, d$y))))
	lines(d$t, d$y,
		  col = "red")
	legend("topright", c("inferred", "data"), lty = 1, col = c("blue", "red"))
}

PLOT(res1, data1)

#xis = res1[["param.maxpost"]]
#xis = xis.extractor(param = xis)
#plot(xis, type = "l")
