# install.packages("timedeppar")
# install.packages("readxl")

# libraries
library(readxl)
library(timedeppar)

# reading data
data   = read_excel("~/Scrivania/MOD B/esercizi/github/Project/data/syntheticData.xlsx", skip = 1)
data   = as.data.frame(data[c(1,2,4,5,7,8)])
colnames(data) = c("t", "y", "t", "y", "t", "y")
data   = na.omit(data)
data1  = data[1:2]
data2  = data[3:4]
data3  = data[5:6]
n      = dim(data)[1]
n_iter = 10000
res    = 200
count  = 0
delta  = n_iter / res


# maker of ts vector
make_ts <- function(as_logic = FALSE, as_numeric = FALSE, param = NA) {
	if (as_numeric) {
		ts        <- paste0("t", 1:n)
		ts        <- as.numeric(param[ts])
	} else {
		if (as_logic) {
			ts        <- rep(TRUE, n)
			names(ts) <- paste0("t", 1:n)
		} else {
			ts        <- as.list(1:n)
			names(ts) <- paste0("t", 1:n)
		}
	}
	return(ts)
}

# maker of vector of initial parameters
param_ini   = list(A        = 10,
				   #omega    = 11,
				   phi      = 2,
				   sigma_y  = 0.05,
				   tau      = 1,
				   mu       = 0,
				   sigma_xi = 0.05)
param_ini   = c(param_ini,   make_ts(as_logic   = FALSE,
								     as_numeric = FALSE))

# turning off working with log-scale inside infer.timedeppar
param_logic = list(A        = FALSE,
				   #omega    = FALSE,
				   phi      = FALSE,
				   sigma_y  = FALSE,
				   tau      = FALSE,
				   mu       = FALSE,
				   sigma_xi = FALSE)
param_logic = c(param_logic, make_ts(as_logic   = TRUE,
								     as_numeric = FALSE))
param_logic = unlist(param_logic)
#################################


# loglikelihood function
loglikeli <- function(param, data) {
	# feedback
	count <<- count + 1
	if (count %% delta == 0) {
		prog = (count * 100) / (delta * res)
		print(paste(prog, "%"))
	}
	#
	#t0        = data$t - 1    # t_{i-1}
	#t1        = data$t        # t_{i}
	#t2        = data$t + 1    # t_{i+1}
	y         = data$y
	ts        = make_ts(as_numeric = TRUE, param = param)
	t0        = ts - 1
	t1        = ts
	t2        = ts + 1
	A         = param$A
	#omega     = param$omega
	omega     = 11
	phi       = param$phi
	sigma_y   = param$sigma_y
	tau       = param$tau
	mu        = param$mu
	sigma_xi  = param$sigma_xi
	#
	# epsilon_{y,  i} = [y_i - (A * cos(omega * t_i + phi))] / sigma_y
	y_likeli  = y - (A * cos(omega * t1 + phi)) / sigma_y
	y_likeli  = sum(dnorm(y_likeli,  0, 1, TRUE))
	#
	# epsilon_{xi, i} = C * sqrt{tau / 2} / sigma_xi
	# where C = xi_{i+1} - xi_i - (mu - xi_i) / tau
	# or    C = t_{i+1} + (1/tau - 2) * t_i + (1 - 1/tau) * t_{i-1} - mu/tau
	xi_lideli = (t2 + ((1 / tau) - 2) * t1 + (1 - (1 / tau)) * t0 - (mu / tau)) * sqrt(tau / 2) / sigma_xi
	xi_lideli = sum(dnorm(xi_lideli, 0, 1, TRUE))
	#
	return(y_likeli + xi_lideli)
}

logprior <- function(param) {
	A         = param_ini$A
	#omega     = param_ini$omega
	#omega     = 11
	phi       = param_ini$phi
	sigma_y   = param_ini$sigma_y
	sigma_xi  = param_ini$sigma_xi
	tau       = param_ini$tau
	mu        = param_ini$mu
	#
	return(dnorm(param["A"],              A,    A / 3, TRUE) +   # A
		   dnorm(param["phi"],           phi, phi / 2, TRUE) +   # phi
		   dnorm(param["sigma_y"],   sigma_y,     0.5, TRUE) +   # sigma_y
		   dnorm(param["sigma_xi"], sigma_xi,     0.5, TRUE) +   # sigma_xi
		   dnorm(param["tau"],           tau,       1, TRUE) +   # tau
		   dnorm(param["mu"],             mu,       1, TRUE))    # mu
}

CALL <- function(data) {
	res <- infer.timedeppar(loglikeli      = loglikeli,
                            param.ini      = param_ini,
                            param.log      = param_logic,
                            #param.ou.fixed = c(omega = 11),
						    param.logprior = logprior,
                            n.iter         = n_iter,
                            data           = data,
						    verbose        = 0)
	return(res)
}


res1 = CALL(data1)
res2 = CALL(data2)
res3 = CALL(data3)
