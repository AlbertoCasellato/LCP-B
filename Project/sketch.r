# install.packages("timedeppar")
# install.packages("readxl")

# libraries
set.seed(42)
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
omega  = 2 * pi / 11
n_iter = 750             # default: 10000
res    = 10              # default: 200
count  = 0
delta  = n_iter / res
work_with_xis = TRUE     # default: TRUE


# maker of ts vector (WORKING WITH ts)
make_ts <- function(as_logic   = FALSE,
					as_range   = FALSE,
					as_numeric = FALSE,
					param      = NA) {
	if (as_numeric) {
		ts            <- paste0("t", 1:n)
		ts            <- as.numeric(param[ts])
	} else {
		if (as_logic) {
			ts        <- rep(FALSE, n)
			names(ts) <- paste0("t", 1:n)
		} else if (as_range) {
			ts        <- mapply(function(i) c(i + 0.2, i + 1.8), 0:(n - 1), SIMPLIFY = FALSE)
			names(ts) <- paste0("t", 1:n)
		} else {
			ts        <- as.list(1:n)
			names(ts) <- paste0("t", 1:n)
		}
	}
	return(ts)
}

# maker of xis vector (WORKING WITH xis)
make_xis <- function(as_logic   = FALSE,
					 as_range   = FALSE,
					 as_numeric = FALSE,
					 param      = NA) {
	if (as_numeric) {
		xis            <- paste0("xi", 2:n)
		xis            <- as.numeric(param[xis])
	} else {
		if (as_logic) {
			xis        <- rep(FALSE, (n - 1))
			names(xis) <- paste0("xi", 2:n)
		} else if (as_range) {
			xis        <- lapply(2:n, function(i) c(0.4, 1.6))
			names(xis) <- paste0("xi", 2:n)
		} else {
			xis        <- as.list(rep(1, (n - 1)))
			names(xis) <- paste0("xi", 2:n)
		}
	}
	return(xis)
}

if (work_with_xis) {
	maker = make_xis
} else {
	maker = make_ts
}


# maker of vector of initial parameters
param_ini   = list(A        = 10,
				   phi      = 2,
				   sigma_y  = 0.05,
				   tau      = 1,
				   mu       = 0,
				   sigma_xi = 0.05)
param_ini   = c(param_ini,   maker(as_logic   = FALSE,
								   as_range   = FALSE,
								   as_numeric = FALSE))

# turning off working with log-scale inside infer.timedeppar
param_logic = list(A        = FALSE,
				   #omega    = FALSE,
				   phi      = FALSE,
				   sigma_y  = FALSE,
				   tau      = FALSE,
				   mu       = FALSE,
				   sigma_xi = FALSE)
param_logic = c(param_logic, maker(as_logic   = TRUE,
								   as_range   = FALSE,
								   as_numeric = FALSE))
param_logic = unlist(param_logic)

# ranges
param_ranges = list(A        = c(5,    15),
				    phi      = c(0,    2 * pi),
				    sigma_y  = c(0.01, 10),
				    tau      = c(0.05, 20),
				    mu       = c(-10,  10),
				    sigma_xi = c(0.01, 10))
param_ranges = c(param_ranges, maker(as_logic   = FALSE,
								     as_range   = TRUE,
								     as_numeric = FALSE))
#########################################################


# loglikelihood function
loglikeli <- function(param) {
	# feedback
	count <<- count + 1
	#print(count)
	if (count %% delta == 0) {
		prog = (count * 100) / (delta * res)
		print(paste(prog, "%"))
	}
	#
	#t0        = data$t - 1    # t_{i-1}
	#t1        = data$t        # t_{i}
	#t2        = data$t + 1    # t_{i+1}
	#y         = data$y
	#ts        = make_ts(as_numeric = TRUE, param = param)
	#t0        = ts - 1
	#t1        = ts
	#t2        = ts + 1
	xis       = make_xis(as_numeric = TRUE, param = param)
	xi0       = xis[1:(n - 2)]
	xi1       = xis[2:(n - 1)]
	A         = param$A
	#omega     = param$omega
	phi       = param$phi
	sigma_y   = param$sigma_y
	tau       = param$tau
	mu        = param$mu
	sigma_xi  = param$sigma_xi
	#
	# WORKING WITH ts:
	# epsilon_{y,  i} = [y_i - (A * cos(omega * t_i + phi))] / sigma_y
	#y_likeli  = (y - A * cos(omega * t1 + phi)) / sigma_y
	#
	# WORKING WITH xis:
	#y_likeli  = y - (A * cos(omega * (sum_{j=i}^{2}(xi_j) + 1) + phi)) / sigma_y
	y_likeli  = (y - A * cos(omega * c(1, cumsum(xis) + 1) + phi)) / sigma_y
	y_likeli  = sum(dnorm(y_likeli,  0, 1, TRUE))
	#
	# epsilon_{xi, i} = C * sqrt{tau / 2} / sigma_xi
	# WORKING WITH ts:
	# C = t_{i+1} + (1/tau - 2) * t_i + (1 - 1/tau) * t_{i-1} - mu/tau
	#xi_likeli = (t2 + ((1 / tau) - 2) * t1 + (1 - (1 / tau)) * t0 - (mu / tau)) * sqrt(tau / 2) / sigma_xi
	#
	# WORKING WITH xis:
	# C = xi_{i+1} - xi_i - (mu - xi_i) / tau
	#xi_likeli = (xi1 - xi0 - (1 / tau) * (mu - xi0)) * sqrt(tau / 2) / sigma_xi
	#xi_likeli = sum(dnorm(xi_likeli, 0, 1, TRUE))
	#
	#print(y_likeli + xi_likeli)
	#return(y_likeli + xi_likeli)
	return(y_likeli)
}

logprior <- function(param) {
	A         = param_ini$A
	#omega     = param_ini$omega
	phi       = param_ini$phi
	sigma_y   = param_ini$sigma_y
	sigma_xi  = param_ini$sigma_xi
	tau       = param_ini$tau
	mu        = param_ini$mu
	#
	prior     = dnorm(param["A"],              A,    A / 2, TRUE) +   # A
				dnorm(param["phi"],           phi,     phi, TRUE) +   # phi
				dnorm(param["sigma_y"],   sigma_y,       1, TRUE) +   # sigma_y
				dnorm(param["sigma_xi"], sigma_xi,       1, TRUE) +   # sigma_xi
				dnorm(param["tau"],           tau,       2, TRUE) +   # tau
				dnorm(param["mu"],             mu,       2, TRUE)# +   # mu
				#sum(dnorm(make_xis(as_numeric = TRUE,
				#			       param      = param), 1, 1, TRUE))  # xis
	#print(prior)
	return(prior)
}

CALL <- function() {
	res <- infer.timedeppar(loglikeli           = loglikeli,
							loglikeli.keepstate = FALSE,
							param.ini           = param_ini,
							param.range         = param_ranges,
							param.log           = param_logic,
							param.logprior      = logprior,
							task                = "start",
							n.iter              = n_iter,
							verbose             = 0)
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
	xis    = make_xis(as_numeric = TRUE, param = params)
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



xis = res2[["param.maxpost"]]
xis = make_xis(as_numeric = TRUE, param = xis)
