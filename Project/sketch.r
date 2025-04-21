# install.packages("readxl")

library(readxl)
library(timedeppar)

data  = read_excel("~/Scrivania/MOD B/esercizi/github/Project/data/syntheticData.xlsx", skip = 1)
data  = as.data.frame(data[c(1,2,4,5,7,8)])
colnames(data) = c("t", "y", "t", "y", "t", "y")
data  = na.omit(data)
data1 = data[1:2]
data2 = data[3:4]
data3 = data[5:6]


n          <- dim(data)[1]
tmp        <- as.list(1:n)
names(tmp) <- paste0("t", 1:n)
param_ini = list(A        = 10,
				 #omega    = 11,
				 phi      = 2,
				 sigma_y  = 0.05,
				 tau      = 0.5,
				 mu       = 0,
				 sigma_xi = 0.05)
param_ini = c(param_ini, tmp)

loglikeli <- function(param, data)
{
	#t0        = data$t - 1    # t_{i-1}
	#t1        = data$t        # t_{i}
	#t2        = data$t + 1    # t_{i+1}
	nomi      <- paste0("t", 1:500)
	richiamo  <- as.list(c(param[nomi]))
	richiamo  = as.numeric(unlist(richiamo))
	t0        = richiamo - 1
	t1        = richiamo
	t2        = richiamo + 1
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
	y_likeli  = log(20 + abs(sum(data$y - (A * cos(omega * t1 + phi)) / sigma_y)))
	#
	# epsilon_{xi, i} = C * sqrt{tau / 2} / sigma_xi
	# where C = xi_{i+1} - xi_i - (mu - xi_i) / tau
	# or    C = t_{i+1} + (1/tau - 2) * t_i + (1 - 1/tau) * t_{i-1} - mu/tau
	xi_lideli = log(20 + abs(sum((t2 + ((1 / tau) - 2) * t1 + (1 - (1 / tau)) * t0 - (mu / tau)) * sqrt(tau / 2) / sigma_xi)))
	#
	return(y_likeli + xi_lideli)
}

logprior <- function(param)
{
	A         = param_ini$A
	#omega     = param$omega
	#omega     = 11
	phi       = param_ini$phi
	#sigma_y   = param_ini$sigma_y
	tau       = param_ini$tau
	mu        = param_ini$mu
	#sigma_xi  = param_ini$sigma_xi
	nomi      <- paste0("t", 1:500)
	richiamo  <- as.list(c(param[nomi]))
	richiamo  = as.numeric(unlist(richiamo))
	return(log(20 + rnorm(1, A, A / 3))     +   # A
	       log(20 + rnorm(1, phi, phi / 3)) +   # phi
	       log(20 + rnorm(1, 0, 0.05))  +   # sigma_y
	       log(20 + rnorm(1, 0, 0.05))  +   # sigma_xi
		   log(20 + rnorm(1, tau, tau / 3)) +   # tau
		   log(20 + rnorm(1,  mu,  mu / 3)) +   # mu
	       sum(20 + log(rnorm(500, richiamo, 0.3))))
}


res1 <- infer.timedeppar(loglikeli      = loglikeli,
                         param.ini      = param_ini,
                         param.log      = c(A       = FALSE,
										   #omega    = FALSE,
										   phi      = FALSE,
										   sigma_y  = FALSE,
										   tau      = FALSE,
										   mu       = FALSE,
										   sigma_xi = FALSE),
                         #param.ou.fixed = c(omega   = 11),
						 param.logprior = logprior,
                         n.iter         = 10000,
                         data           = data1)


res2 <- infer.timedeppar(loglikeli      = loglikeli,
                         param.ini      = param_ini,
                         param.log      = c(A       = FALSE,
										   #omega    = FALSE,
										   phi      = FALSE,
										   sigma_y  = FALSE,
										   tau      = FALSE,
										   mu       = FALSE,
										   sigma_xi = FALSE),
                         #param.ou.fixed = c(omega   = 11),
						 param.logprior = logprior,
                         n.iter         = 10000,
                         data           = data2)


res3 <- infer.timedeppar(loglikeli      = loglikeli,
                         param.ini      = param_ini,
                         param.log      = c(A       = FALSE,
										   #omega    = FALSE,
										   phi      = FALSE,
										   sigma_y  = FALSE,
										   tau      = FALSE,
										   mu       = FALSE,
										   sigma_xi = FALSE),
                         #param.ou.fixed = c(omega   = 11),
						 param.logprior = logprior,
                         n.iter         = 10000,
                         data           = data3)
