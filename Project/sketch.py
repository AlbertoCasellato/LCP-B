import matplotlib.pyplot as plt
import cmdstanpy
import pandas as pd
import numpy  as np
import pickle

SAVE = lambda data, filepath: pickle.dump(data, open(filepath + ".pkl", "wb"))   # save file
LOAD = lambda filepath:       pickle.load(open(filepath       + ".pkl", "rb"))   # read file

# control parameter
isOnCVVM = True
#
if isOnCVVM:
    data_file = "~/data/syntheticData.xlsx"
    stan_file = "~/STAN/model.stan"
    iter_s    =  50000      # OR:  50000
    iter_w    = 150000      # OR: 150000
    cs        = 2         # OR: 4
    #ths_per_c = 2         # OR: 2
    parl_cs   = 2
else:
    data_file = "~/Scrivania/MOD B/esercizi/github/Project/data/syntheticData.xlsx"
    stan_file = "~/Scrivania/MOD B/esercizi/github/Project/STAN/model.stan"
    iter_s    = 1000
    iter_w    = 1000
    cs        = 1
    #ths_per_c = 2
    parl_cs   = 1
#################

# compile model
model        = cmdstanpy.CmdStanModel(stan_file   = stan_file,
                                      cpp_options = {"STAN_THREADS" : True})

# read data
data         = pd.read_excel(data_file, skiprows = 1)
data         = data.iloc[1:, [0, 1, 3, 4, 6, 7]]
data.columns = ["t", "y", "t", "y", "t", "y"]
data         = data.dropna()
data.reset_index(drop = True, inplace = True)
data.index   = data.index + 1
data1        = data.iloc[:, [0, 1]]
data2        = data.iloc[:, [2, 3]]
data3        = data.iloc[:, [4, 5]]
n            = data.shape[0]

# select data
data = {
        "n"     : n,
        "y"     : data3.iloc[:, [1]].values.reshape(n),
        "omega" : 2 * np.pi / 11
}

# starting points
init = {
    'A'        : 10.0,
    'phi'      :  2.0,
    'sigma_y'  :  0.05,
    'tau'      :  1.0,
    'mu'       :  1.0,
    'sigma_xi' :  0.03,
    'xi'       : [1] * (n - 1)
}

# sampling
res = model.sample(data              = data,
                   inits             = init,
                   iter_sampling     = iter_s,
                   iter_warmup       = iter_w,
                   chains            = cs,
                   #threads_per_chain = ths_per_c,
                   parallel_chains   = parl_cs,
                   seed              = 42,
                   adapt_delta       = 0.98,
                   max_treedepth     = 20,
                   metric            = "dense_e",
                   adapt_engaged     = True,
                   show_console      = True)

# results
#A        = res.stan_variable('A')
#phi      = res.stan_variable('phi')
#sigma_y  = res.stan_variable('sigma_y')
#tau      = res.stan_variable('tau')
#mu       = res.stan_variable('mu')
#sigma_xi = res.stan_variable('sigma_xi')
#xi       = res.stan_variable('xi')

# plotting
#res.summary()[1:7]
#plt.plot(res.summary()["Mean"][7:])
#plt.show()

#SAVE(res, "~/res1")
#SAVE(res, "~/res2")
SAVE(res, "~/res3")
