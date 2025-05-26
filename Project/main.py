import matplotlib.pyplot as plt
import cmdstanpy
import pandas as pd
import numpy  as np
import pickle
from   math import pi

SAVE = lambda data, filepath: pickle.dump(data, open(filepath + ".pkl", "wb"))   # save file
LOAD = lambda filepath:       pickle.load(open(filepath       + ".pkl", "rb"))   # read file

# control parameter
isOnCVVM = False
n_sigma  = 4
#
if isOnCVVM:
    data_file = "~/data/solar_data.xlsx"
    harm_file = "~/data/solar_harmonics.xlsx"
    stan_file = "~/STAN/model.stan"
    iter_s    =  50000      # OR:  50000
    iter_w    = 150000      # OR: 150000
    cs        = 2         # OR: 4
    #ths_per_c = 2         # OR: 2
    parl_cs   = 2
else:
    data_file = "~/Scrivania/MOD B/esercizi/github/Project/data/solar_data.xlsx"
    harm_file = "~/Scrivania/MOD B/esercizi/github/Project/data/solar_harmonics.xlsx"
    stan_file = "~/Scrivania/MOD B/esercizi/github/Project/STAN/model.stan"
    iter_s    = 50
    iter_w    = 50
    cs        = 1
    #ths_per_c = 2
    parl_cs   = 1
#################

# compile model
model        = cmdstanpy.CmdStanModel(stan_file   = stan_file,
                                      cpp_options = {"STAN_THREADS" : True})

# read data
data         = pd.read_excel(data_file, skiprows = 0)
data         = data.iloc[4:, [1, 2]]
data.columns = ["age", "y"]
data["age"]  = data["age"] + 55
data         = data[data["age"] < 30000 + 55]
data["y"]    = data["y"] - data["y"].mean()
data.reset_index(drop = True, inplace = True)
data.index   = data.index + 1
n            = data.shape[0]
#
# read harmonics
harm         = pd.read_excel(harm_file, skiprows = 0)
harm         = harm.iloc[:, [1,2,3,4,5,6]]
harm.columns = list(harm.iloc[3])
harm         = harm.iloc[4:]
harm.reset_index(drop = True, inplace = True)
harm.index   = harm.index + 1
harm.loc[1, "sigma_phase"] = harm["Phase"].iloc[0] * (harm["sigma_phase"] / harm["Phase"]).iloc[0:21].mean()
#
# make intervals
harm["omega"]     = 2 * np.pi /  harm["Period"]
harm["omega_min"] = 2 * np.pi / (harm["Period"] + (n_sigma * harm["sigma_period"]))
harm["omega_max"] = 2 * np.pi / (harm["Period"] - (n_sigma * harm["sigma_period"]))
harm              = harm.iloc[:, 2:]
harm["A"]         = harm["Amplitude"]
harm["A_min"]     = harm["Amplitude"]           - (n_sigma * harm["sigma_amplitude"])
harm["A_max"]     = harm["Amplitude"]           + (n_sigma * harm["sigma_amplitude"])
harm              = harm.iloc[:, 2:]
harm["phi"]       = harm["Phase"]
harm["phi_min"]   = harm["Phase"]               - (n_sigma * harm["sigma_phase"])
harm["phi_max"]   = harm["Phase"]               + (n_sigma * harm["sigma_phase"])
harm.loc[harm["phi_max"] > 2 * pi, "phi_max"]   = 2 * pi
harm              = harm.iloc[:, 2:]
harm[harm < 0]    = 0.01
#
# make order
harm              = pd.concat([harm.iloc[8:9], harm.iloc[:8], harm.iloc[9:]])
harm.reset_index(drop = True, inplace = True)
harm.index        = harm.index + 1
l                 = harm.shape[0]
#
# make xis
ts  = data.iloc[:, [0]].values.reshape(n)
xis = ts[1:] - ts[:-1]
# np.insert(np.cumsum(xis) + ts[0], 0, ts[0])


# select data
data = {
        "n"           : n,
        "l"           : l,
        "y"           : data.iloc[:, [1]].values.reshape(n),
        "t0"          : ts[0],
        "omega_fixed" : harm["omega"].iloc[0],
        "omega_min"   : harm["omega_min"].iloc[1:].values.reshape(l - 1),
        "omega_max"   : harm["omega_max"].iloc[1:].values.reshape(l - 1),
        "A_min"       : harm["A_min"].values.reshape(l),
        "A_max"       : harm["A_max"].values.reshape(l),
        "phi_min"     : harm["phi_min"].values.reshape(l),
        "phi_max"     : harm["phi_max"].values.reshape(l)
}

# starting points
init = {
    "A"        : harm["A"].values.reshape(l),
    "phi"      : harm["phi"].values.reshape(l),
    "omega"    : harm["omega"].iloc[1:].values.reshape(l - 1),
    "sigma_y"  :  2.0,
    "tau"      :  1.0,
    "mu"       : 20,
    "sigma_xi" :  3.0,
    "xi"       : xis
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

SAVE(res, "~/res")
