import cmdstanpy
import pandas as pd
import numpy  as np

data_file = "~/Scrivania/MOD B/esercizi/github/Project/data/syntheticData.xlsx"
stan_file = "~/Scrivania/MOD B/esercizi/github/Project/STAN/model.stan"
model     = cmdstanpy.CmdStanModel(stan_file = stan_file)

# read data
data  = pd.read_excel(data_file, skiprows = 1)
data  = data.iloc[1:, [0, 1, 3, 4, 6, 7]]
data.columns = ["t", "y", "t", "y", "t", "y"]
data  = data.dropna()
data.reset_index(drop = True, inplace = True)
data.index   = data.index + 1
data1 = data.iloc[:, [0, 1]]
data2 = data.iloc[:, [2, 3]]
data3 = data.iloc[:, [4, 5]]
n     = data.shape[0]

# select data
data = {
        "n"     : n,
        "y"     : data1.iloc[:, [1]].values.reshape(n),
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
risultati = model.sample(data            = data,
                         inits           = init,
                         chains          = 4,
                         parallel_chains = 1,
                         seed            = 42,
                         show_console    = True)

# Visualizza i risultati
print(risultati)

# Estrai i campioni per i parametri di interesse
A_samples        = risultati.stan_variable('A')
phi_samples      = risultati.stan_variable('phi')
sigma_y_samples  = risultati.stan_variable('sigma_y')
tau_samples      = risultati.stan_variable('tau')
mu_samples       = risultati.stan_variable('mu')
sigma_xi_samples = risultati.stan_variable('sigma_xi')
xi_samples       = risultati.stan_variable('xi')
