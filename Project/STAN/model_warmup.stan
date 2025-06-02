data {
    int<lower=1>                             n;             // number of armonics
    int<lower=1>                             l;             // number of armonics
    vector[n]                                y;             // data points
    real<lower=0>                            t0_lower;
    real<lower=0>                            t0_upper;
    real<lower=0>                            om_f_lower;
    real<lower=0>                            om_f_upper;
    real<lower=t0_lower,   upper=t0_upper>   t0;            // starting time
    real<lower=om_f_lower, upper=om_f_upper> omega_fixed;   // omega fixed (first one for convention)
    vector<lower=0>[l]                       A_ini;         // starting values of A                     // A
    vector<lower=0>[l]                       sigma_A;       // sigmas of A
    vector<lower=0>[l]                       phi_ini;       // starting values of phi                   // phi
    vector<lower=0>[l]                       sigma_phi;     // sigmas of phi
    real<lower=0>                            xi_lower;
    real<lower=0>                            xi_upper;
    real<lower=0>                            s_xi_lower;
    real<lower=0>                            s_xi_upper;
    real<lower=0>                            mu_lower;
    real<lower=0>                            mu_upper;
    real<lower=0>                            A_lower;
    real<lower=0>                            A_upper;
    real<lower=0>                            om_lower;
    real<lower=0>                            om_upper;
    real<lower=0>                            s_y_lower;
    real<lower=0>                            s_y_upper;
    real<lower=0>                            mu_mean;
    real<lower=0>                            mu_std;
    real<lower=0>                            s_xi_mean;
    real<lower=0>                            s_xi_std;
    real<lower=0>                            xi1_mean;
    real<lower=0>                            xi1_std;
    real<lower=0>                            s_y_mean;
    real<lower=0>                            s_y_std;
}

parameters {
    vector<lower=A_lower,  upper=A_upper>[l]      A;          // amplitudes
    vector<lower=0.001,    upper=6.284>[l]        phi;        // phases
    real<lower=s_y_lower,  upper=s_y_upper>       sigma_y;    // sigma of y
    real<lower=0.01,       upper=350>             tau;        // inverse of gamma
    real<lower=mu_lower,   upper=mu_upper>        mu;         // average of xi
    real<lower=s_xi_lower, upper=s_xi_upper>      sigma_xi;   // sigma of xi
    vector<lower=xi_lower, upper=xi_upper>[n - 1] xi;         // xis
}

transformed parameters {
    vector[n]     t;            // absolute timestamps
    vector[n]     y_hat;        // pre-computed y model
    vector[n - 2] xi_hat;       // pre-computed xi model
    real eff_sigma_xi;          // effective sigma xi
    real gamma = 1 / tau;       // gamma parameter
    //
    // make ts from xi
    t[1]            = t0;
    for (i in 2:n) {
        t[i]        = t[i - 1] + xi[i - 1];
    }
    //
    // pre-compute y_hat[n]
    for (j in 1:n) {
        y_hat[j]    = dot_product(A, cos((omega_fixed .* t[j]) + phi));
    }
    //
    // pre-compute xi_hat[n - 2] and eff_sigma_xi
    eff_sigma_xi    = sigma_xi * sqrt(2 * gamma);
    xi_hat          = xi[1:(n - 2)] + mu * gamma - xi[1:(n - 2)] .* gamma;
}

model {
    tau   ~ uniform(0.05, 300);
    //
    A     ~ lognormal(A_ini,     sigma_A);
    phi   ~ lognormal(phi_ini,   sigma_phi);
    //
    // model for xi
    mu            ~ lognormal(mu_mean,   mu_std);
    sigma_xi      ~ lognormal(s_xi_mean, s_xi_std);
    xi[1]         ~ lognormal(xi1_mean,  xi1_std);
    xi[2:(n - 1)] ~ lognormal(xi_hat,    eff_sigma_xi);
    //
    // model for y
    sigma_y  ~ lognormal(s_y_mean, s_y_std);
    y        ~ normal(y_hat,       sigma_y);
}

generated quantities {
    vector[n] y_rep;
    //
    for (k in 1:n) {
        y_rep[k] = normal_rng(y_hat[k], sigma_y);
    }
}
