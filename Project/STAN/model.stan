data {
    int<lower=1>                           n;             // number of armonics
    int<lower=1>                           l;             // number of armonics
    vector[n]                              y;             // data points
    real<lower=0.20,    upper=0.25>        t0;            // starting time
    real<lower=0.005,   upper=1000>        omega_fixed;   // omega fixed (first one for convention)
    vector<lower=0>[l - 1]                 omega_min;     // omega inferior limit
    vector<lower=0>[l - 1]                 omega_max;     // omega superior limit
    vector<lower=0>[l]                     A_min;         // A inferior limit
    vector<lower=0>[l]                     A_max;         // A superior limit
    vector<lower=0>[l]                     phi_min;       // phi inferior limit
    vector<lower=0>[l]                     phi_max;       // phi superior limit
}

parameters {
    vector<lower=0.01,  upper=5>[l]        A;             // amplitudes
    vector<lower=0.01,  upper=6.3>[l]      phi;           // phases
                                                          // comment the following block to turn off
    vector<lower=0.015, upper=1.6>[l - 1]  omega;         // others omega (T from 0.5y to 120ky)
    real<lower=0.01,    upper=8>           sigma_y;       // sigma of ys
    real<lower=0.01,    upper=320>         tau;           // inverse of gamma
    real<lower=0.15,    upper=0.25>        mu;            // average
    real<lower=0.005,   upper=3.5>         sigma_xi;      // sigma of xis
    vector<lower=0.1,   upper=0.3>[n - 1]  xi;            // xis
}

transformed parameters {
    vector[n]     t;                // absolute timestamps
    vector[n]     y_hat;            // mean function at each t
    vector[n - 2] xi_hat;           // average of xi
    vector[l]     omega_full;       // concatenate fixed + free omegas
    real rho = 1 / tau;
    //
    t[1]            = t0;
    for (k in 2:n)
        t[k]        = t[k - 1] + xi[k - 1];
    //
    omega_full[1]   = omega_fixed;
    omega_full[2:l] = omega;
    //
    // pre-compute y_hat[n]
    for (k in 1:n)
        y_hat[k]    = dot_product(A, cos((omega_full .* t[k]) + phi));
    //
    // pre-compute hi_hat[n -2]
    xi_hat          = xi[1:(n - 2)] + mu * rho - xi[1:(n - 2)] .* rho;
}

model {
    A        ~ uniform(A_min,     A_max);
    phi      ~ uniform(phi_min,   phi_max);
    omega    ~ uniform(omega_min, omega_max);
    xi[1]    ~ uniform(0.18, 0.23);
    //
    sigma_y  ~ uniform(0.01,    8);
    tau      ~ uniform(0.1,   300);
    mu       ~ uniform(0.16, 0.24);
    sigma_xi ~ uniform(0.005,   3);
    //
    {
        real xis = sigma_xi * sqrt(2 / tau);
        xi[2:(n - 1)] ~ normal(xi_hat, xis);
    }
    //
    y ~ normal(y_hat, sigma_y);
}

generated quantities {
    vector[n] y_rep;
    //
    for (k in 1:n)
        y_rep[k] = normal_rng(y_hat[k], sigma_y);
}
