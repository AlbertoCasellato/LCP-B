// Optimized harmonic mixture with AR(1) time increments
// Carlo Albert project – v2 (2025‑05‑27)
//
// Key changes vs. reference implementation
//   •  O(n) instead of O(n²): pre‑compute cumulative times once.
//   •  Vectorised priors and likelihood; single y ~ normal() statement.
//   •  Dot‑product formulation removes inner harmonic loop.
//   •  All heavy algebra in transformed parameters, leaving a compact model.
//   •  Added generated quantities block for PPC / LOO.
//
// ────────────────────────────────────────────────────────────────────────────────

data {
    int<lower=1>                  n;             // number of observations
    int<lower=1>                  l;             // number of harmonics
    vector[n]                     y;             // observed values
    real<lower=0.20,  upper=0.25> t0;            // initial time
    real<lower=0.005, upper=1000> omega_fixed;   // first omega (convention)
    vector<lower=0>[l-1]          omega_min;     // lower bounds for free omegas
    vector<lower=0>[l-1]          omega_max;     // upper bounds for free omegas
    vector<lower=0>[l]            A_min;         // lower bounds for amplitudes
    vector<lower=0>[l]            A_max;         // upper bounds for amplitudes
    vector<lower=0>[l]            phi_min;       // lower bounds for phases
    vector<lower=0>[l]            phi_max;       // upper bounds for phases
}

parameters {
    vector<lower=0.01,  upper=5>[l]      A;          // amplitudes
    vector<lower=0.01,  upper=2*pi()>[l] phi;        // phases (wrapped to 2π)
    vector<lower=0.015, upper=1.6>[l-1]  omega;      // additional omegas
    real<lower=0.01,    upper=8>         sigma_y;    // observation SD
    real<lower=0.1,     upper=300>       tau;        // AR(1) time‑scale (1/γ)
    real<lower=0.15,    upper=0.25>      mu;         // xi long‑term mean
    real<lower=0.005,   upper=3.5>       sigma_xi;   // xi innovation SD
    vector<lower=0.1,   upper=0.3>[n-1]  xi;         // time increments
}

transformed parameters {
    vector[n] t;                // absolute timestamps
    vector[n] y_hat;            // mean function at each t
    vector[l] omega_full;       // concatenate fixed + free omegas
    //
    // ── build timeline ──
    t[1] = t0;
    for (k in 2:n)
        t[k] = t[k-1] + xi[k-1];
    //
    // ── merge omega vectors ──
    omega_full[1]   = omega_fixed;
    omega_full[2:l] = omega;
    //
    // ── expected signal ──
    for (k in 1:n)
        y_hat[k] = dot_product(A, cos(omega_full .* t[k] + phi));
}

model {
    // ── priors / constraints supplied via bounds ──
    A     ~ uniform(A_min,     A_max);
    phi   ~ uniform(phi_min,   phi_max);
    omega ~ uniform(omega_min, omega_max);
    xi[1] ~ uniform(0.18, 0.23);
    //
    // AR(1) dynamics for xi
    {
        real sxi = sigma_xi * sqrt(2 / tau);
        real rho = 1 / tau;
        for (z in 1:(n-2))
            xi[z+1] ~ normal(xi[z] + rho * (mu - xi[z]), sxi);
    }
    //
    // likelihood
    y ~ normal(y_hat, sigma_y);
}

generated quantities {
    vector[n] y_rep = normal_rng(y_hat, sigma_y); // posterior predictive draw
}
