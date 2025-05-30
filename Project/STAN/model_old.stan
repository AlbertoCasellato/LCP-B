data {
    int<lower=1>                  n;             // number of armonics
    int<lower=1>                  l;             // number of armonics
    vector[n]                     y;             // data points
    real<lower=0.20,  upper=0.25> t0;            // starting time
    real<lower=0.005, upper=1000> omega_fixed;   // omega fixed (first one for convention)
    vector<lower=0>[l - 1]        omega_min;     // omega inferior limit
    vector<lower=0>[l - 1]        omega_max;     // omega superior limit
    vector<lower=0>[l]            A_min;         // A inferior limit
    vector<lower=0>[l]            A_max;         // A superior limit
    vector<lower=0>[l]            phi_min;       // phi inferior limit
    vector<lower=0>[l]            phi_max;       // phi superior limit
}

parameters {
    vector<lower=0.01,  upper=5>[l]         A;          // amplitudes
    vector<lower=0.01,  upper=6.3>[l]       phi;        // phases
                                                        // comment the following block to turn off
    vector<lower=0.015, upper=1.6>[l - 1]   omega;      // others omega (T from 0.5y to 120ky)
    real<lower=0.01,    upper=8>            sigma_y;    // sigma of ys
    real<lower=0.01,    upper=320>          tau;        // inverse of gamma
    real<lower=0.15,    upper=0.25>         mu;         // average
    real<lower=0.005,   upper=3.5>          sigma_xi;   // sigma of xis
    vector<lower=0.1,   upper=0.3>[n - 1]   xi;         // xis
}

model {
    vector[l] M_1;
    vector[l] M_k;
    real      sxi;
    real      uat;
    //
    // priors [1/2]
    for (i in 1:l) {
        A[i]         ~ uniform(A_min[i],   A_max[i]);
        phi[i]       ~ uniform(phi_min[i], phi_max[i]);
    }
    sigma_y          ~ uniform(0.01,    8);
    tau              ~ uniform(0.1,   300);
    mu               ~ uniform(0.16, 0.24);
    sigma_xi         ~ uniform(0.005,   3);
    xi[1]            ~ uniform(0.18, 0.23);
    //
    // priors [2/2] + distributions for y[n] [1/2]
    M_1[1]           = A[1] * cos(omega_fixed * t0  + phi[1]);
    for (j in 2:l) {                                                  // comment this block to turn off
        omega[j - 1] ~ uniform(omega_min[j - 1], omega_max[j - 1]);   // same
        M_1[j]       = A[j] * cos((omega[j - 1] * t0) + phi[j]);      // same
    }                                                                 // same
    //
    // distributions for y[n] [2/2]
    y[1]             ~ normal(sum(M_1[1:l]), sigma_y);
    for (k in 2:n) {
        M_k[1]       = A[1] * cos(omega_fixed  * (t0 + sum(xi[1:(k - 1)])) + phi[1]);
        for (w in 2:l) {                                                                // comment this block to turn off
            M_k[w]   = A[w] * cos(omega[w - 1] * (t0 + sum(xi[1:(k - 1)])) + phi[w]);   // same
        }                                                                               // same
        y[k]         ~ normal(sum(M_k[1:l]), sigma_y);
    }
    //
    // distributions for xi[n-1]
    sxi              = sigma_xi * sqrt(2 / tau);
    uat              = (1 / tau);
    for (z in 1:(n - 2)) {
        xi[z + 1]    ~ normal((xi[z] + uat * (mu - xi[z])), sxi);
    }
}
