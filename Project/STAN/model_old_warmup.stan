data {
    int<lower=1>     n;
    vector[n]        y;
    real<lower=0.01> omega;
}

parameters {
    real<lower=0.01> A;
    real<lower=0.01> phi;
    real<lower=0.01> sigma_y;
    real<lower=0.01> tau;
    real<lower=0.01> mu;
    real<lower=0.01> sigma_xi;
    vector<lower=0.01, upper=1.9>[n - 1] xi;
}

model {
    // vector[n - 2] xi_eps;
    //
    // priors
    A        ~ uniform(2,      18);
    phi      ~ uniform(0.01,    6.28318);
    sigma_y  ~ uniform(0.01,    8);
    mu       ~ uniform(0.2,     1.8);
    sigma_xi ~ uniform(0.01,    4);
    tau      ~ uniform(0.1,   300);
    //
    // distributions for y[n]
    y[1]     ~ normal((A * cos(omega + phi)), sigma_y);
    for (i in 2:n) {
        y[i] ~ normal((A * cos(omega * (1 + sum(xi[1:(i - 1)])) + phi)), sigma_y);
    }
    //
    // distributions for xi[n-1]
    // xi[1] ~ uniform(0.01, 2);
    for (i in 1:(n - 2)) {
        xi[i + 1] ~ normal((xi[i] + (1 / tau) * (mu - xi[i])), sigma_xi * sqrt(2 / tau));
        // xi[i + 1] ~ uniform(0.01, 2);
        // xi_eps[i] = (xi[i + 1] - xi[i] - (1 / tau) * (mu - xi[i])) * sqrt(tau / 2);
        // xi_eps[i] ~ normal(0, sigma_xi);
    }
}
