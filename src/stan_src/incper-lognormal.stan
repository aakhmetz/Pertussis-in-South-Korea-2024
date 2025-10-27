data {
    int<lower = 1> N; 
    array[N] real exp2onset_R;
    array[N] real<upper = exp2onset_R> exp2onset_L;
}

parameters {
	vector<lower = 0, upper = 1>[N] exposure_raw, onset_raw; // non-transformed exposure and onset times
	real log_mean_incper;
    real<lower = 0> sqr_cv_incper;
}

transformed parameters {
    real mean_incper = exp(log_mean_incper); // mean of the time intervals

    real param1_incper, param2_incper;
    real cv_incper = sqrt(sqr_cv_incper),
        sd_incper = mean_incper .* cv_incper,
        log_sd_incper = log(sd_incper);

    // Lognormal distribution
    {
        real sqr_param2_lognormal = log1p(sqr_cv_incper);
        param2_incper = sqrt(sqr_param2_lognormal);
        param1_incper = log_mean_incper - sqr_param2_lognormal / 2.0;
    }
}

model {
    // priors
	mean_incper ~ normal(7, 7);
	sqr_cv_incper ~ gamma(2, .5);

	onset_raw ~ beta(1, 1); 
	exposure_raw ~ beta(1, 1);
  
    vector[N] onset = to_vector(exp2onset_R) + onset_raw, exposure; 
    for (n in 1 : N) {
        real exposure_R = exp2onset_R[n] - exp2onset_L[n] + 1; 
        exposure_R = onset[n] < exposure_R ? onset[n] : exposure_R;
        exposure[n] = exposure_R * exposure_raw[n];
    }

	/* likelihood */
    target += lognormal_lpdf(onset - exposure | param1_incper, param2_incper);
}

generated quantities {
    real incper_pred  = lognormal_rng(param1_incper, param2_incper);
}