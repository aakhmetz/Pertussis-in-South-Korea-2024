data {
    int<lower = 1> N; 
    array[N] real exp2onset_R;
    array[N] real<upper = exp2onset_R> exp2onset_L;
}

transformed data {
    int D = 3; // number of distributions for incper and TOST (Gamma, Weibull, Lognormal)
}

parameters {
	array[D] vector<lower = 0, upper = 1>[N] exposure_raw, onset_raw; // non-transformed exposure and onset times
	real log_mean_incper;
    real<lower = 0> param1_Weibull_incper; 

    simplex[D] weight; // mixing proportions
}

transformed parameters {
    real mean_incper = exp(log_mean_incper); // mean of the time intervals

    vector[D] param1_incper, param2_incper;
    
    // Weibull distribution
    param1_incper[2] = param1_Weibull_incper;
    real sqr_cv_incper = tgamma(1 + 2 / param1_incper[2]) ./ square(tgamma(1 + 1 / param1_incper[2])) - 1,
        cv_incper = sqrt(sqr_cv_incper),
        sd_incper = mean_incper .* cv_incper,
        log_sd_incper = log(sd_incper);
    param2_incper[2] = mean_incper / tgamma(1.0 + 1.0 / param1_incper[2]);

    // Gamma distribution
    param1_incper[1] = 1 / sqr_cv_incper;
    param2_incper[1] = param1_incper[1] ./ mean_incper;

    // Lognormal distribution
    {
        real sqr_param2_lognormal = log1p(sqr_cv_incper);
        param2_incper[3] = sqrt(sqr_param2_lognormal);
        param1_incper[3] = log_mean_incper - sqr_param2_lognormal / 2.0;
    }


    vector[D] lps = log(weight); // internal component likelihoods
    {
        for (d in 1 : D) {
            vector[N] onset = to_vector(exp2onset_R) + onset_raw[d], exposure; 
            for (n in 1 : N) {
                real exposure_R = exp2onset_R[n] - exp2onset_L[n] + 1; 
                exposure_R = onset[n] < exposure_R ? onset[n] : exposure_R;
                exposure[n] = exposure_R * exposure_raw[d, n];
            }

            if (d == 1) {
                lps[d] += gamma_lpdf(onset - exposure | param1_incper[d], param2_incper[d]);
            } else if (d == 2) {
                lps[d] += weibull_lpdf(onset - exposure | param1_incper[d], param2_incper[d]);
            } else if (d == 3) {
                lps[d] += lognormal_lpdf(onset - exposure | param1_incper[d], param2_incper[d]);
            }
        }
    }
}

model {
    // priors
	mean_incper ~ normal(7, 7);
    param1_Weibull_incper ~ exponential(1.0);

    for (d in 1 : D) {
        onset_raw[d] ~ beta(1, 1); 
        exposure_raw[d] ~ beta(1, 1);
    }
  
    weight ~ dirichlet(rep_vector(1, D));

	/* likelihood */
    target += log_sum_exp(lps);
}

generated quantities {
    vector<lower = 0, upper = 1>[D] q = softmax(lps);
    int comp_incper = categorical_rng(q);

    real incper_pred  = (comp_incper == 3) ? lognormal_rng(param1_incper[comp_incper], param2_incper[comp_incper]) : (comp_incper == 2) ? weibull_rng(param1_incper[comp_incper], param2_incper[comp_incper]) : gamma_rng(param1_incper[comp_incper], param2_incper[comp_incper]);
}