functions {
    //// ODE for the integral
    real mean_gamma_1st_integral_function(real xi, real y, array[] real omega_r, array[] real x_r, array[] int x_i) {
        real ts = omega_r[1], param1 = omega_r[2], param2 = omega_r[3];
        return xi * exp(gamma_lpdf(ts * xi | param1, param2));
    }

    real mean_weibull_1st_integral_function(real xi, real y, array[] real omega_r, array[] real x_r, array[] int x_i) {
        real ts = omega_r[1], param1 = omega_r[2], param2 = omega_r[3];
        return xi * exp(weibull_lpdf(ts * xi | param1, param2));
    }

    real mean_lognormal_1st_integral_function(real xi, real y, array[] real omega_r, array[] real x_r, array[] int x_i) {
        real ts = omega_r[1], param1 = omega_r[2], param2 = omega_r[3];
        return xi * exp(lognormal_lpdf(ts * xi | param1, param2));
    }

    real mean_1st_integral(int d, real ts, real param1, real param2) {
		real value;
        if (d == 1)
			value = integrate_1d(mean_gamma_1st_integral_function, 1e-6, 1 - 1e-6, {ts, param1, param2}, rep_array(0.0, 0), rep_array(0, 0));
		else if (d == 2)
			value = integrate_1d(mean_weibull_1st_integral_function, 1e-6, 1 - 1e-6, {ts, param1, param2}, rep_array(0.0, 0), rep_array(0, 0));
		else if (d == 3)
			value = integrate_1d(mean_lognormal_1st_integral_function, 1e-6, 1 - 1e-6, {ts, param1, param2}, rep_array(0.0, 0), rep_array(0, 0));
        return square(ts) * value;
	}
}

data {
	/* linelist data for unique IDs */
	// all the data are in the same time units in days
    int<lower = 1> N_IDs; // number of unique IDs (onsets are all the same for those IDs)
    /* array of symptom onsets and exposures corresponding to unique IDs */
    array[N_IDs] real onsetR;
    array[N_IDs] real<upper = onsetR> onsetL;

    real t_int; // time of the intervention

    /* array of IDs corresponding to infectors and infectees in the transmission pairs */
    int<lower = 1> N; // number of case records 
    array[N] int<lower = 1, upper = N_IDs> ID_infector, ID_infectee;

    /* incubation period parameters */
    real<lower = 0> mean_incper, sd_incper;

    /* number of days for posterior predictive draws */
    int<lower = 1> T;
    int<lower = 0> t_preceeding; // number of days before the first onset for the posterior predictive draws
}

transformed data {
    int D = 3; // number of distributions for incper and TOST (Gamma, Weibull, Lognormal)
    real exposure_scale = 7.0;

    // here we identify the list of index cases = i.e. cases who serve solely as infectors
	// and do not have any infectee ID in the transmission pairs
	array[N_IDs] int index_case = rep_array(1, N_IDs);
	for (ID in 1 : N_IDs)
		for (n in 1 : N) 
			if (ID == ID_infectee[n]) {
				index_case[ID] = 0; // this ID is an infectee
				break; // no need to check further
			}
	int N_index_cases = sum(index_case);
	array[N_index_cases] int index_case_IDs;
	{
		int idx = 0;
		for (ID in 1 : N_IDs) 
			if (index_case[ID]) {
				idx += 1;
				index_case_IDs[idx] = ID;
			}
	}
	// we identify the maximum number of secondary cases per primary case 
	// this will be used for determining the array size in the tuples
	int N_max_secondary_cases = 0;
	for (ID in 1 : N_IDs) {
		int N_ = 0;
		for (n in 1 : N) {
			if (ID_infector[n] == ID) N_ += 1;
		}
		if (N_ > N_max_secondary_cases) 
			N_max_secondary_cases = N_;
	}
	// we create an array of tuples for the epidemiological links
	// the first element of the tuple is the number of secondary cases,
	// the second element is an array of IDs of secondary cases
	array[N_IDs] tuple(int, array[N_max_secondary_cases] int) epilink;
	for (ID in 1 : N_IDs) {
		epilink[ID].2 = rep_array(0, N_max_secondary_cases);
		int idx = 0;
		for (n in 1 : N) 
			if (ID_infector[n] == ID) {
				idx += 1;  
				epilink[ID].2[idx] = ID_infectee[n];
			}
		epilink[ID].1 = idx;
	}
	
	// Validate network structure
	{
		// Verify each infectee appears only once (tree property)
		array[N_IDs] int infectee_count = rep_array(0, N_IDs);
		for (n in 1 : N) {
			infectee_count[ID_infectee[n]] += 1;
			if (infectee_count[ID_infectee[n]] > 1) {
				reject("Network validation failed: ID ", ID_infectee[n], " appears as infectee multiple times. Expected tree structure.");
			}
		}
		// Verify tree structure consistency
		// Count unique infectees
		array[N_IDs] int is_infectee = rep_array(0, N_IDs);
		for (n in 1 : N) {
			is_infectee[ID_infectee[n]] = 1;
		}
		int N_unique_infectees = sum(is_infectee);
		if (N_index_cases + N_unique_infectees != N_IDs) {
			reject("Network validation failed: Index cases (", N_index_cases, ") + unique infectees (", N_unique_infectees, ") = ", N_index_cases + N_unique_infectees, " but expected ", N_IDs, " total IDs.");
		}
	}

    // parameters of the lognormal distribution for the incubation period
    real sqr_param2_incper = log1p(square(sd_incper / mean_incper)),
        param1_incper = log(mean_incper) - sqr_param2_incper / 2.0,
        param2_incper = sqrt(sqr_param2_incper);
}

parameters {
	vector<lower = 0, upper = 1>[N_IDs] onset_raw; 
    vector[N_index_cases] log_exposure_raw_index_cases; // raw exposures for index cases
    vector<lower = 0, upper = 1>[N_IDs - N_index_cases] exposure_raw_secondary_cases; // raw exposures for secondary cases

	real log_mean_gentime_raw;
    real log_param1_Weibull_gentime; 

    real log_R0;
    real log_phi;
    real logit_epsilon;

    simplex[D] weight; // mixing proportions
}

transformed parameters {
    real epsilon = inv_logit(logit_epsilon);
    vector[N_index_cases] exposure_raw_index_cases = exp(log_exposure_raw_index_cases);
    real phi = exp(log_phi);

    real log_mean_gentime = 2.5 + log_mean_gentime_raw * 5.0;
    real mean_gentime = exp(log_mean_gentime); // mean of the time intervals
    vector[D] param1_gentime, param2_gentime;
    
    // Weibull distribution
    param1_gentime[2] = exp(log_param1_Weibull_gentime);
    real sqr_cv_gentime = tgamma(1 + 2 / param1_gentime[2]) / square(tgamma(1 + 1 / param1_gentime[2])) - 1;
    param2_gentime[2] = mean_gentime / tgamma(1.0 + 1.0 / param1_gentime[2]);

    // Gamma distribution
    param1_gentime[1] = 1 / sqr_cv_gentime;
    param2_gentime[1] = param1_gentime[1] / mean_gentime;

    // Lognormal distribution
    {
        real sqr_param2_lognormal = log1p(sqr_cv_gentime);
        param2_gentime[3] = sqrt(sqr_param2_lognormal);
        param1_gentime[3] = log_mean_gentime - sqr_param2_lognormal / 2.0;
    }


    vector[D] lps = log(weight); // internal component likelihoods
    real lp_incper;
    {
        vector[N_IDs] exposure, onset; 
        // We traverse the network starting from index cases using a queue-like approach
        // First, process all index cases and add them to a processing queue
        array[N_IDs] int processing_queue = rep_array(0, N_IDs);
        int queue_size = 0;
        int index_censored_exposureL_infectors = 0;
        for (n in 1 : N_index_cases) {
            int index_ID = index_case_IDs[n];
            onset[index_ID] = fma(onsetR[index_ID] - onsetL[index_ID], onset_raw[index_ID], onsetL[index_ID]);
            exposure[index_ID] = onset[index_ID] - exposure_scale * exposure_raw_index_cases[n];
            // Add this index case to the processing queue
            queue_size += 1;
            processing_queue[queue_size] = index_ID;
        }
        // Process the queue - each processed node adds its secondary cases to the queue
        int queue_position = 1;
        int idx_secondary = 1;
        while (queue_position <= queue_size) {
            int primary_ID = processing_queue[queue_position];
            // Process all secondary cases for the current ID
            int N_secondary_cases = epilink[primary_ID].1;
            for (i in 1 : N_secondary_cases) {
                int secondary_ID = epilink[primary_ID].2[i];
                // Process this secondary case
                real onsetL_ = onsetL[secondary_ID] < exposure[primary_ID] ? exposure[primary_ID] : onsetL[secondary_ID];
                onset[secondary_ID] = fma(onsetR[secondary_ID] - onsetL_, onset_raw[secondary_ID], onsetL_);
                exposure[secondary_ID] = fma(onset[secondary_ID] - exposure[primary_ID], exposure_raw_secondary_cases[idx_secondary], exposure[primary_ID]);
                // Add this secondary case to the queue for further processing
                queue_size += 1;
                processing_queue[queue_size] = secondary_ID;
                idx_secondary += 1;
            }
            queue_position += 1;
        }

        // incubation periods
        lp_incper = lognormal_lpdf(onset - exposure + 1e-6 | param1_incper, param2_incper);

        real log1m_epsilon = log1m(epsilon);

        for (d in 1 : D) {
            real log_R; // case reproduction numbers
            for (ID in 1 : N_IDs) {
                if (exposure[ID] > t_int) {
                    // this ID is a post-intervention case
                    log_R = log_R0 + log1m(epsilon);
                } else {
                    // this ID is a pre-intervention case
                    if (d == 1) 
                        log_R = log_R0 + log_sum_exp(gamma_lcdf(t_int - exposure[ID] | param1_gentime[d], param2_gentime[d]), log1m_epsilon + gamma_lccdf(t_int - exposure[ID] | param1_gentime[d], param2_gentime[d]));
                    else if (d == 2)
                        log_R = log_R0 + log_sum_exp(weibull_lcdf(t_int - exposure[ID] | param1_gentime[d], param2_gentime[d]), log1m_epsilon + weibull_lccdf(t_int - exposure[ID] | param1_gentime[d], param2_gentime[d]));
                    else if (d == 3)
                        log_R = log_R0 + log_sum_exp(lognormal_lcdf(t_int - exposure[ID] | param1_gentime[d], param2_gentime[d]), log1m_epsilon + lognormal_lccdf(t_int - exposure[ID] | param1_gentime[d], param2_gentime[d]));
                }
                lps[d] += neg_binomial_2_log_lpmf(epilink[ID].1 | log_R, phi);
            }
            lps[d] += neg_binomial_2_log_lpmf(0 | log_R0 + log1m_epsilon, phi); // case No. 42 has no secondary cases, so we add a zero to the likelihood

            if (d == 1) {
                lps[d] += gamma_lpdf(exposure[ID_infectee] - exposure[ID_infector] | param1_gentime[d], param2_gentime[d]);
                for (n in 1 : N)
                    if (t_int > exposure[ID_infector[n]]) {
                        if (exposure[ID_infectee[n]] > t_int)
                            lps[d] += log1m_epsilon;
                        lps[d] -= log_sum_exp(log1m_epsilon + gamma_lccdf(t_int - exposure[ID_infector[n]] | param1_gentime[d], param2_gentime[d]), gamma_lcdf(t_int - exposure[ID_infector[n]] | param1_gentime[d], param2_gentime[d]));
                    }
            } else if (d == 2) {
                lps[d] += weibull_lpdf(exposure[ID_infectee] - exposure[ID_infector] | param1_gentime[d], param2_gentime[d]);
                for (n in 1 : N)
                    if (t_int > exposure[ID_infector[n]]) {
                        if (exposure[ID_infectee[n]] > t_int)
                            lps[d] += log1m_epsilon;
                        lps[d] -= log_sum_exp(log1m_epsilon + weibull_lccdf(t_int - exposure[ID_infector[n]] | param1_gentime[d], param2_gentime[d]), weibull_lcdf(t_int - exposure[ID_infector[n]] | param1_gentime[d], param2_gentime[d]));
                    } 
            } else if (d == 3) {
                lps[d] += lognormal_lpdf(exposure[ID_infectee] - exposure[ID_infector] | param1_gentime[d], param2_gentime[d]);
                for (n in 1 : N)
                    if (t_int > exposure[ID_infector[n]]) {
                        if (exposure[ID_infectee[n]] > t_int)
                            lps[d] += log1m_epsilon;
                        lps[d] -= log_sum_exp(log1m_epsilon + lognormal_lccdf(t_int - exposure[ID_infector[n]] | param1_gentime[d], param2_gentime[d]), lognormal_lcdf(t_int - exposure[ID_infector[n]] | param1_gentime[d], param2_gentime[d]));
                    }
            }
        }
    }
}

model {
    // priors
    log_mean_gentime_raw ~ std_normal();
    log_param1_Weibull_gentime ~ std_normal();

    onset_raw ~ beta(1, 1);
    log_exposure_raw_index_cases ~ std_normal();
    exposure_raw_secondary_cases ~ beta(1, 1);

    log_R0 ~ std_normal();
    log_phi ~ std_normal();
    logit_epsilon ~ std_normal();

    weight ~ dirichlet(rep_vector(1, D));

	/* likelihood */
    target += lp_incper + log_sum_exp(lps);
}

generated quantities {
    real R0 = exp(log_R0);

    real cv_gentime = sqrt(sqr_cv_gentime), 
        sd_gentime = mean_gentime * cv_gentime;

    vector<lower = 0, upper = 1>[D] q = softmax(lps);
    int comp_gentime = categorical_rng(q);

    real gentime_pred = (comp_gentime == 3) ? lognormal_rng(param1_gentime[comp_gentime], param2_gentime[comp_gentime]) : (comp_gentime == 2) ? weibull_rng(param1_gentime[comp_gentime], param2_gentime[comp_gentime]) : gamma_rng(param1_gentime[comp_gentime], param2_gentime[comp_gentime]);
    
    // two other epi intervals
    real TOST_pred, serint_pred;
    {
        vector[2] incper_pred;
        for (j in 1 : 2)
            incper_pred[j] = lognormal_rng(param1_incper, param2_incper);
        TOST_pred = gentime_pred - incper_pred[1];
        serint_pred = gentime_pred + incper_pred[2] - incper_pred[1];
    }

    vector[T + t_preceeding] R_t, gentime_mean_t;
    for (t in 1 : T + t_preceeding) {
        real t_ = t + 0.5 - t_preceeding;
        if (comp_gentime == 1) {
            R_t[t] = (t_ < t_int) ? 
                R0 * (gamma_cdf(t_int - t_ | param1_gentime[comp_gentime], param2_gentime[comp_gentime]) + (1 - epsilon) * (1.0 - gamma_cdf(t_int - t_ | param1_gentime[comp_gentime], param2_gentime[comp_gentime]))) : 
                R0 * (1 - epsilon);
            gentime_mean_t[t] = (t_ < t_int) ? 
                ((1 - epsilon) * mean_gentime + epsilon * mean_1st_integral(comp_gentime, t_int - t_, param1_gentime[comp_gentime], param2_gentime[comp_gentime])) / ((1 - epsilon) * (1 - gamma_cdf(t_int - t_ | param1_gentime[comp_gentime], param2_gentime[comp_gentime])) + gamma_cdf(t_int - t_ | param1_gentime[comp_gentime], param2_gentime[comp_gentime])) :
                mean_gentime;
        } else if (comp_gentime == 2) {
            R_t[t] = (t_ < t_int) ? 
                R0 * (weibull_cdf(t_int - t_ | param1_gentime[comp_gentime], param2_gentime[comp_gentime]) + (1 - epsilon) * (1.0 - weibull_cdf(t_int - t_ | param1_gentime[comp_gentime], param2_gentime[comp_gentime]))) : 
                R0 * (1 - epsilon);
            gentime_mean_t[t] = (t_ < t_int) ?
                ((1 - epsilon) * mean_gentime + epsilon * mean_1st_integral(comp_gentime, t_int - t_, param1_gentime[comp_gentime], param2_gentime[comp_gentime])) / ((1 - epsilon) * (1 - weibull_cdf(t_int - t_ | param1_gentime[comp_gentime], param2_gentime[comp_gentime])) + weibull_cdf(t_int - t_ | param1_gentime[comp_gentime], param2_gentime[comp_gentime])) :
                mean_gentime;
        } else {
            R_t[t] = (t_ < t_int) ? 
                R0 * (lognormal_cdf(t_int - t_ | param1_gentime[comp_gentime], param2_gentime[comp_gentime]) + (1 - epsilon) * (1.0 - lognormal_cdf(t_int - t_ | param1_gentime[comp_gentime], param2_gentime[comp_gentime]))) : 
                R0 * (1 - epsilon);
            gentime_mean_t[t] = (t_ < t_int) ? 
                ((1 - epsilon) * mean_gentime + epsilon * mean_1st_integral(comp_gentime, t_int - t_, param1_gentime[comp_gentime], param2_gentime[comp_gentime])) / ((1 - epsilon) * (1 - lognormal_cdf(t_int - t_ | param1_gentime[comp_gentime], param2_gentime[comp_gentime])) + lognormal_cdf(t_int - t_ | param1_gentime[comp_gentime], param2_gentime[comp_gentime])) :
                mean_gentime;
        }
    }
}