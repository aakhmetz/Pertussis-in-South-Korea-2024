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
}

transformed data {
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
    //// to note the first element provides the case reproduction number
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

    // parameters of the lognormal distribution for incper
    real sqr_param2_incper = log1p(square(sd_incper / mean_incper)),
        param1_incper = log(mean_incper) - sqr_param2_incper / 2.0,
        param2_incper = sqrt(sqr_param2_incper);
}

parameters {
	vector<lower = 0, upper = 1>[N_IDs] onset_raw; 
    vector[N_index_cases] log_exposure_raw_index_cases; // raw exposures for index cases
    vector<lower = 0, upper = 1>[N_IDs - N_index_cases] exposure_raw_secondary_cases; // raw exposures for secondary cases

	real log_mean_gentime_raw;;
    real<lower = 0> sqr_cv_gentime; 

    real log_R0;
    real log_phi; 
    real logit_epsilon;
}

transformed parameters {
    real epsilon = inv_logit(logit_epsilon);
    vector[N_index_cases] exposure_raw_index_cases = exp(log_exposure_raw_index_cases);
    real phi = exp(log_phi);

    real log_mean_gentime = 2.5 + log_mean_gentime_raw * 5.0;
    real mean_gentime = exp(log_mean_gentime),
        cv_gentime = sqrt(sqr_cv_gentime),
        sd_gentime = mean_gentime .* cv_gentime,
        log_sd_gentime = log(sd_gentime);

    // Gamma distribution
    real param1_gentime = 1 / sqr_cv_gentime,
        param2_gentime = param1_gentime / mean_gentime;

    vector[N_IDs] exposure, onset; 
    {
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
    }
}

model {
    // priors
    log_mean_gentime_raw ~ std_normal();
	sqr_cv_gentime ~ gamma(2, .5);

	onset_raw ~ beta(1, 1); 
	log_exposure_raw_index_cases ~ std_normal();
	exposure_raw_secondary_cases ~ beta(1, 1);

	// negative_censored_exposureL_index_cases ~ normal(20, 10);

    log_R0  ~ std_normal();
    log_phi ~ std_normal();
    logit_epsilon ~ std_normal();

    // incubation period
    target += lognormal_lupdf(onset - exposure | param1_incper, param2_incper);

    // offspring distribution
    vector[N_IDs] log_Rt; // case reproduction numbers
    for (ID in 1 : N_IDs) {
        if (exposure[ID] > t_int) {
            // this ID is a post-intervention case
            log_Rt[ID] = log_R0 + log1m(epsilon);
        } else {
            // this ID is a pre-intervention case
            log_Rt[ID] = log_R0 + log_sum_exp(gamma_lcdf(t_int - exposure[ID] | param1_gentime, param2_gentime), log1m(epsilon) + gamma_lccdf(t_int - exposure[ID] | param1_gentime, param2_gentime));
        }
        target += neg_binomial_2_log_lupmf(epilink[ID].1 | log_Rt[ID], phi);
    }
    target += neg_binomial_2_log_lupmf(0 | log_R0 + log1m(epsilon), phi); // case No. 42 has no secondary cases, so we add a zero to the likelihood

    target += gamma_lpdf(exposure[ID_infectee] - exposure[ID_infector] | param1_gentime, param2_gentime);
    for (n in 1 : N) {
        if (t_int > exposure[ID_infector[n]]) {
            if (exposure[ID_infectee[n]] > t_int)
                target += log1m(epsilon);
            target += -log_sum_exp(log1m(epsilon) + gamma_lccdf(t_int - exposure[ID_infector[n]] | param1_gentime, param2_gentime), gamma_lcdf(t_int - exposure[ID_infector[n]] | param1_gentime, param2_gentime));
        }
    }
}

generated quantities {
    real R0 = exp(log_R0);

    vector[2] incper_pred;
    for (j in 1 : 2)
        incper_pred[j] = lognormal_rng(param1_incper, param2_incper);
    real gentime_pred = gamma_rng(param1_gentime, param2_gentime);
    // two other epi intervals
    real TOST_pred = gentime_pred - incper_pred[1];
    real serint_pred = gentime_pred + incper_pred[2] - incper_pred[1];
}