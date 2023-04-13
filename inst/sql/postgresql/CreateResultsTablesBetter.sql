
-- Drop old tables if exists

DROP TABLE IF EXISTS all_ipcs;
DROP TABLE IF EXISTS analysis;
DROP TABLE IF EXISTS database;
DROP TABLE IF EXISTS exposure;
DROP TABLE IF EXISTS mses;
DROP TABLE IF EXISTS powers;
DROP TABLE IF EXISTS priors;
DROP TABLE IF EXISTS summary;
DROP TABLE IF EXISTS time_to_signal;
DROP TABLE IF EXISTS type1s;

-- Create tables

CREATE TABLE all_ipcs (
     outcome_id INTEGER NOT NULL,
     exposure_id INTEGER NOT NULL,
     negative_control_id INTEGER NOT NULL,
     effect_size NUMERIC NOT NULL,
     PRIMARY KEY(outcome_id, exposure_id, negative_control_id, effect_size)
);

CREATE TABLE analysis (
     method VARCHAR(255) NOT NULL,
     analysis_id INTEGER NOT NULL,
     description VARCHAR(255) ,
     time_at_risk VARCHAR(255) ,
     PRIMARY KEY(method, analysis_id)
);

CREATE TABLE database (
     database_id VARCHAR(255) NOT NULL,
     database_name VARCHAR(255) ,
     description TEXT ,
     vocabulary_version VARCHAR(255) ,
     min_obs_period_date DATE ,
     max_obs_period_date DATE ,
     study_package_version VARCHAR(255) ,
     is_meta_analysis INTEGER ,
     PRIMARY KEY(database_id)
);

CREATE TABLE exposure (
     exposure_id INTEGER NOT NULL,
     exposure_name VARCHAR(255) ,
     total_shots INTEGER ,
     base_exposure_id INTEGER NOT NULL,
     base_exposure_name VARCHAR(255) ,
     shot VARCHAR(255) ,
     start_date DATE ,
     end_date DATE ,
     history_start_date DATE ,
     history_end_date DATE ,
     PRIMARY KEY(exposure_id, base_exposure_id)
);

CREATE TABLE mses (
     database_id VARCHAR(255) NOT NULL,
     method VARCHAR(255) NOT NULL,
     analysis_id INTEGER NOT NULL,
     exposure_id INTEGER NOT NULL,
     effect_size NUMERIC NOT NULL,
     maxSPRT_mse NUMERIC ,
     calibrated_mse NUMERIC ,
     maxSPRT_coverage NUMERIC ,
     calibrated_coverage NUMERIC ,
     Bayesian_mse NUMERIC ,
     adjusted_mse NUMERIC ,
     Bayesian_coverage NUMERIC ,
     adjusted_covearge NUMERIC ,
     PRIMARY KEY(database_id, method, analysis_id, exposure_id, effect_size)
);

CREATE TABLE powers (
     database_id VARCHAR(255) NOT NULL,
     method VARCHAR(255) NOT NULL,
     exposure_id INTEGER NOT NULL,
     analysis_id INTEGER NOT NULL,
     period_id INTEGER NOT NULL,
     y NUMERIC ,
     effect_size NUMERIC NOT NULL,
     stats VARCHAR(255) ,
     approach VARCHAR(255) NOT NULL,
     stage NUMERIC ,
     power NUMERIC ,
     true_rr VARCHAR(255) ,
     PRIMARY KEY(database_id, method, exposure_id, analysis_id, period_id, effect_size, approach)
);

CREATE TABLE priors (
     Mean INTEGER ,
     Sd INTEGER ,
     prior_id INTEGER NOT NULL,
     PRIMARY KEY(prior_id)
);

CREATE TABLE summary (
     database_id VARCHAR(255) NOT NULL,
     method VARCHAR(255) NOT NULL,
     analysis_id INTEGER NOT NULL,
     exposure_id INTEGER NOT NULL,
     outcome_id INTEGER NOT NULL,
     period_id INTEGER NOT NULL,
     prior_id INTEGER NOT NULL,
     post_mean NUMERIC ,
     post_map NUMERIC ,
     post_median NUMERIC ,
     adjusted_post_mean NUMERIC ,
     adjusted_post_map NUMERIC ,
     adjusted_post_median NUMERIC ,
     ci_95_lb NUMERIC ,
     ci_95_ub NUMERIC ,
     p_1 NUMERIC ,
     p_0 NUMERIC ,
     adjusted_ci_95_lb NUMERIC ,
     adjusted_ci_95_ub NUMERIC ,
     adjusted_p_1 NUMERIC ,
     adjusted_p_0 NUMERIC ,
     PRIMARY KEY(database_id, method, analysis_id, exposure_id, outcome_id, period_id, prior_id)
);

CREATE TABLE time_to_signal (
     database_id VARCHAR(255) NOT NULL,
     exposure_id INTEGER NOT NULL,
     method VARCHAR(255) NOT NULL,
     analysis_id INTEGER NOT NULL,
     effect_size NUMERIC NOT NULL,
     time_to_sens INTEGER ,
     approach VARCHAR(255) NOT NULL,
     time_at_risk VARCHAR(255) ,
     sensitivity NUMERIC NOT NULL,
     PRIMARY KEY(database_id, exposure_id, method, analysis_id, effect_size, approach, sensitivity)
);

CREATE TABLE type1s (
     database_id VARCHAR(255) NOT NULL,
     method VARCHAR(255) NOT NULL,
     exposure_id INTEGER NOT NULL,
     analysis_id INTEGER NOT NULL,
     period_id INTEGER NOT NULL,
     y NUMERIC ,
     effect_size NUMERIC NOT NULL,
     stats VARCHAR(255) ,
     approach VARCHAR(255) NOT NULL,
     stage NUMERIC ,
     analysis VARCHAR(255) ,
     time_at_risk VARCHAR(255) ,
     PRIMARY KEY(database_id, method, exposure_id, analysis_id, period_id, effect_size, approach)
);

