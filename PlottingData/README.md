# Mapping results to be plotted

Each file is a dataframe containing vectors of arm-level aggregate data alongside other details including the trial name, intervention name, LMS reference chart used etc.

## Percentile results
- z_for_perc: reported zBMI values corresponding to the mapped percentile data (i.e. the target values for mapped percentile data)
- perc_analytic: estimated zBMI values mapped from percentile using the analytic method
- perc_sampling: estimated zBMI values mapped from percentile using the sampling method 
- perc_opt_dist: estimated zBMI values mapped from percentile using the optimization method (distribution estimates)
- perc_opt_samp: estimated zBMI values mapped from percentile using the optimization method (sample estimates)

## BMI results
- z_for_bmi: reported zBMI values corresponding to the mapped BMI data (i.e. the target values for mapped BMI data)
- bmi_sampling_aNorm: estimated zBMI values mapped from BMI using the sampling method & a normal distribution for age
- bmi_sampling_aUnif: estimated zBMI values mapped from BMI using the sampling method & a uniform distribution for age
- bmi_opt_dist: estimated zBMI values mapped from BMI using the optimization method (distribution estimates)
- bmi_opt_samp: estimated zBMI values mapped from BMI using the optimization method (sample estimates)
