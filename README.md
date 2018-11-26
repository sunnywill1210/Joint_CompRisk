# Joint_CompRisk
R coding for examples 

Example 1: Hodgkin disease
hodgkin.txt: data file
CSHandCIFcomparison: produce figures for CSH and CIF for different age groups.
hodgkin_twosample_CSHACH: Conduct two sample joint tests for CSH and ACH, and produce confidence region
hodgkin_twosample_CSHCIF: Conduct Two sample joint tests for CSH and CIF
hodgkin_twosample_CSHOCH: Conduct Two sample joint tests for CSH and OCH, and produce confidence region

Example 2: BMT
bmt_data.csv: data file
bmt_dataCleaning: cleaning data and create necessary variables
BMT_twosample_CSHACH: Conduct two sample joint tests for CSH and ACH, and produce confidence region
BMT_twosample_CSHCIF: Conduct Two sample joint tests for CSH and CIF
BMT_twosample_CSHOCH: Conduct Two sample joint tests for CSH and OCH, and produce confidence region
BMT_regression_CSHACH: Conduct joint tests for CSH and ACH based on regression, and produce confidence region

Example 3: Follicular
follicular: data file
follicular_regression_CSHACH: Conduct joint tests for CSH and ACH based on regression, and produce confidence region

correlation_CSHCIF: program to calculate correlation between CSH and CIF test statistics. Called by "hodgkin_twosample_CSHCIF" and "BMT_twosample_CSHCIF"
