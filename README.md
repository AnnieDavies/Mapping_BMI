# Mapping_BMI
Source codes and data for the paper "Mapping between aggregate level BMI data on different measurement scales in a meta-analysis of childhood obesity prevention interventions" Davies, Ades and Higgins (2024).

## Data
A folder containing the raw data for each age group.

## Charts
A folder containing the LMS reference charts per sex. In each chart, age is tabulated against the corresponding LMS values.

## functions
A folder containing the functions used to perform the mapping methods and format the data.

## PlottingData
A folder containing the results of the mapping methods plotted in the paper.

## 1_GenerateMappedData.R
An R file to map the data (read in from folder 'Data') from percentile and BMI onto zBMI using each of the mapping methods. 

## 2_CreatePlottingVectors.R
An R file to read in the results from 1_GenerateMappedData.R and create vectors of results for plotting. Results of this code are saved in the folder PlottingData.

