# Mapping_BMI
Source codes and data for the paper "Mapping between measurement scales in meta-analysis, with application to measures of body mass index in children" Davies, Ades and Higgins (2024).

## Charts
A folder containing the LMS reference charts per sex. In each chart, age is tabulated against the corresponding LMS values.

## Data
A folder containing the raw obesity prevention trial data for each age group.

## functions
A folder containing the functions used to perform the mapping methods on the obesity data.

## PlottingData
A folder containing the results of the application to the obesity data.

## simulaton_funcs
A folder containing the functions used to perform the simulation studies.

## 1_GenerateMappedData.R
An R file to map the data (read in from folder 'Data') from percentile and BMI onto zBMI using each of the mapping methods. 

## 2_CreatePlottingVectors.R
An R file to read in the results from 1_GenerateMappedData.R and create vectors of results for plotting. Results of this code are saved in the folder PlottingData.

## 3_SimulationStudy.R
An example of an R code to perform a simulation study of the mapping methods (with SD=1 for zBMI and SD=3.5 for BMI).
