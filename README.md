# Mapping_BMI
Source codes and data for the paper "Mapping between measurement scales in meta-analysis, with application to measures of body mass index in children" Davies, Ades and Higgins (2024).

## Charts
A folder containing the LMS reference charts per sex. In each chart, age is tabulated against the corresponding LMS values.

## Data
A folder containing the raw obesity prevention trial data for each age group.

## PlottingData
A folder containing the results of the application to the obesity data. Files in this folder are read into the Python script PlotApplication.ipynb to create the plots in the paper (main and supplement).

## SimulationData
A folder containing the results of the simulation study. Files in this folder are read into the Python script PlotSimulationResults.ipynb to create the plots in the paper (main and supplement).

Files FromB_IOTF_SDmid.csv and FromZ_IOTF_SDmid.csv are created from the code 3_SimulationStudy.R. The other files are created by running the same code but changing the standard deviations (see paper). 

## functions
A folder containing the functions used to perform the mapping methods on the obesity data.

## simulaton_funcs
A folder containing the functions used to perform the simulation studies.

## 1_GenerateMappedData.R
An R file to map the data (read in from folder 'Data') from percentile and BMI onto zBMI using each of the mapping methods. Outputs the RData file MappedData.RData.

## 2_CreatePlottingVectors.R
An R file to read in MappedData.RData (the results from 1_GenerateMappedData.R) and create vectors of results for plotting. Results of this code are saved in the folder PlottingData and in the RData file PlottingVecs.RData.

## 3_SimulationStudy.R
An example of an R code to perform a simulation study of the mapping methods (with SD=1 for zBMI and SD=3.5 for BMI). Results are saved in the folder SimulationData (FromB_IOTF_SDmid.csv and FromZ_IOTF_SDmid.csv).

## MappedData.RData
Data saved after running 1_GenerateMappedData.R.

## PlotApplication.ipynb
A Python script to plot the results of the obesity data analysis. It reads in data from the folder PlotttingData.

## PlotSimulationResults.ipynb
A Python script to plot the results of the simulations. It reads in data from the folder SimulationData.

## PlottingVecs.RData
Data saved after running 2_CreatePlottingVectors.R. Individual files (dataframes) saved in the folder PlottingData.
