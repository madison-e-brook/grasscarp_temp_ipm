# grasscarp_temp_ipm

# File name explanations

dd refers to degree days

appendix2-pm-parameter-dd-figures.R-code for figures of the pm(a) relationship for appendix 2dd-for-comparison-locations.R-get DD estimates for three locations to compare model with; Upper Mississippi, Lake Erie, Lake Superiordd-lambda-figure.R-code that makes the figure for the relationship between dd and lambda valuesdd-regression-figure-appendix2.R-code to get regression plot for appendix 2default-model-code.R-code for the IPM that is used in many other scripts

diagnostic-plot-bin-number-max-size-appendix1.R-code for finding appropriate value for upper limit of IPM and bin number, and associated plotsdownload-cpc-temp.R-downloads cpc data from internetelasticity-figure.R-code that makes the elasticity analysis figureelasticity-XX-survival-XXXX-degree-days.R-codes that calculate elasticity values and saves csvs at provided survival value (XX) and at given degree days (XXXX)

fecundity-data-appendix2.R-fecundity figures for appendix 2

get-DD-fecundity-locations.R-script used to get the dd for the locations with relative fecunditiesgrowth-figures-appendix2.R-figures for von bert, growth increment and growth distributions for appendix 2

lambdas-across-survival-range.R-code calculates the population growth rates across a range of maximum survival values for 3 different degree days, and creates the csv for the survival ranges.latin-hypercube-sampling.R-code that gets the predicted lambda values at the two maximum survival values, and calculates the confidence intervals. 

main-function-graphs.R-code to make the figures for pm(a), s(a), ginc(z), b(z)make-raster-lambda-map.R-code that makes the raster files (0.6 survival North america lambda map.grd and 0.9 survival North america lambda map.grd) used in making-lambda-map-plot.Rmaking-lambda-map-plot.R-makes the map of north America values of lambda for both survival maximums.max-age-diagnostics.R-code to determine the appropriate maximum age of the population

pm-age-sd-from-regression.R-script that generates temperature dependent age at maturity relationship.pm-parameter-dd-relationship.R-uses estimated ages at maturity and the variance to generate the relationship between the pm(a) parameters and DD

survival-figure-multiple-dd.R-code that makes the figure of maximum survival affect at multiple degree daystemp-dependent-growth-and-survival.R-code for adding temperature dependent growth rate k, Linf and maximum survival to model, and associated figures