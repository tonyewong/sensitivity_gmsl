##==============================================================================
## todo
##
## *  Sea level data and model values are in millimeters
##
## Note:  The file HadCRUT.4.4.0.0.annual_ns_avg_unc.txt was created by adding
##        the uncertainty estimate from one of the ensemble members, from:
## HadCRUT.4.4.0.0.annual_ns_avg_realisations/HadCRUT.4.4.0.0.annual_ns_avg.1.txt
##        (all members' uncertainty estimates are the same, if you check the files)
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


## Clear workspace
rm(list=ls())


##==============================================================================
## Read sea level (calibration) and temperature (forcing) data
##    * normalize relative to 1961-1990 mean for each
##    * this means the initial sea level parameter should also be relative to this

year_norm <- c(1961,1990)

## sea level data has years as half-years, so take floor and recognize that the
## values are for the year *beginning* with the value of 'year'
data_sealevel <- read.table('GMSL_ChurchWhite2011_yr_2015.txt')
colnames(data_sealevel) <- c('year','sealevel','uncertainty')
data_sealevel$year <- floor(data_sealevel$year)

## normalize sea level data
ind_norm <- which(data_sealevel$year==year_norm[1]):which(data_sealevel$year==year_norm[2])
data_sealevel$sealevel <- data_sealevel$sealevel - mean(data_sealevel$sealevel[ind_norm])

## temperature data
data_temperature <- read.csv('HadCRUT.4.4.0.0.annual_ns_avg_unc.csv')

## normalize temperature data
ind_norm <- which(data_temperature$year==year_norm[1]):which(data_temperature$year==year_norm[2])
data_temperature$temperature <- data_temperature$temperature - mean(data_temperature$temperature[ind_norm])
##==============================================================================


##==============================================================================
## Set up parameters to sample (not calibrating any statistical parameters...?)

parnames    <- c('sl_temp_sens', 'temp_equil', 'sl0')    # parameter names
parameters0 <- c(     3.4      ,     -0.3    ,  0   )    # initial parameter estimates
bound_lower <- c(      0       ,     -1.5    , -2*max(data_sealevel$uncertainty))
bound_upper <- c(     3.5      ,      1.5    ,  2*max(data_sealevel$uncertainty))
n_parameters <- length(parnames)
##==============================================================================


##==============================================================================
## Set up for optimization

source('gmsl_r07.R')
source('likelihood.R')

forcing <- data.frame(cbind(data_temperature$year, data_temperature$temperature))
names(forcing) <- c('year','temperature')

indices <- vector('list', 2); names(indices) <- c('normalize', 'mod2obs')
indices$normalize <- which(forcing$year==year_norm[1]):which(forcing$year==year_norm[2])
indices$mod2obs   <- which(forcing$year==min(data_sealevel$year)):which(forcing$year==max(data_sealevel$year))
##==============================================================================


##==============================================================================
## Preliminary DE optimization
NP_deoptim <- 50
niter_deoptim <- 50
F_deoptim <- 0.8
CR_deoptim <- 0.9

out.deoptim <- DEoptim(neg_log_likelihood, lower=bound_lower, upper=bound_upper,
                       DEoptim.control(NP=NP_deoptim,itermax=niter_deoptim,F=F_deoptim,CR=CR_deoptim,trace=FALSE),
                       parnames=parnames, data_sealevel=data_sealevel,
                       indices=indices, temperature_forcing=forcing)
##==============================================================================


#install.packages('sensitivity')
library(sensitivity)


##==============================================================================
## Method of Morris

mom <- morris(model=neg_log_likelihood, factors=parnames, r = 20, binf=bound_lower, bsup=bound_upper,
              scale=TRUE, design=list(type="oat", levels=5, grid.jump=3),
              parnames=parnames, data_sealevel=data_sealevel,
              indices=indices, temperature_forcing=forcing)
##==============================================================================


##==============================================================================
## Sobol' Method

#install.packages('lhs')
library(lhs)

n_sample <- 1000
n_bootstrap <- 100

## Need preliminary function to map ranges from [0,1] and back
map_range <- function(X, lbin, ubin, lbout, ubout){
    Y <- lbout + (ubout-lbout)*( (X-lbin)/(ubin-lbin) )
    return(Y)
}

## Sample parameters (need 2 data frames)
parameters_lhs1 <- randomLHS(n_sample, n_parameters)
parameters_lhs2 <- randomLHS(n_sample, n_parameters)

## Sobol' wrapper - assumes uniform distributions on parameters
gmsl_sobol <- function(parameters_lhs, parnames, data_sealevel, indices,
                       temperature_forcing, bound_lower, bound_upper) {
  parameters <- sapply(1:length(parnames), function(pp) {
                  map_range(parameters_lhs[,pp], 0, 1, bound_lower[pp], bound_upper[pp])})
  nll_output <- neg_log_likelihood(parameters, parnames, data_sealevel, indices, temperature_forcing)
  nll_output_centered <- nll_output - mean(nll_output)
  return(nll_output_centered)
}

## Actually run the Sobol'
t.out <- system.time(s.out <- sobolSalt(model=gmsl_sobol,
                             parameters_lhs1,
                             parameters_lhs2,
                             scheme='B',
                             nboot=n_bootstrap))


##==============================================================================


##==============================================================================
## Latin Hypercube

todo
##==============================================================================


##==============================================================================
## End
##==============================================================================
