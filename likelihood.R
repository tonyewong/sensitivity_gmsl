##==============================================================================
## likelihood.R
##
## Compute (log) likelihood of the observational data, given the values of the
## model parameters from which the model simulation was generated.
##
## Assumption is that the observations are independent and identically
## distributed (IID). This is almost certainly not true.  Really ought to be
## using an AR1 (for example) residual model to 'whiten' the residuals.
##
## Input:
##
## Output:
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


##==============================================================================
## Log-likelihood
log_likelihood = function(parameters,
                          parnames,
                          data_sealevel,
                          indices,
                          temperature_forcing)
{
  n_data <- length(indices$mod2obs)
  n_simulations <- round( length(parameters)/length(parnames) )
  if(n_simulations==1) {
    sl_temp_sens <- parameters[match("sl_temp_sens",parnames)]
    temp_equil   <- parameters[match("temp_equil",parnames)]
    sl0          <- parameters[match("sl0",parnames)]
  } else {
    sl_temp_sens <- parameters[,match("sl_temp_sens",parnames)]
    temp_equil   <- parameters[,match("temp_equil",parnames)]
    sl0          <- parameters[,match("sl0",parnames)]
  }

  gmsl_model <- sapply(1:n_simulations, function(ss) {
                  gmsl_r07(sl_temp_sens=sl_temp_sens[ss], temp_equil=temp_equil[ss],
                           sl0=sl0[ss], temperature_forcing=temperature_forcing$temperature)})

  gmsl_model_norm <- sapply(1:n_simulations, function(ss) {
                       gmsl_model[,ss] - mean(gmsl_model[indices$normalize,ss])})

  # Initialize
  loglike <- rep(0, n_simulations)

  # Calculate the residuals
  residuals <- sapply(1:n_simulations, function(ss) {
                 data_sealevel$sealevel - gmsl_model_norm[indices$mod2obs,ss]})

  # Calculate the likelihood. The observations are not correlated. They are independent
  loglike <- sapply(1:n_simulations, function(ss) {
               sum(dnorm(residuals[,ss], mean=rep(0,n_data), sd=data_sealevel$uncertainty, log=TRUE))})

  return(loglike)
}
##==============================================================================


##==============================================================================
## Negative log-likelihood
neg_log_likelihood = function(parameters,
                              parnames,
                              data_sealevel,
                              indices,
                              temperature_forcing)
{return(-log_likelihood(parameters, parnames, data_sealevel, indices, temperature_forcing))}
##==============================================================================


##==============================================================================
## End
##==============================================================================
