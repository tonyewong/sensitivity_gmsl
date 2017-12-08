##==============================================================================
## Simple, mechanistically-motivated model of global mean sea level (GMSL)
## Based on Rahmstorf 2007, Equation 1
## Integrates using first-order explicit (Euler) scheme.  Tests suggest that the
## scheme is stable (by using smaller step-sizes)
##==============================================================================
##  Requires (input variables):
##  - Tg        global temperature [degC]
##
##  Simulates (output variables):
##  - gmsl      global mean sea level (mm sle)]
#
##  Parameters:
##  - sl_temp_sens  sensitivity of equilibrium sea level [mm sle/year/degC]
##  - temp_equil    equilibrium temperature (at which there is no change in sea
##                  level) [deg C]
##  - sl0           initial sea-level (take equal to 0 and normalize to match obs)
##                  (consistent with normalizing the forcing temperatures to 0 at
##                  the beginning of the simulation, which is done in the driver)
##  - tstep         time step (in years)
##==============================================================================

gmsl_r07 <- function( sl_temp_sens = 3.4,
                      temp_equil = 0,
                      sl0 = 0,
                      tstep = 1,
                      temperature_forcing)
{
  np      <- length(temperature_forcing)
  gmsl    <- rep(NA, np)
  gmsl[1] <- sl0
  for (i in 2:np) {
    gmsl[i] <- gmsl[i-1] + tstep*sl_temp_sens*(temperature_forcing[i-1] - temp_equil)
  }
  return(gmsl)
}

##==============================================================================
## End
##==============================================================================
