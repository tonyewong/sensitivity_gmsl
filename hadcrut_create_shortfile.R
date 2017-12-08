##==============================================================================
## hadcrut_create_shortfile.R
##
## Note:  The file HadCRUT.4.4.0.0.annual_ns_avg_unc.txt was created by adding
##        the uncertainty estimate from one of the ensemble members, from:
## HadCRUT.4.4.0.0.annual_ns_avg_realisations/HadCRUT.4.4.0.0.annual_ns_avg.1.txt
##        (all members' uncertainty estimates are the same, if you check the files)
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

data_temperature <- read.table('HadCRUT.4.4.0.0.annual_ns_avg.txt')
times <- data_temperature[,1]
temperatures <- data_temperature[,2]

data_uncertainty <- read.table('HadCRUT.4.4.0.0.annual_ns_avg.1.txt')
uncertainties <- data_uncertainty[,3]

temperature_to_write <- cbind(times, temperatures, uncertainties)
colnames(temperature_to_write) <- c('year','temperature','uncertainty')

write.csv(file='HadCRUT.4.4.0.0.annual_ns_avg_unc.csv', x=temperature_to_write, row.names=FALSE)

##==============================================================================
## End
##==============================================================================
