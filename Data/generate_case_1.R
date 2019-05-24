

rm(list = ls())

### Code to generate an example  of normal background and anomlous densities
###  this can be used to replicate the results in Table 2
###  It requires to select the size of the type of radition source (celsium or cobalt),
### and its distance to the spetometer
## once these parameters have been selected, the code will produce the respective pre an post change densities
working_directory = "C:/Users/oscar/Desktop/Anomaly_detection/supplementary_materials/supplementary_materials/Experiments"


setwd(paste(working_directory,"/Code_J",sep=""))
source('utils.R')
setwd(paste(working_directory,"/Data_J",sep=""))

# How strong is the background rate, and how large is the source in milliCuries?
background_cps = 39  
source_size = 100#650  #mCi
# Distance to the source
dist_grid = c(50, 100,  150)
i = 1## index of distance used, if for example i = 1, then distance to source is 50m

##################################3


K =  length(dist_grid)
m = 2048

f0 =    matrix(0, K+1,m)


#f0[2,] =  sim_spectrum
#f0[3,] =  sim_spectrum2


###########################

for(i in 1:K)
{
  
  # Read in background
  background_train = as.numeric(read.csv('background_train_mean.csv', header=FALSE))
  plot(background_train, type='l')
  n_bins = length(background_train)
  
  # choose anomaly
  celsium  = 1 # if celsium = 1 then the anomaly is celsium is celsium-137. If set to zero, we use cobalt
  
  # Read in anomaly
  
  if(celsium==1)
  {
    cs137_spectrum = as.numeric(read.csv('2013-cs137-5cm.csv', header=FALSE))
  }
  if(celsium ==0)
  {
    cs137_spectrum = as.numeric(read.csv('co60-4min-320cps-c7-20cm.csv', header=FALSE))
  }
  
  
  # Winsorize: assign all counts after bin 2048 to bin 2048, and then truncate
  cs137_spectrum[2048] = sum(cs137_spectrum[2048:4096])
  cs137_spectrum = head(cs137_spectrum, 2048)
  
  
  # Normalize and plot
  cs137_spectrum = cs137_spectrum/sum(cs137_spectrum)
  plot(cs137_spectrum, type='l')
  
  
  # Create composite spectrum from background + anomaly
  this_dist =  dist_grid[i]
  sim_spectrum = inject_source(background_train, background_cps, cs137_spectrum, source_size, this_dist)
  sim_spectrum = sim_spectrum/sum(sim_spectrum)
  
  f0[i+1,] =  sim_spectrum
}
f0[1,] =  background_train


