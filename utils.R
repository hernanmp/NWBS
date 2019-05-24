sizeToCPS = function(mCi, dist) {
# function adapted from Alex Reinhart's code
# Convert a source size in mCi (miliCuries) to an approximate number of counts
# at dist (meters), using calibrations with Cs-137 sources of known sizes.

    knownSize = 0.000844 # mCi
    knownDist = 0.05 # meters
    knownCounts = 630 # cps
    mu = 0.0100029 # attenuation coefficient, in units of m^{-1}, for 660 keV in air

   (mCi / knownSize) * knownCounts * (knownDist / dist)^2 * exp(-mu * (knownDist + dist))
}


inject_source = function(background_spectrum, background_cps, source_spectrum, source_size, distance) {
# Creates a new spectrum based on a combination of the background and source spectra
# background_spectrum and source_spectrum: spectral densities, i.e. normalized to 1
# background_cps: how many counts per second (on average) does the background contribute to the detector?
# source_size: strength of source in mCi
# distance: distance to source in meters
	source_cps = sizeToCPS(source_size, distance)
	# take linear combination of background and source densities, weighted by CPS
	new_spectrum = background_cps*background_spectrum + source_cps*source_spectrum
  #plot(source_spectrum)
return(	new_spectrum)
}


cumulative_frac = function(x) cumsum(x)/sum(x)
