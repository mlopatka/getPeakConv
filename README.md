# getPeakConv
Probabilistic peak detection for first order chromatographic data 

  getPeaksCov - Probabilistic peak detection of first order chromatographic
  or time-series data. For a provided pair of vectors (time and intensity)
  the function uses probabilistic reasoning to evaluate an exhaustive set 
  of (realistic) peak models in a windowed region. Prior probabilities are
  evaluated by the method published by Davis and Giddings 1983.
  Probabilistic reasoning is used to estimate the posterior probability
  that a given point is affected by a chromatographic peak.
 
  Input:
  p = getPeaksConv(x, y, bandWidth, ResSigma, ALPHA, MAX_PEAKS_PER_ROW)
  
  x : Time index of each measurement
  y : Intensity value of each measurement
  bandWidth: The standard deviation of a prototypical Gaussian peak for the
  chromatographic system expressed in terms of retention time
  ResSigma: Standard deviation of a blank measurement (baseline noise)
  MAX_PEAKS_PER_ROW (optional): The boound on number of peaks allowed to overlap in a
  chromatographic space define by 2 plates in the chromatographic system.
  Default = 2
  ALPHA (optional): The saturation of the chromatogram [1 > ALPHA > 0]
  Default = 0.25
 
  Output:
  p : Estiamted posterior probibility that the point at p is affected by a
  chromatographic peak
 
  Example usage:
  p = getPeaksConv(retention_time_vector, intenstiy_vector, 2.5, 1.2e-5, 0.3, 3);
 
  If this software is useful to your academic work, please cite our
  publication in lieu of thanks:
 
 %% Lopatka M., Vivó-Truyols G., Sjerps M.J., "Probabilistic peak detection
 %% for first-order chromatographic data." Analitica Chimica Acta. 2014.
 %% DOI: 10.1016/j.aca.2014.02.015
 
  Author: Martin Lopatka <martin.lopatka@gmail.com> Created: 29th August, 2013
  Gabriel Vivó-Truyols <g.vivotruyols@uva.nl> Revised: 23rd April, 2015
