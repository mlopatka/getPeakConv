# getPeakConv
Probabilistic peak detection for first order chromatographic data 

This fast implemetation uses MATLAB's conv function to perform the sliding window operation along the chromatogram's length. The code uses Bayes' Theorem to evaluate the exhaustive set of realistic peak configuration models within a given window size (expresed as data points in a chromatogram).
