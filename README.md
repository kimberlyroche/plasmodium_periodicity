# plasmodium_periodicity

These are scripts associated with the P. falciparum oscillator analysis of Smith et al. (2020) that perform interpolation, scoring, and "wrapping" of transcriptome and microscopy data.

### Folder: data
Expression and microscopy data and periodic gene lists associated with 4 strains of P. falciparum

### Folder: working
Intermediate output, including interpolated data and lists of identified expression minima and maxima

### Folder: output
Figures S6 through S10 in .svg format

### File: score_expression_wraps.py
Performs interpolation of expression data and exhausting wrapping and wrap error-calculation

### File: score_microscopy_wraps.py
Performs interpolation of microscopy data and exhausting wrapping and wrap error-calculation

### File: perform_expression_wrap.py
Wraps expression data around desired start, end time points

### File: perform_microscopy_wrap.py
Wraps microscopy data around desired start, end time points

### File: plot_extrema_distribution.py
Calculates distances between peak-peak and trough-trough pairs renders histograms of these distributions across all genes
