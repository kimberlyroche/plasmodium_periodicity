# this script performs interpolation of microscopy data from 3-hr to 1-hr time points
# using linear interpolation, then calculates error (RMSE) associated with all possible
# single wraps
# 
# output are:
#   plots of per-wrap RMSE in output/S7

import sys

from include import get_microscopy_interpolated
from include import get_min_error_wrap

# usage: python score_microscopy_wraps.py [strain idx]
# e.g.: python score_microscopy_wraps.py 0

strain_idx = -1
interpolation_method = 'PCHIP'

if len(sys.argv) >= 2:
	strain_idx = int(sys.argv[1])
else:
	print("No input for strain_idx")
	exit()

microscopy_interpolated = get_microscopy_interpolated(strain_idx, method=interpolation_method)

# this will calculate the error and export the RMSE plots as .svg files
get_min_error_wrap(microscopy_interpolated, strain_idx, microscopy=True, verbose=True, output_dir="output/S7/")
