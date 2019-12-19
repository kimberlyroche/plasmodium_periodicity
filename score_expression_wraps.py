# this script performs interpolation of transcriptome data from 3-hr to 1-hr time points
# using a PCHIP spline, then calculates error (RMSE) associated with all possible single
# wraps
# 
# output are:
#   plots of per-wrap RMSE in output/S6
#   interpolated expression data in working/{strain}_expression_interpolated_{interpolation}.csv

import sys

from include import get_expression_interpolated
from include import scale_expression_data
from include import get_min_error_wrap

# usage: python score_expression_wraps.py [strain idx]
# e.g.: python score_expression_wraps.py 0

strain_idx = -1
interpolation_method = 'PCHIP'

if len(sys.argv) >= 2:
	strain_idx = int(sys.argv[1])
else:
	print("No input for strain_idx")
	exit()

expression_interpolated = get_expression_interpolated(strain_idx, file_suffix='_p25_genelist',
	genelist_dir='data/p25/', method=interpolation_method)

scaled_expression_interpolated = scale_expression_data(expression_interpolated)

# this will calculate the error and export the RMSE plots as .svg files
get_min_error_wrap(scaled_expression_interpolated, strain_idx, microscopy=False, output_dir="output/S6/")
