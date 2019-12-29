# this script performs wraps previously interpolated transcriptome data around
# the specified indices
# 
# output are:
#   interpolated, wrapped, and average expression data in working/{strain}_expression_wrapped_{interpolation}.csv

import sys

from include import get_expression_interpolated
from include import get_concensus_wrap_indices
from include import get_strain_label
from include import apply_wrap_expression
#from include import scale_expression_data
#from include import get_min_error_wrap

# usage: python perform_expression_wrap.py [strain idx]
# e.g.: python perform_expression_wrap.py 0

strain_idx = -1
interpolation_method = 'PCHIP'

if len(sys.argv) >= 2:
	strain_idx = int(sys.argv[1])
else:
	print("No input for strain_idx")
	exit()

expression_interpolated = get_expression_interpolated(strain_idx, file_suffix='_p25_genelist',
	genelist_dir='data/p25/', method=interpolation_method)

# these indices must be set in include.py manually!
(start, end) = get_concensus_wrap_indices(strain_idx)

expression_wrapped = apply_wrap_expression(expression_interpolated, start, end, 0)

expression_wrapped.to_csv("working/" + get_strain_label(strain_idx) + "_expression_wrapped_" + interpolation_method + ".csv")

#scaled_expression_interpolated = scale_expression_data(expression_interpolated)

#scaled_expression_interpolated.to_csv("working/"+get_strain_label(strain_idx)+"_scaled_expression_wrapped_" + interpolation_method + ".csv")

# this will calculate the error and export the RMSE plots as .svg files
#get_min_error_wrap(scaled_expression_interpolated, strain_idx, microscopy=False, output_dir="output/S6/")
