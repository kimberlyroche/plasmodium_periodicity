# this script performs wraps previously interpolated transcriptome data around
# the specified indices
# 
# output are:
#   interpolated, wrapped, and average expression data in working/{strain}_microscopy_wrapped_{interpolation}.csv
#   plots of original, wrapped, and wrapped/averaged microscopy measurements

import sys

from include import get_microscopy_interpolated
from include import get_mapping
from include import get_series_length
from include import get_best_microscopy_wrap_indices
from include import get_concensus_wrap_indices
from include import apply_wrap_microscopy
from include import plot_microscopy_original
from include import plot_microscopy_overlapped
from include import plot_microscopy_averaged

# usage: python perform_microscopy_wrap.py [strain idx]
# e.g.: python perform_microscopy_wrap.py 0

strain_idx = -1
output_dir = 'output/S7'
interpolation_method = 'PCHIP'

if len(sys.argv) >= 2:
	strain_idx = int(sys.argv[1])
else:
	print("No input for strain_idx")
	exit()

microscopy_interpolated = get_microscopy_interpolated(strain_idx, method=interpolation_method)

# these indices must be set in include.py manually!
# these are the "best" (lowest-error) wraps by microscopy
(i, j) = get_best_microscopy_wrap_indices(strain_idx)
mapping = get_mapping(0, i, j, get_series_length(strain_idx)-1)

microscopy_wrapped = apply_wrap_microscopy(microscopy_interpolated, mapping, i, j, strain_idx)

plot_microscopy_original(microscopy_interpolated, strain_idx, i, j, output_dir="output/S7/")
plot_microscopy_overlapped(microscopy_interpolated, mapping, i, j, strain_idx, output_dir="output/S7/")
plot_microscopy_averaged(microscopy_interpolated, mapping, i, j, strain_idx, output_dir="output/S7/")

# these are the consensus "best" wraps
(i, j) = get_concensus_wrap_indices(strain_idx)
mapping = get_mapping(0, i, j, get_series_length(strain_idx)-1)

microscopy_wrapped = apply_wrap_microscopy(microscopy_interpolated, mapping, i, j, strain_idx)

plot_microscopy_original(microscopy_interpolated, strain_idx, i, j, output_dir="output/S10/")
plot_microscopy_overlapped(microscopy_interpolated, mapping, i, j, strain_idx, output_dir="output/S10/")
plot_microscopy_averaged(microscopy_interpolated, mapping, i, j, strain_idx, output_dir="output/S10/")
