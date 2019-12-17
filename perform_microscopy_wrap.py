import sys

from include import get_microscopy_interpolated
from include import get_mapping
from include import get_series_length
from include import get_wrap_indices
from include import get_min_error_wrap
from include import apply_wrap_microscopy
from include import plot_microscopy_original
from include import plot_microscopy_overlapped
from include import plot_microscopy_averaged

# usage: python perform_microscopy_wrap.py [strain idx]
# e.g.: python perform_microscopy_wrap.py 0

strain_idx = -1
file_suffix = '_p25_genelist'
genelist_dir = 'data/p25/'
output_dir = 'output/'
interpolation_method = 'linear'

if len(sys.argv) >= 2:
	strain_idx = int(sys.argv[1])
else:
	print("No input for strain_idx")
	exit()

microscopy_interpolated = get_microscopy_interpolated(strain_idx, method=interpolation_method)

(i, j) = get_wrap_indices(strain_idx)
mapping = get_mapping(0, i, j, get_series_length(strain_idx)-1)

microscopy_wrapped = apply_wrap_microscopy(microscopy_interpolated, mapping, i, j, strain_idx)

plot_microscopy_original(microscopy_interpolated, strain_idx, i, j)
plot_microscopy_overlapped(microscopy_interpolated, mapping, i, j, strain_idx)
plot_microscopy_averaged(microscopy_interpolated, mapping, i, j, strain_idx)
