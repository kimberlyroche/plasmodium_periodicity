import sys
import random

from include import get_strain_label
from include import get_expression_interpolated
from include import get_expression
from include import get_pchip_optima
from include import get_splined_peak
from include import get_troughs_for_peak
from include import get_peak_width_markers

from include import plot_expression_splined

# usage: python calculate_peak_widths.py [strain idx]
# e.g.: python calculate_peak_widths.py 0

strain = -1
file_suffix = '_p25_genelist'
genelist_dir = 'data/p25/'
working_dir = 'working/'
interpolation_method = 'linear'

if len(sys.argv) >= 2:
	strain = int(sys.argv[1])
else:
	print("No input for strain")
	exit()

fig_sz = 10

strain_label = get_strain_label(strain)

expression = get_expression(strain, '_p25_genelist', genelist_dir=genelist_dir)
expression_interpolated = get_expression_interpolated(strain, file_suffix, genelist_dir=genelist_dir, method=interpolation_method)

genes = list(expression_interpolated.index)
no_genes = len(expression_interpolated.index)

# calculate and write out peak widths as determined by critical points of PCHIP spline on each gene's expression data
outfile = open(working_dir+strain_label+'_peak_markers.txt', 'w')
for i in range(no_genes):
	if i % 1000 == 0:
		print("Iteration %d / %d..." % (i, no_genes))
	g = genes[i]
	(t1, half1, peak, half2, t2) = get_peak_width_markers(expression, expression_interpolated, g, strain)
	outfile.write("%s\t%d\t%d\t%d\t%d\t%d\n" % (g, t1, half1, peak, half2, t2))

outfile.close()
