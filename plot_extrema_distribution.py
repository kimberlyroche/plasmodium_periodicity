# this script plots histograms for peak-to-peak and trough-to-trough distances
# given all critical points (mins, maxes) identified for a given strain and recorded
# in working/{strain}_maxes_{interpolation}.csv and working/{strain}_mins_{interpolation}.csv
# 
# output are:
#   paired peak list in working/{strain}_maxes_{interpolation}.csv
#   paired trough list in working/{strain}_mins_{interpolation}.csv
#   histrogram of peak-peak distances in output/S8/{strain}_peakpeak_maxes_{interpolation}.svg
#   histrogram of trough-trough distances in output/S8/{strain}_peakpeak_mins_{interpolation}.svg

import numpy as np
import sys

from include import get_strain_label
from include import get_expression_interpolated
from include import get_expression
from include import get_optima_pair_distance

import matplotlib
# Force matplotlib to not use any Xwindows backend.
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

# usage: python plot_extrema_distribution.py [strain idx]
# e.g.: python plot_extrema_distribution.py 0

strain = -1
genelist_dir = "data/p25/"
working_dir = "working/"
output_dir = "output/"
interpolation_method = "PCHIP"

if len(sys.argv) >= 2:
	strain = int(sys.argv[1])
else:
	print("No input for strain")
	exit()

strain_label = get_strain_label(strain)

expression = get_expression(strain, '_p25_genelist', genelist_dir=genelist_dir)
expression_interpolated = get_expression_interpolated(strain, file_suffix='_p25_genelist',
	genelist_dir=genelist_dir, method=interpolation_method)

file_maxes = open(working_dir + strain_label + "_maxes_" + interpolation_method + ".csv", "w") # peaks
file_maxes.write("gene,expr,timepoint\n")
file_mins = open(working_dir + strain_label + "_mins_" + interpolation_method + ".csv", "w") # troughs
file_mins.write("gene,expr,timepoint\n")
max_deltas = []
min_deltas = []
for i in range(len(expression_interpolated.index)):
	if(i % 1000 == 0):
		print("Iteration " + str(i))
	gene = expression_interpolated.index[i]
	(max_d, min_d, all_maxes, all_mins) = get_optima_pair_distance(expression_interpolated, gene, return_all_optima=True)
	if(max_d >= 0):
		max_deltas.append(max_d)
	if(min_d >= 0):
		min_deltas.append(min_d)
	for m in all_maxes:
		file_maxes.write("%s,%f,%d\n" % (gene, m, all_maxes[m]))
	for m in all_mins:
		file_mins.write("%s,%f,%d\n" % (gene, m, all_mins[m]))

file_maxes.close()
file_mins.close()

fig = plt.figure(1, figsize=(15, 10))
ax = plt.subplot(111)
plt.hist(max_deltas, bins=70, rwidth=1)
plt.xlim(xmin=0., xmax=70.)
plt.xticks(np.arange(0, 70, step=2))
plt.ylim(ymin=0., ymax=1200.)
plt.savefig(output_dir + "S8/" + strain_label + "_peakpeak_maxes_" + interpolation_method + ".svg", bbox_inches="tight")
plt.close(fig)

fig = plt.figure(1, figsize=(15, 10))
ax = plt.subplot(111)
plt.hist(min_deltas, bins=70, rwidth=1)
plt.xlim(xmin=0., xmax=70.)
plt.xticks(np.arange(0, 70, step=2))
plt.ylim(ymin=0., ymax=1200.)
plt.savefig(output_dir + "S9/" + strain_label + "_peakpeak_mins_" + interpolation_method + ".svg", bbox_inches="tight")
plt.close(fig)
