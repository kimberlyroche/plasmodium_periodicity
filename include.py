import math
import pandas as pd
import numpy as np
import sys
import os
from sklearn import preprocessing
from scipy import interpolate

import matplotlib
# for cluster usage:
#   force matplotlib to not use any Xwindows backend
#   use this if running on a cluster; remove this line for desktop use
# matplotlib.use('Agg')

import matplotlib.pyplot as plt

import svgutils.transform as sg

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# ============================================================================================
#   HOUSEKEEPING
# ============================================================================================

# we'll adhere to these (arbitrary) indices for each P. falciparum strain throughout
def get_strain_label(idx):
	# strain indices
	#	0: 3D7
	#	1: D6
	#	2: SA250
	#	3: FVO
	if idx < 0 or idx > 3:
		return None
	strain_labels = ['3D7', 'D6', 'SA250', 'FVO']
	return strain_labels[idx]

# normalize expression data
def scale_expression_data(df):
	expression_T = df.transpose().values
	min_max_scaler = preprocessing.MinMaxScaler()
	expression_T_scaled = min_max_scaler.fit_transform(expression_T)
	return pd.DataFrame(expression_T_scaled, columns=df.index, index=df.columns).transpose()

# get indices to of best wrap determined previously by microscopy (given a desired cycle length)
# having evaluated these already (via the functions included), we're fixing them here
# cycle lengths fixed as:
#    3D7 = 39
#    D6 = 38
#    SA250 = 54
#    FVO = 45
def get_wrap_indices(idx):
	spans = [(0, 38), (0, 35), (0, 53), (0, 44)]
	return spans[idx]

# get hour label for a given strain's data file
def get_hour_label(idx):
	if idx < 0 or idx > 4:
		return None
	hour_labels = ['Time_hours', 'Hour', 'Hour', 'Time', 'Hour']
	return hour_labels[idx]

# get phase labels for a given strain's data file
def get_phase_labels(idx):
	if idx < 0 or idx > 4:
		return None
	phase_labels = [
		['ring', 'troph', 'schizont'],
		['Ring', 'Troph', 'Schizont'],
		['Ring', 'Troph', 'Schizont'],
		['Percent_ring_sum', 'Percent_troph_sum', 'Percent_schizont_sum'],
		['Ring', 'Troph', 'Schizont']
	]
	return phase_labels[idx]

# ============================================================================================
#   INTERPOLATION
# ============================================================================================

# interpolate expression measurements to 3-hr frequency to 1-hr
# method can be 'linear' or 'PCHIP'
def interpolate_expression_measurements(filtered_df, interpolated_df, method = 'linear'):
	if method == 'linear':
		for i in range(len(filtered_df.columns)):
			if i < len(filtered_df.columns)-1:
				lower = filtered_df.columns[i]
				upper = filtered_df.columns[i+1]
				# print("\tInterpolating between %d and %d" % (int(lower), int(upper)))
				dx = int(upper) - int(lower)
				for gene_id in filtered_df.index:
					dy = filtered_df.loc[gene_id][upper] - filtered_df.loc[gene_id][lower]
					slope = dy/dx
					lower_value = filtered_df.loc[gene_id][lower]
					for j in range(int(lower), int(upper)+1):
						multiplier = j-int(lower)
						if j == int(lower):
							interpolated_df.set_value(gene_id, j, lower_value)
						elif j == int(upper):
							interpolated_df.set_value(gene_id, j, filtered_df.loc[gene_id][upper])
						else:
							interpolated_df.set_value(gene_id, j, lower_value + slope*multiplier)
	elif method == 'PCHIP':
		count = 1
		for gene_id in filtered_df.index:
			if count % 1000 == 0:
				print("Gene #%d" % count)
			gene_df = filtered_df.loc[gene_id,:]
			orig_x = np.asarray([int(i) for i in filtered_df.columns.tolist()])
			orig_y = np.asarray(gene_df.tolist()) # returns a Series before tolist()
			pchip_spline = interpolate.PchipInterpolator(orig_x, orig_y)
			interp_x = np.arange(0, int(max(orig_x))+1, 1)
			interp_y = pchip_spline(interp_x)
			for i in range(len(interp_x)):
				interpolated_df.set_value(gene_id, interp_x[i], interp_y[i])
			count = count + 1
	else:
		print("Can't interpolate expression by method '%s'" % method)
		exit()

# get interpolated expression measurements
def get_expression_interpolated(strain_idx, file_suffix=None, working_dir = 'working/', data_dir = 'data/', genelist_dir = '', method = 'linear'):
	# min_max_scaler will throw a type error about the row column, headers unless read in with read_csv
	# this should be debugged (TBD) but as a workaround if expression needs to be interpolated, perform
	# the interpolation, write the results out, then read it back in
	interpolated_filepath = working_dir + get_strain_label(strain_idx) + "_expression_interpolated.csv"
	print("Interpolating genes : " + interpolated_filepath)
	if not os.path.isfile(interpolated_filepath):
		if file_suffix == None:
			# no period genelist to filter by
			print("Interpolating all genes...")
			expression = pd.read_csv(data_dir + "strain_" + get_strain_label(strain_idx) + "_fpkm_all_timepoints.txt", delimiter='\t', header=0, index_col=0)
			max_col = int(expression.columns[len(expression.columns)-1])
			expr_interpolated = pd.DataFrame(index=expression.index, columns=range(max_col+1))
			interpolate_expression_measurements(expression, expr_interpolated, method)
			expr_interpolated.to_csv(interpolated_filepath)
		else:
			periodic_genes = pd.read_csv(genelist_dir + get_strain_label(strain_idx) + file_suffix, header=None, index_col = 0)
			print("Interpolating %d genes..." % len(list(periodic_genes.index)))
			expression = pd.read_csv(data_dir + "strain_" + get_strain_label(strain_idx) + "_fpkm_all_timepoints.txt", delimiter='\t', header=0, index_col=0)
			expr_filtered = expression.loc[list(periodic_genes.index)]
			max_col = int(expr_filtered.columns[len(expr_filtered.columns)-1])
			expr_interpolated = pd.DataFrame(index=expr_filtered.index, columns=range(max_col+1))
			interpolate_expression_measurements(expr_filtered, expr_interpolated, method)
			expr_interpolated.to_csv(interpolated_filepath)
	expr_interpolated = pd.read_csv(interpolated_filepath, index_col=0)
	return expr_interpolated

# ============================================================================================
#   WRAPPING
# ============================================================================================

# calculate a mapping of original time points to overlapped time points
#
# setup is like this with indices (k, i, j, l)
#
# . . x x x x                                         * * * * * * * . . . . . .
# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#     ^       ^                                     ^             ^
#     k       i                                     j             l
#
# where the final wrap will be
#
#            (l)                 (k)
# (j+1)       |                   | (i-1)
# |           |                   |     |
# * * * * * * *                   x x x x
# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
# ^                                     ^
# i                                     j
#
def get_mapping(k, i, j, l):
	mapping = {}
	for idx in range(i, j+1):
		mapping[idx] = -1
	for idx in range(k, i):
		overlap_idx = idx - i + j + 1
		mapping[overlap_idx] = idx
	for idx in range(j+1, l+1):
		overlap_idx = idx - j + i - 1
		mapping[overlap_idx] = idx
	return mapping

# get the indices of time points overlapping given the wrap defined by (k, i, j, l)
def get_overlapped_pairs(k, i, j, l):
	overlapped_pairs = []
	for idx in range(k, i):
		overlap_idx = idx - i + j + 1
		overlapped_pairs.append((idx, overlap_idx))
	for idx in range(j+1, l+1):
		overlap_idx = idx - j + i - 1
		overlapped_pairs.append((idx, overlap_idx))
	return overlapped_pairs

# calculate microscopy overlap error
def calc_microscopy_overlap_error(df_micro, overlapped_pairs, hour_label, phase_labels):
	microscopy_error = 0.
	for idx in range(len(overlapped_pairs)):
		t1 = overlapped_pairs[idx][0]
		t2 = overlapped_pairs[idx][1]
		# calculate microscopy error
		for phase in range(3):
			# ring = 0, troph = 1, schizont = 2
			microscopy_error = microscopy_error + (df_micro.at[df_micro.index[t1], phase_labels[phase]] - df_micro.at[df_micro.index[t2], phase_labels[phase]])**2
	microscopy_error = np.sqrt(microscopy_error / len(overlapped_pairs))
	return microscopy_error

# calculate expression overlap error
def calc_expr_overlap_error(df_expr, overlapped_pairs):
	use_normalized_error = False
	expression_error = 0.
	for idx in range(len(overlapped_pairs)):
		t1 = overlapped_pairs[idx][0]
		t2 = overlapped_pairs[idx][1]
		# calculate expression error
		for gene in df_expr.index:
			expression_error = expression_error + (df_expr.at[gene, df_expr.columns[t1]] - df_expr.at[gene, df_expr.columns[t2]])**2
	if use_normalized_error:
		gene_maxes = df_expr.max(axis=1).to_frame()
		gene_mins = df_expr.min(axis=1).to_frame()
		gene_ranges = gene_maxes - gene_mins
		gene_range_sum = 0.
		for gene in gene_ranges.index:
			gene_range_sum = gene_range_sum + (gene_ranges.at[gene,0])**2
		expression_error = np.sqrt(expression_error) / (len(overlapped_pairs) * np.sqrt(gene_range_sum))
	else:
		expression_error = np.sqrt(expression_error / len(overlapped_pairs))
	return expression_error

# calculate microscopy and/or expression overlap error
def calc_overlap_error(df_micro, df_expr, k, i, j, l, hour_label, phase_labels, microscopy_only=False, expression_only=False):
	overlapped_pairs = get_overlapped_pairs(k, i, j, l)
	if len(overlapped_pairs) > 0:
		m_error = 0.
		e_error = 0.
		if not expression_only:
			m_error = calc_microscopy_overlap_error(df_micro, overlapped_pairs, hour_label, phase_labels)
		if not microscopy_only:
			e_error = calc_expr_overlap_error(df_expr, overlapped_pairs)
		return (m_error, e_error)
	else:
		return (0., 0.)

# calculates a minimum error (RMSE) wrapping for a microscopy/expression time series passed as
# a parameter this function is wrapped by get_min_error_wrap()
def get_min_error_wrap_sub(df, strain_idx, microscopy, verbose=False):
	# this assumes no truncation: k = 0, l = span-1
	span = 0
	if microscopy:
		span = len(df.index)
	else:
		span = len(df.columns)
	if verbose:
		print("i\tj\tlength\terror")
	# assume maximum number of wraps == 2 for now (reasonable from microscopy)
	window_min = math.ceil(span/2.)
	# i (starting index) can take values in range 0 ... span-window_min+1
	scores_matrix = [{} for i in range(0, span-window_min+1)]
	min_error = np.inf
	min_error_core = None
	min_error_mapping = None
	hour_label = get_hour_label(strain_idx)
	phase_labels = get_phase_labels(strain_idx)
	for i in range(span):
		for j in range(i + window_min - 1, span):
			# force at least one point to wrap
			test_span_len = j - i + 1
			if test_span_len < span:
				error = 0.
				if microscopy:
					(error, discard) = calc_overlap_error(df, None, 0, i, j, span-1, hour_label, phase_labels, microscopy_only=True)
				else:
					# use expression scores
					(discard, error) = calc_overlap_error(None, df, 0, i, j, span-1, hour_label, phase_labels, expression_only=True)
				if error < min_error:
					min_error = error
					min_error_core = [i, j]
					min_error_mapping = get_mapping(0, i, j, (span-1))
				scores_matrix[i][test_span_len] = error
				if verbose:
					print("%d\t%d\t%d\t%.4f" % (i, j, (j-i+1), error))
	return (min_error, min_error_core[0], min_error_core[1], min_error_mapping, scores_matrix)

# use get_min_error_wrap_sub() to find minimum error (RMSE) wrap for microscopy or expression data
def get_min_error_wrap(df, strain_idx, microscopy=True, plot=True, verbose=False, output_dir='./'):
	strain_label = get_strain_label(strain_idx)
	(min_error, min_error_start, min_error_end, min_error_mapping, all_scores) = get_min_error_wrap_sub(df, strain_idx, microscopy, verbose)
	if plot:
		plot_all_scores(all_scores, strain_label, microscopy, output_dir)
	print("Min " + strain_label + " error achieved with (%d, %d): %.3f" % (min_error_start, min_error_end, min_error))
	# rationale behind the extra 1hr: a cycle's boundaries aren't at time points [x, y] but halfway between x's predecessor and x and y's successor and y
	#	so inflate out halfway to the next time point on either end: 2 x 30 min = 1 hr
	period = min_error_end - min_error_start + 1
	print("Period length: %d" % period)
	return (min_error_start, min_error_end, min_error_mapping)

# apply a wrapping defined by (trucate, start, end) to expression time series passed as a parameter ('original_df')
# 'truncate' is the index (0 if not applicable) from before which measurements are dropped from the series
# 'offset' is the zero-hour measurement to which to circularly shift all measurements (0 if not applicable)
def apply_wrap_expression(original_df, start, end, truncate, offset=0):
	# [start, end] are inclusive and indicate a single cell cycle
	# wrap edges and average values
	period = end - start + 1
	wrapped_df = pd.DataFrame(np.zeros((len(original_df.index), period)), index=original_df.index, columns=list(range(period)))
	counts = [0 for j in range(period)]
	for j in original_df.columns:
		for i in original_df.index:
			if int(j) < start:
				if int(j) >= truncate:
					# map j to (j + period - start)
					adj_j = int(j) + period - start
					wrapped_df.set_value(i, adj_j, wrapped_df.at[i, adj_j] + original_df.at[i, j])
					if i == original_df.index[0]:
						counts[adj_j] = counts[adj_j] + 1.
			elif int(j) > end:
				# map j to (j - period - start)
				adj_j = int(j) - period - start
				wrapped_df.set_value(i, adj_j, wrapped_df.at[i, adj_j] + original_df.at[i, j])
				if i == original_df.index[0]:
					counts[adj_j] = counts[adj_j] + 1.
			else:
				# map j to (j - start)
				adj_j = int(j) - start
				wrapped_df.set_value(i, adj_j, wrapped_df.at[i, adj_j] + original_df.at[i, j])
				if i == original_df.index[0]:
					counts[adj_j] = counts[adj_j] + 1.
	for j in range(len(counts)):
		if counts[j] > 1.:
			# average over
			for i in original_df.index:
				wrapped_df.set_value(i, j, wrapped_df.at[i,j] / counts[j])
	if offset > 0:
		# reorder expression so the "zero hour" is the first time point
		first_pt = [str(val) for val in wrapped_df.columns[offset:]]
		second_pt = [str(val) for val in wrapped_df.columns[:offset]]
		new_columns = first_pt + second_pt
		wrapped_df.columns = [str(val) for val in wrapped_df.columns]
		wrapped_df = wrapped_df[new_columns]
		wrapped_df.columns = [i for i in range(period)]
	return wrapped_df

# ============================================================================================
#   VISUALIZATION
# ============================================================================================

# plot all wrap error scores in 2D (error as a function of length (end - start))
def plot_all_scores(all_scores, strain_label, microscopy=True, output_dir='./'):
	fig = plt.figure(figsize=(10, 10))
	markersize = 100
	for i in range(len(all_scores)):
		xs = []
		ys = []
		# commenting out random color generation
		# c = np.random.rand(3,)
		c = 'black'
		for key in all_scores[i]:
			xs.append(key)
			ys.append(all_scores[0][key])
		plt.scatter(xs, ys, linewidth=0, s=markersize, c=[c])
	if microscopy:
		plt.savefig(output_dir + strain_label + '_microscopy_error.svg', bbox_inches='tight')
	else:
		plt.savefig(output_dir + strain_label + '_expression_error.svg', bbox_inches='tight')
	plt.close(fig)

# composite figure code (future)

# #create new SVG figure
# fig = sg.SVGFigure("16cm", "6.5cm")

# # load matpotlib-generated figures
# fig1 = sg.fromfile('sigmoid_fit.svg')
# fig2 = sg.fromfile('anscombe.svg')

# # get the plot objects
# plot1 = fig1.getroot()
# plot2 = fig2.getroot()
# plot2.moveto(280, 0, scale=0.5)

# # add text labels
# txt1 = sg.TextElement(25,20, "A", size=12, weight="bold")
# txt2 = sg.TextElement(305,20, "B", size=12, weight="bold")

# # append plots and labels to figure
# fig.append([plot1, plot2])
# fig.append([txt1, txt2])

# # save generated SVG files
# fig.save("fig_final.svg")
























