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
matplotlib.use('Agg')

import matplotlib.pyplot as plt

#import svgutils.transform as sg

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

# get time series length for a given strain's data file
def get_series_length(idx):
	lengths = [61, 70, 64, 70]
	return lengths[idx]

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
#    D6 = 36
#    SA250 = 54
#    FVO = 45
def get_best_microscopy_wrap_indices(idx):
	spans = [(0, 38), (0, 35), (0, 54), (0, 42)]
	return spans[idx]

# get consensus best wraps, fixed here as:
#    3D7 = 39
#    D6 = 36
#    SA250 = 54
#    FVO = 45
def get_concensus_wrap_indices(idx):
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

# get 0-hour time point for a given strain
# these are chosen by eye to be the time point closest to the troph-to-schizont transition in the
#	wrapped data for all strains, giving a closest-to-common starting point for all strains
def get_zero_hour(idx):
	zero_hours = [23, 20, 38, 36]
	return zero_hours[idx]

# get unaltered/original expression measurements
def get_expression(strain_idx, file_suffix, working_dir = 'working/', data_dir = 'data/', genelist_dir = ''):
	# min_max_scaler will throws a type error about the row column, headers unless read in with read_csv
	# didn't feel like debugging; if expression needs to be interpolated, do that, write it out, and read it back in
	return pd.read_csv(data_dir + "strain_" + get_strain_label(strain_idx) + "_fpkm_all_timepoints.txt", delimiter='\t', header=0, index_col=0)

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

# interpolate microscopy measurements from 3-hr frequency to 1-hr
# method can be 'linear' or 'PCHIP'
def interpolate_microscopy_measurements(original_df, interpolated_df, strain_idx, method='linear'):
	hour_label = get_hour_label(strain_idx)
	phase_labels = get_phase_labels(strain_idx)
	if method == 'linear':
		for i in range(len(original_df.index)):
			if i < len(original_df.index)-1:
				lower = original_df.at[i,hour_label]
				upper = original_df.at[i+1,hour_label]
				dx = int(upper) - int(lower)
				for label in phase_labels:
					lower_value = original_df.at[i,label]
					upper_value = original_df.at[i+1,label]
					dy = upper_value - lower_value
					slope = dy/dx
					for j in range(int(lower), int(upper)+1):
						multiplier = j-int(lower)
						interpolated_df.set_value(j, hour_label, j)
						if j == int(lower):
							interpolated_df.set_value(j, label, lower_value)
						elif j == int(upper):
							interpolated_df.set_value(j, label, upper_value)
						else:
							interpolated_df.set_value(j, label, lower_value + slope*multiplier)
	elif method == 'PCHIP':
		for label in phase_labels:
			phase_df = original_df.loc[:,[hour_label, label]]
			orig_x = phase_df[hour_label].tolist()
			orig_y = phase_df[label].tolist()
			pchip_spline = interpolate.PchipInterpolator(orig_x, orig_y)
			interp_x = np.arange(0, int(max(orig_x))+1, 1)
			interp_y = pchip_spline(interp_x)
			interpolated_df[hour_label] = interp_x
			interpolated_df[label] = interp_y
	else:
		print("Can't interpolate microscopy by method '%s'" % method)
		exit()

# get interpolated expression measurements
def get_expression_interpolated(strain_idx, file_suffix=None, working_dir = 'working/', data_dir = 'data/', genelist_dir = '', method = 'linear'):
	# min_max_scaler will throw a type error about the row column, headers unless read in with read_csv
	# this should be debugged (TBD) but as a workaround if expression needs to be interpolated, perform
	# the interpolation, write the results out, then read it back in
	interpolated_filepath = working_dir + get_strain_label(strain_idx) + "_expression_interpolated_" + method + ".csv"
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

# get interpolated microscopy measurements
def get_microscopy_interpolated(strain_idx, working_dir = 'working/', data_dir = 'data/', method = 'linear'):
	path = working_dir + get_strain_label(strain_idx) + "_microscopy_interpolated_" + method + ".csv"
	print("Path is: " + path)
	if os.path.isfile(path):
		microscopy_interpolated = pd.read_csv(path, index_col=0)
		return microscopy_interpolated
	else:
		microscopy = pd.read_table(data_dir + get_strain_label(strain_idx) + "_timeseries_full.txt")
		microscopy_interpolated = pd.DataFrame(np.zeros((len(microscopy.index)*3-2,4)), columns=[get_hour_label(strain_idx)]+get_phase_labels(strain_idx))
		interpolate_microscopy_measurements(microscopy, microscopy_interpolated, strain_idx, method)
		microscopy_interpolated.to_csv(path)
		return microscopy_interpolated

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
def get_min_error_wrap(df, strain_idx, microscopy=True, plot=True, verbose=False, output_dir='output/'):
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

# apply a wrapping defined by (start, end) to microscopy time series passed as a parameter ('original_df')
def apply_wrap_microscopy(original_df, row_mapping, start, end, strain_idx):
	hour_label = get_hour_label(strain_idx)
	phase_labels = get_phase_labels(strain_idx)
	# 4 incicates the columns we'll retain: hours, ring %, troph %, schizont %
	wrapped_df = pd.DataFrame(np.zeros(((end-start+1), 4)), columns=[hour_label]+phase_labels)
	idx = 0
	for key in row_mapping:
		wrapped_df.set_value(idx, hour_label, original_df.at[key, hour_label])
		if row_mapping[key] >= 0:
			wrapped_df.set_value(idx, phase_labels[0], (original_df.at[key, phase_labels[0]]+original_df.at[row_mapping[key], phase_labels[0]])/2)
			wrapped_df.set_value(idx, phase_labels[1], (original_df.at[key, phase_labels[1]]+original_df.at[row_mapping[key], phase_labels[1]])/2)
			wrapped_df.set_value(idx, phase_labels[2], (original_df.at[key, phase_labels[2]]+original_df.at[row_mapping[key], phase_labels[2]])/2)
		else:
			wrapped_df.set_value(idx, phase_labels[0], original_df.at[key, phase_labels[0]])
			wrapped_df.set_value(idx, phase_labels[1], original_df.at[key, phase_labels[1]])
			wrapped_df.set_value(idx, phase_labels[2], original_df.at[key, phase_labels[2]])
		idx = idx + 1
	return wrapped_df

# ============================================================================================
#   PEAK-TO-PEAK / TROUGH-TO-TROUGH CALCULATIONS
# ============================================================================================

# get maxes and mins from PCHIP spline of a given gene's expression
def get_pchip_optima(df, gene):
	x = [int(x) for x in df.columns]
	min_x = 0.
	max_x = float(x[len(x)-1])
	y = list(df.loc[gene])
	pchip_spline = interpolate.PchipInterpolator(x, y)
	# get interpolated points
	interp_x = np.arange(0, int(max_x+1), 0.05)
	# derivative
	der = pchip_spline.derivative(1)
	interp_der_y = der(interp_x, 0)
	# collect optima
	maxes = {}
	mins = {}
	for r in der.roots():
		if r >= min_x and r <= max_x:
			# get 2nd derivative (max vs. min)
			der2 = pchip_spline.derivative(2)
			interp_der2_y = der2([r], 0)
			if interp_der2_y[0] < 0:
				maxes[pchip_spline([r], 0)[0]] = r
			else:
				mins[pchip_spline([r], 0)[0]] = r
	return (maxes, mins)

def get_optima_pair_distance(df, gene, return_all_optima=False):
	(maxes, mins) = get_pchip_optima(df, gene)
	# format of maxes: {90.373199999999997: 54.0, 99.218999999999994: 15.0, 14.877700000000001: 33.0}
	# format of mins:  {1.8947619999999998: 3.0, 3.4686900000000001: 42.0, 30.000900000000001: 30.0}
	max_delta = -1
	min_delta = -1
	if(len(maxes) > 2):
		dkeys = list(maxes.keys())
		dkeys.sort()
		dkeys.reverse()
		max1 = maxes[dkeys[0]]
		max2 = maxes[dkeys[1]]
		max_delta = int(abs(max2-max1))
	if(len(mins) > 2):
		dkeys = list(mins.keys())
		dkeys.sort()
		min1 = mins[dkeys[0]]
		min2 = mins[dkeys[1]]
		min_delta = int(abs(min2-min1))
	return (max_delta, min_delta, maxes, mins)

# filter a single absolute maximum from list of local maxima
def get_splined_peak(maxes):
	# format of maxes: {1152.2: 0.0, 1144.25: 48.0, 1286.21: 9.0, 1623.8900000000001: 36.0}
	dkeys = list(maxes.keys())
	dkeys.sort()
	dkeys.reverse()
	# return peak (x, y)
	return (maxes[dkeys[0]], dkeys[0])

# get nearest troughs for a given maximum
def get_troughs_for_peak(peak_x, peak_y, mins):
	# format of mins: {736.51199999999994: 42.0, 77.524100000000004: 27.0, 895.822: 3.0}
	mins_sorted_by_x = sorted(mins.items(), key=lambda kv: kv[1])
	# format of mins_sorted_by_x: [(895.822, 3.0), (77.524100000000004, 27.0), (736.51199999999994, 42.0)]
	t1_x = t1_y = -1
	t2_x = t2_y = -1
	for t in mins_sorted_by_x:
		if t[1] < peak_x:
			t1_x = t[1]
			t1_y = t[0]
		if t[1] > peak_x and t2_x < 0:
			t2_x = t[1]
			t2_y = t[0]
	return ((t1_x, t1_y), (t2_x, t2_y))

# get half peak width markers from peak and trouph locations determined by derivative of spline for a
# given gene's expression
def get_peak_width_markers(original_df, interpolated_df, gene, strain, verbose=False):
	(maxes, mins) = get_pchip_optima(original_df, gene)
	(peak_x, peak_y) = get_splined_peak(maxes)
	((trough1_x, trough1_y), (trough2_x, trough2_y)) = get_troughs_for_peak(peak_x, peak_y, mins)

	halfway1 = -1
	closest_x_left = -1
	halfway2 = -1
	closest_x_right = -1
	if peak_x >= 0: # no reason it shouldn't be but
		if trough1_x >= 0:
			halfway1 = (peak_y - trough1_y)/2. + trough1_y
			closest_x_left = -1
			closest_error = np.inf
			for t in range(int(trough1_x), int(peak_x)+1):
				this_error = abs(interpolated_df.loc[gene][t] - halfway1)
				if this_error < closest_error:
					closest_error = this_error
					closest_x_left = t
		if trough2_x >= 0:
			halfway2 = (peak_y - trough2_y)/2. + trough2_y
			closest_x_right = -1
			closest_error = np.inf
			for t in range(int(peak_x), int(trough2_x)):
				this_error = abs(interpolated_df.loc[gene][t] - halfway2)
				if this_error < closest_error:
					closest_error = this_error
					closest_x_right = t

	if verbose:
		print("Positions: (%d -- %d -- %d -- %d -- %d)" % (trough1_x, closest_x_left, peak_x, closest_x_right, trough2_x))

	return (trough1_x, closest_x_left, peak_x, closest_x_right, trough2_x)

# ============================================================================================
#   VISUALIZATION
# ============================================================================================

# plot all wrap error scores in 2D (error as a function of length (end - start))
def plot_all_scores(all_scores, strain_label, microscopy=True, output_dir='output/'):
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

# plot PCHIP splined expression for a given gene
def plot_expression_splined(df, gene, strain_idx, fig_sz):
	strain_label = get_strain_label(strain_idx)
	(spline, x, y, interp_x, interp_y) = get_pchip_spline(df, gene)
	fig = plt.figure()
	ax = plt.subplot(111)
	plt.plot(x, y, 'x', interp_x, interp_y)
	plt.title(gene + " (PCHIP spline)")
	plt.savefig(strain_label+'_'+gene+'_PCHIPsplined.png', bbox_inches='tight')
	plt.close(fig)

def choose_ylim(strain_idx):
	y_lim=100.
	if(strain_idx == 0):
		y_lim=1.
	return(y_lim)

# plot microscopy data linearly
def plot_microscopy_original(original_df, strain_idx, start=-1, end=-1, output_dir='output/'):
	y_lim = choose_ylim(strain_idx)
	strain_label = get_strain_label(strain_idx)
	hour_label = get_hour_label(strain_idx)
	phase_labels = get_phase_labels(strain_idx)
	fig = plt.figure(1, figsize=(15, 5))
	ax = plt.subplot(111)
	if start >= 0:
		ax.fill_between([0,start], 0, y_lim, facecolor='gray', alpha=0.3)
	if end >= 0:
		ax.fill_between([(end+1),original_df.at[len(original_df.index)-1,hour_label]], 0, y_lim, facecolor='gray', alpha=0.3)
	ax.plot(original_df[[hour_label]], original_df[[phase_labels[0]]], marker="o", color="red")
	ax.plot(original_df[[hour_label]], original_df[[phase_labels[1]]], marker="o", color="green")
	ax.plot(original_df[[hour_label]], original_df[[phase_labels[2]]], marker="o", color="blue")
	ax.set_xlim(0, get_series_length(strain_idx)-1)
	ax.set_ylim(0, y_lim)
	plt.savefig(output_dir + strain_label+'_microscopy_original.svg', bbox_inches='tight')
	plt.close(fig)

# plot microscopy data with specified overlap
def plot_microscopy_overlapped(original_df, min_error_mapping, start, end, strain_idx, output_dir='output/'):
	y_lim = choose_ylim(strain_idx)
	strain_label = get_strain_label(strain_idx)
	hour_label = get_hour_label(strain_idx)
	phase_labels = get_phase_labels(strain_idx)
	alpha = 0.25
	fig = plt.figure(1, figsize=(15, 5))
	ax = plt.subplot(111)
	color = ["red", "green", "blue"]
	for c in range(3):
		# --------------------------------------------------------------------------------
		# plot core points
		# --------------------------------------------------------------------------------
		xs = list(range(start,end+1))
		ys = []
		for x in xs:
			ys.append(original_df.at[x, phase_labels[c]])
		ax.plot(xs, ys, marker="o", color=color[c], alpha=alpha)
		# --------------------------------------------------------------------------------
		# plot wrapped points
		# --------------------------------------------------------------------------------
		xs = []
		ys = []
		for m in min_error_mapping:
			if min_error_mapping[m] < 0:
				break
			xs.append(m)
			ys.append(original_df.at[min_error_mapping[m], phase_labels[c]])
		ax.plot(xs, ys, marker="o", color=color[c])
		xs = []
		ys = []
		for m2 in min_error_mapping:
			if m2 < m or min_error_mapping[m2] < 0:
				continue
			xs.append(m2)
			ys.append(original_df.at[min_error_mapping[m2], phase_labels[c]])
		ax.plot(xs, ys, marker="o", color=color[c])
	ax.set_xlim(start, end)
	ax.set_ylim(0, y_lim)
	plt.savefig(output_dir + strain_label+'_microscopy_overlap_'+str(start)+'_'+str(end)+'.svg', bbox_inches='tight')
	plt.close(fig)

# plot microscopy data with overlapped points averaged
def plot_microscopy_averaged(original_df, min_error_mapping, start, end, strain_idx, offset=0, output_dir='output/'):
	y_lim = choose_ylim(strain_idx)
	wrapped_df = apply_wrap_microscopy(original_df, min_error_mapping, start, end, strain_idx)
	strain_label = get_strain_label(strain_idx)
	hour_label = get_hour_label(strain_idx)
	phase_labels = get_phase_labels(strain_idx)
	fig = plt.figure(1, figsize=(15, 5))
	ax = plt.subplot(111)
	if offset > 0:
		xs = [0. for i in range(len(wrapped_df.index))]
		ring_ys = [0. for i in range(len(wrapped_df.index))]
		troph_ys = [0. for i in range(len(wrapped_df.index))]
		schizont_ys = [0. for i in range(len(wrapped_df.index))]
		for i in wrapped_df.index:
			hour =  int(wrapped_df.loc[i][hour_label])
			shifted = get_shifted_hour(hour, offset, (end - start + 1))
			xs[shifted] = shifted
			ring_ys[shifted] = wrapped_df.loc[i][phase_labels[0]]
			troph_ys[shifted] = wrapped_df.loc[i][phase_labels[1]]
			schizont_ys[shifted] = wrapped_df.loc[i][phase_labels[2]]
		ax.plot(xs, ring_ys, marker="o", color="red")
		ax.plot(xs, troph_ys, marker="o", color="green")
		ax.plot(xs, schizont_ys, marker="o", color="blue")
	else:
		ax.plot(wrapped_df[[hour_label]], wrapped_df[[phase_labels[0]]], marker="o", color="red")
		ax.plot(wrapped_df[[hour_label]], wrapped_df[[phase_labels[1]]], marker="o", color="green")
		ax.plot(wrapped_df[[hour_label]], wrapped_df[[phase_labels[2]]], marker="o", color="blue")
	ax.set_xlim(start, end)
	ax.set_ylim(0, y_lim)
	plt.savefig(output_dir + strain_label+'_microscopy_average_'+str(start)+'_'+str(end)+'.svg', bbox_inches='tight')
	plt.close(fig)

# plot PCHIP splined expression for a given gene
def plot_expression_splined(df, gene, strain_idx, fig_sz):
	strain_label = get_strain_label(strain_idx)
	(spline, x, y, interp_x, interp_y) = get_pchip_spline(df, gene)
	fig = plt.figure()
	ax = plt.subplot(111)
	plt.plot(x, y, 'x', interp_x, interp_y)
	plt.title(gene + " (PCHIP spline)")
	plt.savefig(strain_label+'_'+gene+'_PCHIPsplined.png', bbox_inches='tight')
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
























