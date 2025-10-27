import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_neb(images, splines, file_path, filename, filename_prefix):
	"""
	This method is designed to create the NEB energy profile plots.

	Parameters
	----------
	images : list of (float, float)
		These are the images of each NEB iteration, given as (arc length, energy).
	splines : list of (float, float)
		These are the splines of each NEB iteration, given as (arc length, energy).
	file_path : str.
		This is the path to the folder containing the NEB calculation.
	filename : str.
		This is the filename of the .interp file.
	filename_prefix : str.
		This is the prefix of the .interp file.
	"""

	# First, obtain the number of iterations that ORCA performed. 
	no_of_iters = len(images)

	# Second, create a subplot plot to plot the figure
	if len(images) > 1:
		fig, (ax1, ax2) = plt.subplots(1,2, gridspec_kw={'width_ratios': [30, 1]})
	else:
		fig, ax1, = plt.subplots(1,1)

	# Third, obtain the list of colours to use for plotting
	if len(images) > 1:
		colours = mpl.colormaps['rainbow'](np.linspace(0, 1, len(images)))[::-1]
	else:
		colours = ['tab:blue']

	# Fourth, add the plots for each NEB iteration to the figure
	zorder = 1000
	for index, ((arcS, Eimg), (arcS2, Eimg2), colour) in enumerate(zip(images, splines, colours)):
		ax1.plot   (arcS2, Eimg2, color=colour,  linestyle='-', zorder=zorder+1)
		ax1.scatter(arcS,  Eimg,  color='black', s=2, zorder=zorder)
		zorder += 2

	# Fifth, add a colourbar to the figure
	if len(images) > 1:
		no_of_ticks = 10
		ticklabels  = [int(index) for index in np.linspace(0, len(images)-1, no_of_ticks, endpoint=True)]
		ticks       = [index/(len(images)-1) for index in ticklabels]
		cb  = mpl.colorbar.ColorbarBase(ax2, cmap=mpl.colormaps['rainbow_r'], orientation='vertical', ticks=ticks)
		cb.set_ticklabels(ticklabels=ticklabels)

	# Sixth, save all the NEB optimisation iterations together in one plot
	ax1.set_xlabel("Displacement (Å)", fontsize=11)
	ax1.set_ylabel("Energy (eV)", fontsize=11)
	if len(images) == 1:
		if 'final' in filename.lower():
			ax1.set_title("Final Iteration")
		else:
			ax1.set_title(f"Only one iteration provided in {filename}")
	else:
		ax1.set_title( "Iterations: %i to %i" % (0, no_of_iters-1) )
		ax2.set_ylabel("Iterations")
	plt.savefig(file_path+'/'+filename_prefix+'_NEB_optimization.png', dpi=400)

	# Seventh, close the figure and clear matplotlib
	plt.close(fig)
	plt.clf()
	plt.cla()

	# Eighth, save the energy profile of the last NEB iteration.
	plt.plot   (arcS2, Eimg2, color='tab:blue', linestyle='-', zorder=zorder)
	plt.scatter(arcS,  Eimg,  color='black', s=2, zorder=zorder+1)
	plt.xlabel("Displacement (Å)", fontsize=11)
	plt.ylabel("Energy (eV)", fontsize=11)
	plt.savefig(file_path+'/'+filename_prefix+'_NEB_last_iteration.png', dpi=400)

# --------------------------------------------------------------------------------------------------------------------


