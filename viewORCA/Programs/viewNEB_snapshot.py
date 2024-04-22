import os
import numpy as np
from ase import Atoms, Atom
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator
from ase.visualize import view
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from copy import deepcopy

class CLICommand:
	"""This module allows the user to view the progress of an ORCA Nudged Elastic Band (NEB) job (as it is running or after it has run) by viewing the energy profiles of all the iterations of the NEB optimisation job.
	"""

	@staticmethod
	def add_arguments(parser):
		parser.add_argument('-p', '--path',     nargs=1, help='This is the path to the folder containing the ORCA NEB job.', default=['.'])
		parser.add_argument('-f', '--filename', nargs=1, help='This is the name of the interpolation file. If this is not given, viewORCA will use the ".interp" file in the NEB folder.', default=[None])

	@staticmethod
	def run(args):
		"""
		Run this program. This will allow the user to view the optimisation process for an ORCA NEB job. 

		Parameters
		----------
		args : argparse.Namespace
			This is the namespace that give all the information about the arguments.
		"""

		# First, obtain the local path to the ORCA NEB job.
		path_to_images = args.path
		if len(path_to_images) != 1:
			raise Exception('Error: input "path" must only have one input: '+str(path_to_images))
		path_to_images = path_to_images[0]

		# Second, obtain the filename of the interpolation file you want to use
		filename = args.filename
		if len(path_to_images) != 1:
			raise Exception('Error: input "filename" must only have one input: '+str(filename))
		filename = filename[0]

		# Third, if filename is None, locate the iterpolation file by locating the file contain the ".interp" suffix.
		if filename is None:
			filename = get_interpolation_filename(path_to_images)

		# Fourth, run the method.
		Run_Method(path_to_images, filename=filename)

def Run_Method(path_to_images, filename):
	"""
	This will allow the user to view the optimisation process for an ORCA NEB job. 

	Parmeters
	---------
	path_to_images : str.
		This is the path to the ORCA NEB jobs. 
	filename : str.
		This is the name of the file to check
	"""

	# First, obtain the ORCA images and splines from the interpolation file. 
	images, splines = read_orca_interpolation_file(path_to_images+'/'+filename)

	# Second, obtain the number of iterations that ORCA performed. 
	no_of_iters = len(images)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
	# Third, perform several checks.

	# 3.1: Make sure that the number of images equal the number of splines obtained from the interpolation file. 
	if len(images) != len(splines):
		to_string  = 'Error: images list is not the same as the splines list.\n'
		to_string += f'len(images) = {len(images)}\n'
		to_string += f'len(splines) = {len(splines)}\n'
		to_string += f'This may indicate the {filename} file has been cut short, and missing some spline data.'
		to_string += 'Check this.'
		raise RuntimeError(to_string)

	# 3.2: Give a note to the user if the interpolation file only contain one iteration. 
	if no_of_iters == 1:
		if 'final' in filename.lower():
			print('*** Note that %s contains only the last iteration of a NEB/CI-NEB run  ***' % (path_to_images+'/'+filename))
		else:
			print('%s contains only one iteration?   ***' % filename)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
	# Fourth, create a plot giving the optimisation information from the ORCA NEB job. 

	# 4.1: Create a subplot plot to plot the figure
	if len(images) > 1:
		fig, (ax1, ax2) = plt.subplots(1,2, gridspec_kw={'width_ratios': [30, 1]})
	else:
		fig, ax1, = plt.subplots(1,1)

	# 4.2: Obtain the list of colours to use for plotting
	if len(images) > 1:
		colours = mpl.colormaps['rainbow'](np.linspace(0, 1, len(images)))[::-1]
	else:
		colours = ['tab:blue']

	# 4.3: Add the plots for each NEB iteration to the figure
	zorder = 1000
	for index, ((arcS, Eimg), (arcS2, Eimg2), colour) in enumerate(zip(images, splines, colours)):
		ax1.plot   (arcS2, Eimg2, color=colour,  linestyle='-', zorder=zorder+1)
		ax1.scatter(arcS,  Eimg,  color='black', s=2, zorder=zorder)
		zorder += 2

	# 4.4: Add a colourbar to the figure
	if len(images) > 1:
		no_of_ticks = 10
		ticklabels  = [int(index) for index in np.linspace(0, len(images)-1, no_of_ticks, endpoint=True)]
		ticks       = [index/(len(images)-1) for index in ticklabels]
		cb  = mpl.colorbar.ColorbarBase(ax2, cmap=mpl.colormaps['rainbow_r'], orientation='vertical', ticks=ticks)
		cb.set_ticklabels(ticklabels=ticklabels)

	# 4.5: Save all the NEB optimisation iterations together in one plot
	ax1.set_xlabel("Displacement (Å)", fontsize=11)
	ax1.set_ylabel("Energy (eV)", fontsize=11)
	ax1.set_title( "Iterations: %i to %i" % (0, no_of_iters-1) )
	ax2.set_ylabel("Iterations")
	plt.savefig('neb_optimization.png', dpi=400)

	# 4.6: Close the figure and clear matplotlib
	plt.close(fig)
	plt.clf()
	plt.cla()

	# 4.7: Save the energy profile of the last NEB iteration.
	plt.plot(  arcS2, Eimg2, color='tab:blue', linestyle='-', zorder=zorder)
	plt.scatter(arcS, Eimg,  color='black', s=2, zorder=zorder+1)
	plt.xlabel("Displacement (Å)", fontsize=11)
	plt.ylabel("Energy (eV)", fontsize=11)
	plt.savefig('neb_last_iteration.png', dpi=400)

# ================================================================================================

def get_interpolation_filename(path_to_images):
	"""
	This method is designed to locate the ".interp" file from the path_to_images folder.

	Parmeters
	---------
	path_to_images : str.
		This is the path to the ORCA NEB jobs. 
	"""

	# First, check that path_to_images exists.
	if not os.path.exists(path_to_images):
		raise Exception('Error: Can not find folder --> '+str(path_to_images))

	# Second, look in the path_to_images folder for any files that ends with the ".interp" suffix.
	interpolation_files = [file for file in os.listdir(path_to_images) if file.endswith('.interp')]

	# Third, exclude any interpolation file that contain final in their name.
	non_final_interpolation_files = [file for file in interpolation_files if not file.endswith('.final.interp')]

	# Fourth, if there are no interpolation files, indicate this is the case
	if len(non_final_interpolation_files) == 0:
		to_string  = f'Error: There are no interpolation files in folder --> {path_to_images}\n'
		to_string += f'Interpolation files are those that end with the ".interp" suffix.\n'
		to_string += 'Check this out.'
		raise Exception(to_string)

	# Fifth, make sure there is only one .interp file in the non_final_interpolation_files list.
	if not len(non_final_interpolation_files) == 1:
		to_string  = f'Error: There is more than 1 interpolation file in folder --> {path_to_images}\n'
		to_string += f'Choose one of these files and run it by including --filename ...\n'
		to_string += f'Interpolation files: {interpolation_files}\n'
		to_string += f'For example --> viewORCA neb_opt --filename {interpolation_files[0]}\n'
		to_string += 'Check this out.'
		raise Exception(to_string)

	# Sixth, return non_final_interpolation_files[0]
	return non_final_interpolation_files[0]

# ================================================================================================

def read_orca_interpolation_file(path_to_orca_interpolation_file):
	"""
	This method is designed to read the ORCA interpolation file and obtain energy and arc_length information from it.

	Parmeters
	---------
	path_to_orca_interpolation_file : str.
		This is the path to the ORCA interpolation file.
	"""

	# First, initialise the list for holding image and spline data.
	images = []
	spline = []

	# Second, initialise the temporary lists for recording the current iteration being recorded.
	current_image_arc_length  = []
	current_image_energy      = []
	current_spline_arc_length = []
	current_spline_energy     = []

	# Third, initialise booleans to indicate what is being recorded:
	#         * 0: Nothing
	#         * 1: Images
	#         * 2: Spline
	record_type  = 0

	# Fourth, open the interpolation file for the NEB.
	with open(path_to_orca_interpolation_file) as interpTXT:

		# Fifth, for each line in interpTXT.
		for line in interpTXT:

			# Sixth, determine what is being recorded, based on beginning of line.
			if   line.startswith('Iteration'):
				continue
			elif len(line.strip()) == 0:
				if not (record_type == 0):
					record_data(record_type, images, spline, current_image_arc_length, current_image_energy, current_spline_arc_length, current_spline_energy)
					record_type = 0
				continue
			elif line.startswith('Images:'):
				record_type = 1
				continue
			elif line.startswith('Interp.:'):
				record_type = 2
				continue

			# Seventh, separate the lines into there components
			image_number, arc_length, energy = line.rstrip().split()

			# Eighth, convert the arc_length and energy to floats
			arc_length = float(arc_length) * 0.529177249 # Bohr to angstrom 
			energy     = float(energy) *  27.2114079527 # Hartrees to electron-volts

			# Ninth, add the arc_length and energy to either the image or spline list.
			if   record_type == 1:
				current_image_arc_length.append(arc_length)
				current_image_energy    .append(energy)
			elif record_type == 2:
				current_spline_arc_length.append(arc_length)
				current_spline_energy    .append(energy)

	# Tenth, if you have got to the end of file but there is still data to record, record this.
	if record_type in [1, 2]:
		record_data(record_type, images, spline, current_image_arc_length, current_image_energy, current_spline_arc_length, current_spline_energy)
		record_type = 0

	# Tenth, return images and spline
	return images, spline

def record_data(record_type, images, spline, current_image_arc_length, current_image_energy, current_spline_arc_length, current_spline_energy):
	"""
	This method is designed to record the data to from the current lists to the images/spline list. 

	Parmeters
	---------
	record_type : int
		This integer indicates if the data is being recorded to the images list (record_type=1) or the spline list (record_type=2).
	images : list
		This list contains the arc length and energies for all the images for the NEB optimisation iterations.
	spline : list
		This list contains the arc length and energies for all the splines for the NEB optimisation iterations.
	current_image_arc_length : list
		This list contains the arc lengths of the current image being obtained.
	current_image_energy : list
		This list contains the energies of the current image being obtained.
	current_spline_arc_length : list
		This list contains the arc lengths of the current spline being obtained.
	current_spline_energy : list
		This list contains the energies of the current spline being obtained.
	"""

	# First, record iteration of image to main list.
	if record_type == 1:

		# 1.1: Check that current_image_arc_length and current_image_energy are the same length.
		if not (len(current_image_arc_length) == len(current_image_energy)):
			raise Exception('Error: no of arc_lengths and energies given are not the same')

		# 1.2: Append the current arc length and energy to the images list.
		images.append((tuple(current_image_arc_length), tuple(current_image_energy)))

		# 1.3: Reset the current_image_arc_length and current_image_energy lists.
		for _ in range(len(current_image_arc_length)):
			current_image_arc_length.pop()
		for _ in range(len(current_image_energy)):
			current_image_energy.pop()

	# Second, record iteration of spline to main list.
	if record_type == 2:

		# 2.1: Check that current_image_arc_length and current_image_energy are the same length.
		if not (len(current_spline_arc_length) == len(current_spline_energy)):
			raise Exception('Error: no of arc_lengths and energies given are not the same')

		# 2.2: Append the current arc length and energy to the spline list.
		spline.append((tuple(current_spline_arc_length), tuple(current_spline_energy)))

		# 2.3: Reset the current_spline_arc_length and current_spline_energy lists.
		for _ in range(len(current_spline_arc_length)):
			current_spline_arc_length.pop()
		for _ in range(len(current_spline_energy)):
			current_spline_energy.pop()

# ================================================================================================
