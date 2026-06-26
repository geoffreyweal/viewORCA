"""
read_orca_interpolation_file.py, Geoffrey Weal, 23/4/24

This method is designed to read the ORCA interpolation file and obtain energy and arc_length information from it.
"""

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

# ----------------------------------------------------------------------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------------------------------------------------------------------
