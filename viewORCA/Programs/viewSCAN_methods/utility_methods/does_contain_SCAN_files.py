"""
does_contain_SCAN_files.py, Geoffrey Weal, 18/4/24

This method is designed to determine if path_to_images contains SCAN files.
"""
import os

def does_contain_SCAN_files(path_to_images):
	"""
	This method determines if this path contains the SCAN files.

	Parmeters
	---------
	path_to_images : str.
		This is the path to the ORCA SCAN jobs. 

	Returns
	-------
	This will return True if SCAN files are found. False if not.
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
	# First, determine if this folder contains orca.XXX.xyz files

	# 1.1: Initialise a list to hold the names of xyz files. 
	xyz_filenames = []
	xyz_numbers   = []

	# 1.2: Go though the files in path_to_images and obtain all the image xyz file for the SCAN Job.
	for filename in os.listdir(path_to_images):

		# 1.2.1: Check that filename is a file.
		if not os.path.isfile(path_to_images+'/'+filename):
			continue

		# 1.2.2: Check that the file is a xyz file.
		if not filename.endswith('.xyz'):
			continue

		# 1.2.3: If the second to last component is a three digit number (can start with zeros), this is a ORCA SCAN image.
		filename_components = filename.split('.')
		if (not filename_components[-2].isdigit()) or (not len(filename_components) >= 3):
			continue

		# 1.2.4: Append filename_components to xyz_filenames
		xyz_filenames.append(filename)
		xyz_numbers.append(int(filename_components[-2]))

	# 1.3: Check there are no duplicate numbers in the set of xyz files.
	if not len(xyz_numbers) == len(set(xyz_numbers)):
		to_string  = f'Error: There seems to be duplicate SCAN numbers.\n'
		to_string += f'xyz_filenames = {xyz_filenames}.\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# 1.4: Check that the xyz numbers are consecutive from 1 to len(xyz_numbers)+1
	if not sorted(xyz_numbers) == list(range(1,len(xyz_numbers)+1)):
		to_string  = f'Error: Expected that the xyz files from the SCAN calculation to be labeled from 1 to {len(xyz_numbers)+1}.\n'
		to_string += f'However, there seems to be some missing xyz files.\n'
		missing_numbers = sorted(set(range(1,len(xyz_numbers)+1)) - set(xyz_numbers))
		to_string += f'Missing xyz numbers: {missing_numbers}.\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	# Second, indicate if there are SCAN xyz files.
	have_xyz_files = len(xyz_filenames) > 0

	# Third, determine if the 'orca_trj.xyz' is found in path_to_images.
	have_orca_trj = os.path.exists('orca_trj.xyz')

	# Fourth, return if the files needed to check SCAN calculations are found
	return have_xyz_files and have_orca_trj

# ------------------------------------------------------------------------------

