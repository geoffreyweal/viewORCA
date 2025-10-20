"""
does_contain_SCAN_files.py, Geoffrey Weal, 18/4/24

This method is designed to determine if path_to_ORCA_scan_folder contains SCAN files.
"""
import os
from viewORCA.Programs.viewSCAN_methods.utility_methods.get_trajectory_filename import get_trajectory_filename

def does_contain_SCAN_files(path_to_ORCA_scan_folder):
	"""
	This method determines if this path contains the SCAN files.

	Parmeters
	---------
	path_to_ORCA_scan_folder : str.
		This is the path to the ORCA SCAN job. 

	Returns
	-------
	This will return True if SCAN files are found. False if not.
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
	# First, obtain the name of the trajectory file

	# 1.1: Obtain the filename of the _trj.xyz file. 
	have_orca_trj, trj_filename = get_trajectory_filename(path_to_ORCA_scan_folder)

	# 1.2: If the _trj.xyz file is not found, return False. 
	if not have_orca_trj:
		return False

	# 1.3: Obtain the prefix of the filename of the _trj.xyz file.
	#      * This should be the name of the ORCA input file.
	filename_prefix = trj_filename.replace('_trj.xyz','')

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
	# Second, determine if this folder contains filename_prefix.XXX.xyz files

	# 2.1: Initialise a list to hold the names of xyz files. 
	xyz_filenames  = []
	xyz_numbers    = []
	wrong_prefixes = []

	# 2.2: Go though the files in path_to_ORCA_scan_folder and obtain all the image xyz file for the SCAN Job.
	for filename in os.listdir(path_to_ORCA_scan_folder):

		# 2.2.1: Check that filename is a file.
		if not os.path.isfile(path_to_ORCA_scan_folder+'/'+filename):
			continue

		# 2.2.2: Check that the file is a xyz file.
		if not filename.endswith('.xyz'):
			continue

		# 2.2.3: If the second to last component is a three digit number (can start with zeros), this is a ORCA SCAN image.
		filename_components = filename.split('.')
		if (not filename_components[-2].isdigit()) or (not len(filename_components) >= 3):
			continue

		if not filename.startswith(filename_prefix):
			wrong_prefixes.append(filename)

		# 2.2.4: Append filename_components to xyz_filenames
		xyz_filenames.append(filename)
		xyz_numbers.append(int(filename_components[-2]))

	# 2.3: Make sure that there are no .XXX.xyz file that have a different prefix to filename.
	if len(wrong_prefixes) > 0:
		to_string  = f'Error: There are some .XXX.xyz files (where XXX is a number) that have a different prefix compared to the _trj.xyz file.\n'
		to_string += f'Expected Prefix = {filename_prefix}.\n'
		to_string += f'trj.xyz filename = {trj_filename}.\n'
		to_string += f'.XXX.xyz files with different prefix = {sorted(wrong_prefixes)}.\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# 2.4: Check there are no duplicate numbers in the set of xyz files.
	if not len(xyz_numbers) == len(set(xyz_numbers)):
		to_string  = f'Error: There seems to be duplicate SCAN numbers.\n'
		to_string += f'xyz_filenames = {xyz_filenames}.\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# 2.5: Check that the xyz numbers are consecutive from 1 to len(xyz_numbers)+1
	if not sorted(xyz_numbers) == list(range(1,len(xyz_numbers)+1)):
		to_string  = f'Error: Expected that the xyz files from the SCAN calculation to be labeled from 1 to {len(xyz_numbers)+1}.\n'
		to_string += f'However, there seems to be some missing xyz files.\n'
		missing_numbers = sorted(set(range(1,len(xyz_numbers)+1)) - set(xyz_numbers))
		to_string += f'Missing xyz numbers: {missing_numbers}.\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Determine if there are any files that include the suffix "_trj.xyz".  

	# Third, indicate if there are SCAN xyz files.
	have_xyz_files = len(xyz_filenames) > 0

	# Fourth, return if the files needed to check SCAN calculations are found
	return have_xyz_files

# ------------------------------------------------------------------------------

