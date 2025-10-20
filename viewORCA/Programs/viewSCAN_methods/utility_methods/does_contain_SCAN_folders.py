"""
does_contain_SCAN_folders.py, Geoffrey Weal, 18/4/24

This method is designed to determine if path_to_images contains SCAN folders.
"""
import os

def does_contain_SCAN_folders(path_to_images):
	"""
	This method determines if this path contains folders of multiple SCAN files. 

	Parmeters
	---------
	path_to_images : str.
		This is the path to the ORCA SCAN jobs. 

	Returns
	-------
	This will return True if SCAN files are found. False if not.
	"""

	# First, initialise a list to hold the names of xyz files. 
	SCAN_foldernames    = []
	SCAN_folder_numbers = []

	# Second, go though the files in path_to_images and obtain all the image xyz file for the SCAN Job.
	for foldername in os.listdir(path_to_images):

		# 2.1: If this is not a folder, move on.
		if not os.path.isdir(path_to_images+'/'+foldername):
			continue

		# 2.2: If the foldername does not start with "SCAN_", continue.
		if not foldername.startswith('SCAN_'):
			continue

		# 2.3: If the foldername does not have the format SCAN_XXX, where XXX is a number, move on.
		if not foldername.replace('SCAN_','').isdigit():
			continue

		# 2.4: Get the folder number for this SCAN job.
		SCAN_folder_number = int(foldername.replace('SCAN_',''))

		# 2.5: Add foldername to SCAN_foldernames.
		SCAN_foldernames.append(foldername)
		SCAN_folder_numbers.append(SCAN_folder_number)

	# Third, check there are no duplicate numbers in the set of xyz files.
	if not len(SCAN_folder_numbers) == len(set(SCAN_folder_numbers)):
		to_string  = f'Error: There seems to be duplicate SCAN numbers.\n'
		to_string += f'SCAN folders = {SCAN_foldernames}.\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# 1.4: Check that the xyz numbers are consecutive from 1 to len(xyz_numbers)+1
	if not sorted(SCAN_folder_numbers) == list(range(1,len(SCAN_folder_numbers)+1)):
		to_string  = f'Error: Expected that the SCAN folders to be labeled from 1 to {len(SCAN_folder_numbers)+1}.\n'
		to_string += f'However, there seems to be some missing SCAN folders.\n'
		missing_numbers = sorted(set(range(1,len(SCAN_folder_numbers)+1)) - set(SCAN_folder_numbers))
		to_string += f'Missing xyz numbers: {missing_numbers}.\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Second, indicate if there are SCAN folders.
	has_SCAN_folders = len(SCAN_foldernames) > 0

	# Third, return if there are SCAN folders.
	return has_SCAN_folders

# ------------------------------------------------------------------------------

