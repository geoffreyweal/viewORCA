"""
get_trajectory_filename.py, Geoffrey Weal, 23/4/24

This method is designed to locate the "_trj.xyz" file from the path_to_ORCA_scan_folder folder.
"""
import os

def get_trajectory_filename(path_to_ORCA_scan_folder):
	"""
	This method is designed to locate the "_trj.xyz" file from the path_to_ORCA_scan_folder folder.

	Parmeters
	---------
	path_to_ORCA_scan_folder : str.
		This is the path to the ORCA scan job. 
	"""

	# First, check that path_to_ORCA_scan_folder exists.
	if not os.path.exists(path_to_ORCA_scan_folder):
		#raise Exception('Error: Can not find folder --> '+str(path_to_ORCA_scan_folder))
		return False, None

	# Second, look in the path_to_ORCA_scan_folder folder for any files that ends with the "_trj.xyz" suffix.
	trajectory_files = [file for file in os.listdir(path_to_ORCA_scan_folder) if file.endswith('_trj.xyz')]

	# Third, if there are no trajectory files, indicate this is the case.
	if len(trajectory_files) == 0:
		'''
		to_string  = f'Error: There are no ORCA trajectory files in folder --> {path_to_ORCA_scan_folder}\n'
		to_string += f'"viewORCA scan" requires a "_trj.xyz" file to read the SCAN calculation from.\n'
		to_string += f'ORCA trajectory files are those that end with the "_trj.xyz" suffix.\n'
		to_string += 'Check this out.'
		raise Exception(to_string)
		'''
		return False, None

	# Fourth, make sure there is only one .interp file in the trajectory_files list.
	if not len(trajectory_files) == 1:
		to_string  = f'Error: There is more than 1 ORCA trajectory file in folder --> {path_to_ORCA_scan_folder}\n'
		to_string += f'"viewORCA scan" requires a "_trj.xyz" file to read the SCAN calculation from.\n'
		to_string += f'Choose one of these files and run it by including --trj_name ...\n'
		to_string += f'ORCA trajectory files: {trajectory_files}\n'
		to_string += f'For example --> viewORCA scan --trj_name {trajectory_files[0]}\n'
		to_string += 'Check this out.'
		raise Exception(to_string)

	# Sixth, return trajectory_files[0]
	return True, trajectory_files[0]