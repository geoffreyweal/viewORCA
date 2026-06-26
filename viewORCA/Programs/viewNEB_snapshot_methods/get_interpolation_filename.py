"""
get_interpolation_filename.py, Geoffrey Weal, 23/4/24

This method is designed to locate the ".interp" file from the path_to_NEB_job folder.
"""
import os

def get_interpolation_filename(path_to_NEB_job):
	"""
	This method is designed to locate the ".interp" file from the path_to_NEB_job folder.

	Parmeters
	---------
	path_to_NEB_job : str.
		This is the path to the ORCA NEB job. 
	"""

	# First, check that path_to_NEB_job exists.
	if not os.path.exists(path_to_NEB_job):
		raise Exception('Error: Can not find folder --> '+str(path_to_NEB_job))

	# Second, look in the path_to_NEB_job folder for any files that ends with the ".interp" suffix.
	interpolation_files = [file for file in os.listdir(path_to_NEB_job) if file.endswith('.interp')]

	# Third, exclude any interpolation file that contain final in their name.
	non_final_interpolation_files = [file for file in interpolation_files if not file.endswith('.final.interp')]

	# Fourth, if there are no interpolation files, indicate this is the case
	if len(non_final_interpolation_files) == 0:
		to_string  = f'Error: There are no interpolation files in folder --> {path_to_NEB_job}\n'
		to_string += f'"viewORCA neb_snap" requires a ".interp" file to read the NEB optimisation process from.\n'
		to_string += f'Interpolation files are those that end with the ".interp" suffix.\n'
		to_string += 'Check this out.'
		raise Exception(to_string)

	# Fifth, make sure there is only one .interp file in the non_final_interpolation_files list.
	if not len(non_final_interpolation_files) == 1:
		to_string  = f'Error: There is more than 1 interpolation file in folder --> {path_to_NEB_job}\n'
		to_string += f'"viewORCA neb_snap" requires a ".interp" file to read the NEB optimisation process from.\n'
		to_string += f'Choose one of these files and run it by including --path ...\n'
		to_string += f'Interpolation files: {interpolation_files}\n'
		to_string += f'For example --> viewORCA neb_snap --path {interpolation_files[0]}\n'
		to_string += 'Check this out.'
		raise Exception(to_string)

	# Sixth, return non_final_interpolation_files[0]
	return non_final_interpolation_files[0]