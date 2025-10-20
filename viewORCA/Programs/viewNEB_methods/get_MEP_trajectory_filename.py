"""
get_MEP_trajectory_filename.py, Geoffrey Weal, 23/4/24

This method is designed to locate the "_MEP_trj.xyz" file from the path_to_ORCA_neb_folder folder.
"""
import os

def get_MEP_trajectory_filename(path_to_ORCA_neb_folder):
	"""
	This method is designed to locate the "_MEP_trj.xyz" file from the path_to_ORCA_neb_folder folder.

	Parmeters
	---------
	path_to_ORCA_neb_folder : str.
		This is the path to the ORCA geometric optimisation job. 
	"""

	# First, check that path_to_ORCA_neb_folder exists.
	if not os.path.exists(path_to_ORCA_neb_folder):
		raise Exception('Error: Can not find folder --> '+str(path_to_ORCA_neb_folder))

	# Second, look in the path_to_ORCA_neb_folder folder for any files that ends with the "_MEP_trj.xyz" suffix.
	NEB_trajectory_files = [file for file in os.listdir(path_to_ORCA_neb_folder) if file.endswith('_MEP_trj.xyz')]

	# Third, if there are no NEB trajectory files, indicate this is the case.
	if len(NEB_trajectory_files) == 0:
		to_string  = f'Error: There are no ORCA NEB trajectory files in folder --> {path_to_ORCA_neb_folder}\n'
		to_string += f'"viewORCA neb" requires a "_MEP_trj.xyz" file to read the NEB reaction pathway from.\n'
		to_string += f'ORCA NEB trajectory files are those that end with the "_MEP_trj.xyz" suffix.\n'
		to_string += 'Check this out.'
		raise Exception(to_string)

	# Fourth, make sure there is only one _MEP_trj.xyz file in the NEB_trajectory_files list.
	if not len(NEB_trajectory_files) == 1:
		to_string  = f'Error: There is more than 1 ORCA trajectory file in folder --> {path_to_ORCA_neb_folder}\n'
		to_string += f'"viewORCA neb" requires a "_MEP_trj.xyz" file to read the NEB reaction pathway from.\n'
		to_string += f'Choose one of these files and run it by including --path ...\n'
		to_string += f'ORCA trajectory files: {NEB_trajectory_files}\n'
		to_string += f'For example --> viewORCA opt --path {NEB_trajectory_files[0]}\n'
		to_string += 'Check this out.'
		raise Exception(to_string)

	# Sixth, return NEB_trajectory_files[0]
	return NEB_trajectory_files[0]