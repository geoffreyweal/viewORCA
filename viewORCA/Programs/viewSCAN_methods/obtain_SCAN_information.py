"""
obtain_SCAN_information.py, Geoffrey Weal, 18/4/24

This method is designed to obtain the images and energies of SCAN jobs.
"""
import os
from viewORCA.Programs.viewOPT_methods.get_trajectory_filename                                        import get_trajectory_filename
from viewORCA.Programs.viewSCAN_methods.obtain_SCAN_information_method.get_ORCA_SCAN_images           import get_ORCA_SCAN_images
from viewORCA.Programs.viewSCAN_methods.obtain_SCAN_information_method.assign_energies_to_SCAN_images import assign_energies_to_SCAN_images

def obtain_SCAN_information(path_to_ORCA_scan_folder):
	"""
	This method is designed to obtain the images and energies of SCAN jobs.

	Parameters
	----------
	path_to_ORCA_scan_folder : str.
		This is the path to the ORCA SCAN jobs. 

	Returns
	-------
	all_SCAN_images : list of ase.Atoms
		These are the images of each SCAN step.
	all_scan_images_energies : list of floats
		These are the energies of the images for each SCAN step.
	"""

	# Prelimiary Step 1: Determine if path_to_ORCA_scan_folder is a folder or a file.
	if   not os.path.exists(path_to_ORCA_scan_folder):
		# Path does not exist.
		to_string  = 'Error: Path given does not exist.\n'
		to_string += f'Given path: {path_to_ORCA_scan_folder}\n'
		to_string += 'Check this.'
		raise Exception(to_string)
	elif os.path.isdir(path_to_ORCA_scan_folder):
		# Path is to a folder: Get the "_trj.xyz" file from path_to_ORCA_scan_folder
		filename = get_trajectory_filename(path_to_ORCA_scan_folder)
		path_to_trj_filepath = path_to_ORCA_scan_folder+'/'+filename
	elif os.path.isfile(path_to_ORCA_scan_folder):
		# Path is to a file: Check it is a trajectory file
		if not path_to_ORCA_scan_folder.endswith('_trj.xyz'):
			to_string  = 'Error: You file is not an ORCA trajectory file.\n'
			to_string += 'The file must end with the suffix "_trj.xyz".\n'
			to_string += f'Given filepath: {path_to_ORCA_scan_folder}\n'
			to_string += 'Check this.'
			raise Exception(to_string)
		path_to_trj_filepath = path_to_ORCA_scan_folder
	else:
		# Path is to neither a file or directory.
		to_string  = 'Error: Path given is to neither a file or directory.\n'
		to_string += f'Given filepath: {path_to_ORCA_scan_folder}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 2: Check that path_to_trj_filepath is a file.
	if not os.path.isfile(path_to_trj_filepath):
		to_string  = 'Error: Run_Method requires a file for the "path_to_trj_filepath" variable".\n'
		to_string += f'path_to_trj_filepath =  {path_to_trj_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 3: Check that path_to_trj_filepath ends with the '.interp' suffix
	#                    * This indicates it is a ORCA interpolation file.
	if not path_to_trj_filepath.endswith('_trj.xyz'):
		to_string  = 'Error: You file is not an ORCA trajectory file.\n'
		to_string += 'The file must end with the suffix "_trj.xyz".\n'
		to_string += f'Given filepath: {path_to_trj_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 4: Get the directory name that the _trj.xyz file is in.
	file_path = os.path.dirname(path_to_trj_filepath)

	# Prelimiary Step 5: Obtain the filename of the _trj.xyz file.
	filename = os.path.basename(path_to_trj_filepath)

	# Prelimiary Step 6: Obtain the prefix of the filename of the _trj.xyz file.
	#                    * This should be the name of the ORCA input file.
	filename_prefix = filename.replace('_trj.xyz','')

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# First, obtain the SCAN images from the folder you are in
	SCAN_images = get_ORCA_SCAN_images(path_to_ORCA_scan_folder)

	# Second, obtain the energies from the whole ORCA SCAN process (including geometry optimisations) for the SCAN images of interest,
	scan_images_energies = assign_energies_to_SCAN_images(file_path+'/'+filename_prefix+'_trj.xyz', scan_images=SCAN_images)

	# Third, return SCAN_images and scan_images_energies
	return SCAN_images, scan_images_energies

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
