"""
obtain_SCAN_information.py, Geoffrey Weal, 18/4/24

This method is designed to obtain the images and energies of SCAN jobs.
"""
from viewORCA.Programs.viewSCAN_methods.obtain_SCAN_information_method.get_ORCA_SCAN_images           import get_ORCA_SCAN_images
from viewORCA.Programs.viewSCAN_methods.obtain_SCAN_information_method.assign_energies_to_SCAN_images import assign_energies_to_SCAN_images

def obtain_SCAN_information(path_to_images):
	"""
	This method is designed to obtain the images and energies of SCAN jobs.

	Parameters
	----------
	path_to_images : str.
		This is the path to the ORCA SCAN jobs. 

	Returns
	-------
	all_SCAN_images : list of ase.Atoms
		These are the images of each SCAN step.
	all_scan_images_energies : list of floats
		These are the energies of the images for each SCAN step.
	"""

	# First, obtain the SCAN images from the folder you are in
	SCAN_images = get_ORCA_SCAN_images(path_to_images)

	# Second, obtain the energies from the whole ORCA SCAN process (including geometry optimisations) for the SCAN images of interest,
	scan_images_energies = assign_energies_to_SCAN_images('orca_trj.xyz', scan_images=SCAN_images)

	# Third, return SCAN_images and scan_images_energies
	return SCAN_images, scan_images_energies