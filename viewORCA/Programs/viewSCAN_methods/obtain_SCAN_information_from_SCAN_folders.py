"""
obtain_SCAN_information_from_SCAN_folders.py, Geoffrey Weal, 18/4/24

This method is designed to obtain the SCAN jobs from multiple folders and append them together.

This method is designed to be used in cases where a SCAN procedure needed to be restarted multiple times due to ORCA failing for some reason. The folders should be named SCAN_1, SCAN_2, SCAN_3, .... These folders will be pasted together as SCAN_1->SCAN_2->SCAN_3->...
"""
import os
from viewORCA.Programs.viewSCAN_methods.obtain_SCAN_information import obtain_SCAN_information

def obtain_SCAN_information_from_SCAN_folders(path_to_images):
	"""
	This method is designed to obtain the SCAN jobs from multiple folders and append them together.

	This method is designed to be used in cases where a SCAN procedure needed to be restarted multiple times due to ORCA failing for some reason. The folders should be named SCAN_1, SCAN_2, SCAN_3, .... These folders will be pasted together as SCAN_1->SCAN_2->SCAN_3->...
	
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

	# First, initialise a list to hold the names of xyz files. 
	SCAN_foldernames    = []

	# Second, go though the files in path_to_images and obtain all the image xyz file for the SCAN Job.
	for foldername in os.listdir(path_to_images):

		# 2.1: If this is not a folder, move on.
		if not os.path.isdir(path_to_images+'/'+foldername):
			continue

		# 2.2: If the foldername does not start with "SCAN_", continue.
		if foldername.startswith('SCAN_'):
			continue

		# 2.3: If the foldername does not have the format SCAN_XXX, where XXX is a number, move on.
		if foldername.replace('SCAN_','').is_digit():
			continue

		# 2.4: Add foldername to SCAN_foldernames.
		SCAN_foldernames.append(foldername)

	# Third, sort the SCAN folders from SCAN_1, SCAN_2, SCAN_3, ...
	SCAN_foldernames.sort(key=lambda foldername: int(foldername.replace('SCAN_','')))

	# Fourth, initialise list for storing the images and energies for each SCAN step. 
	all_SCAN_images = []
	all_scan_images_energies = []

	# Fifth, store a value that represents the energy to add to each image from each SCAN folder
	additional_energy = 0.0 # may not be needed

	# Sixth, for each SCAN folder, in order of SCAN_1, SCAN_2, SCAN_3, ...
	for SCAN_foldername in SCAN_foldernames:

		# 6.1: Obtain the number of the SCAN folder
		SCAN_folder_number = int(SCAN_foldername.replace('SCAN_',''))

		# 6.2: Obtain the images and the energies of the SCAN images in the SCAN folder. 
		SCAN_images, scan_images_energies = obtain_SCAN_information(path_to_images+'/'+SCAN_foldername)

		# 6.3: If we are looking at a restarted SCAN job, do the following: 
		if not SCAN_folder_number == 1:
			# We are looking at restarted SCAN jobs

			# 6.3.1: Remove the first image in SCAN_images, as it will be the same as the lst image in all_scan_images_energies
			del SCAN_images[0]

			# 6.3.2: Remove the a
			del scan_images_energies[0]

		# 6.4: Add the SCAN images to the all_SCAN_images folder. 
		all_SCAN_images += SCAN_images

		# 6.5: Add the SCAN energies to the all_scan_images_energies folder. 
		all_scan_images_energies += scan_images_energies

	# Seventh, return all_SCAN_images and all_scan_images_energies
	return all_SCAN_images, all_scan_images_energies

