"""
get_ORCA_SCAN_images.py, Geoffrey Weal, 18/4/24

This method is designed to obtain all the images for the SCAN process.
"""
import os
from ase.io import read

def get_ORCA_SCAN_images(path_to_images):
	"""
	This method is designed to obtain all the images for the SCAN process.

	Parameters
	----------
	path_to_images : str.
		This is the path to the folder containing the SCAN image files. 

	Returns
	-------
	Return the images of the SCAN
	"""

	# First, initialise a list to hold the SCAN images
	SCAN_images = []

	# Second, go though the files in path_to_images and obtain all the image xyz file for the SCAN Job.
	for filename in os.listdir(path_to_images):

		# 2.1: Check that filename is a file
		if not os.path.isfile(path_to_images+'/'+filename):
			continue

		# 2.2: Check that the file is a xyz file
		if not filename.endswith('.xyz'):
			continue

		# 2.3: If the second to last component is a three digit number (can start with zeros), this is a ORCA SCAN image.
		filename_components = filename.split('.')
		if (not filename_components[-2].isdigit()) or (not len(filename_components) >= 3):
			continue

		# 2.4: Obtain the component number as an integer
		component_number = int(filename_components[-2])

		# 2.5: Read the ORCA SCAN image xyz file
		system = read(path_to_images+'/'+filename)

		# 2.6: Append (component_number, system) to SCAN_images
		SCAN_images.append( (component_number, system) )

	# Third, sort SCAN_images by image number
	SCAN_images.sort()

	# Fourth, check that the component numbers in SCAN_images are consecutive
	component_numbers = [component_number for (component_number, system) in SCAN_images]
	if not (component_numbers == list(range(1,len(component_numbers)+1))):
		raise Exception('Error')

	# Fifth, return SCAN_images, but just the ase.Atoms components
	return [system for (component_number, system) in SCAN_images]

# ------------------------------------------------------------------------------
