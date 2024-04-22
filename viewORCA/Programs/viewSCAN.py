import os
from ase.io import write
from ase.calculators.singlepoint import SinglePointCalculator
from ase.visualize import view

from viewORCA.Programs.viewSCAN_methods.utility_methods.does_contain_SCAN_files   import does_contain_SCAN_files
from viewORCA.Programs.viewSCAN_methods.utility_methods.does_contain_SCAN_folders import does_contain_SCAN_folders
from viewORCA.Programs.viewSCAN_methods.obtain_SCAN_information                   import obtain_SCAN_information
from viewORCA.Programs.viewSCAN_methods.obtain_SCAN_information_from_SCAN_folders import obtain_SCAN_information_from_SCAN_folders

class CLICommand:
	"""This module allows the user to view the output of an ORCA SCAN job visually. 
	"""

	@staticmethod
	def add_arguments(parser):
		parser.add_argument('-p', '--path', nargs=1, help='This is the path to the folder containing the ORCA SCAN job.', default=['.'])
		parser.add_argument('-v', '--view', nargs=1, help='If you want to view the ASE GUI for the SCAN job immediate after this program has run, set this to True.', default=['True'])

	@staticmethod
	def run(args):
		"""
		Run this program. This will allow the user to view a ORCA SCAN job. 

		Parameters
		----------
		args : argparse.Namespace
			This is the namespace that give all the information about the arguments.
		"""

		# First, obtain the local path to the ORCA SCAN job.
		path_to_images = args.path
		if len(path_to_images) != 1:
			raise Exception('Error: input "path" must only have one input: '+str(path_to_images))
		path_to_images = path_to_images[0]

		# Second, obtain the tag to indicate if the user want to view the ASE GUI of the ORCA SCAN job after running this program. 
		view_SCAN = args.view
		if len(view_SCAN) != 1:
			raise Exception('Error: input "view" must only have one input: '+str(view_SCAN))
		view_SCAN = view_SCAN[0]
		if   view_SCAN.lower() in ['true', 't']:
			view_SCAN = True
		elif view_SCAN.lower() in ['false', 'f']:
			view_SCAN = False
		else:
			raise Exception('Error: Input for viewOpt should be either True or False')

		# Third, run the method.
		Run_Method(path_to_images, view_SCAN)

def Run_Method(path_to_images, view_SCAN):
	"""
	This method is designed to allow the user to view the output of an ORCA SCAN job visually. 

	Parmeters
	---------
	path_to_images : str.
		This is the path to the ORCA SCAN jobs. 
	view_SCAN : bool.
		This indicates if the user wants to view the ORCA SCAN job after running this program.
	"""

	# First, determine if there are is a SCAN 'orca_trj.xyz' file.
	contains_SCAN_files = does_contain_SCAN_files(path_to_images)

	# Second, determine if there are SCAN folders to put SCAN data from.
	contains_SCAN_folders = does_contain_SCAN_folders(path_to_images)

	# Third, obtain the SCAN details.
	if   (not contains_SCAN_files) and (not contains_SCAN_folders):
		to_string  = f'Error: SCAN files were not found in path: {path_to_images}\n'
		to_string += f'Absolute path: {os.path.abspath(path_to_images)}\n'
		to_string += 'Check this'
		raise Exception(to_string)
	elif      contains_SCAN_files  and (not contains_SCAN_folders):
		SCAN_images, scan_images_energies = obtain_SCAN_information(path_to_images)
	elif (not contains_SCAN_files) and      contains_SCAN_folders :
		SCAN_images, scan_images_energies = obtain_SCAN_information_from_SCAN_folders(path_to_images)
	elif      contains_SCAN_files  and      contains_SCAN_folders :
		to_string  = f'Error: SCAN files and folders were both found in path: {path_to_images}\n'
		to_string += f'Either have only:\n'
		to_string += f'  -> name.XXX.xyz files and name_traj.xyz (SCAN files from a single SCAN job)\n'
		to_string += f'or\n'
		to_string += f'  -> SCAN_XXX folders (SCAN process that needed to be restarted several time. Scan files are found in each SCAN folder)\n'
		to_string += f'(where XXX is a number)\n'
		to_string += f'\n'
		to_string += f'Absolute path: {os.path.abspath(path_to_images)}\n'
		to_string += 'Check this'
		raise Exception(to_string)
	else:
		raise Exception('Error: Got here, how? Indicates a programming error.')

	# Third, assign energies to the appropriate SCAN image
	for index, energy in enumerate(scan_images_energies):
		SCAN_images[index].calc = SinglePointCalculator(atoms=SCAN_images[index], energy=energy)

	# Fourth, save the SCAN process as an XYZ file
	write('SCAN_images.xyz', SCAN_images)

	# Fifth, view the SCAN images with energies. 
	if view_SCAN:
		view(SCAN_images)

