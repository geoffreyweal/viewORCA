"""
viewIRC.py, Geoffrey Weal, 23/4/24

This module allows the user to view the output of an ORCA Intrinsic Reaction Coordinate (IRC) job visually. 
"""
import os

from ase.io                      import write
from ase.calculators.singlepoint import SinglePointCalculator
from ase.visualize               import view

from viewORCA.Programs.viewIRC_methods.get_trajectory_filename import get_trajectory_filename
from viewORCA.Programs.viewIRC_methods.get_ORCA_IRC_images     import get_ORCA_IRC_images

class CLICommand:
	"""This module allows the user to view the output of an ORCA Intrinsic Reaction Coordinate (IRC) job visually. 
	"""

	@staticmethod
	def add_arguments(parser):
		parser.add_argument('-p', '--path', nargs=1, help='This is the path to the folder containing the ORCA IRC job. You can also give the direct path to the "_IRC_Full_trj.xyz" file here (which is the file that is used to show the ORCA IRC process).', default=['.'])
		parser.add_argument('-v', '--view', nargs=1, help='If you want to view the ASE GUI for the IRC job immediate after this program has run, set this to True.', default=['True'])

	@staticmethod
	def run(args):
		"""
		Run this program. This will allow the user to view a ORCA IRC job. 

		Parameters
		----------
		args : argparse.Namespace
			This is the namespace that give all the information about the arguments.
		"""

		# First, obtain the local path to the ORCA IRC job.
		path_to_ORCA_irc_folder = args.path
		if len(path_to_ORCA_irc_folder) != 1:
			raise Exception('Error: input "path" must only have one input: '+str(path_to_ORCA_irc_folder))
		path_to_ORCA_irc_folder = path_to_ORCA_irc_folder[0]

		# Second, obtain the tag to indicate if the user want to view the ASE GUI of the ORCA IRC job after running this program. 
		view_IRC = args.view
		if len(view_IRC) != 1:
			raise Exception('Error: input "view" must only have one input: '+str(view_IRC))
		view_IRC = view_IRC[0]
		if   view_IRC.lower() in ['true', 't']:
			view_IRC = True
		elif view_IRC.lower() in ['false', 'f']:
			view_IRC = False
		else:
			raise Exception('Error: Input for "viewORCA irc" should be either True or False')

		# Third, run the method.
		Run_Method(path_to_ORCA_irc_folder, view_IRC)

def Run_Method(path_to_ORCA_irc_folder, view_IRC):
	"""
	This method is designed to allow the user to view the output of an ORCA geometric optimisation job visually. 

	Parmeters
	---------
	path_to_ORCA_irc_folder : str.
		This is the path to the ORCA geometric optimisation jobs. 
	view_IRC : bool.
		This indicates if the user wants to view the ORCA IRC job after running this program.
	"""

	# Prelimiary Step 1: Determine if path_to_ORCA_irc_folder is a folder or a file.
	if   not os.path.exists(path_to_ORCA_irc_folder):
		# Path does not exist.
		to_string  = 'Error: Path given does not exist.\n'
		to_string += f'Given path: {path_to_ORCA_irc_folder}\n'
		to_string += 'Check this.'
		raise Exception(to_string)
	elif os.path.isdir(path_to_ORCA_irc_folder):
		# Path is to a folder: Get the "_IRC_Full_trj.xyz" file from path_to_ORCA_irc_folder
		filename = get_trajectory_filename(path_to_ORCA_irc_folder)
		path_to_IRC_Full_trj_filepath = path_to_ORCA_irc_folder+'/'+filename
	elif os.path.isfile(path_to_ORCA_irc_folder):
		# Path is to a file: Check it is a trajectory file
		if not path_to_ORCA_irc_folder.endswith('_IRC_Full_trj.xyz'):
			to_string  = 'Error: You file is not an ORCA IRC trajectory file.\n'
			to_string += 'The file must end with the suffix "_IRC_Full_trj.xyz".\n'
			to_string += f'Given filepath: {path_to_ORCA_irc_folder}\n'
			to_string += 'Check this.'
			raise Exception(to_string)
		path_to_IRC_Full_trj_filepath = path_to_ORCA_irc_folder
	else:
		# Path is to neither a file or directory.
		to_string  = 'Error: Path given is to neither a file or directory.\n'
		to_string += f'Given filepath: {path_to_ORCA_irc_folder}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 2: Check that path_to_IRC_Full_trj_filepath is a file.
	if not os.path.isfile(path_to_IRC_Full_trj_filepath):
		to_string  = 'Error: Run_Method requires a file for the "path_to_IRC_Full_trj_filepath" variable".\n'
		to_string += f'path_to_IRC_Full_trj_filepath =  {path_to_IRC_Full_trj_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 3: Check that path_to_IRC_Full_trj_filepath ends with the '.interp' suffix
	#                    * This indicates it is a ORCA interpolation file.
	if not path_to_IRC_Full_trj_filepath.endswith('_IRC_Full_trj.xyz'):
		to_string  = 'Error: You file is not an ORCA IRC trajectory file.\n'
		to_string += 'The file must end with the suffix "_IRC_Full_trj.xyz".\n'
		to_string += f'Given filepath: {path_to_IRC_Full_trj_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 4: Get the directory name that _IRC_Full_trj.xyz file is in.
	file_path = os.path.dirname(path_to_IRC_Full_trj_filepath)

	# Prelimiary Step 5: Obtain the filename of the _IRC_Full_trj.xyz file.
	filename = os.path.basename(path_to_IRC_Full_trj_filepath)

	# Prelimiary Step 6: Obtain the prefix of the filename of the _IRC_Full_trj.xyz file.
	#                    * This should be the name of the ORCA input file.
	filename_prefix = filename.replace('_IRC_Full_trj.xyz','')

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

	# First, obtain the images and the energies for the whole ORCA NEB process,
	irc_images, irc_images_energies = get_ORCA_IRC_images(file_path+'/'+filename_prefix+'_IRC_Full_trj.xyz')

	# Second, assign energies to the appropriate SCAN image
	for index, (irc_image, irc_images_energy) in enumerate(zip(irc_images, irc_images_energies)):
		irc_images[index].calc = SinglePointCalculator(atoms=irc_images[index], energy=irc_images_energy)

	# Third, save the SCAN process as an XYZ file
	write(file_path+'/'+filename_prefix+'_IRC_images.xyz', irc_images)

	# Fourth, view the SCAN images with energies. 
	if view_IRC:
		view(irc_images)

# ------------------------------------------------------------------------------
