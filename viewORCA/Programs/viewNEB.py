"""
viewNEB.py, Geoffrey Weal, 23/4/24

This module allows the user to view the output of an ORCA Nudged Elastic Band (NEB) job visually.  
"""
import os

from ase.io                      import write
from ase.calculators.singlepoint import SinglePointCalculator
from ase.visualize               import view

from viewORCA.Programs.viewNEB_methods.get_MEP_trajectory_filename import get_MEP_trajectory_filename
from viewORCA.Programs.viewNEB_methods.get_ORCA_NEB_images         import get_ORCA_NEB_images

class CLICommand:
	"""This module allows the user to view the output of an ORCA Nudged Elastic Band (NEB) job visually. 
	"""

	@staticmethod
	def add_arguments(parser):
		parser.add_argument('-p', '--path', nargs=1, help='This is the path to the folder containing the ORCA NEB job. You can also give the direct path to the "_MEP_trj.xyz" file here (which is the file that is read to obtain information about the ORCA NEB process).', default=['.'])
		parser.add_argument('-v', '--view', nargs=1, help='If you want to view the ASE GUI for the NEB job immediate after this program has run, set this to True.', default=['True'])

	@staticmethod
	def run(args):
		"""
		Run this program. This will allow the user to view a ORCA NEB job. 

		Parameters
		----------
		args : argparse.Namespace
			This is the namespace that give all the information about the arguments.
		"""

		# First, obtain the local path to the ORCA NEB job.
		path_to_ORCA_neb_folder = args.path
		if len(path_to_ORCA_neb_folder) != 1:
			raise Exception('Error: input "path" must only have one input: '+str(path_to_ORCA_neb_folder))
		path_to_ORCA_neb_folder = path_to_ORCA_neb_folder[0]

		# Second, obtain the tag to indicate if the user want to view the ASE GUI of the ORCA NEB job after running this program. 
		view_NEB = args.view
		if len(view_NEB) != 1:
			raise Exception('Error: input "view" must only have one input: '+str(view_NEB))
		view_NEB = view_NEB[0]
		if   view_NEB.lower() in ['true', 't']:
			view_NEB = True
		elif view_NEB.lower() in ['false', 'f']:
			view_NEB = False
		else:
			raise Exception('Error: Input for viewOpt should be either True or False')

		# Third, run the method.
		Run_Method(path_to_ORCA_neb_folder, view_NEB)

def Run_Method(path_to_ORCA_neb_folder, view_NEB):
	"""
	This method is designed to allow the user to view the output of an ORCA NEB job visually. 

	Parmeters
	---------
	path_to_ORCA_neb_folder : str.
		This is the path to the ORCA NEB jobs. 
	view_NEB : bool.
		This indicates if the user wants to view the ORCA NEB job after running this program.
	"""

	# Prelimiary Step 1: Determine if path_to_ORCA_neb_folder is a folder or a file.
	if   not os.path.exists(path_to_ORCA_neb_folder):
		# Path does not exist.
		to_string  = 'Error: Path given does not exist.\n'
		to_string += f'Given path: {path_to_ORCA_neb_folder}\n'
		to_string += 'Check this.'
		raise Exception(to_string)
	elif os.path.isdir(path_to_ORCA_neb_folder):
		# Path is to a folder: Get the "_MEP_trj.xyz" file from path_to_ORCA_neb_folder
		filename = get_MEP_trajectory_filename(path_to_ORCA_neb_folder)
		path_to_MEP_trj_filepath = path_to_ORCA_neb_folder+'/'+filename
	elif os.path.isfile(path_to_ORCA_neb_folder):
		# Path is to a file: Check it is a trajectory file
		if not path_to_ORCA_neb_folder.endswith('_MEP_trj.xyz'):
			to_string  = 'Error: You file is not an ORCA trajectory file.\n'
			to_string += 'The file must end with the suffix "_MEP_trj.xyz".\n'
			to_string += f'Given filepath: {path_to_ORCA_neb_folder}\n'
			to_string += 'Check this.'
			raise Exception(to_string)
		path_to_MEP_trj_filepath = path_to_ORCA_neb_folder
	else:
		# Path is to neither a file or directory.
		to_string  = 'Error: Path given is to neither a file or directory.\n'
		to_string += f'Given filepath: {path_to_ORCA_neb_folder}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 2: Check that path_to_MEP_trj_filepath is a file.
	if not os.path.isfile(path_to_MEP_trj_filepath):
		to_string  = 'Error: Run_Method requires a file for the "path_to_MEP_trj_filepath" variable".\n'
		to_string += f'path_to_MEP_trj_filepath =  {path_to_MEP_trj_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 3: Check that path_to_MEP_trj_filepath ends with the '._MEP_trj.xyz' suffix
	#                    * This indicates it is a ORCA NEB trajectory file.
	if not path_to_MEP_trj_filepath.endswith('_MEP_trj.xyz'):
		to_string  = 'Error: You file is not an ORCA trajectory file.\n'
		to_string += 'The file must end with the suffix "_MEP_trj.xyz".\n'
		to_string += f'Given filepath: {path_to_MEP_trj_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 4: Get the directory name that the _MEP_trj.xyz file is in.
	file_path = os.path.dirname(path_to_MEP_trj_filepath)

	# Prelimiary Step 5: Obtain the filename of the _MEP_trj.xyz file.
	filename = os.path.basename(path_to_MEP_trj_filepath)

	# Prelimiary Step 6: Obtain the prefix of the filename of the _MEP_trj.xyz file.
	#                    * This should be the name of the ORCA input file.
	filename_prefix = filename.replace('_MEP_trj.xyz','')

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

	# First, obtain the images and the energies for the whole ORCA NEB process,
	neb_images, neb_images_energies = get_ORCA_NEB_images(file_path+'/'+filename_prefix+'_MEP_trj.xyz')

	# Second, assign energies to the appropriate NEB image.
	for index, (neb_image, neb_images_energy) in enumerate(zip(neb_images, neb_images_energies)):
		neb_images[index].calc = SinglePointCalculator(atoms=neb_images[index], energy=neb_images_energy)

	# Third, save the NEB reaction pathway as an .xyz file.
	write(file_path+'/'+filename_prefix+'_NEB_images.xyz', neb_images)

	# Fourth, view the SCAN images with energies. 
	if view_NEB:
		view(neb_images)

# ------------------------------------------------------------------------------
