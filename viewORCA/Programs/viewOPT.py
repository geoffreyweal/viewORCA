"""
viewOpt.py, Geoffrey Weal, 23/4/24

This module allows the user to view the output of an ORCA geometric optimisation job visually. 
"""
import os

from ase.io                      import write
from ase.calculators.singlepoint import SinglePointCalculator
from ase.visualize               import view

from viewORCA.Programs.viewOPT_methods.get_trajectory_filename import get_trajectory_filename
from viewORCA.Programs.viewOPT_methods.get_ORCA_OPT_images     import get_ORCA_OPT_images

class CLICommand:
	"""This module allows the user to view the output of an ORCA geometric optimisation job visually. 
	"""

	@staticmethod
	def add_arguments(parser):
		parser.add_argument('-p', '--path', nargs=1, help='This is the path to the folder containing the ORCA geometric optimisation job. You can also give the direct path to the "_trj.xyz" file here (which is the file that is used to show the ORCA optimisation process).', default=['.'])
		parser.add_argument('-v', '--view', nargs=1, help='If you want to view the ASE GUI for the optimisation job immediate after this program has run, set this to True.', default=['True'])

	@staticmethod
	def run(args):
		"""
		Run this program. This will allow the user to view a ORCA geometric optimisation job. 

		Parameters
		----------
		args : argparse.Namespace
			This is the namespace that give all the information about the arguments.
		"""

		# First, obtain the local path to the ORCA geometric optimisation job.
		path_to_ORCA_opt_folder = args.path
		if len(path_to_ORCA_opt_folder) != 1:
			raise Exception('Error: input "path" must only have one input: '+str(path_to_ORCA_opt_folder))
		path_to_ORCA_opt_folder = path_to_ORCA_opt_folder[0]

		# Second, obtain the tag to indicate if the user want to view the ASE GUI of the ORCA geometric optimisation job after running this program. 
		view_Opt = args.view
		if len(view_Opt) != 1:
			raise Exception('Error: input "view" must only have one input: '+str(view_Opt))
		view_Opt = view_Opt[0]
		if   view_Opt.lower() in ['true', 't']:
			view_Opt = True
		elif view_Opt.lower() in ['false', 'f']:
			view_Opt = False
		else:
			raise Exception('Error: Input for "viewORCA opt" should be either True or False')

		# Third, run the method.
		Run_Method(path_to_ORCA_opt_folder, view_Opt)

def Run_Method(path_to_ORCA_opt_folder, view_Opt):
	"""
	This method is designed to allow the user to view the output of an ORCA geometric optimisation job visually. 

	Parmeters
	---------
	path_to_ORCA_opt_folder : str.
		This is the path to the direcotry containing the ORCA geometric optimisation calculation. 
	view_Opt : bool.
		This indicates if the user wants to view the ORCA geometric optimisation job after running this program.
	"""

	# Prelimiary Step 1: Determine if path_to_ORCA_opt_folder is a folder or a file.
	if   not os.path.exists(path_to_ORCA_opt_folder):
		# Path does not exist.
		to_string  = 'Error: Path given does not exist.\n'
		to_string += f'Given path: {path_to_ORCA_opt_folder}\n'
		to_string += 'Check this.'
		raise Exception(to_string)
	elif os.path.isdir(path_to_ORCA_opt_folder):
		# Path is to a folder: Get the "_trj.xyz" file from path_to_ORCA_opt_folder
		filename = get_trajectory_filename(path_to_ORCA_opt_folder)
		path_to_trj_filepath = path_to_ORCA_opt_folder+'/'+filename
	elif os.path.isfile(path_to_ORCA_opt_folder):
		# Path is to a file: Check it is a trajectory file
		if not path_to_ORCA_opt_folder.endswith('_trj.xyz'):
			to_string  = 'Error: You file is not an ORCA trajectory file.\n'
			to_string += 'The file must end with the suffix "_trj.xyz".\n'
			to_string += f'Given filepath: {path_to_ORCA_opt_folder}\n'
			to_string += 'Check this.'
			raise Exception(to_string)
		path_to_trj_filepath = path_to_ORCA_opt_folder
	else:
		# Path is to neither a file or directory.
		to_string  = 'Error: Path given is to neither a file or directory.\n'
		to_string += f'Given filepath: {path_to_ORCA_opt_folder}\n'
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

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

	# First, obtain the images and the energies for the whole ORCA NEB process,
	opt_images, opt_images_energies = get_ORCA_OPT_images(file_path+'/'+filename_prefix+'_trj.xyz')

	# Second, assign energies to the appropriate SCAN image
	for index, (neb_image, opt_images_energy) in enumerate(zip(opt_images, opt_images_energies)):
		opt_images[index].calc = SinglePointCalculator(atoms=opt_images[index], energy=opt_images_energy)

	# Third, save the SCAN process as an XYZ file
	write(file_path+'/'+filename_prefix+'_OPT_images.xyz', opt_images)

	# Fourth, view the SCAN images with energies. 
	if view_Opt:
		view(opt_images)

# ------------------------------------------------------------------------------
