"""
viewNEB_interpolation.py, Geoffrey Weal, 19/5/24

This module allows the user to check the interpolation pathway before running an NEB in ORCA. 
"""
import os

from ase.io        import read, write
from ase.neb       import NEB
from ase.visualize import view

class CLICommand:
	"""This module allows the user to check the interpolation pathway before running an NEB in ORCA. 
	"""

	@staticmethod
	def add_arguments(parser):
		parser.add_argument('reactants_filepath',      nargs=1, help='This is the path to the xyz file holding the reactants.') # ' (Default: reactants.xyz)', default=['reactants.xyz'])
		parser.add_argument('products_filepath',       nargs=1, help='This is the path to the xyz file holding the products.') # ' (Default: products.xyz)',  default=['products.xyz'])
		parser.add_argument('no_of_steps',             nargs=1, help='This is the number of steps you would like to use in your NEB calculation.')
		parser.add_argument('-i', '--use_idpp_method', nargs=1, help='This tag indicates if you want to use the IDPP method (Default: True).', default=['True'])
		parser.add_argument('-v', '--view',            nargs=1, help='If you want to view the ASE GUI for the NEB job immediate after this program has run, set this to True (Default: True).', default=['True'])

	@staticmethod
	def run(args):
		"""
		Run this program. This module allows the user to check the interpolation pathway before running an NEB in ORCA. 

		Parameters
		----------
		args : argparse.Namespace
			This is the namespace that give all the information about the arguments.
		"""

		# First, obtain the path to the reactants file. 
		reactants_filepath = args.reactants_filepath
		if len(reactants_filepath) != 1:
			raise Exception('Error: input "reactants_filepath" must only have one input: '+str(reactants_filepath))
		reactants_filepath = reactants_filepath[0]

		# Second, obtain the path to the products file. 
		products_filepath = args.products_filepath
		if len(products_filepath) != 1:
			raise Exception('Error: input "products_filepath" must only have one input: '+str(products_filepath))
		products_filepath = products_filepath[0]

		# Third, obtain the number of steps you want to perform in your NEB calculation. 
		no_of_steps = args.no_of_steps
		if len(no_of_steps) != 1:
			raise Exception('Error: input "no_of_steps" must only have one input: '+str(no_of_steps))
		no_of_steps = no_of_steps[0]

		# First, determine if you want to use the IDPP interpolation method
		use_idpp_method = args.use_idpp_method
		if len(use_idpp_method) != 1:
			raise Exception('Error: input "use_idpp_method" must only have one input: '+str(use_idpp_method))
		use_idpp_method = use_idpp_method[0]
		if   use_idpp_method.lower() in ['true', 't']:
			use_idpp_method = True
		elif use_idpp_method.lower() in ['false', 'f']:
			use_idpp_method = False
		else:
			raise Exception('Error: Input for viewOpt should be either True or False')

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
		Run_Method(reactants_filepath, products_filepath, no_of_steps, use_idpp_method, view_NEB)

def Run_Method(reactants_filepath, products_filepath, no_of_steps, use_idpp_method, view_NEB):
	"""
	This module allows the user to check the interpolation pathway before running an NEB in ORCA. 

	Parmeters
	---------
	reactants_filepath : str.
		This is the path to the xyz file holding the reactants.
	products_filepath : str.
		This is the path to the xyz file holding the products.
	no_of_steps : str.
		This is the number of steps you would like to use in your NEB calculation.
	use_idpp_method : bool.
		This tag indicates if you want to use the IDPP method.
	view_NEB : bool.
		This indicates if the user wants to view the ORCA NEB job after running this program.
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Prelimiary Step 1: Determine if reactants_filepath is a folder or a file.
	if   not os.path.exists(reactants_filepath):
		# Path does not exist.
		to_string  = 'Error: Reactants file given does not exist.\n'
		to_string += f'Given reactants path: {reactants_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)
	elif os.path.isdir(reactants_filepath):
		# Path is a folder. 
		to_string  = 'Error: Reactants file given is a folder, not a file.\n'
		to_string += f'Given reactants path: {reactants_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 2: Check that reactants_filepath is a file.
	if not os.path.isfile(reactants_filepath):
		to_string  = 'Error: Run_Method requires a file for the "reactants_filepath" variable".\n'
		to_string += f'reactants_filepath =  {reactants_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 3: Check that reactants_filepath ends with the '.xyz' suffix
	#                    * This indicates it is a ORCA NEB trajectory file.
	if not reactants_filepath.endswith('.xyz'):
		to_string  = 'Error: You reactants file is not an xyz.\n'
		to_string += 'The file must end with the suffix ".xyz".\n'
		to_string += f'Given filepath: {reactants_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Prelimiary Step 4: Determine if products_filepath is a folder or a file.
	if   not os.path.exists(products_filepath):
		# Path does not exist.
		to_string  = 'Error: Products file given does not exist.\n'
		to_string += f'Given products path: {products_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)
	elif os.path.isdir(products_filepath):
		# Path is a folder. 
		to_string  = 'Error: Products file given is a folder, not a file.\n'
		to_string += f'Given products path: {products_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 5: Check that products_filepath is a file.
	if not os.path.isfile(products_filepath):
		to_string  = 'Error: Run_Method requires a file for the "products_filepath" variable".\n'
		to_string += f'products_filepath =  {products_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 6: Check that products_filepath ends with the '.xyz' suffix
	#                    * This indicates it is a ORCA NEB trajectory file.
	if not products_filepath.endswith('.xyz'):
		to_string  = 'Error: You products file is not an xyz.\n'
		to_string += 'The file must end with the suffix ".xyz".\n'
		to_string += f'Given filepath: {products_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Prelimiary Step 7: Make sure that no_of_steps is an integer.
	if not no_of_steps.isdigit():
		to_string  = 'Error: no_of_steps needs to be an integer.\n'
		to_string += f'no_of_steps given: {no_of_steps}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Prelimiary Step 8: Get the directory name that the _MEP_trj.xyz file is in.
	file_path = os.path.dirname(reactants_filepath)
	if file_path == '':
		file_path = '.'

	# Prelimiary Step 9: Obtain the filename of the _MEP_trj.xyz file.
	filename = 'interpolation.xyz'

	# Prelimiary Step 10: Obtain the prefix of the filename of the _MEP_trj.xyz file.
	#                     * This should be the name of the ORCA input file.
	filename_prefix = filename.replace('.xyz','')

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

	# First, read in the reactants file.
	reactants = read(reactants_filepath)

	# Second, read in the products file.
	products = read(products_filepath)

	# Third, convert the number of steps into a int (it will be a string current when read from the terminal).
	no_of_steps = int(no_of_steps)

	# Fourth, create all the images for the NEB.
	interpolation_images = [reactants] + [reactants.copy() for i in range(no_of_steps)] +[products]

	# Fifth, create the NEB object
	neb = NEB(interpolation_images)

	# Sixth, interpolate the images
	if use_idpp_method:
		neb.interpolate(method='idpp')
	else:
		neb.interpolate()

	# Seventh, save the NEB interpolation pathway as an .xyz file.
	write(file_path+'/NEB_interpolation_check.xyz', interpolation_images)

	# Eighth, view the interpolation images. 
	if view_NEB:
		view(interpolation_images)

# ------------------------------------------------------------------------------
