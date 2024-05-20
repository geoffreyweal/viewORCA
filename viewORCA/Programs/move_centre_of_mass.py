"""
move_centre_of_mass.py, Geoffrey Weal, 14/5/24

This module allows the user to move the centre of mass of a system. 
"""
import os
import numpy as np
from ase.io import read, write

class CLICommand:
	"""This module allows the user to move the centre of mass of a system. 
	"""

	@staticmethod
	def add_arguments(parser):
		parser.add_argument('filepath',                       nargs=1,   help='This is the path to the xyz file you want to move the centre of mass of.')
		parser.add_argument('-c', '--move_centre_of_mass_to', nargs=1, help='This is the (x,y,z) position that you want to move the centre of mass of the system to. If "centre" is given, centre the chemical system. Default=centre', default=['centre'])

	@staticmethod
	def run(args):
		"""
		Run this program. This method will allow the user to obtain the number of electrons in the system of interest, and therefore obtain the possible multiplicities that the system could be in.

		Parameters
		----------
		args : argparse.Namespace
			This is the namespace that give all the information about the arguments.
		"""

		# First, obtain the local filepath for the molecule you want to obtain the number of electrons for
		filepath = args.filepath
		if len(filepath) != 1:
			raise Exception('Error: input "filepath" must only have one input: '+str(filepath))
		filepath = filepath[0]

		# Second, obtain the charge of the system
		move_centre_of_mass_to = args.move_centre_of_mass_to
		if (len(move_centre_of_mass_to) == 1) and (move_centre_of_mass_to[0] in ['centre', 'center']):
			move_centre_of_mass_to = 'centre'
		elif (len(move_centre_of_mass_to) != 1):
			toString  = f'Error, you must give three integers/float values for the centre of mass. These represent the (x,y,z) position that you want to move the centre of mass of your chemical system too.'
			toString += f'Input you gave for move_centre_of_mass_to = {move_centre_of_mass_to}'
			exit(toString)
		else:
			move_centre_of_mass_to = move_centre_of_mass_to[0]

		# Third, run the method.
		Run_Method(filepath, move_centre_of_mass_to)

def Run_Method(filepath, move_centre_of_mass_to='centre'):
	"""
	This module allows the user to move the centre of mass of a system. 

	Parmeters
	---------
	filepath : str.
		This is the path to the system file (preferably an xyz file). 
	move_centre_of_mass_to : str. or list/tuple
		This is the (x,y,z) position that you want to move the centre of mass of the system to. If 'centre' is given, centre the chemical system. Default='centre'
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# PART I: Perform checks on filepath

	# First, determine if filepath exists.
	if not os.path.exists(filepath):
		raise Exception(f'Error: Filepath does not exist. filepath: {filepath}')

	# Second, determine if filepath is a file or not.
	if not os.path.isfile(filepath):
		raise Exception(f'Error: Filepath exists, but is not a file. filepath: {filepath}')

	# Third, check move_centre_of_mass_to and obtain the value to move the centre of mass to. 
	if move_centre_of_mass_to == 'centre':
		move_centre_of_mass_to = np.array([0,0,0])
	else:
		try:
			move_centre_of_mass_to = np.array([float(pos) for pos in move_centre_of_mass_to.split(',')])
		except Exception as exception:
			toString  = f'Error, you must give three integers/float values for the centre of mass. These represent the (x,y,z) position that you want to move the centre of mass of your chemical system too.'
			toString += f'Input you gave for move_centre_of_mass_to = {move_centre_of_mass_to}'
			exit(toString)
	if (len(move_centre_of_mass_to) != 3):
		toString  = f'Error, you must give three integers/float values for the centre of mass. These represent the (x,y,z) position that you want to move the centre of mass of your chemical system too.'
		toString += f'Input you gave for move_centre_of_mass_to = {move_centre_of_mass_to}'
		exit(toString)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Part II: Move the centre of mass of the chemical system to the new position. 

	# Fourth, obtain the system from the filepath.
	system = read(filepath).copy()

	# Fifth, obtain the current centre of mass.
	current_centre_of_mass = np.array(system.get_center_of_mass())

	# Sixth, get the displacement vector to move the centre of mass - from current_centre_of_mass --to--> move_centre_of_mass_to
	com_displacement = move_centre_of_mass_to - current_centre_of_mass

	# Seventh, get the current positions of the atoms in the chemical system
	current_positions = system.get_positions()

	# Seventh, move the atoms 
	new_positions = current_positions + com_displacement

	# Eighth, set the new positions of the atoms in the system
	system.set_positions(new_positions)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Part III: Save the chemical system as a new file. 

	# Ninth, get pathname and basename
	filename = os.path.basename(filepath)
	dirname  = os.path.dirname (filepath)

	# Tenth, obtain the suffix of the filepath.
	if '.' not in filename:
		suffix = ''
	else:
		_, suffix = os.path.splitext(filename)

	# Eleventh, get the new filepath to save the file to.
	if len(dirname.strip()) == 0:
		new_filepath = ''
	else:
		new_filepath = dirname + '/'
	new_filepath += filename.replace(suffix,'') + '_com_' + str(move_centre_of_mass_to[0]) + '_' + str(move_centre_of_mass_to[1]) + '_' + str(move_centre_of_mass_to[2]) + suffix

	# Twelfth, save the molecule to file. 
	write(new_filepath, system)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

