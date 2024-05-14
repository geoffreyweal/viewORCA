"""
get_multiplicity.py, Geoffrey Weal, 14/5/24

This module allows the user to obtain the number of electrons in the system of interest, and therefore obtain the possible multiplicities that the system could be in.
"""
import os
from ase.io import read

class CLICommand:
	"""This module allows the user to obtain the number of electrons in the system of interest, and therefore obtain the possible multiplicities that the system could be in.
	"""

	@staticmethod
	def add_arguments(parser):
		parser.add_argument('filepath',       nargs=1, help='This is the path to the xyz file you want to obtain the number of electrons of.')
		parser.add_argument('-c', '--charge', nargs=1, help='This is the overall charge of the molecule', default=['0'])

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
		charge = args.charge
		if len(charge) != 1:
			raise Exception('Error: input "charge" must only have one input: '+str(charge))
		charge = charge[0]

		# Third, run the method.
		Run_Method(filepath, charge)

def Run_Method(filepath, charge):
	"""
	This method will allow the user to obtain the number of electrons in the system of interest, and therefore obtain the possible multiplicities that the system could be in.

	Parmeters
	---------
	filepath : str.
		This is the path to the system file (preferably an xyz file). 
	charge : bool.
		This is the overall charge of the system. 
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# PART I: Perform checks

	# First, determine if filepath exists.
	if not os.path.exists(filepath):
		raise Exception(f'Error: Filepath does not exist. filepath: {filepath}')

	# Second, determine if filepath is a file or not.
	if not os.path.isfile(filepath):
		raise Exception(f'Error: Filepath exists, but is not a file. filepath: {filepath}')

	# Third, make sure that charge is a digit
	try:
		charge = int(charge)
	except Exception as exception:
		raise Exception(f'Error: charge must be an integer. charge = {charge}')

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Part II: Obtan the number of electrons in the molecule. 

	# Fourth, obtain the system from the filepath.
	system = read(filepath)

	# Fifth, get the atomic numbers of the atoms in the system. This is equivalent to the number of electrons in the neutrally charged system. 
	total_number_of_electrons_in_neutral_system = int(sum(system.get_atomic_numbers()))

	# Sixth, obtain the number of extra electrons in the system.
	no_of_extra_electrons = -int(charge)

	# Seventh, obtain the total number of electrons in the system. 
	total_number_of_electrons = total_number_of_electrons_in_neutral_system + no_of_extra_electrons

	# Eighth, make sure that the total number of electrons in the system is not a negative number, which makes no sense.
	if total_number_of_electrons < 0:
		raise Exception('Error: Number of electrons in the system is a negative number?: total_number_of_electrons = '+str(total_number_of_electrons))

	# Ninth, print the number of electrons in the system. 
	print(f'Number of electrons in {os.path.basename(filepath)} is: {total_number_of_electrons}')

	# Tenth, print if there are an odd or even number of electrons
	if total_number_of_electrons % 2 == 0:
		number_type = 'even'
		multiplcity_number_type = 'odd'
		lowest_multiplcity_number = 1
		example_multiplcity_numbers = '1,3,5,7,9,11,..'
	elif total_number_of_electrons % 2 == 1:
		number_type = 'odd'
		multiplcity_number_type = 'even'
		lowest_multiplcity_number = 2
		example_multiplcity_numbers = '2,4,6,8,10,12,..'
	else:
		raise Exception(f'Error: The total number of electrons in the system is neither even or odd? total_number_of_electrons = {total_number_of_electrons}')
	print(f'The total number of electrons in the system is {number_type}.')

	# Eleventh, Indicate what the multiplicity could be, and give examples of possibly multiplicities.
	print(f'The multiplicity can be an {multiplcity_number_type} number greater than or equal to {lowest_multiplcity_number} based on the M=2S+1 rule. Examples: {example_multiplcity_numbers}')


	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

