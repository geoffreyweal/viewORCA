"""
combine_systems_together.py, Geoffrey Weal, 14/5/24

This module allows the user to combine a set of chemical systems together. 
"""
import os
import numpy as np
from ase import Atoms
from ase.io import read, write

class CLICommand:
	"""This module allows the user to combine a set of chemical systems together. 
	"""

	@staticmethod
	def add_arguments(parser):
		parser.add_argument('filepaths', nargs='*', help='These are the filepath of all the molecules you want to combine together into one molecule.')

	@staticmethod
	def run(args):
		"""
		Run this program. This module allows the user to combine a set of chemical systems together. 

		Parameters
		----------
		args : argparse.Namespace
			This is the namespace that give all the information about the arguments.
		"""

		# First, obtain the local filepaths for the chemical systems you want to combine together. 
		filepaths = args.filepaths

		# Second, run the method.
		Run_Method(filepaths)

def Run_Method(filepaths):
	"""
	This module allows the user to combine a set of chemical systems together. 

	Parmeters
	---------
	filepaths : list of str.
		These are the paths to the chemical systems you want to combine together, in order of first entry to last entry.
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# PART I: Perform checks on all the filepaths in the filepath list. 

	# First, perform several checks on the chemical system files you want to combine together. 
	for filepath in filepaths:

		# 1.1: Determine if filepath exists.
		if not os.path.exists(filepath):
			raise Exception(f'Error: Filepath does not exist. filepath: {filepath}')

		# 1.2: Determine if filepath is a file or not.
		if not os.path.isfile(filepath):
			raise Exception(f'Error: Filepath exists, but is not a file. filepath: {filepath}')

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Part II: Combine all the chemical systems together

	# Second, read in all the chemical systems:
	systems = [read(filepath).copy() for filepath in filepaths] 

	# Third, check that the molecules are at some distance from each other. 
	shortest_distances = {}
	for index1 in range(len(systems)):
		filepath1 = filepaths[index1]
		system1   = systems[index1]
		for index2 in range(index1+1,len(systems)):
			filepath2 = filepaths[index2]
			system2   = systems[index2]
			shortest_distance = get_shortest_distance(system1, system2)
			shortest_distances[(index1,index2)] = shortest_distance

	# Fourth, print shortest_distances
	for (index1, index2), shortest_distances in sorted(shortest_distances.items(), key=lambda x:x[0]):
		filepath1 = filepaths[index1]
		filepath2 = filepaths[index2]
		print(f'{filepath1} and {filepath2}: {round(shortest_distances,4)} A')

	overall_system = Atoms()

	for system in systems:
		overall_system += system.copy()



	# Twelfth, save the molecule to file. 
	write('overall_chemical_system.xyz', overall_system)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def get_shortest_distance(system1, system2):
	"""

	"""

	shortest_distance = float('inf')

	for position_system1 in system1.get_positions():
		for position_system2 in system2.get_positions():
			distance = np.linalg.norm(position_system2-position_system1)
			if distance < shortest_distance:
				shortest_distance = distance

	return shortest_distance














