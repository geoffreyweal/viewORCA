import os
from ase import Atoms, Atom
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator
from ase.visualize import view

class CLICommand:
	"""This module allows the user to view the output of an ORCA Intrinsic Reaction Coordinate (IRC) job visually. 
	"""

	@staticmethod
	def add_arguments(parser):
		parser.add_argument('-p', '--path', nargs=1, help='This is the path to the folder containing the ORCA IRC job.', default=['.'])
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
		path_to_images = args.path
		if len(path_to_images) != 1:
			raise Exception('Error: input "path" must only have one input: '+str(path_to_images))
		path_to_images = path_to_images[0]

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
		Run_Method(path_to_images, view_IRC)

def Run_Method(path_to_images, view_IRC):
	"""
	This method is designed to allow the user to view the output of an ORCA geometric optimisation job visually. 

	Parmeters
	---------
	path_to_images : str.
		This is the path to the ORCA geometric optimisation jobs. 
	view_IRC : bool.
		This indicates if the user wants to view the ORCA IRC job after running this program.
	"""

	# First, obtain the images and the energies for the whole ORCA NEB process,
	irc_images, irc_images_energies = get_ORCA_IRC_images(path_to_images+'/'+'orca_IRC_Full_trj.xyz')

	# Second, assign energies to the appropriate SCAN image
	for index, (irc_image, irc_images_energy) in enumerate(zip(irc_images, irc_images_energies)):
		irc_images[index].calc = SinglePointCalculator(atoms=irc_images[index], energy=irc_images_energy)

	# Third, save the SCAN process as an XYZ file
	write('IRC_images.xyz', irc_images)

	# Fourth, view the SCAN images with energies. 
	if view_IRC:
		view(irc_images)

# ------------------------------------------------------------------------------

def get_ORCA_IRC_images(path_to_orca_IRC_trj_xyz):
	"""
	This method is designed to obtain the images and energies from the ORCA IRC file.

	Parameters
	----------
	path_to_orca_IRC_trj_xyz : str.
		This is the path to the orca IRC trajectory file.
	"""

	neb_images          = []
	neb_images_energies = []

	# First: Read the trajectory file
	with open(path_to_orca_IRC_trj_xyz) as trajectoryFILE:

		# Second: Reset system files
		new_image = read_title = True
		system = Atoms()

		# Third: For each line in the ORCA trajectory file
		for line in trajectoryFILE:

			if new_image:
				# 3.1: If new image, read the number of atoms expect in the file
				number_of_atoms = int(line.rstrip())
				new_image = False
				read_title = True

			elif read_title:
				# 3.2: Read the energy from the title
				energy = float(line.rstrip().split()[-1])
				read_title = False

			else:
				# 3.3: Obtain the atom information.
				symbol, xx, yy, zz = line.rstrip().split()
				xx = float(xx); yy = float(yy); zz = float(zz); 
				system.append(Atom(symbol=symbol,position=(xx,yy,zz)))

				if len(system) == number_of_atoms:
					# 3.4: If you have obtained all the atoms for this traj step:

					# 3.4.1: Save the image and its energy to 
					neb_images.append(system.copy())
					neb_images_energies.append(energy)

					# 3.4.2: Reset system files
					new_image = read_title = True
					system = Atoms()

				elif len(system) > number_of_atoms:
					# 3.5: Something weird happened to get to this point.
					raise Exception('Error: More atoms than should have.')

	# Fourth: Check that all the entries in neb_images_energies are energies
	if any([(energy is None) for energy in neb_images_energies]):
		raise Exception('Error: Energy not given')

	return neb_images, neb_images_energies

# ------------------------------------------------------------------------------

