import os
from ase import Atoms, Atom
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator
from ase.visualize import view

class CLICommand:
	"""This module allows the user to view the output of an ORCA geometric optimisation job visually. 
	"""

	@staticmethod
	def add_arguments(parser):
		parser.add_argument('-p', '--path', nargs=1, help='This is the path to the folder containing the ORCA geometric optimisation job.', default=['.'])
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
		path_to_images = args.path
		if len(path_to_images) != 1:
			raise Exception('Error: input "path" must only have one input: '+str(path_to_images))
		path_to_images = path_to_images[0]

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
			raise Exception('Error: Input for viewOpt should be either True or False')

		# Third, run the method.
		Run_Method(path_to_images, view_Opt)

def Run_Method(path_to_images, view_Opt):
	"""
	This method is designed to allow the user to view the output of an ORCA geometric optimisation job visually. 

	Parmeters
	---------
	path_to_images : str.
		This is the path to the ORCA geometric optimisation jobs. 
	view_Opt : bool.
		This indicates if the user wants to view the ORCA geometric optimisation job after running this program.
	"""

	# First, obtain the images and the energies for the whole ORCA NEB process,
	opt_images, opt_images_energies = get_ORCA_OPT_images(path_to_images+'/'+'orca_trj.xyz')

	# Second, assign energies to the appropriate SCAN image
	for index, (neb_image, opt_images_energy) in enumerate(zip(opt_images, opt_images_energies)):
		opt_images[index].calc = SinglePointCalculator(atoms=opt_images[index], energy=opt_images_energy)

	# Third, save the SCAN process as an XYZ file
	write('OPT_images.xyz', opt_images)

	# Fourth, view the SCAN images with energies. 
	if view_Opt:
		view(opt_images)

# ------------------------------------------------------------------------------

def get_ORCA_OPT_images(path_to_orca_OPT_trj_xyz):
	"""
	This method is designed to obtain the images and energies from the ORCA geometry optimisation file.

	Parameters
	----------
	path_to_orca_OPT_trj_xyz : str.
		This is the path to the orca NEB trajectory file.
	"""

	# First, initialise the lists for recording the images and their energies for the optimisation calculation. 
	opt_images          = []
	opt_images_energies = []

	# Second, read the trajectory file
	with open(path_to_orca_OPT_trj_xyz) as trajectoryFILE:

		# Third, reset system files
		new_image = read_title = True
		system = Atoms()

		# Fourth, for each line in the ORCA trajectory file
		for line in trajectoryFILE:

			if new_image:
				# 4.1: If new image, read the number of atoms expect in the file
				number_of_atoms = int(line.rstrip())
				new_image = False
				read_title = True

			elif read_title:
				# 4.2: Read the energy from the title
				energy = float(line.rstrip().split()[-1])
				read_title = False

			else:
				# 4.3: Obtain the atom information.
				symbol, xx, yy, zz = line.rstrip().split()
				xx = float(xx); yy = float(yy); zz = float(zz); 
				system.append(Atom(symbol=symbol,position=(xx,yy,zz)))

				if len(system) == number_of_atoms:
					# 4.4: If you have obtained all the atoms for this traj step:

					# 4.4.1: Save the image and its energy to 
					opt_images.append(system.copy())
					opt_images_energies.append(energy)

					# 4.4.2: Reset system files
					new_image = read_title = True
					system = Atoms()

				elif len(system) > number_of_atoms:
					# 4.5: Something weird happened to get to this point.
					raise Exception('Error: More atoms than should have.')

	# Fifth, check that all the entries in opt_images_energies are energies
	if any([(energy is None) for energy in opt_images_energies]):
		raise Exception('Error: Energy not given')

	# Sixth, return opt_images and opt_images_energies
	return opt_images, opt_images_energies

# ------------------------------------------------------------------------------

