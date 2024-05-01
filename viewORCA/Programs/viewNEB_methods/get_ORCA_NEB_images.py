"""
get_ORCA_NEB_images.py, Geoffrey Weal, 23/4/24

This method is designed to obtain the images and energies from the ORCA NEB file.
"""
from ase                         import Atoms, Atom
from ase.io                      import read

def get_ORCA_NEB_images(path_to_orca_NEB_trj_xyz):
	"""
	This method is designed to obtain the images and energies from the ORCA NEB file.

	Parameters
	----------
	path_to_orca_NEB_trj_xyz : str.
		This is the path to the orca NEB trajectory file.
	"""

	# First, initialise the lists for recording the images and their energies for the NEB calculation. 
	neb_images          = []
	neb_images_energies = []

	# Second, read the trajectory file
	with open(path_to_orca_NEB_trj_xyz) as trajectoryFILE:

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
					neb_images.append(system.copy())
					neb_images_energies.append(energy)

					# 4.4.2: Reset system files
					new_image = read_title = True
					system = Atoms()

				elif len(system) > number_of_atoms:
					# 4.5: Something weird happened to get to this point.
					raise Exception('Error: More atoms than should have.')

	# Fifth, check that all the entries in neb_images_energies are energies
	if any([(energy is None) for energy in neb_images_energies]):
		raise Exception('Error: Energy not given')

	# Sixth, return neb_images and neb_images_energies
	return neb_images, neb_images_energies