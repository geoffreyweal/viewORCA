"""
get_ORCA_IRC_images.py, Geoffrey Weal, 23/4/24

This method is designed to obtain the images and energies from the ORCA IRC file.
"""
import warnings
from ase    import Atoms, Atom
from ase.io import read

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
		warnings.warn('Warning: Some of your molecules may not have an energy value assigned to it. This may be because your optimisation is still running or did not complete fully. Will continue running viewORCA, but be warned.')

	return neb_images, neb_images_energies
