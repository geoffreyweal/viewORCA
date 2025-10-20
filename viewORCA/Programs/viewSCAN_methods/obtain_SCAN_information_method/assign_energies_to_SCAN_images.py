"""
assign_energies_to_SCAN_images.py, Geoffrey Weal, 18/4/24

This method is designed to assign energies from the orca trajectory file to the xyz files for the SCAN images.
"""
import warnings
from ase import Atoms, Atom

def assign_energies_to_SCAN_images(path_to_orca_trj_xyz, scan_images=[]):
	"""
	This method is designed to assign energies from the orca trajectory file to the xyz files for the SCAN images.

	Parameters
	----------
	path_to_orca_trj_xyz : str.
		This is the path to the orca trajectory file.
	scan_images : list of ase.Atoms
		These are all the steps in the SCAN trajectory.
	"""

	# First, set the energies for the SCAN job to None initally. 
	scan_images_energies = [None for _ in range(len(scan_images))]

	# Second, read the trajectory file
	with open(path_to_orca_trj_xyz) as trajectoryFILE:

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

					# 4.4.1: Compare trajectory step to images in scan_images
					for index, image in enumerate(scan_images):
						if compare(image, system):
							scan_images_energies[index] = energy
							break

					# 4.4.2: Reset system files
					new_image = read_title = True
					system = Atoms()

				elif len(system) > number_of_atoms:
					# 4.5: Something weird happened to get to this point.
					raise Exception('Error: More atoms than should have.')

	# Fifth, check that all the entries in scan_images_energies are energies
	if any([(energy is None) for energy in scan_images_energies]):
		warnings.warn('Warning: Some of your molecules may not have an energy value assigned to it. This may be because your optimisation is still running or did not complete fully. Will continue running viewORCA, but be warned.')

	# Sixth, return scan_images_energies
	return scan_images_energies

# ------------------------------------------------------------------------------

def get_distance(position1, position2):
	"""
	This method gives the bond length between two atoms in the molecule. Does not use Minimum Image Convention. 

	Parameters
	----------
	position1 : np.array
		A 1x3 row matrix containing the x, y, and z coordinates of the atom1. Given in Å. 
	position2 : np.array
		A 1x3 row matrix containing the x, y, and z coordinates of the atom2. Given in Å. 

	Returns
	-------
	float
		Returns the distance between those atom1 and atom2. Given in Å. 
		
	"""
	return (sum([(p1-p2)**2.0 for p1, p2 in zip(position1,position2)]))**0.5

def compare(image, system):
	"""
	This method is designed to check if image and system are the same

	Parameters
	----------
	image : ase.Atoms
		This is the SCAN image. 
	system : ase.Atoms
		This is the image in the ORCA trajectory to compare to. This contains the energy of the trajectory step.
	
	Returns
	-------
		True if image and system are the same. False if not
	"""

	for image_atom, system_atom in zip(image, system):
		if not (image_atom.symbol == system_atom.symbol):
			raise Exception('Error')
		if not get_distance(image_atom.position, system_atom.position) < 0.00001:
			return False

	return True

# ------------------------------------------------------------------------------
