"""
viewNEB_snapshot.py, Geoffrey Weal, 23/4/24

This module allows the user to view the progress of an ORCA Nudged Elastic Band (NEB) job (as it is running or after it has run) by viewing the energy profiles of all the iterations of the NEB optimisation job.
"""
import os

from viewORCA.Programs.viewNEB_snapshot_methods.get_interpolation_filename   import get_interpolation_filename
from viewORCA.Programs.viewNEB_snapshot_methods.read_orca_interpolation_file import read_orca_interpolation_file
from viewORCA.Programs.viewNEB_snapshot_methods.plot_neb                     import plot_neb

class CLICommand:
	"""This module allows the user to view the progress of an ORCA Nudged Elastic Band (NEB) job (as it is running or after it has run) by viewing the energy profiles of all the iterations of the NEB optimisation job.
	"""

	@staticmethod
	def add_arguments(parser):
		parser.add_argument('-p', '--path', nargs=1, help='This is the path to the folder containing the ORCA NEB job. You can also give the path directly to the ".interp" file you want to process.', default=['.'])

	@staticmethod
	def run(args):
		"""
		Run this program. This will allow the user to view the optimisation process for an ORCA NEB job. 

		Parameters
		----------
		args : argparse.Namespace
			This is the namespace that give all the information about the arguments.
		"""

		# First, obtain the local path to the ORCA NEB job.
		path_to_NEB_job = args.path
		if len(path_to_NEB_job) != 1:
			raise Exception('Error: input "path" must only have one input: '+str(path_to_NEB_job))
		path_to_NEB_job = path_to_NEB_job[0]

		# Third, run the method.
		Run_Method(path_to_NEB_job)

def Run_Method(path_to_NEB_job):
	"""
	This will allow the user to view the optimisation process for an ORCA NEB job. 

	Parmeters
	---------
	path_to_NEB_job : str.
		This is the directory containing the ORCA NEB calculation. 
	"""

	# Prelimiary Step 1: Determine if path_to_NEB_job is a folder or a file.
	if   not os.path.exists(path_to_NEB_job):
		# Path does not exist.
		to_string  = 'Error: Path given does not exist.\n'
		to_string += f'Given path: {path_to_NEB_job}\n'
		to_string += 'Check this.'
		raise Exception(to_string)
	elif os.path.isdir(path_to_NEB_job):
		# Path is to a folder: Get the ".interp" file from path_to_NEB_job
		filename = get_interpolation_filename(path_to_NEB_job)
		NEB_interpolation_filepath = path_to_NEB_job+'/'+filename
	elif os.path.isfile(path_to_NEB_job):
		# Path is to a file: Check it is a interpolation file
		if not path_to_NEB_job.endswith('.interp'):
			to_string  = 'Error: You file is not a iterpolation file.\n'
			to_string += 'The file must end with the suffix ".interp".\n'
			to_string += f'Given filepath: {path_to_NEB_job}\n'
			to_string += 'Check this.'
			raise Exception(to_string)
		NEB_interpolation_filepath = path_to_NEB_job
	else:
		# Path is to neither a file or directory.
		to_string  = 'Error: Path given is to neither a file or directory.\n'
		to_string += f'Given filepath: {path_to_NEB_job}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 2: Check that NEB_interpolation_filepath is a file.
	if not os.path.isfile(NEB_interpolation_filepath):
		to_string  = 'Error: Run_Method requires a file for the "NEB_interpolation_filepath" variable".\n'
		to_string += f'NEB_interpolation_filepath =  {NEB_interpolation_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 3: Check that NEB_interpolation_filepath ends with the '.interp' suffix
	#                    * This indicates it is a ORCA interpolation file.
	if not NEB_interpolation_filepath.endswith('.interp'):
		to_string  = 'Error: You file is not a iterpolation file.\n'
		to_string += 'The file must end with the suffix ".interp".\n'
		to_string += f'Given filepath: {NEB_interpolation_filepath}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Prelimiary Step 4: Get the directory name that the interpolation file is in.
	file_path = os.path.dirname(NEB_interpolation_filepath)

	# Prelimiary Step 5: Obtain the filename of the interpolation file.
	filename = os.path.basename(NEB_interpolation_filepath)

	# Prelimiary Step 6: Obtain the prefix of the filename of the interpolation file.
	#                    * This should be the name of the ORCA input file.
	filename_prefix = filename.replace('.final.interp','').replace('.interp','')

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

	# First, obtain the ORCA images and splines from the interpolation file. 
	images, splines = read_orca_interpolation_file(NEB_interpolation_filepath)

	# Second, obtain the number of iterations that ORCA performed. 
	no_of_iters = len(images)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
	# Third, perform several checks.

	# 3.1: Make sure that the number of images equal the number of splines obtained from the interpolation file. 
	if len(images) != len(splines):
		to_string  = 'Error: images list is not the same as the splines list.\n'
		to_string += f'len(images) = {len(images)}\n'
		to_string += f'len(splines) = {len(splines)}\n'
		to_string += f'This may indicate the {filename} file has been cut short, and missing some spline data.'
		to_string += 'Check this.'
		raise RuntimeError(to_string)

	# 3.2: Give a note to the user if the interpolation file only contain one iteration. 
	if no_of_iters == 1:
		if 'final' in filename.lower():
			print('*** Note that %s contains only the last iteration of a NEB/CI-NEB run  ***' % (NEB_interpolation_filepath+'/'+filename))
		else:
			print('%s contains only one iteration?   ***' % filename)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
	
	# Fourth, create a plot giving the optimisation information from the ORCA NEB job. 
	plot_neb(images, splines, file_path, filename, filename_prefix)

# ================================================================================================
