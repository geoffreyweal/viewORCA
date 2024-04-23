import os
from setuptools import setup

def get_version_number():
  path_to_written_version = 'viewORCA/__init__.py'
  with open(path_to_written_version) as initPY:
    for line in initPY:
      if line.startswith('__version__'):
        version = eval(line.rstrip().replace('__version__ = ',''))
        break
  return version

def get_long_description():
  this_directory = os.path.abspath(os.path.dirname(__file__))
  readme_file = os.path.join(this_directory, 'README.md')
  with open(readme_file, encoding='utf-8') as f:
    long_description = f.read()
  return long_description

def find_packages(root):
  packages = []
  for root, dirs, files in os.walk(root, topdown=False):
    for name in files:
      if name.endswith('.py'):
        file = os.path.relpath(os.path.join(root, name))
        dirname = os.path.dirname(file)
        rel_dirname = os.path.relpath(dirname)
        if not rel_dirname in packages:
          packages.append(rel_dirname)
  return sorted(packages)

setup(name='viewORCA',
      packages=find_packages(root='viewORCA'),
      scripts=['bin/viewORCA'],
      version=get_version_number(),
      description="This program is designed to allow the user to inspect calculations from ORCA. This is to be used along side the ORCA Reaction Mechanism Tutorial.",
      long_description=get_long_description(),
      long_description_content_type='text/markdown',
      author='Dr. Geoffrey R. Weal',
      author_email='geoffrey.weal@gmail.ac.nz',
      url = 'https://geoffreyweal.github.io/ORCA_Mechanism_Procedure/',
      download_url = 'https://github.com/geoffreyweal/viewORCA/archive/v'+str(get_version_number())+'.tar.gz',
      license='GNU AFFERO GENERAL PUBLIC LICENSE',
      zip_safe=False,
      keywords = ['computational-chemistry', 'DFT', 'atomic-simulation-environment', 'ORCA'],
      install_requires=['packaging', 'numpy', 'scipy', 'matplotlib', 'ase>=3.19.0'],
      classifiers=[
        'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Science/Research',      # Define that your audience are developers
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Unix',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Natural Language :: English',
        'License :: OSI Approved :: GNU Affero General Public License v3',   # Again, pick a license
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        ],
      )

