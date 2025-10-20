"""
check_prerequisites.py, Geoffrey Weal, 19/4/24

This method is designed to check that you have all the pre-requisite programs for running the viewORCA program. 
"""
import sys
import importlib
import importlib.util

def check_prerequisites(version):
    """
    This method is designed to check that you have all the pre-requisite programs for running the viewORCA program. 

    Parameters
    ----------
    version : str.
        This is the version of viewORCA
    """

    # First, check that you are using a value version of python
    if sys.version_info[0] == 2:
        toString = ''
        toString += '\n'
        toString += '================================================'+'\n'
        toString += 'This is the viewORCA Program'+'\n'
        toString += 'Version: '+str(version)+'\n'
        toString += '\n'
        toString += 'The viewORCA program requires Python3. You are attempting to execute this program in Python2.'+'\n'
        toString += 'Make sure you are running the viewORCA program in Python3 and try again'+'\n'
        toString += 'This program will exit before beginning'+'\n'
        toString += '================================================'+'\n'
        raise ImportError(toString)
    if sys.version_info[1] < 4:
        toString = ''
        toString += '\n'
        toString += '================================================'+'\n'
        toString += 'This is the viewORCA Program'+'\n'
        toString += 'Version: '+str(version)+'\n'
        toString += '\n'
        toString += 'The viewORCA program requires Python 3.4 or greater.'+'\n'
        toString += 'You are using Python '+str('.'.join(sys.version_info))
        toString += '\n'
        toString += 'Use a version of Python 3 that is greater or equal to Python 3.4.\n'
        toString += 'This program will exit before beginning'+'\n'
        toString += '================================================'+'\n'
        raise ImportError(toString)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    '''
    # Second, make sure you have the packaging package
    packaging_spec = importlib.util.find_spec("packaging")
    packaging_found = (packaging_spec is not None)
    if not packaging_found:
        toString = ''
        toString += '\n'
        toString += '================================================'+'\n'
        toString += 'This is the viewORCA Program'+'\n'
        toString += 'Version: '+str(__version__)+'\n'
        toString += '\n'
        toString += 'The viewORCA program requires the "packaging" program.'+'\n'
        toString += '\n'
        toString += 'Install packaging through pip by following the instruction in https://github.com/geoffreyweal/viewORCA'+'\n'
        toString += 'These instructions will ask you to install packaging by typing the following into your terminal\n'
        toString += '\n'
        toString += 'pip install --user --upgrade packaging\n'
        toString += '\n'
        toString += 'This program will exit before beginning'+'\n'
        toString += '================================================'+'\n'
        raise ImportError(toString)
    '''

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Third, make sure you have the ase package
    ase_spec = importlib.util.find_spec("ase")
    ase_found = (ase_spec is not None)
    if not ase_found:
        toString = ''
        toString += '\n'
        toString += '================================================'+'\n'
        toString += 'This is the viewORCA Program'+'\n'
        toString += 'Version: '+str(version)+'\n'
        toString += '\n'
        toString += 'The viewORCA program requires ASE.'+'\n'
        toString += '\n'
        toString += 'Install ASE through pip by following the instruction in https://github.com/geoffreyweal/viewORCA'+'\n'
        toString += 'These instructions will ask you to install ase by typing the following into your terminal\n'
        toString += 'pip3 install --user --upgrade ase\n'
        toString += '\n'
        toString += 'This program will exit before beginning'+'\n'
        toString += '================================================'+'\n'
        raise ImportError(toString) 

    # Fourth, make sure that the version of ase is valid for the viewORCA program.
    import ase
    ase_version_minimum = '3.19.0'
    from packaging import version
    #from distutils.version import StrictVersion
    #if StrictVersion(ase.__version__) < StrictVersion(ase_version_minimum):
    if version.parse(ase.__version__) < version.parse(ase_version_minimum):
        toString = ''
        toString += '\n'
        toString += '================================================'+'\n'
        toString += 'This is the viewORCA Program'+'\n'
        toString += 'Version: '+str(version)+'\n'
        toString += '\n'
        toString += 'The viewORCA program requires ASE greater than or equal to '+str(ase_version_minimum)+'.'+'\n'
        toString += 'The current version of ASE you are using is '+str(ase.__version__)+'.'+'\n'
        toString += '\n'
        toString += 'Install ASE through pip by following the instruction in https://github.com/geoffreyweal/viewORCA'+'\n'
        toString += 'These instructions will ask you to install ase by typing the following into your terminal\n'
        toString += 'pip3 install --user --upgrade ase\n'
        toString += '\n'
        toString += 'This program will exit before beginning'+'\n'
        toString += '================================================'+'\n'
        raise ImportError(toString)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Fifth, make sure you have the numpy package
    packaging_spec = importlib.util.find_spec("numpy")
    packaging_found = (packaging_spec is not None)
    if not packaging_found:
        toString = ''
        toString += '\n'
        toString += '================================================'+'\n'
        toString += 'This is the viewORCA Program'+'\n'
        toString += 'Version: '+str(__version__)+'\n'
        toString += '\n'
        toString += 'The viewORCA program requires the "numpy" program.'+'\n'
        toString += '\n'
        toString += 'Install numpy through pip by following the instruction in https://github.com/geoffreyweal/viewORCA'+'\n'
        toString += 'These instructions will ask you to install numpy by typing the following into your terminal\n'
        toString += '\n'
        toString += 'pip install --user --upgrade numpy\n'
        toString += '\n'
        toString += 'This program will exit before beginning'+'\n'
        toString += '================================================'+'\n'
        raise ImportError(toString)

    # Sixth, make sure you have the matplotlib package
    packaging_spec = importlib.util.find_spec("scipy")
    packaging_found = (packaging_spec is not None)
    if not packaging_found:
        toString = ''
        toString += '\n'
        toString += '================================================'+'\n'
        toString += 'This is the viewORCA Program'+'\n'
        toString += 'Version: '+str(__version__)+'\n'
        toString += '\n'
        toString += 'The viewORCA program requires the "scipy" program.'+'\n'
        toString += '\n'
        toString += 'Install scipy through pip by following the instruction in https://github.com/geoffreyweal/viewORCA'+'\n'
        toString += 'These instructions will ask you to install scipy by typing the following into your terminal\n'
        toString += '\n'
        toString += 'pip install --user --upgrade scipy\n'
        toString += '\n'
        toString += 'This program will exit before beginning'+'\n'
        toString += '================================================'+'\n'
        raise ImportError(toString)

    # Seventh, make sure you have the matplotlib package
    packaging_spec = importlib.util.find_spec("matplotlib")
    packaging_found = (packaging_spec is not None)
    if not packaging_found:
        toString = ''
        toString += '\n'
        toString += '================================================'+'\n'
        toString += 'This is the viewORCA Program'+'\n'
        toString += 'Version: '+str(__version__)+'\n'
        toString += '\n'
        toString += 'The viewORCA program requires the "matplotlib" program.'+'\n'
        toString += '\n'
        toString += 'Install matplotlib through pip by following the instruction in https://github.com/geoffreyweal/viewORCA'+'\n'
        toString += 'These instructions will ask you to install matplotlib by typing the following into your terminal\n'
        toString += '\n'
        toString += 'pip install --user --upgrade matplotlib\n'
        toString += '\n'
        toString += 'This program will exit before beginning'+'\n'
        toString += '================================================'+'\n'
        raise ImportError(toString)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
