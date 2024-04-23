# ``viewORCA`` Manual

There are several modules available for viewing ORCA calculations. This webpage describes all the modules available in ``viewORCA``:

???+ info

	You can view all the module available in ``viewORCA`` by typing ``viewORCA --help`` into the terminal:

	```bash
	username@computer name Desktop % viewORCA --help
	usage: viewORCA [-h] [-T] {help,opt,scan,neb,neb_snap,irc} ...

	Program for viewing jobs from ORCA

	optional arguments:
	  -h, --help            show this help message and exit
	  -T, --traceback

	Sub-commands:
	  {help,opt,scan,neb,neb_snap,irc}
	    help                Help for sub-command.
	    opt                 This module allows the user to view the output of an ORCA geometric optimisation job visually.
	    scan                This module allows the user to view the output of an ORCA SCAN job visually.
	    neb                 This module allows the user to view the output of an ORCA Nudged Elastic Band (NEB) job visually.
	    neb_snap            This module allows the user to view the progress of an ORCA Nudged Elastic Band (NEB) job (as it is running or after it has run) by viewing the energy profiles of all the iterations of the NEB optimisation job.
	    irc                 This module allows the user to view the output of an ORCA Intrinsic Reaction Coordinate (IRC) job visually.
	```


## ``viewORCA opt``: View Geometric Optimisation Calculation

This module is designed to allow the user to view the output of an ORCA geometric optimisation job visually. 

To use this module:

1. In the terminal, change directory into the directory containing your geometric optimisation job, and 
2. type ``viewORCA opt`` into the terminal. 

```bash
# First, change directory into the folder containing your ORCA geometric optimisation job.
cd path_to_geometric_optimisation_job

# Second, run "viewORCA opt".
viewORCA opt
```

``viewORCA opt`` will extract information from the ``ORCA_INPUT_FILENAME_trj.xyz`` file and show the images of the optimisation, along with the energy profile of the optimisation. 

* For this reason, you **MUST** have the ``ORCA_INPUT_FILENAME_trj.xyz`` file in your geometric optimisation folder for ``viewORCA opt`` to work.
* The prefix ``ORCA_INPUT_FILENAME`` is the same prefix you gave to your ORCA input file (being ``ORCA_INPUT_FILENAME.inp``). 

This program will also create an xyz file called ``ORCA_INPUT_FILENAME_OPT_images.xyz`` that contains the images and the energy profile from the optimisation. 

* You can view this by typing ``ase gui ORCA_INPUT_FILENAME_OPT_images.xyz`` into the terminal.

	```bash
	ase gui ORCA_INPUT_FILENAME_OPT_images.xyz
	```

Optional inputs for ``viewORCA opt`` include:

* ``--path``: This is the path to the folder containing the ORCA optimisation calculation. If not given, the default path is the current directory you are in. 
* ``--view``: If tag indicates if you want to view the ASE GUI for the optimisation job immediate after this program has run. By default this is ``True``. If you only want to create the ``ORCA_INPUT_FILENAME_OPT_images.xyz`` file and not view it immediately, set this to ``False``.
	* Setting this to ``False`` is ideal if you are running ``viewORCA`` on an HPC, where ``ase gui`` windows may be more responsive if ``ase gui`` is run from your own computer rather than directly on the HPC. 
	* You can view ``ORCA_INPUT_FILENAME_OPT_images.xyz`` on your computer by typing the following into the terminal: ``ase gui ORCA_INPUT_FILENAME_OPT_images.xyz``.

???+ example

	An example of a geometric optimisation viewed using ``viewORCA opt`` is shown below, along with the energy profile for this optimisation. This example is given in [Examples/viewORCA_opt](https://github.com/geoffreyweal/viewORCA/tree/main/Examples/viewORCA_opt)

	<figure markdown="span">
	<img src="Figures/view_ORCA_Manual/opt/optimisation.gif?raw=true" alt="NEB Images from a geometric optimisation" width="600"/>
	<figcaption>NEB Images from a geometric optimisation</figcaption>
	</figure>

	The energy profile for this example is given below:

	<figure markdown="span">
	    <img src="Figures/view_ORCA_Manual/opt/energy_profile.png?raw=true" alt="Energy Profile from a geometric optimisation" width="600"/>
	    <figcaption>Energy Profile from a geometric optimisation</figcaption>
	</figure>


## ``viewORCA scan``: View SCAN Calculation

This module is designed to allow the user to view the output of an ORCA SCAN job visually. 

To use this module:

1. In the terminal, change directory into the directory containing your SCAN job, and 
2. type ``viewORCA scan`` into the terminal. 

```bash
# First, change directory into the folder containing your ORCA SCAN job.
cd path_to_SCAN_job

# Second, run "viewORCA scan".
viewORCA scan
```

``viewORCA scan`` will show the images of the SCAN calculation, along with the energy profile of the SCAN calculation. 

* This program will also create an xyz file called ``SCAN_images.xyz`` that contains the images and the energy profile from the SCAN calculation. You can view this by typing ``ase gui SCAN_images.xyz`` into the terminal.

	```bash
	ase gui SCAN_images.xyz
	```

Optional inputs for ``viewORCA scan`` include:

* ``--path``: This is the path to the folder containing the ORCA SCAN calculation. If not given, the default path is the current directory you are in. 
* ``--view``: If tag indicates if you want to view the ASE GUI for the SCAN job immediate after this program has run. By default this is ``True``. If you only want to create the ``SCAN_images.xyz`` file and not view it immediately, set this to ``False``.
	* Setting this to ``False`` is ideal if you are running ``viewORCA`` on an HPC, where ``ase gui`` windows may be more responsive if ``ase gui`` is run from your own computer rather than directly on the HPC. 
	* You can view ``SCAN_images.xyz`` on your computer by typing the following into the terminal: ``ase gui SCAN_images.xyz``.

???+ example

	An example of a SCAN calculation viewed using ``viewORCA scan`` is shown below, along with the energy profile for this SCAN calculation. This example is given in [Examples/viewORCA_scan](https://github.com/geoffreyweal/viewORCA/tree/main/Examples/viewORCA_scan)

	<figure markdown="span">
	<img src="Figures/view_ORCA_Manual/scan/SCAN_example.gif?raw=true" alt="NEB Images from a SCAN calculation" width="600"/>
	<figcaption>NEB Images from a SCAN calculation</figcaption>
	</figure>

	The energy profile for this example is given below:

	<figure markdown="span">
	    <img src="Figures/view_ORCA_Manual/scan/SCAN_energy.png?raw=true" alt="Energy Profile from a SCAN calculation" width="600"/>
	    <figcaption>Energy Profile from a SCAN calculation</figcaption>
	</figure>


## ``viewORCA neb``: View Nudged Elastic Band (NEB) Calculation

This module is designed to allow the user to view the output of an ORCA NEB job visually. 

To use this module:

1. In the terminal, change directory into the directory containing your NEB job, and 
2. type ``viewORCA neb`` into the terminal. 

```bash
# First, change directory into the folder containing your ORCA NEB job.
cd path_to_NEB_job

# Second, run "viewORCA neb".
viewORCA neb
```

``viewORCA neb`` will extract information from the ``ORCA_INPUT_FILENAME_MEP_trj.xyz`` file and show the images of the NEB calculation, along with the energy profile of the NEB calculation. 

* For this reason, you **MUST** have the ``ORCA_INPUT_FILENAME_MEP_trj.xyz`` file in your NEB calculation folder for ``viewORCA neb`` to work.
* The prefix ``ORCA_INPUT_FILENAME`` is the same prefix you gave to your ORCA input file (being ``ORCA_INPUT_FILENAME.inp``). 

This program will also create an xyz file called ``ORCA_INPUT_FILENAME_NEB_images.xyz`` that contains the images and the energy profile from the NEB calculation. 

* You can view this by typing ``ase gui ORCA_INPUT_FILENAME_NEB_images.xyz`` into the terminal.

	```bash
	ase gui ORCA_INPUT_FILENAME_NEB_images.xyz
	```

Optional inputs for ``viewORCA neb`` include:

* ``--path``: This is the path to the folder containing the ORCA NEB calculation. If not given, the default path is the current directory you are in. 
* ``--view``: If tag indicates if you want to view the ASE GUI for the NEB job immediate after this program has run. By default this is ``True``. If you only want to create the ``ORCA_INPUT_FILENAME_NEB_images.xyz`` file and not view it immediately, set this to ``False``.
	* Setting this to ``False`` is ideal if you are running ``viewORCA`` on an HPC, where ``ase gui`` windows may be more responsive if ``ase gui`` is run from your own computer rather than directly on the HPC. 
	* You can view ``ORCA_INPUT_FILENAME_NEB_images.xyz`` on your computer by typing the following into the terminal: ``ase gui ORCA_INPUT_FILENAME_NEB_images.xyz``.

???+ example

	An example of a NEB calculation viewed using ``viewORCA neb`` is shown below, along with the energy profile for this NEB calculation. This example is given in [Examples/viewORCA_neb_and_neb_snap](https://github.com/geoffreyweal/viewORCA/tree/main/Examples/viewORCA_neb_and_neb_snap)

	<figure markdown="span">
	<img src="Figures/view_ORCA_Manual/neb/NEB_example.gif?raw=true" alt="NEB Images from a NEB calculation" width="500"/>
	<figcaption>NEB Images from a NEB calculation</figcaption>
	</figure>

	The energy profile for this example is given below:

	<figure markdown="span">
	    <img src="Figures/view_ORCA_Manual/neb/NEB_energy.png?raw=true" alt="Energy Profile from a NEB calculation" width="600"/>
	    <figcaption>Energy Profile from a NEB calculation</figcaption>
	</figure>


## ``viewORCA neb_snap``: View NEB Optimisation Process

This module is designed to allow the user to view the energy profile of all the iterations performed during the NEB optimisation process. This module is an adaptation of [neb_snapshots.py](https://github.com/via9a/neb_visualize.py) by Vilhjálmur Ásgeirsson and Benedikt Orri Birgirsson of the Universiy of Iceland. 

* It is helpful to view this as your NEB calculation is running to check if there are any issues occurring as ORCA us running the NEB calculation. 

To use this module:

1. In the terminal, change directory into the directory containing your NEB job, and 
2. type ``viewORCA neb_snap`` into the terminal. 

```bash
# First, change directory into the folder containing your ORCA NEB job.
cd path_to_NEB_job

# Second, run "viewORCA neb_snap".
viewORCA neb_snap
```

``viewORCA neb_snap`` will extract information from the ``ORCA_INPUT_FILENAME.interp`` file and create plots of the energy profiles of all the iterations of the NEB optimisation process. 

* For this reason, you **MUST** have the ``ORCA_INPUT_FILENAME.interp`` file in your NEB calculation folder for ``viewORCA neb_snap`` to work.
* The prefix ``ORCA_INPUT_FILENAME`` is the same prefix you gave to your ORCA input file (being ``ORCA_INPUT_FILENAME.inp``). 

``viewORCA neb_snap`` will create two png files containing the energy profiles of the NEB calculation:

* ``ORCA_INPUT_FILENAME_NEB_optimization.png``: This plot displays the NEB energy profiles of all iterations performed during the NEB optimisation process. 
* ``ORCA_INPUT_FILENAME_NEB_last_iteration.png``: This plot shows the energy profile of the last iteration performed from the NEB optimisation process. 

Optional inputs for ``viewORCA neb_snap`` include:

* ``--path``: This is the path to the folder containing the ORCA NEB calculation. If not given, the default path is the current directory you are in. 
	* You can also provide the path to the ``.interp`` file you want to process if you want to process a specific ``.interp`` file of interest. This file must end with the ``.interp`` suffix. 

???+ example

	An example of a NEB optimisation calculation viewed using ``viewORCA neb_snap`` is shown below, along with the energy profile of the last iteration for this NEB calculation. This example is given in [Examples/viewORCA_neb_and_neb_snap](https://github.com/geoffreyweal/viewORCA/tree/main/Examples/viewORCA_neb_and_neb_snap)

	<figure markdown="span">
	    <img src="Figures/view_ORCA_Manual/neb_snap/orca_NEB_optimization.png?raw=true" alt="ORCA_INPUT_FILENAME_NEB_optimization.png" width="600"/>
	    <figcaption>Energy Profile of each NEB optimisation iteration performed by ORCA. This is given by the 'ORCA_INPUT_FILENAME_NEB_optimization.png' file.</figcaption>
	</figure>

	<figure markdown="span">
	    <img src="Figures/view_ORCA_Manual/neb_snap/orca_NEB_last_iteration.png?raw=true" alt="ORCA_INPUT_FILENAME_NEB_last_iteration.png" width="600"/>
	    <figcaption>Energy Profile of each NEB optimisation iteration performed by ORCA. This is given by the 'ORCA_INPUT_FILENAME_NEB_last_iteration.png' file.</figcaption>
	</figure>


## ``viewORCA irc``: View Intrinsic Reaction Coordinate (IRC) Calculation

This module is designed to allow the user to view the output of an ORCA IRC job visually. 

To use this module:

1. In the terminal, change directory into the directory containing your IRC job, and 
2. type ``viewORCA irc`` into the terminal. 

```bash
# First, change directory into the folder containing your ORCA IRC job.
cd path_to_geometric_optimisation_job

# Second, run "viewORCA irc".
viewORCA irc
```

``viewORCA irc`` will extract information from the ``ORCA_INPUT_FILENAME_IRC_Full_trj.xyz`` file and show the images of the IRC calculation, along with the energy profile of the IRC calculation. 

* For this reason, you **MUST** have the ``ORCA_INPUT_FILENAME_IRC_Full_trj.xyz`` file in your NEB calculation folder for ``viewORCA neb`` to work.
* The prefix ``ORCA_INPUT_FILENAME`` is the same prefix you gave to your ORCA input file (being ``ORCA_INPUT_FILENAME.inp``). 

This program will also create an xyz file called ``ORCA_INPUT_FILENAME_IRC_images.xyz`` that contains the images and the energy profile from the IRC calculation. 

* You can view this by typing ``ase gui ORCA_INPUT_FILENAME_IRC_images.xyz`` into the terminal

	```bash
	ase gui ORCA_INPUT_FILENAME_IRC_images.xyz
	```

Optional inputs for ``viewORCA irc`` include:

* ``--path``: This is the path to the folder containing the ORCA IRC calculation. If not given, the default path is the current directory you are in. 
* ``--view``: If tag indicates if you want to view the ASE GUI for the IRC job immediate after this program has run. By default this is ``True``. If you only want to create the ``SCAN_images.xyz`` file and not view it immediately, set this to ``False``.
	* Setting this to ``False`` is ideal if you are running ``viewORCA`` on an HPC, where ``ase gui`` windows may be more responsive if ``ase gui`` is run from your own computer rather than directly on the HPC. 
	* You can view ``ORCA_INPUT_FILENAME_IRC_images.xyz`` on your computer by typing the following into the terminal: ``ase gui ORCA_INPUT_FILENAME_IRC_images.xyz``.

???+ example

	An example of a IRC calculation viewed using ``viewORCA irc`` is shown below, along with the energy profile for this IRC calculation. This example is given in [Examples/viewORCA_irc](https://github.com/geoffreyweal/viewORCA/tree/main/Examples/viewORCA_irc)

	<figure markdown="span">
	<img src="Figures/view_ORCA_Manual/irc/IRC_example.gif?raw=true" alt="NEB Images from a IRC calculation" width="500"/>
	<figcaption>NEB Images from a IRC calculation</figcaption>
	</figure>

	The energy profile for this example is given below:

	<figure markdown="span">
	    <img src="Figures/view_ORCA_Manual/irc/IRC_energy_original.png?raw=true" alt="Energy Profile from a IRC calculation" width="600"/>
	    <figcaption>Energy Profile from a IRC calculation</figcaption>
	</figure>




