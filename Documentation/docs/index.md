# ``viewORCA``: A Program for Viewing ORCA Calculations

This program is designed to allow the user to view ORCA calculations. Currently, viewORCA allows the user to view:

* Geometric Optimisation Calculations: View the local optimisation of an ORCA calculation, along with associated energy profile of the optimisation. 
* SCAN Calculations: View the geometrically converged states of the images of a SCAN calculation, along with associated energy profile of the SCAN calculation. 
* Nudged Elastic Band (NEB) Calculations: View the last NEB profile, along with the energy profile of the NEB calculation.
* Intrinsic Reaction Coordinate (IRC) Calculations: View the last NEB profile, along with the energy profile of the NEB calculation.


## ORCA

I have written this procedure for: 

* ORCA 5.0.3
* ORCA 5.0.4

This method should be valid for future versions of ORCA, but just in case any problems occur it may be due to ORCA version issues. 


## Installing the ``viewORCA`` Program

See [Installation: Setting Up the ``viewORCA`` Program and Pre-Requisites Packages](Installation.md) for more information. 

* Note that you can install the ``viewORCA`` program through ``pip3`` and ``conda``. 


## Instructions for using ``viewORCA``

This program is fundamental to the [Procedure for Investigating Reaction Mechanisms in ORCA](https://geoffreyweal.github.io/ORCA_Mechanism_Procedure) tutorial. Click on this tutorial to see how ``viewORCA`` can be used to view ORCA calculations. 

* Also, [click here to access the ``viewORCA`` manual](viewORCA_Manual.md).


## Examples for using ``viewORCA``

Examples of the various modules as used for visualising the various ORCA calculation types can be found in [the Examples folder](https://github.com/geoffreyweal/viewORCA/tree/main/Examples) in Github. 


## Questions and Feedback

If you have any questions about using ``viewORCA``, or would like to provide feedback, click on the [``Issues`` page here](https://github.com/geoffreyweal/viewORCA/issues), then click on the ``New issue`` button, and write you question/give you feedback. 

Any feedback would be amazing and appreciated! 

Thanks!


## Acknowledgements

The ``viewORCA neb_snap`` module is an adaptation of [neb_snapshots.py](https://github.com/via9a/neb_visualize.py) by Vilhjálmur Ásgeirsson and Benedikt Orri Birgirsson of the Universiy of Iceland. 

