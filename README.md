# DoubleRQ
Data for the paper "An Efficient Approximation to the Pull Policy for Hybrid Manufacturing and Remanufacturing Systems with Setup Costs".

Since this has not yet been published, I will not put who the co-authors are.



###############################################################################################################
File Name: Test for Validity of Approximation Results_Github Upload.xlsx

Outlines the results of the numerical studies and approximations presented in the paper

Section Headers:

Problem Parameters - Simulation and Approximation Parameters

Obtained from Modified Federgruen and Zheng	 - Double RQ parameters obtained when using the Modified Federgruen and Zheng Algorithm

Obtained from Simulation-Optimization	 - pull policy parameters obtained when using the Simulation-Optimization Method

vdl and Teunter Approximated Parameters	- pull policy parameters obtained when using van der Laan and Teunter's (2006) heuristic/approximation

Results and Comparison - Percent Cost Difference between the described approaches

###############################################################################################################
File Name: Test for Validity of Approximation Results_Github Upload

Contains the results of the Statistical Comparison between the Double (r, Q) Approximation and Procedure 1 in the paper


###############################################################################################################
CODE: Files can be found in Double RQ Code.zip

Input Files:

RandomNumbers.csv and SampleParameters.csv

These file names can be changed, but it also has to be changed in the code

Initialization:

There are two initialization files:

1.) Simulation_Initilization.py

This file runs the numerical study which calculates the Double RQ approximation and the Simulation-optimization method.

2.) Initialization.py

This file runs the numerical study which calculates the van der Laan and Teunter (2006) estimate and compares it with the 2RQ and v2RQ policies

