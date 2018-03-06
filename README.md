# Aileron Structural Simulation
Code for the structure assignment of AE3212-II SVV 2018 (TU Delft) of Group 1.
The assignment consists in analyzing the aileron of a Boeing 737 loaded with an actuator and aerodynamic loads. 

The code solves the static undetermined problem and computes the distribution of forces, moments and stresses in the structure.
Furthermore it is able to calculate the deflection, the slope and the curvature of the aileron.

## Overview
main.py

Boom.py  #It finds the position and the area of the booms

force.py  #It defines coefficients of the force for different uses, especially in the calculation of the internal forces

inertia.py  #It calculates the inertia from teh booms and it rotates it between different coordinate frames, body frame and the general one

Input_variables.py  #INPUT VARIABLES given for the Boeing 737

internal_forces.py  #Calculates internal forces and moments. It computes deflection, curvature, slope and normal stress of the structure

plots.py  #Generates plot of forces, moments, stresses, displacements, curvatures, slopes (Also in 3D) and it calculates Von Mises stresses

reaction_forces.py  #Calculates the reaction forces of the statically indeterminate structure using equation of motion and displacement relations

shear_force.py  #Calculates the shear flow in the structure, the shear center and the  aileron's twist

## Data

This folder contains the final report for this project developed by group 1 and the results.
The upward and downward rotation of the different validation cases is always 28 deg.

Main source of reference: Megson, T. H. G. (2013). Aircraft Structures for Engineering Students (5th ed.). Oxford: Elsevier Aerospace Engineering Series.

## Requirements
Python  2.7.13

matplotlib  2.1.2

plotly 2.4.1

numpy 1.13.1
