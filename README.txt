README FILE

                 Lennard Jones Particle Interaction Simulation
  ------------Authors: Alexander Kintrea(B105147) & Gonzalo Gil(B101878)------------

Welcome, dearest marker! This Code simulates a Lennard Jones (LJ) Particle Potential using
Python and VMD for Argon in Solid, Liquid, Gas states - relevent 14/03/2019.

------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

The main code used is VV_Particle.py, which allows for usage from the command
line using the format:

python3 VV_Particle.py input_file.txt trajectory.xyz

------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

VV_Particle: The Main Code, containing functions from a Particle3D.py class,
MDUtilities.py for initial particle conditions in a lattice, and an
observables.py module with observable data from the LJ interaction.

This code also prints the energy data, Mean Squared Displacement (MSD) data, and
radial Distribution Function (RDF) data into libre office/excel files for the user to
analyse plots of Energy, MSD, and RDF are plotted with matplotlib. We have included a
folder of examples of our data files.

We use Minimum Image Convention and Periodic Boundary Conditions to make the computing
process smoother, while keeping the simulation accurate.

------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

input_file: Input file that contains simulation parameters, containing
reduced density, number of particles, and reduced temperature. Usage is
reduced density, number of particles, reduced Temperature and cutoff distance.

We have already set these parameters in files named:

argon_gas.txt,
argon_liquid.txt,
argon_solid.txt

...for the three different states of matter.

Feel free to adjust these files at your leisure.

------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

trajectory.xyz: Output file to be read by VMD, a simulation visualizer.

------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

Particle3D: This is the same as the previous checkpoint - we didn't change it at all.

------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

Acknowledgments:


Dr M Martinez-Canales and Dr A Hermann for organising and lecturing the course
Humphrey, W., Dalke, A. and Schulten, K., ``VMD - Visual Molecular Dynamics'' J. Molec. Graphics 1996, 14.1, 33-38
Every tutor and colleague who helped us along the way !

------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

Have a nice day - and hope you have fun with our code.
