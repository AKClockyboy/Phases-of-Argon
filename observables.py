"""
    External code from VV_Particle.py with observable functions, includes:
    Radial distribution
    Mean Squared Displacement (MSD) and MSD Grapher
    Energy Grapher
    Minimum Image Convention

Authors: Alexander Kintrea(B105147) & Gonzalo Gil(B101878)
"""

#importing modules and packages
import sys
import math
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as pyplot
import MDU
from Particle3D import Particle3D

def rdf(particles,box_size,rho,n_bins):
    """
    Radial Distribution function. Obtains particle position data to Calculate
    the probability of particle positions in relation to a reference particle

    -Inputs are: List of Particle3D objects
                 Box size as an array
                 Reduced Box density
                 Number of bins to use for RDF plot
    -Returns: List to use for RDF histogram
              RDF Values
    """
    sep = []
    for i in range(len(particles)):
        for j in range(i+1, len(particles)):
            sep.append(LA.norm(MIC((particles[i].position - particles[j].position),box_size)))
    dr = (np.amax(sep)-np.amin(sep))/n_bins # Calculating bin width
    n = len(particles)
    rdf, r_list = np.histogram(sep, bins=n_bins)
    r_listC  = r_list[:-1]+dr/2 # Removing final element and adding half bin width
    rdf = rdf/(4*(math.pi)*rho*n*dr*(r_listC**2)) # applying normalising factor

    return r_listC, rdf

def radial_D_grapher(radial_list,numstep):
    """
    Graphs the radial distribution function of particles in the system

    -Inputs are: List of RDF Values
                 Number of timesteps in simulation
    -Returns: Matplotlib plot of RDF
    """

    pyplot.plot(radial_list, numstep, "r-")
    pyplot.xlabel("Distance")
    pyplot.ylabel("rdf(r)")
    pyplot.title("Radial Distribution Function")
    pyplot.show()

def msd(particles, particles_0, box_size, N):
    """
    Calculates the mean squared displacement of particles in relation to their
    initial positions. Uses MIC and appends to a list to then plot.

    -Inputs are: List of Particle3D objects
                 List of Particle3D objects in initial positions
                 Box size as an array
                 Number of particles
    -Returns: Mean Squared displacement Values for all particles
    """
    msd_sum = 0
    for line in range(len(particles)):
        msd_sum += LA.norm(MIC(particles[line].position - particles_0[line].position,box_size))
    msd = msd_sum/N

    return msd

def msd_grapher(time_list,msd_list):
    """
    Graphs mean squared displacement of particles over time

    -Inputs are: List of timesteps
                 List of MSD values
    -Returns: Matplotlib plot of MSD
    """
    pyplot.title('MSD')
    pyplot.xlabel('Time')
    pyplot.ylabel('Particle Separation')
    pyplot.plot(time_list, msd_list, "r-")
    pyplot.show()

def energy_grapher(time_list, poteng, kineng, energy):
    """
    Graphs energies over time

    -Inputs are: List of timesteps
                 List of Potential energy Values
                 List of Kinetic Energy Values
                 List of Total Energy values
    -Returns: Matplotlib plot of Energy over time
    """

    pyplot.title('Energies of the system')
    pyplot.xlabel('Time')
    pyplot.ylabel('Energy')
    pyplot.plot(time_list, energy, "g-", label = "Total Energy")
    pyplot.plot(time_list, kineng, "r-", label = "Kinetic Energy")
    pyplot.plot(time_list, poteng, "b-", label = "Potential Energy")
    pyplot.legend()
    pyplot.show()

def MIC(r, L):
    """
    Use the Minimum Image Convention method on particle positions

    -Inputs are: Particle3D Particle Positions
                 Box Size as an array

    -Returns: Particle position adhering to the Minimum Image Convention
    """
    position = np.mod(r,L)
    O = np.zeros(len(position))
    for i in range(len(position)):
        for j in range(len(position)):
            if i == j:
                if (position[i] > L[j]/2): #check particle position
                    O[i] = np.mod(position[i],L[j]/2)-L[j]/2 #project image of particle
                else:
                    O[i] = position[i] #project image of particle
    return(O)
