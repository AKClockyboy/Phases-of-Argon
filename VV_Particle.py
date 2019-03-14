"""
Computer Modelling Project:
Velocity Verlet time integration of a user defined
number of particles interacting via a Lennard-Jones Potential.
Produces plots of the radial distribution of the system, particle Mean Squared
Displacement and its energy fluctuations.
Also saves information obtained to extrnal files for analysis in the report.

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
import observables as obs

def force_LJP(particles, L, rc):
    """
    Function to return force on a particle given by the Lennard Jones Potential:
    F = 48[(r^-12) - (r^6)](ri-rj) using reduced units and a cutoff distance

    -Inputs are: List of Particle3D objects
                 Box size as an array
                 cutoff distance
    -Returns: Force on particles from other particles in system

    """
    natoms = len(particles)
    LJF = np.zeros((natoms,natoms,3))
    for i in range(natoms):
        for j in range(i+1, natoms):
            if i != j:
                rij = particles[i].position - particles[j].position
                r_mic = (obs.MIC(rij, L))
                r_mag = LA.norm(r_mic)

                if r_mag < rc:
                    LJF[i,j] += -48*(((1/r_mag**14)-(1/(2*r_mag**8)))*r_mic)
                    LJF[j,i] -= LJF[i,j]

    return LJF

def pot_energy_LJP(particles, L, rc):
    """
    Function to return potential between N particles:
    U = 4(r^-12 - r^-6))

    -Inputs are: List of Particle3D objects
                 Box size as an array
                 cutoff distance
    -Returns: Potential energy on particles from other particles in system
    """
    natoms = len(particles)
    LJP = 0
    for i in range(natoms):
        for j in range(i+1, natoms):
            r = particles[i].position - particles[j].position
            r_mic = (obs.MIC(r, L))
            r_mag = LA.norm(r_mic)
            if r_mag < rc:
                LJP += 4*(1/(r_mag**12) - 1/(r_mag**6))
    return LJP

def initial_conditions(rho, particles, particles_0,temp):
    """
    Sets initial conditions of particle lists using MDU. Also, makes initial
    list of particles that remains constant for the entire simulation, used
    for the Mean Squared Displacement calculations.

    -Inputs are: Reduced box density
                 List of Particle3D objects
                 List of Particle3D objects in initial positions
                 Reduced Temperature
    -Returns: Initial conditions using MDUtilities
    """
    MDU_positions = MDU.set_initial_positions(rho,particles)
    MDU_positions_0 = MDU.set_initial_positions(rho,particles_0)
    MDU_velocities = MDU.set_initial_velocities(temp,particles)

    return MDU_positions, MDU_positions_0, MDU_velocities

def main():
    # Read name of output file from command line
    if len(sys.argv)!=3:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + "param.firstinput" + " traj.xyz")
        quit()
    else:
        infile_name = sys.argv[1]
        output_file = sys.argv[2]

    #Read simulation parameters and constants (rho, No. particles, T, r_cutoff)
    infile = open(infile_name,"r")
    infile = np.loadtxt(infile)
    print("The simulation parameters, rho, number of particles, temperature and cut-off distance, are: ", infile)
    rho = float(infile[0])
    N = int(infile[1])
    temp = float(infile[2])
    rc = float(infile[3])
    n_bins = 50 #number of bins to plot radial distribution histogram later

    #open empty particle lists
    particles = []
    particles_0 = []

    #make initial particle objects to be used by molecular dynamics module
    for i in range(N):
        particles.append(Particle3D(np.array([0,0,0]), np.array([0,0,0]), 1, i))
        particles_0.append(Particle3D(np.array([0,0,0]), np.array([0,0,0]), 1, i))


    #Use MDUtilities to set box size, intial positions and velocities
    #also make inital position particle list to be used for MSD
    box_size = MDU.set_initial_positions(rho,particles)
    initial_conditions(rho, particles, particles_0,temp)
    # Set up numstep time parameters, using reduced units - you can change these freely!
    dt = 0.001
    numstep = 5000
    time = 0

    #Make initial energy parameters, to be updated in time loop
    init_kineng = [particles[i].kinetic_energy() for i in range(len(particles))]
    init_poteng = (pot_energy_LJP(particles,box_size,rc))
    init_total_kineng = sum(init_kineng)
    init_energy = init_total_kineng + init_poteng

    #initial force value
    force = force_LJP(particles,box_size,rc)

    #initialise lists that will be appended to in time loop
    time_list = [0.0]
    msd_list = []
    msd_timelist = []
    poteng = [init_poteng]
    kineng = [init_total_kineng]
    energy = [init_energy]

    #Open files for writing
    f = open(output_file,'w') #VMD File
    g = open("energy_file.csv",'w') #Pot_energy, Kin_Energy and Tot_Energy File
    h = open("RDF.csv",'w')
    e = open("MSD.csv",'w')

    #Titles and initial Values for Files
    g.write("time, KEng, PEng, TotEng \n" )
    g.write("%s, %s, %s, %s \n" % ('%.5f'% time, kineng[-1], poteng[-1], energy[-1]))
    h.write(str("separation, rdf \n"))
    e.write(str("time, msd \n"))


    # Start the time integration loop
    for t in range(numstep):

        # Update particle position. This uses force from the previous timestep
        #second line in for loop uses Periodic Boundary Conditions
        for line in range(len(particles)):
            particles[line].leap_pos2nd(dt, np.sum((force),axis=0)[line])
            particles[line].position = np.mod(particles[line].position, box_size)

        #New force value
        force_new = force_LJP(particles,box_size,rc)

        #Update particle velocity with Velocity Verlet Method
        for line in range(len(particles)):
            (particles[line].leap_velocity(dt, 0.5*np.sum((force+force_new),axis=0)[line]))

        #Update time and force
        time += dt
        time_list.append(time)
        force = force_new

        #Cute percentge gauge that I saw someone else do - thought it was a great addition!
        print("\033[95m Completion : " + str(round(100*time/(numstep*dt),5)) + '%\033[0m', end= '\r')

        #Update particle energies and append to lists. Use Particle3D's kinetic
        #energy function and the internal potential energy function
        kineng1 = [particles[i].kinetic_energy() for i in range(len(particles))]
        poteng1 = (pot_energy_LJP(particles, box_size,rc))
        kineng_sum = sum(kineng1)
        poteng_sum = poteng1
        kineng.append(kineng_sum)
        poteng.append(poteng_sum)
        energy.append(kineng_sum + poteng_sum)

        #Write energy Values to file
        g.write("%s, %s, %s, %s \n" % ('%.5f'% time, kineng[-1], poteng[-1], energy[-1]))

        #Calculating radial distribution and write to file
        if t%10:
            r_listC, rdf = obs.rdf(particles,box_size,rho,n_bins)
            for i in range(len(particles)):
			             h.write("%s, %s, \n" % (r_listC[i], rdf[i]))#... and writing to file

        #Calculating values for mean squared distribution
        if t%10:
            msd = obs.msd(particles, particles_0, box_size, N)
            msd = msd*100/numstep
            msd_list.append(msd)
            msd_timelist.append(time)
            e.write("%s, %s, \n" % ('%.5f'% time, msd_list[-1])) # ...and writing to file

        #Writing to a .xyz file for VMD to read
        f.write(str(N)+"\n")
        f.write("Point = " + str(t+1) +"\n")
        for i in range(len(particles)):
            f.write("s"+ str(i) + " " + str(particles[i].position[0]) + " " + str(particles[i].position[1]) + " " + str(particles[i].position[2]) + "\n")

    #Sending all info to be graphed
    obs.energy_grapher(time_list,poteng,kineng,energy)
    obs.msd_grapher(msd_timelist,msd_list)
    obs.radial_D_grapher(r_listC,rdf)
main()
