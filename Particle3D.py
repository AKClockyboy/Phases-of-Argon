"""
Particle3D, a class to describe the mass, position and velocity of 3D particles.
This code contains functions to define and arrange particles in a usable format
"""
import numpy as np
import scipy
import math

class Particle3D(object):


    def __init__(self, pos, vel, mass, label):
        """
        Initialise a Particle3D instance
        :param pos: position as np array
        :param vel: velocity as magnitude of np array
        :param mass: mass as float
        """

        self.position = pos
        self.velocity = vel
        self.mass = mass
        self.label = label



    def __str__(self):

        """
        Defining an output format.
        For particle p=(2.0, 0.5, 1.0) this will print as
        "x = 2.0, v = 0.5, m = 1.0"

        I don't actually use this, but thought I'd keep it in case I want to use it in the project next semester!
        """

        return str(self.position) + " " + str(self.velocity) + " " + str(self.mass) + " " + str(self.label)


    def kinetic_energy(self):
        """
        Function to return kinetic energy as
        1/2*mass*vel^2
        """

        return 0.5*self.mass*(np.linalg.norm(self.velocity))**2

    # Time integration methods

    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)/m

        :param dt: timestep as float
        :param force: force on particle as float
        """

        self.velocity = self.velocity + dt*force/self.mass


    def leap_pos1st(self, dt):

        """
        Defining the leap position for the symplectic Euler time integration methods
        => x(t+dt) = x(t) + dt*v(t)
        """

        self.position = self.position + dt*(self.velocity)


    def leap_pos2nd(self, dt, force):

        """
        Defining the leap position for the Velocity Verlet time integration methods
        => x(t+dt) = x(t) + dt*v(t) + 0.5*dt^2*F(t)/m
        """
        self.position = self.position + dt*(self.velocity) + 0.5*dt**2*force/self.mass


    @staticmethod
    def separation(p1,p2):
        """
        Static method for returning separation of two particles given initial parametres
        Defining the difference in position between two particles as:
        (r12) = r1-r2
        """
        r1 = p1.position
        r2 = p2.position
        sep = r1 - r2

        return sep


    @staticmethod
    def posdif(p1, p2):
        """
        Static method for returning the magnitude of separation between two particles given initial parametres
        """
        r1 = p1.position
        r2 = p2.position
        r = np.linalg.norm((r1 - r2))

        return r


    @staticmethod

    def from_file(file_handle, m):
        """
        Static Method to obtain information about particles position, velocity, and mass from a given file_handle
        """
        file_handle = file_handle.readline()
        tokens = file_handle.rsplit(" ")
        pos = np.array([tokens[0],tokens[1],tokens[2]],float)
        vel = np.array([tokens[3],tokens[4],tokens[5]],float)
        mass = m
        return Particle3D(pos, vel, mass)
