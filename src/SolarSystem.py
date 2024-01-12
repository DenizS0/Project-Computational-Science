#
# Names: Dagmar Salomons & Deniz Saglam
# Student IDs: ... & 13827189
#
# 3D Simulation of the solar system

import matplotlib.pyplot as plt
import numpy as np

G_CONSTANT = 6.6743e-11

# object representing solar system
# contains all other objects
class System:
    def __init__(self):
        self.suns = []
        self.planets = []

    def make_sun(self, new_sun):
        self.suns.append(new_sun)

    def make_planet(self, new_planet):
        self.planets.append(new_planet)

"""
Creates sun object

Args:
    mass (float) = mass of the sun
    velocity (float) = starting velocity of the sun
    position (np.ndarray) = starting position of the sun
    system (class System) = system of the sun

"""
class Sun:
    def __init__(self, mass, velocity, position, system):
        self.mass = mass
        self.velocity = velocity
        self.position = position
        self.system = system

"""
Creates planet object

Args:
    mass (float) = mass of the planet
    velocity (float) = starting velocity of the planet
    position (np.ndarray) = starting position of the planet
    system (class System) = system of the planet

"""
class Planet:
    def __init__(self, mass, velocity, position, system):
        self.mass = mass
        self.velocity = velocity
        self.position = position
        self.system = system