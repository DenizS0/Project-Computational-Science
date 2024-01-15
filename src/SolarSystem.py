#
# Names: Dagmar Salomons & Deniz Saglam
# Student IDs: ... & 13827189
#
# 3D Simulation of the solar system

import matplotlib.pyplot as plt
import numpy as np

G_CONSTANT = 6.6743e-11

"""
object representing solar system

contains all body objects

"""
class System:
    def __init__(self):
        self.bodies = []

        self.axis = plt.axes(projection='3d')

    def make_body(self, new_body):
        self.bodies.append(new_body)

    def calculate_acceleration(self):
        for n, body in enumerate(len(self.bodies)):
            body.acceleration(self.bodies, n)

"""
Creates body object

Args:
    __init__():
        mass (float): mass of the body
        velocity (ndarray): starting velocity of the body
        position (ndarray): starting position of the body
        system (class System): system of the body
    acceleration():
        bodies (list): all bodies in the system
        n (int): which body is being calculated

"""
class Body:
    def __init__(self, mass, velocity, position, system):
        self.mass = mass
        self.velocity = velocity
        self.position = position
        self.system = system

        self.system.make_body(self)

    # calculate acceleration according to a = gm/r^2
    def acceleration(self, bodies, n):
        self.accel = np.array([0, 0, 0])

        for i in range(len(bodies)):
            if i != n:
                force_vector = bodies[i].position - self.position
                distance = np.linalg.norm(force_vector)

                self.accel += (G_CONSTANT * bodies[i].mass * force_vector) / distance**3



# System()
# plt.show()