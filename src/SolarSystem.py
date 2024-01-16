#
# Names: Dagmar Salomons & Deniz Saglam
# Student IDs: 13981714 & 13827189
#
# 3D Simulation of the solar system

import matplotlib.pyplot as plt
import numpy as np

G_CONSTANT = 4.0 * np.pi**2
dt = 1.0 / 365

"""
object representing solar system

contains all body objects

"""
class System:
    def __init__(self):
        self.bodies = []
        self.figure = plt.figure()

    def make_body(self, new_body):
        self.bodies.append(new_body)

    def run_sim(self):
        self.calculate_position()
        self.calculate_acceleration()
        self.calculate_velocity()
        self.draw_plot()

    def draw_plot(self):
        self.axes = plt.axes(projection = '3d')
        self.axes.set_xlim(-6, 6)
        self.axes.set_ylim(-6, 6)
        for body in self.bodies:
            body.draw_body()

    def calculate_position(self):
        for body in self.bodies:
            body.new_pos()

    def calculate_acceleration(self):
        for n, body in enumerate(self.bodies):
            body.acceleration(self.bodies, n)

    def calculate_velocity(self):
        for body in self.bodies:
            body.new_vel()

"""
Creates body object

Args:
    __init__():
        mass (float): mass of the body
        velocity (ndarray): starting velocity of the body
        accel (ndarray): starting acceleration of the body
        position (ndarray): starting position of the body
        system (class System): system of the body
    acceleration():
        bodies (list): all bodies in the system
        n (int): which body is being calculated

"""
class Body:
    def __init__(self, mass, velocity, accel, position, system):
        self.mass = mass
        self.velocity = velocity
        self.position = position
        self.system = system
        self.accel = accel
        self.system.make_body(self)

    # calculate new position acc. to velocity verlet
    def new_pos(self):
        self.position = self.position + dt * self.velocity + .5 * self.accel * dt**2

    # calculate new acceleration acc. to a = gm/r^2
    def acceleration(self, bodies, n):
        self.old_accel = self.accel
        self.accel = np.array([0, 0, 0])

        for i in range(len(bodies)):
            if i != n:
                force_vector = bodies[i].position - self.position
                distance = np.linalg.norm(force_vector)

                self.accel = self.accel + (G_CONSTANT * bodies[i].mass * force_vector) / distance**3

    # calculate new velocity acc. to velocity verlet
    def new_vel(self):
        self.velocity = self.velocity + .5 * dt * (self.old_accel + self.accel)

    def draw_body(self):
        self.system.axes.plot(self.position[0], self.position[1], self.position[2], marker="o")


solarsys = System()
earth = Body(3.0e-6, np.array([0, 6.279, 0]), np.array([0,0,0]), np.array([1,0,0]), solarsys)
sun = Body(1, np.array([0,0,0]), np.array([0,0,0]), np.array([0,0,0]), solarsys)

for t in range(365):
    solarsys.run_sim()
    plt.pause(.001)
    plt.clf()
