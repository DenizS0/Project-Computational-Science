#
# Names: Dagmar Salomons & Deniz Saglam
# Student IDs: 13981714 & 13827189
#
# 3D Simulation of the solar system

import matplotlib.pyplot as plt
import numpy as np

G_CONSTANT = 6.67e-11 * (3.16e7)**2
dt = 1.0 / 365

plt.style.use('dark_background')

"""
object representing solar system

contains all body objects

"""
class System:
    def __init__(self):
        self.bodies = []
        self.energy = []

        self.figure = plt.figure()
        self.axes = self.figure.add_subplot(111, projection='3d')
        self.axes.set_xlim(-550.0e10, 550.0e10)
        self.axes.set_ylim(-550.0e10, 550.0e10)
        self.axes.set_zlim(-550.0e10, 550.0e10)

    def make_body(self, new_body):
        self.bodies.append(new_body)

    def run_sim(self):
        self.calculate_position()
        self.calculate_acceleration()
        self.calculate_velocity()
        self.draw_plot()

        self.total_energy()

    def total_energy(self):
        total = 0

        for body in self.bodies:
            total += .5 * body.mass * np.linalg.norm(body.velocity)**2

            for i in range(len(self.bodies)):
                if self.bodies[i] != body:
                    total += (-0.5 * G_CONSTANT * body.mass * self.bodies[i].mass) / np.linalg.norm(self.bodies[i].position - body.position)

        self.energy.append(total)


    def draw_plot(self):
        self.bodies.sort(key=lambda item: item.position[0])
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
    def __init__(self, mass, velocity, accel, position, size, color, system):
        self.mass = mass
        self.velocity = velocity
        self.accel = accel
        self.position = position
        self.size = size
        self.color = color
        self.system = system
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
        self.system.axes.plot(self.position[0], self.position[1], self.position[2], marker="o", color=self.color, markersize = self.size)


solarsys = System()
sun = Body(2e30, np.array([0,0,0]), np.array([0,0,0]) * 3.16e7, np.array([0,0,0]), 25, "yellow", solarsys)
mercury = Body(0.330e24, np.array([0, np.cos(7 * np.pi/180) * 47.4, np.sin(7 * np.pi/180) * 47.4]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(7 * np.pi/180) * 57.9e9, 0, np.sin(7 * np.pi/180) * 57.9e9]), 3, "burlywood", solarsys)
venus = Body(4.87e24, np.array([0, np.cos(3.4 * np.pi/180) * 35.0, np.sin(3.4 * np.pi/180) * 35.0]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(3.4 * np.pi/180) * 108.2e9, 0, np.sin(3.4 * np.pi/180) * 108.2e9]), 5, "blanchedalmond", solarsys)
earth = Body(5.97e24, np.array([0, 29.8, 0]) * 3.16e10, np.array([0,0,0]), np.array([149.6e9, 0, 0]), 5, "deepskyblue", solarsys)
#moon = Body(0.073e24, np.array([0, np.cos(5.1 * np.pi/180) * 29.8, np.sin(5.1 * np.pi/180) * 29.8]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(5.1 * np.pi/180) * 150e9, 0, np.sin(5.1 * np.pi/180) * 150e9]), 1, "darkgray", solarsys)
mars = Body(0.642e24, np.array([0, np.cos(1.8 * np.pi/180) * 24.1, np.sin(1.8 * np.pi/180) * 24.1]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(1.8 * np.pi/180) * 228.0e9, 0, np.sin(1.8 * np.pi/180) * 228.0e9]), 4, "goldenrod", solarsys)
jupiter = Body(1898e24, np.array([0, np.cos(1.3 * np.pi/180) * 13.1, np.sin(1.3 * np.pi/180) * 13.1]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(1.3 * np.pi/180) * 778.5e9, 0, np.sin(1.3 * np.pi/180) * 778.5e9]), 12, "peru", solarsys)
saturn = Body(568e24, np.array([0, np.cos(2.5 * np.pi/180) * 9.7, np.sin(2.5 * np.pi/180) * 9.7]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(2.5 * np.pi/180) * 1432.0e9, 0, np.sin(2.5 * np.pi/180) * 1432.0e9]), 10, "khaki", solarsys)
uranus = Body(86.6e24, np.array([0, np.cos(0.8 * np.pi/180) * 6.8, np.sin(0.8 * np.pi/180) * 6.8]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(0.8 * np.pi/180) * 2867.0e9, 0, np.sin(0.8 * np.pi/180) * 2867.0e9]), 8, "paleturquoise", solarsys)
neptune = Body(102e24, np.array([0, np.cos(1.8 * np.pi/180) * 5.4, np.sin(1.8 * np.pi/180) * 5.4]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(1.8 * np.pi/180) * 4515.0e9, 0, np.sin(1.8 * np.pi/180) * 4515.0e9]), 7, "cyan", solarsys)

for t in range(365):
    solarsys.run_sim()
    plt.pause(0.001)

    for artist in plt.gca().get_lines():
        artist.remove()

days = [i for i in range(365)]

energy = [((solarsys.energy[i] - solarsys.energy[0]) / solarsys.energy[0]) for i in range(len(solarsys.energy))]
plt.clf()
plt.plot(days, energy, "b-")
plt.show()
