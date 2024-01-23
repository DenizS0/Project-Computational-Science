#
# Names: Dagmar Salomons & Deniz Saglam
# Student IDs: 13981714 & 13827189
#
# 3D Simulation of the solar system

import matplotlib.pyplot as plt
import numpy as np

G_CONSTANT = 6.67e-11 * (3.16e7)**2

fraction = 100000
dt = 1.0 / fraction

years = 1

hermite = False
animation = False



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
        self.axes.set_xlim(-150e9, 150e9)
        self.axes.set_ylim(-150e9, 150e9)
        self.axes.set_zlim(-150e9, 150e9)

    def make_body(self, new_body):
        self.bodies.append(new_body)

    def run_sim(self):
        if hermite == True:
            self.calculate_acceleration()
            self.calculate_velocity()
            self.calculate_position()
        else:
            self.calculate_position()
            self.calculate_acceleration()
            self.calculate_velocity()

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
            if hermite == True:
                body.hermite_pos()
            else:
                body.new_pos()

    def calculate_acceleration(self):
        for n, body in enumerate(self.bodies):
            body.acceleration(self.bodies, n)

    def calculate_velocity(self):
        for body in self.bodies:
            if hermite == True:
                body.hermite_vel()
            else:
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

        self.jerk = np.array([0, 0, 0])

    # calculate new position acc. to velocity verlet
    def new_pos(self):
        self.position = self.position + dt * self.velocity + .5 * self.accel * dt**2

    # calculate new velocity acc. to velocity verlet
    def new_vel(self):
        self.velocity = self.velocity + .5 * dt * (self.old_accel + self.accel)

    # calculate new acceleration acc. to a = gm/r^2
    def acceleration(self, bodies, n):
        self.old_accel = self.accel
        self.old_jerk = self.jerk
        self.accel = np.array([0, 0, 0])
        self.jerk = np.array([0, 0, 0])

        for i in range(len(bodies)):
            if i != n:
                force_vector = bodies[i].position - self.position
                distance = np.linalg.norm(force_vector)

                self.accel = self.accel + (G_CONSTANT * bodies[i].mass * force_vector) / distance**3

                if hermite == True:
                    vel_vector = bodies[i].velocity - self.velocity
                    self.jerk = self.jerk + G_CONSTANT * bodies[i].mass * (vel_vector / distance**3 - 3 * ((force_vector * np.dot(force_vector, vel_vector)) / distance**5))

    # calculate new position acc. to hermite
    def hermite_vel(self):
        self.old_velocity = self.velocity
        self.velocity = self.velocity + .5 * (self.old_accel + self.accel) * dt + 1/12 * (self.old_jerk - self.jerk) * dt**2
        #print(self.jerk)

    def hermite_pos(self):
        self.position = self.position + .5 * (self.velocity + self.old_velocity) * dt + 1/12 * (self.old_accel - self.accel) * dt**2

    def draw_body(self):
        self.system.axes.plot(self.position[0], self.position[1], self.position[2], marker="o", color=self.color, markersize = self.size)

#starting values taken from nasa horizon data set
def simulation_solar_system():
    solarsys = System()
    sun = Body(1.9885e30, np.array([2.0981e-3, -1.5479e-2, 7.8087e-5]) * 3.16e10, np.array([0,0,0]), np.array([-1.3512e9, -1.3976e7, 3.1592e7]), 25, "yellow", solarsys)
    mercury = Body(3.302e23, np.array([-1.0355e1, -4.6637e1, -2.8598e0]) * 3.16e10, np.array([0,0,0]), np.array([-5.9051e10, 2.5058e8, 5.3457e9]), 3, "burlywood", solarsys)
    venus = Body(48.685e24, np.array([3.6162e0, 3.464e1, 2.6736e-1]) * 3.16e10, np.array([0,0,0]), np.array([1.0657e11, -1.1721e10, -6.3566e9]), 5, "blanchedalmond", solarsys)
    earth = Body(5.9722e24, np.array([-2.5929e1, -1.5639e1, 1.1483e-3]) * 3.16e10, np.array([0,0,0]), np.array([-7.7974e10, 1.2570e11, 2.5858e7]), 5, "deepskyblue", solarsys)
    #moon = Body(0.073e24, np.array([0, np.cos(5.1 * np.pi/180) * 29.8, np.sin(5.1 * np.pi/180) * 29.8]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(5.1 * np.pi/180) * 150e9, 0, np.sin(5.1 * np.pi/180) * 150e9]), 1, "darkgray", solarsys)
    #mars = Body(0.642e24, np.array([0, np.cos(1.8 * np.pi/180) * 24.1, np.sin(1.8 * np.pi/180) * 24.1]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(1.8 * np.pi/180) * 228.0e9, 0, np.sin(1.8 * np.pi/180) * 228.0e9]), 4, "goldenrod", solarsys)
    #jupiter = Body(1898e24, np.array([0, np.cos(1.3 * np.pi/180) * 13.1, np.sin(1.3 * np.pi/180) * 13.1]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(1.3 * np.pi/180) * 778.5e9, 0, np.sin(1.3 * np.pi/180) * 778.5e9]), 12, "peru", solarsys)
    #saturn = Body(568e24, np.array([0, np.cos(2.5 * np.pi/180) * 9.7, np.sin(2.5 * np.pi/180) * 9.7]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(2.5 * np.pi/180) * 1432.0e9, 0, np.sin(2.5 * np.pi/180) * 1432.0e9]), 10, "khaki", solarsys)
    #uranus = Body(86.6e24, np.array([0, np.cos(0.8 * np.pi/180) * 6.8, np.sin(0.8 * np.pi/180) * 6.8]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(0.8 * np.pi/180) * 2867.0e9, 0, np.sin(0.8 * np.pi/180) * 2867.0e9]), 8, "paleturquoise", solarsys)
    #neptune = Body(102e24, np.array([0, np.cos(1.8 * np.pi/180) * 5.4, np.sin(1.8 * np.pi/180) * 5.4]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(1.8 * np.pi/180) * 4515.0e9, 0, np.sin(1.8 * np.pi/180) * 4515.0e9]), 7, "cyan", solarsys)

    for t in range(years * fraction):
        solarsys.run_sim()

        if animation == True:
            solarsys.draw_plot()
            plt.pause(0.001)

            for artist in plt.gca().get_lines():
                artist.remove()

    print(mercury.position)

    days = [i for i in range(years * fraction)]

    energy = [((solarsys.energy[i] - solarsys.energy[0]) / solarsys.energy[0]) for i in range(len(solarsys.energy))]
    plt.clf()
    plt.plot(days, energy, "b-")
    plt.title("Total Energy Difference")
    plt.show()

# n = number of stars, n = 1: simulation_solar_system
def simulation_multiple_star_system(n):
    mass_sun = 2e30

    if n < 1 or n > 3:
        print("Please choose n > 0 and n < 3")
        return None
    if n == 1:
        simulation_solar_system()
        return None
    if n == 2:
        starsys = System()
        distance = 40e9
        mass = mass_sun / 2

        v1 = np.sqrt((G_CONSTANT * mass * (distance/2)) / (distance)**2)
        v2 = -v1
        x = [distance/2, -(distance/2)]
        v = [v1, v2]

        for i in range(n):
            Body(mass, np.array([0,v[i],0]), np.array([0,0,0]), np.array([x[i],0,0]), 10, "yellow", starsys)

        mercury = Body(0.330e24, np.array([0, np.cos(7 * np.pi/180) * 47.4, np.sin(7 * np.pi/180) * 47.4]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(7 * np.pi/180) * 57.9e9, 0, np.sin(7 * np.pi/180) * 57.9e9]), 3, "burlywood", starsys)
        venus = Body(4.87e24, np.array([0, np.cos(3.4 * np.pi/180) * 35.0, np.sin(3.4 * np.pi/180) * 35.0]) * 3.16e10, np.array([0,0,0]), np.array([np.cos(3.4 * np.pi/180) * 108.2e9, 0, np.sin(3.4 * np.pi/180) * 108.2e9]), 5, "blanchedalmond", starsys)
        earth = Body(5.97e24, np.array([0, 29.8, 0]) * 3.16e10, np.array([0,0,0]), np.array([149.6e9, 0, 0]), 5, "deepskyblue", starsys)

        for t in range(1000000):
            starsys.run_sim()
            plt.pause(0.001)

            for artist in plt.gca().get_lines():
                artist.remove()
    if n == 3:
        starsys = System()

        mass = 1
        x = [-1, 1, 0]
        vx = [0.39295, 0.39295, -2 * 0.39295]
        vy = [0.09758, 0.09758, -2 * 0.09758]

        for i in range(n):
            Body(mass, np.array([vx[i], vy[i], 0]), np.array([0,0,0]), np.array([x[i], 0, 0]), 10, "yellow", starsys)

        for t in range(365):
            starsys.run_sim()
            plt.pause(0.001)

            for artist in plt.gca().get_lines():
                artist.remove()

simulation_multiple_star_system(1)