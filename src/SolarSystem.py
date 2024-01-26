#
# Names: Dagmar Salomons & Deniz Saglam
# Student IDs: 13981714 & 13827189
#
# 3D Simulation of the solar system

import matplotlib.pyplot as plt
import numpy as np
from lmfit import models

# gravitational constant adjusted for years
G_CONSTANT = 6.67e-11 * (3.16e7)**2

# 1/fraction decides the size of dt
fraction = 1000
dt = 1.0 / fraction

# amount of years the simulation simulates
years = 1

# size of one side in the box in which the animation plays
axis_size = 10e13

# True = hermite integration; False = velocity verlet integration
hermite = True
# show animation of system (very slow for small dt)
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

        if animation == True:
            self.figure = plt.figure()
            self.axes = self.figure.add_subplot(111, projection='3d')
            self.axes.set_xlim(-axis_size, axis_size)
            self.axes.set_ylim(-axis_size, axis_size)
            self.axes.set_zlim(-axis_size, axis_size)

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
        self.predict_acc = self.old_accel + self.old_jerk * dt
        self.velocity = self.velocity + .5 * (self.old_accel + self.predict_acc) * dt + 1/12 * (self.old_jerk - self.jerk) * dt**2
        #print(self.jerk)

    def hermite_pos(self):
        self.predict_v = self.velocity + self.accel * dt + 0.5 * self.old_jerk * dt**2
        self.position = self.position + .5 * (self.velocity + self.predict_v) * dt + 1/12 * (self.old_accel - self.predict_acc) * dt**2

    def draw_body(self):
        self.system.axes.plot(self.position[0], self.position[1], self.position[2], marker="o", color=self.color, markersize = self.size)

# starting values taken from nasa horizon data set @ 22 jan 2023 00:00
def simulation_solar_system():
    solarsys = System()
    sun = Body(1.9885e30, np.array([2.0981e-3, -1.5479e-2, 7.8087e-5]) * 3.16e10, np.array([0,0,0]), np.array([-1.3512e9, -1.3976e7, 3.1592e7]), 25, "yellow", solarsys)
    mercury = Body(3.302e23, np.array([-1.0355e1, -4.6637e1, -2.8598e0]) * 3.16e10, np.array([0,0,0]), np.array([-5.9051e10, 2.5058e8, 5.3457e9]), 3, "burlywood", solarsys)
    venus = Body(48.685e24, np.array([3.6162e0, 3.464e1, 2.6736e-1]) * 3.16e10, np.array([0,0,0]), np.array([1.0657e11, -1.1721e10, -6.3566e9]), 5, "blanchedalmond", solarsys)
    earth = Body(5.9722e24, np.array([-2.5929e1, -1.5639e1, 1.1483e-3]) * 3.16e10, np.array([0,0,0]), np.array([-7.7974e10, 1.2570e11, 2.5858e7]), 5, "deepskyblue", solarsys)
    moon = Body(7.349e22, np.array([-2.5005e1, -1.5036e1, -9.1161e-3]) * 3.16e10, np.array([0,0,0]), np.array([-7.7780e10, 1.2541e11, -5.0425e6]), 1, "darkgray", solarsys)
    mars = Body(6.4171e23, np.array([-2.3074e1, -1.3303e0, 5.3856e-1]) * 3.16e10, np.array([0,0,0]), np.array([-3.4426e10, 2.3527e11, 5.7741e9]), 4, "goldenrod", solarsys)
    jupiter = Body(1.8982e27, np.array([-3.3404e0, 1.3283e1, 1.9609e-2]) * 3.16e10, np.array([0,0,0]), np.array([7.1676e11, 1.8057e11, -1.6785e10]), 12, "peru", solarsys)
    saturn = Body(5.6834e26, np.array([4.7786e0, 8.0455e0, 3.3070e-1]) * 3.16e10, np.array([0,0,0]), np.array([1.2262e12, -8.0900e11, -3.4756e10]), 10, "khaki", solarsys)
    uranus = Body(86.81e24, np.array([-5.0654e0, 4.2894e0, 8.1455e-2]) * 3.16e10, np.array([0,0,0]), np.array([1.9896e12, 2.1661e12,-1.7730e10]), 8, "paleturquoise", solarsys)
    neptune = Body(102.4e24, np.array([4.8688e-1, 5.4420e0,-1.2344e-1]) * 3.16e10, np.array([0,0,0]), np.array([4.4517e12, -4.3035e11, -9.3732e10]), 7, "cyan", solarsys)

    for t in range(years * fraction):
        solarsys.run_sim()

        if animation == True:
            solarsys.draw_plot()
            plt.pause(0.001)

            for artist in plt.gca().get_lines():
                artist.remove()

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

        for t in range(years * fraction):
            starsys.run_sim()

            if animation == True:
                starsys.draw_plot()
                plt.pause(0.001)

                for artist in plt.gca().get_lines():
                    artist.remove()
        days = [i for i in range(years * fraction)]

        energy = [((starsys.energy[i] - starsys.energy[0]) / starsys.energy[0]) for i in range(len(starsys.energy))]
        plt.clf()
        plt.plot(days, energy, "b-")
        plt.title("Total Energy Difference")
        plt.show()

    if n == 3:
        starsys = System()

        mass = [2.4e30, 1.6e30, 0.246e30] # alpha centauri A, B en proxima
        pm_ra = np.array([-3.710, -3.6, -3.781]) * np.pi/180 * 1/3600 # proper motion right ascension
        pm_dec = np.array([0.482, 0.952, 0.770]) * np.pi/180 * 1/3600 # proper motion declination
        distance = np.array([1.347, 1.347, 1.301]) * 3.0857e16
        ra = np.array([14/24 + 39/1440 + 40.90/86400, 14/24 + 39/1440 + 39.39/86400, 14/24 + 29/1440 + 34.43/86400]) * np.pi/180 # (14h39m40.90s) 
        dec = np.array([-(60/24 + 50/1440 + 6.53/86400), -(60/24 + 50/1440 + 22.10/86400), -(62/24 + 40/1440 + 34.26/86400)]) * np.pi/180
        radial_v = np.array([-21.4, -18.6, -22.204]) * 3.16e10 # radial velocity

        v_alpha = [(pm_ra[i] * distance[i]) for i in range(n)]
        v_delta = [(pm_dec[i] * distance[i]) for i in range(n)]

        vx = [((v_alpha[i] + radial_v[i]) * np.cos(dec[i]) * np.cos(ra[i])) for i in range(n)]
        vy = [((v_alpha[i] + radial_v[i]) * np.cos(dec[i]) * np.sin(ra[i])) for i in range(n)]
        vz = [(v_delta[i] + radial_v[i] * np.sin(dec[i])) for i in range(n)]

        x = [(distance[i] * np.cos(dec[i]) * np.cos(ra[i])) for i in range(n)]
        y = [(distance[i] * np.cos(dec[i]) * np.sin(ra[i])) for i in range(n)]
        z = [(distance[i] * np.sin(dec[i])) for i in range(n)]

        cm_x = np.dot(mass, x) / np.sum(mass)
        cm_y = np.dot(mass, y) / np.sum(mass)
        cm_z = np.dot(mass, z) / np.sum(mass)

        for i in range(n):
            Body(mass[i], np.array([vx[i], vy[i], vz[i]]), np.array([0,0,0]), np.array([x[i] - cm_x, y[i] - cm_y, z[i] - cm_z]), 5, "yellow", starsys)
            print(f"{x[i] - cm_x, y[i] - cm_y, z[i] - cm_z}")

        for t in range(years * fraction):
            starsys.run_sim()

            if animation == True:
                starsys.draw_plot()
                plt.pause(0.001)

                for artist in plt.gca().get_lines():
                    artist.remove()
        #days = [i for i in range(years * fraction)]

        #energy = [((starsys.energy[i] - starsys.energy[0]) / starsys.energy[0]) for i in range(len(starsys.energy))]
        #plt.clf()
        #plt.plot(days, energy, "b-")
        #plt.title("Total Energy Difference")
        #plt.show()


# starting values taken from nasa horizon data set @ 22 jan 2023 00:00
# mass sun = 1.9885e30
def make_solar_system(solarsys, x_value, mass_sun, v):
    sun = Body(mass_sun, np.array([2.0981e-3, -1.5479e-2, 7.8087e-5]) * 3.16e10, np.array([0,0,0]), np.array([-1.3512e9, -1.3976e7, 3.1592e7]), 25, "yellow", solarsys)
    mercury = Body(3.302e23, np.array([-1.0355e1, -4.6637e1, -2.8598e0]) * 3.16e10 * v, np.array([0,0,0]), np.array([-5.9051e10, 2.5058e8, 5.3457e9]), 3, "burlywood", solarsys)
    venus = Body(48.685e24, np.array([3.6162e0, 3.464e1, 2.6736e-1]) * 3.16e10 * v, np.array([0,0,0]), np.array([1.0657e11, -1.1721e10, -6.3566e9]), 5, "blanchedalmond", solarsys)
    earth = Body(5.9722e24, np.array([-2.5929e1, -1.5639e1, 1.1483e-3]) * 3.16e10 * v, np.array([0,0,0]), np.array([x_value, 1.2570e11, 2.5858e7]), 5, "deepskyblue", solarsys)
    #moon = Body(7.349e22, np.array([-2.5005e1, -1.5036e1, -9.1161e-3]) * 3.16e10, np.array([0,0,0]), np.array([-7.7780e10, 1.2541e11, -5.0425e6]), 1, "darkgray", solarsys)
    mars = Body(6.4171e23, np.array([-2.3074e1, -1.3303e0, 5.3856e-1]) * 3.16e10 * v, np.array([0,0,0]), np.array([-3.4426e10, 2.3527e11, 5.7741e9]), 4, "goldenrod", solarsys)
    jupiter = Body(1.8982e27, np.array([-3.3404e0, 1.3283e1, 1.9609e-2]) * 3.16e10 * v, np.array([0,0,0]), np.array([7.1676e11, 1.8057e11, -1.6785e10]), 12, "peru", solarsys)
    saturn = Body(5.6834e26, np.array([4.7786e0, 8.0455e0, 3.3070e-1]) * 3.16e10 * v, np.array([0,0,0]), np.array([1.2262e12, -8.0900e11, -3.4756e10]), 10, "khaki", solarsys)
    uranus = Body(86.81e24, np.array([-5.0654e0, 4.2894e0, 8.1455e-2]) * 3.16e10 * v, np.array([0,0,0]), np.array([1.9896e12, 2.1661e12,-1.7730e10]), 8, "paleturquoise", solarsys)
    neptune = Body(102.4e24, np.array([4.8688e-1, 5.4420e0,-1.2344e-1]) * 3.16e10 * v, np.array([0,0,0]), np.array([4.4517e12, -4.3035e11, -9.3732e10]), 7, "cyan", solarsys)

# create function to describe data
def fit_funtion(t, exponent, offset):
    delta = exponent * t + offset
    return delta

# calculate lyapunov exponent to determine the chaos in the system
def calculate_exponent():
    mass_list = np.array([1/13, 1/10, 1/7, 1/3, 1, 1.2, 1.4, 1.6, 1.8, 2, 4, 6, 8, 10, 15, 20, 25, 30, 40, 50])
    exp_list = []
    exp_err_list = []

    for mass in mass_list:
        solarsys1 = System()
        solarsys2 = System()

        make_solar_system(solarsys1, -7.7974e10, mass * 1.9885e30, np.sqrt(mass))
        make_solar_system(solarsys2, -7.7974001e10, mass * 1.9885e30, np.sqrt(mass))

        delta = []
        time = []
        errors = []

        number = fraction/10 - 1
        for t in range(years * fraction):
            solarsys1.run_sim()
            solarsys2.run_sim()

            if t == number:
                dx = 0
                dv = 0
                for i in range(len(solarsys1.bodies)):
                    dx += np.linalg.norm(solarsys2.bodies[i].position - solarsys1.bodies[i].position)**2
                    dv += np.linalg.norm(solarsys2.bodies[i].velocity - solarsys1.bodies[i].velocity)**2
                
                delta.append(0.5 * np.log(dx + dv))
                time.append(t)
                errors.append(0.01)
                number += fraction/10

        our_model = models.Model(fit_funtion)
        result    = our_model.fit(delta, t=time, weights =errors, exponent=0.15, offset = 21.5)
        exp = result.params['exponent'].value
        exp_err = result.params['exponent'].stderr

        exp_list.append(1/exp)
        exp_err_list.append(exp_err)
    
    plt.clf()
    plt.plot(mass_list, exp_list, 'o')
    plt.show()


#simulation_multiple_star_system(1)
calculate_exponent()