#
# Names: Dagmar Salomons & Deniz Saglam
# Student IDs: 13981714 & 13827189
#
# 3D Simulation of the solar system
# The solar mass is altered to determine the
# lyapunov time dependence of solar mass.

import matplotlib.pyplot as plt
import numpy as np
from lmfit import models

# File containing true positions and velocities of the planets taken from nasa horizon data set
# Used to calculate errors on positions en velocities 
import data

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

    # add body to list of bodies in the solar system
    def make_body(self, new_body):
        self.bodies.append(new_body)

    # run simulation using either hermite of velocity verlet integration
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

    # calculate total energy of the system by adding kinetic and gravitational energy
    # kinetic energy = 0.5 * mass * v^2, gravitational energy = -(M1 * M2 * G) / r
    def total_energy(self):
        total = 0

        # energies of all bodies are added together
        for body in self.bodies:
            total += .5 * body.mass * np.linalg.norm(body.velocity)**2

            # gravitational energy for 1 body is calculated from interaction between all bodies
            for i in range(len(self.bodies)):
                if self.bodies[i] != body:
                    total += (-0.5 * G_CONSTANT * body.mass * self.bodies[i].mass) / np.linalg.norm(self.bodies[i].position - body.position)

        self.energy.append(total)

    # draw every body in the plot
    def draw_plot(self):
        # sort the planets based on x position to make sure planets behind the sun are not visible
        self.bodies.sort(key=lambda item: item.position[0])
        for body in self.bodies:
            body.draw_body()

    # calculate new position based on chosen integration method
    def calculate_position(self):
        for body in self.bodies:
            if hermite == True:
                body.hermite_pos()
            else:
                body.new_pos()

    # calculate new acceleration, differentiated between integration method in called function
    def calculate_acceleration(self):
        for n, body in enumerate(self.bodies):
            body.acceleration(self.bodies, n)

    # calculate new velocity based on chosen integration method
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

        # acceleration of 1 body has to be calculated from every other body in the solar system
        for i in range(len(bodies)):
            if i != n:
                force_vector = bodies[i].position - self.position
                distance = np.linalg.norm(force_vector)

                self.accel = self.accel + (G_CONSTANT * bodies[i].mass * force_vector) / distance**3

                if hermite == True:
                    vel_vector = bodies[i].velocity - self.velocity
                    self.jerk = self.jerk + G_CONSTANT * bodies[i].mass * (vel_vector / distance**3 - 3 * ((force_vector * np.dot(force_vector, vel_vector)) / distance**5))

    # calculate new velocity acc. to hermite
    def hermite_vel(self):
        self.old_velocity = self.velocity
        self.predict_acc = self.old_accel + self.old_jerk * dt
        self.velocity = self.velocity + .5 * (self.old_accel + self.predict_acc) * dt + 1/12 * (self.old_jerk - self.jerk) * dt**2

    # calculate new position acc. to hermite
    def hermite_pos(self):
        self.predict_v = self.velocity + self.accel * dt + 0.5 * self.old_jerk * dt**2
        self.position = self.position + .5 * (self.velocity + self.predict_v) * dt + 1/12 * (self.old_accel - self.predict_acc) * dt**2

    # draw new position in the simulation
    def draw_body(self):
        self.system.axes.plot(self.position[0], self.position[1], self.position[2], marker="o", color=self.color, markersize = self.size)


# Function to create the sun all the planets for a system.
# starting values taken from nasa horizon data set @ 22 jan 2023 00:00
# Arguments:
# solarsys = system to add the sun and planets to.
# x_value = x_value earth
# mass_sun = mass of the sun
# v = extra factor to change the speed, sqrt(factor of solar mass)
def make_solar_system(solarsys, x_value, mass_sun, v):
    sun = Body(mass_sun, np.array([2.0981e-3, -1.5479e-2, 7.8087e-5]) * 3.16e10, np.array([0,0,0]), np.array([-1.3512e9, -1.3976e7, 3.1592e7]), 25, "yellow", solarsys)
    mercury = Body(3.302e23, np.array([-1.0355e1, -4.6637e1, -2.8598e0]) * 3.16e10 * v, np.array([0,0,0]), np.array([-5.9051e10, 2.5058e8, 5.3457e9]), 3, "burlywood", solarsys)
    venus = Body(48.685e24, np.array([3.6162e0, 3.464e1, 2.6736e-1]) * 3.16e10 * v, np.array([0,0,0]), np.array([1.0657e11, -1.1721e10, -6.3566e9]), 5, "blanchedalmond", solarsys)
    earth = Body(5.9722e24, np.array([-2.5929e1, -1.5639e1, 1.1483e-3]) * 3.16e10 * v, np.array([0,0,0]), np.array([x_value, 1.2570e11, 2.5858e7]), 5, "deepskyblue", solarsys)
    mars = Body(6.4171e23, np.array([-2.3074e1, -1.3303e0, 5.3856e-1]) * 3.16e10 * v, np.array([0,0,0]), np.array([-3.4426e10, 2.3527e11, 5.7741e9]), 4, "goldenrod", solarsys)
    jupiter = Body(1.8982e27, np.array([-3.3404e0, 1.3283e1, 1.9609e-2]) * 3.16e10 * v, np.array([0,0,0]), np.array([7.1676e11, 1.8057e11, -1.6785e10]), 12, "peru", solarsys)
    saturn = Body(5.6834e26, np.array([4.7786e0, 8.0455e0, 3.3070e-1]) * 3.16e10 * v, np.array([0,0,0]), np.array([1.2262e12, -8.0900e11, -3.4756e10]), 10, "khaki", solarsys)
    uranus = Body(86.81e24, np.array([-5.0654e0, 4.2894e0, 8.1455e-2]) * 3.16e10 * v, np.array([0,0,0]), np.array([1.9896e12, 2.1661e12,-1.7730e10]), 8, "paleturquoise", solarsys)
    neptune = Body(102.4e24, np.array([4.8688e-1, 5.4420e0,-1.2344e-1]) * 3.16e10 * v, np.array([0,0,0]), np.array([4.4517e12, -4.3035e11, -9.3732e10]), 7, "cyan", solarsys)


# Function to calculate the total energy difference in percentage for every timestep
# compared to the total energy at t = 0. This can be used to check energy conservation.
#
# Argument solarsys = system for which to calculate the difference in total energy.
# Side effect: plot is created to show difference in energy over time.
def calculate_energy_difference(solarsys):
    days = [i for i in range(years * fraction)]

    energy = [((solarsys.energy[i] - solarsys.energy[0]) * 100 / solarsys.energy[0]) for i in range(len(solarsys.energy))]

    plt.clf()
    plt.plot(days, energy, "o-")
    plt.title("Total Energy Difference")
    plt.ylabel("Percentage error (%)")
    plt.xlabel("Number of time steps")
    plt.show()


# Calculate error on positions and velocities for solar system
# Average error of all the planets is calculated and used as error.
# Arguments:
# system = solar system with regular solar mass to calculate errrors
#
# Return: 2 lists with errors on position and velocity
# for 10 timesteps equally distributed over 1 year
def error_solar_system(system):
    counter = 0
    number = fraction/10 - 1

    err_vel_list = []
    err_pos_list = []

    for t in range(years * fraction):
        system.run_sim()

        if t == number:
            err_list_x = []
            err_list_v = []
            for i in range(len(system.bodies)):
                # calculate error on position and velocity
                true_x = np.array([data.true_pos[30 * i + counter * 3], data.true_pos[30 * i + counter * 3 + 1], data.true_pos[30 * i + counter * 3 + 2]])
                true_v = np.array([data.true_vel[30 * i + counter * 3], data.true_vel[30 * i + counter * 3 + 1], data.true_vel[30 * i + counter * 3 + 2]])

                err_x = np.linalg.norm(true_x * 1e3 - system.bodies[i].position)
                err_v = np.linalg.norm(true_v * 3.16e10 - system.bodies[i].velocity)

                err_list_x.append(err_x * 100 / np.linalg.norm(true_x * 1e3))
                err_list_v.append(err_v * 100 / np.linalg.norm(true_v * 3.16e10))

            counter += 1
            number += fraction/10

            err_pos_list.append(np.mean(err_list_x))
            err_vel_list.append(np.mean(err_list_v))

    return err_pos_list, err_vel_list


# Function to run the simulation of the solar system.
# Argument mass = mass of the sun expressed in solar masses.
#   So, a mass of 2 = 2 * solar mass
# The mass of the sun is altered, and the speed of the planets is altered to create stable orbits.
#   The speed is multiplied by sqrt of the mass, this is based on v = sqrt(G * M / r)
# starting values taken from nasa horizon data set @ 22 jan 2023 00:00
def simulation_solar_system(mass):
    m = np.sqrt(mass)

    solarsys = System()
    sun = Body(1.9885e30 * m**2, np.array([2.0981e-3, -1.5479e-2, 7.8087e-5]) * 3.16e10, np.array([0,0,0]), np.array([-1.3512e9, -1.3976e7, 3.1592e7]), 25, "yellow", solarsys)
    mercury = Body(3.302e23, np.array([-1.0355e1, -4.6637e1, -2.8598e0]) * 3.16e10 * m, np.array([0,0,0]), np.array([-5.9051e10, 2.5058e8, 5.3457e9]), 3, "burlywood", solarsys)
    venus = Body(48.685e24, np.array([3.6162e0, 3.464e1, 2.6736e-1]) * 3.16e10 * m, np.array([0,0,0]), np.array([1.0657e11, -1.1721e10, -6.3566e9]), 5, "blanchedalmond", solarsys)
    earth = Body(5.9722e24, np.array([-2.5929e1, -1.5639e1, 1.1483e-3]) * 3.16e10 * m, np.array([0,0,0]), np.array([-7.7974e10, 1.2570e11, 2.5858e7]), 5, "deepskyblue", solarsys)
    moon = Body(7.349e22, np.array([-2.5005e1, -1.5036e1, -9.1161e-3]) * 3.16e10 * m, np.array([0,0,0]), np.array([-7.7780e10, 1.2541e11, -5.0425e6]), 1, "darkgray", solarsys)
    mars = Body(6.4171e23, np.array([-2.3074e1, -1.3303e0, 5.3856e-1]) * 3.16e10 * m, np.array([0,0,0]), np.array([-3.4426e10, 2.3527e11, 5.7741e9]), 4, "goldenrod", solarsys)
    jupiter = Body(1.8982e27, np.array([-3.3404e0, 1.3283e1, 1.9609e-2]) * 3.16e10 * m, np.array([0,0,0]), np.array([7.1676e11, 1.8057e11, -1.6785e10]), 12, "peru", solarsys)
    saturn = Body(5.6834e26, np.array([4.7786e0, 8.0455e0, 3.3070e-1]) * 3.16e10 * m, np.array([0,0,0]), np.array([1.2262e12, -8.0900e11, -3.4756e10]), 10, "khaki", solarsys)
    uranus = Body(86.81e24, np.array([-5.0654e0, 4.2894e0, 8.1455e-2]) * 3.16e10 * m, np.array([0,0,0]), np.array([1.9896e12, 2.1661e12,-1.7730e10]), 8, "paleturquoise", solarsys)
    neptune = Body(102.4e24, np.array([4.8688e-1, 5.4420e0,-1.2344e-1]) * 3.16e10 * m, np.array([0,0,0]), np.array([4.4517e12, -4.3035e11, -9.3732e10]), 7, "cyan", solarsys)

    for t in range(years * fraction):
        solarsys.run_sim()

        if animation == True:
            solarsys.draw_plot()
            plt.pause(0.001)

            # remove every planet, so in the next time step the new positions can be drawn
            for artist in plt.gca().get_lines():
                artist.remove()

def plot_errors():
    # make solar system
    test_solarsystem = System()
    make_solar_system(test_solarsystem, -7.7974e10, 1.9885e30, np.sqrt(1))

    # create plot with error on positions 10 timesteps for 1 year
    x_err, v_err = error_solar_system(test_solarsystem)
    time = [(i * fraction/10 + 99) for i in range(10)]

    plt.plot(time, x_err, 'o', label = 'position error')
    plt.plot(time, v_err, 'o', label = 'velocity error')
    plt.legend()
    plt.title('Percentage errors on position and velocity')
    plt.xlabel('Number of timesteps')
    plt.ylabel('Percentage errors (%)')
    plt.show()

    # create plot with energy difference for every timestep
    calculate_energy_difference(test_solarsystem)

# create function to describe data
def fit_funtion(t, exponent, offset):
    delta = exponent * t + offset
    return delta


# calculate lyapunov exponent to determine the chaos in the system
def calculate_exponent():
    #mass_list = np.array([1/13, 1/10, 1/7, 1/3, 1, 1.2, 1.4, 1.6, 1.8, 2, 4, 6, 8, 10, 15, 20, 25, 30, 40, 50])
    mass_list = np.array([1, 2, 5, 10])
    exp_list = []
    exp_err_list = []

    test_solarsystem = System()
    make_solar_system(test_solarsystem, -7.7974e10, 1.9885e30, np.sqrt(1))
    delta_err_x_list, delta_err_v_list = error_solar_system(test_solarsystem)

    #time_planets = []
    #number_of_planets_list = []

    #distance_pluto = 5.906e12

    for mass in mass_list:
        solarsys1 = System()
        solarsys2 = System()

        #append_time = True

        make_solar_system(solarsys1, -7.7974e10, mass * 1.9885e30, np.sqrt(mass))
        make_solar_system(solarsys2, -7.7974001e10, mass * 1.9885e30, np.sqrt(mass))

        number_of_bodies = len(solarsys1.bodies)

        delta = []
        time = []
        errors = []
        counter = 0

        #index_planets = [i for i in range(number_of_bodies)]

        number = fraction/10 - 1
        #number_of_planets = 0
    
        for t in range(years * fraction):
            #index_remove = []
            solarsys1.run_sim()
            solarsys2.run_sim()

            if t == number:
                dx = 0
                dv = 0
                err_x = 0
                err_v = 0

                for i in range(len(solarsys1.bodies)):
                    # calculate distance between the same particles in different simulation
                    delta_x = np.linalg.norm(solarsys2.bodies[i].position - solarsys1.bodies[i].position)
                    delta_v = np.linalg.norm(solarsys2.bodies[i].velocity - solarsys1.bodies[i].velocity)

                    # calculate phase space distance
                    dx += delta_x**2
                    dv += delta_v**2

                    # calculate error on position and velocity
                    delta_err_x = (delta_err_x_list[counter] / 100) * np.linalg.norm(solarsys1.bodies[i].position)
                    delta_err_v = (delta_err_v_list[counter] * np.sqrt(mass) / 100) * np.linalg.norm(solarsys1.bodies[i].velocity)

                    # calculate error on distance and velocity bewteen the same particles
                    err_dx = np.sqrt(2) * delta_err_x
                    err_dv = np.sqrt(2) * delta_err_v

                    # calculate total error on distance and velocity
                    err_x += 2 * err_dx * delta_x
                    err_v += 2 * err_dv * delta_v

                delta.append(0.5 * np.log(dx + dv))
                time.append(t)
                # append final error for phase space distance on logaritmic scale
                errors.append(0.5 * (1/(dx + dv)) * (err_x + err_v))
                counter += 1
                number += fraction/10

                # for i in index_planets:
                #     position = solarsys1.bodies[i].position
                #     distance = np.linalg.norm(position)

                #     if distance > distance_pluto:
                #         number_of_planets += 1
                #         index_remove.append(i)

                #         if number_of_planets == 1:
                #             time_planets.append(t)
                #             append_time = False
                
                # for i in index_remove:
                #     index_planets.remove(i)
        
        # number_of_planets_list.append(number_of_planets)
        # if append_time == True:
        #     time_planets.append(years * fraction)

        errors = np.array(errors)
        our_model = models.Model(fit_funtion)
        result    = our_model.fit(delta, t=time, weights=1/errors, exponent=0.15, offset = 21.5)
        exp = result.params['exponent'].value
        exp_err = result.params['exponent'].stderr

        exp_list.append(1/exp)
        exp_err_list.append(exp_err / exp**2)

        # plt.plot(time, delta, 'o')
        # plt.plot(time, result.best_fit, 'o-')
        # plt.show()

    plt.clf()
    plt.errorbar(mass_list, exp_list, yerr=exp_err_list, fmt=".w")
    plt.show()

    # plt.plot(mass_list, number_of_planets_list, 'o')
    # plt.show()

    # plt.plot(mass_list, time_planets, 'o')
    # plt.show()


#simulation_solar_system(1)
#calculate_exponent()

plot_errors()
