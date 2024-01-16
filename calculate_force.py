import numpy as np
import matplotlib.pyplot as plt

# Simulation works with multiple planets, moons not included

# constants
# constants changed to make them all as close to 1 as possible, to prevent
# floating point overflow or underflow at any point in the simulation
dt = 1/365 # time in years
G = 4 * np.pi**2    #6.674 * 10**-11

# radius to sun in AU (1 = 149.6 * 10**6 m), x, y, z, vx, vy (6.278 in au/year = 29.78 km/s), vz, ax, ay, az, mass in solar masses
earth = [1, 1, 0, 0, 0, 6.278, 0, 0, 0, 0, 3 * 10**-6]
venus = [0.72, np.cos(0.059) * 0.72, 0, np.sin(0.059) * 0.72, 0, np.cos(0.059) * 7.383, np.sin(0.059) * 7.383, 0, 0, 0, 2.445 * 10**-6]
mercury = [0.4, 0.4, 0, 0, 0, 10.09, 0, 0, 0, 0, 0.0553 * 3 * 10**-6]
moon = [1 - 0.0027, 1 + 0.0027, 0, 0, 0, -6.278, 0, 0, 0, 0, 0.012 * 3 * 10**-6]
sun = [1, 0, 0, 0] # mass = 1.989 * 10**30 kg, x, y, z
planets = [venus, earth, mercury, moon]
inclination = 0.017 # inclination venus

# calculating acceleration
def gravitation_acceleration(M, planet1, planet2):
    ax = (G * M * (planet1[1] - planet2[1])) / ((planet2[1] - planet1[1])**2 + (planet2[2] - planet1[2])**2 + (planet2[3] - planet1[3])**2)**(3/2)
    ay = (G * M * (planet1[2] - planet2[2])) / ((planet2[1] - planet1[1])**2 + (planet2[2] - planet1[2])**2 + (planet2[3] - planet1[3])**2)**(3/2)
    az = (G * M * (planet1[3] - planet2[3])) / ((planet2[1] - planet1[1])**2 + (planet2[2] - planet1[2])**2 + (planet2[3] - planet1[3])**2)**(3/2)

    return ax, ay, az

# first calculate position at full time interval.
# using the new position, calculate acceleration at full time interval.
# Last, calculate new velocity at full time interval.
def update(object, index_object):
    ax = 0
    ay = 0
    az = 0

    # first step, update position
    object[1] = object[1] + object[4] * dt + object[7] * dt * dt * 0.5
    object[2] = object[2] + object[5] * dt + object[8] * dt * dt * 0.5
    object[3] = object[3] + object[6] * dt + object[9] * dt * dt * 0.5

    # second step, update acceleration
    # gravitation of the sun
    ax, ay, az = gravitation_acceleration(sun[0], sun, object)
    # if index_object == 0:
    #     ax = ax * np.cos(inclination)
    #     az = az * np.sin(inclination)

    # gravitation of the planets
    for i in range(len(planets)):
        if i == index_object:
            continue

        ax_temp, ay_temp, az_temp = gravitation_acceleration(planets[i][10], planets[i], object)

        # if index_object == 0:
        #     ax_temp = ax_temp * np.cos(inclination)
        #     az_temp = az_temp * np.sin(inclination)\

        ax += ax_temp
        ay += ay_temp
        az += az_temp

        #if index_object == 3:
        #    print(f"ax_temp, ax moon: {ax_temp, ax}")

    # third step, update velocity
    object[4] = object[4] + (object[7] + ax) * 0.5 * dt
    object[5] = object[5] + (object[8] + ay) * 0.5 * dt
    object[6] = object[6] + (object[9] + az) * 0.5 * dt

    object[7] = ax
    object[8] = ay
    object[9] = az

x_earth = []
y_earth = []
z_earth = []
x_venus = []
y_venus = []
z_venus = []
x_mercury = []
y_mercury = []
z_mercury = []
x_moon = []
y_moon = []

for t in range(365):
    for i in range(len(planets)):
        update(planets[i], i)
        #print(f"ronde: {t}")
    
        #print(f"planet {i} x,y: {planets[i][1], planets[i][2]}")

    x_earth.append(earth[1])
    y_earth.append(earth[2])
    z_earth.append(earth[3])
    x_venus.append(venus[1])
    y_venus.append(venus[2])
    z_venus.append(venus[3])
    x_mercury.append(mercury[1])
    y_mercury.append(mercury[2])
    z_mercury.append(mercury[3])
    x_moon.append(moon[1])
    y_moon.append(moon[2])

ax = plt.figure().add_subplot(projection='3d')
ax.plot(x_earth, y_earth, z_earth)
ax.plot(x_venus, y_venus, z_venus)
ax.plot(x_mercury, y_mercury, z_mercury)
ax.plot(0, 0, 0, 'bo')

# plt.plot(x_earth, y_earth)
# plt.plot(x_venus, y_venus)
# plt.plot(x_mercury, y_mercury)
# plt.plot(x_moon, y_moon)
# plt.plot(0, 0, 0, 'bo')
plt.show()
