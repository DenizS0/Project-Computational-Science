Project Computational Science

3D simulation of the solar system.
Contains: Sol, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune

Acceleration calculated using Newton's second and gravitational laws. Position and speed determined using Hermite or Velocity Verlet integration scheme.

Global variables give some control over the simulation:
	fraction: the amount of timesteps a single year gets split in
	years: the amount of years simulated. 
	axis_size = the size in meters of a single axis in the 3D plot 
	hermite = boolean for hermite or velocity verlet integration
	animation = boolean for showing the animation of the solar system (very slow)
	
Bottom of the code has three function calls.
	simulation_solar_system(mass): simulation of the solar system. takes argument mass which is expressed in solar masses.
	calculate_exponent(): calculates the lyapunov exponent for several different solar masses
	plot_errors(): compares the simulation to the values in data.py to estimate errors.
	
The plots used in the poster come from calculate_exponent() and plot_errors() with fraction=10000, years=1 and hermite = True.

Dependencies:
	matplotlib.pyplot
	numpy
	lmfit
using the latest versions as of 22-01-'24
