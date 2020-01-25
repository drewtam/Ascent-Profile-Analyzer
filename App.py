from math import *

gravity_const = 6.6743 * (10^-11)  # m^3/kg/s^2, universal gravitational constant
gas_const = 0.287053  # kJ/kg/K, universal gas constant

thrust = 2000  # kN, engine thrust, user set parameter from the rocket design
drag_coefficient = 0.05  # m^2, drag coefficient "Cd * Area", user set parameter from the rocket design
mass_vehicle = 100  # kg, user set parameter from the rocket design
exhaust_flow = 0.1  # kg/s, user set parameter from the rocket design

radius_planet = 600000  # meters, user set parameter from the planet of origin (Kerbal is 600km)
mass_planet = 5.2915158 * (10^22)  # kg, user set parameter from the planet of origin (Kerbal is 5.2915158e22 kg)
initial_pressure = 101.3  # kPa, atmospheric initial pressure, user set parameter from the planet of origin
scale_height = 5600  # meters, atmospheric pressure constant of KSP's atmosphere model, user set parameter from the planet of origin
temperature = 293  # Kelvin, atmospheric temperature, user set parameter from the planet of origin (Kerbal is 20C at ground)
top_of_atmosphere = 70000 # meters, highest altitude that atmosphere affects flight & drag, user set parameter from the planet of origin (Kerbal is 70km)

delta_time = 0.01  # seconds, time step, ???maybe a user set parameter??? Controls simulation accuracy vs computational speed. ???maybe dynamically generated based on velocity or some other error calc???

flight_plan_altitude = [3000, 12000, 20000, 30000, 80000]  # meters, flight plan altitude points, user set values for their launch plan
flight_plan_attitude = [85, 72, 40, 30, 5, 0]  # degrees, flight plan attitude, user set values for their launch plan

current_time = 0.00  # seconds, simulation time, initialization

x_pos = 0  # meters, cartesian position of the vehicle relative to planet center (origin), initialization
y_pos = radius_planet  # meters, cartesian position of the vehicle relative to planet center (origin), initialization

alpha_pos = pi / 2  #radians, angular position of the vehicle relative to planet center (origin), initialization
radius_pos = radius_planet  # meters, radial position of the vehicle from planet center (origin), initialization

x_velocity = 0 # m/s, cartesian velocity of the vehicle, initialization
y_velocity = 0 # m/s, cartesian velocity of the vehicle, initialization

def tangential_velocity(Vx, Vy, alpha):
    # m/s, calculates velocity perpendicular to the position vector, assumes planar motion
    # Vx = x velocity, Vy = y velocity, alpha = alpha_pos
    return Vx * sin(alpha) + Vy * cos(alpha)

def radial_velocity(Vx, Vy, alpha):
    # m/s, calculates velocity co-linear to the position vector, assumes planar motion
    # Vx = x velocity, Vy = y velocity, alpha = alpha_pos
    return Vx * cos(alpha) + Vy * sin(alpha)

def altitude():
    # meters, calculates current altitude from the radial position
    return radius_pos - radius_planet


def attitude_deg(altitude2):
    # degrees, angle, read the flight plan, find the desired attitude as a function of altitude
    a = 0
    for i in flight_plan_altitude:
        if altitude2 >= flight_plan_altitude[i]:
            a = i
    return flight_plan_attitude[a]


def attitude_rad():
    #radians, angle, calculate the attitude in radians from the flight plan in degrees
    a = altitude() #calculate current altitude, to pass to flight plan function
    return attitude_deg(a) * pi / 180


def air_pressure():
    # kPa, calculates air pressure from current altitude
    return initial_pressure * exp(-(altitude()) / scale_height)


def air_density():
    # kg/m^3, calculates air density
    return air_pressure() / (gas_const * temperature)


#def angular_position(alpha_pos):
 #   alpha_pos = atan(y_pos / x_pos)  # update radial position from cartesian coordinates


force_x_N = thrust * cos(attitude_rad - alpha_pos + pi / 2) - force_drag * cos(acos()
velocity_ms = 0


def vehicle_mass(mass_vehicle):
    mass_vehicle = mass_vehicle - exhaust_flow * delta_time


orbital_energy_kgm2s2 = 1 / 2 * mass_vehicle * velocity_ms ^ 2 - gravity_const * mass_planet /
angular_momentum_kgm2s = 0
reduced_mass_kg = mass_planet * mass_vehicle / (mass_planet + mass_vehicle)
central_force_kgm3s2 = gravity_const * mass_planet * mass_vehicle
eccentricity = sqrt(
    1 + 2 * orbital_energy_kgm2s2 * angular_momentum_kgm2s ^ 2 / (reduced_mass_kg * central_force_kgm3s2 ^ 2)))


# record initial conditions

# Begin main loop
while apoapsis <= target_ap and altitude() <= top_of_atmosphere:
    # increment time step
    # increment new states
    # store new states
    #loop

