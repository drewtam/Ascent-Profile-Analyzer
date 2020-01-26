import math
# import json
from math import *

gravity_const = 6.6743 * (10**-11)  # m^3/kg/s^2, universal gravitational constant
gas_const = 0.287053  # kJ/kg/K, universal gas constant
g_accel = 9.81  # m/s^2, acceleration due to gravity at the surface (should only be used for normalization)

engine_thrust = 2000  # N, engine thrust, user set parameter from the rocket design
drag_coefficient = 0.05  # m^2, drag coefficient "Cd * Area", user set parameter from the rocket design
mass_vehicle = 100  # kg, user set parameter from the rocket design
exhaust_flow = 0.1  # kg/s, user set parameter from the rocket design
initial_fuel_mass = 60  # kg, user set parameter from the rocket design
specific_impulse = 300  # seconds, user set parameter from the rocket design

radius_planet = 600000.00  # meters, user set parameter from the planet of origin (Kerbal is 600km)
mass_planet = 5.2915158 * (10**22)  # kg, user set parameter from the planet of origin (Kerbal is 5.2915158e22 kg)
initial_pressure = 101.3  # kPa, atmospheric initial pressure, user set parameter from the planet of origin
scale_height = 5600  # meters, atmospheric pressure constant of KSP's atmosphere model, user set parameter from the planet of origin
temperature = 293  # Kelvin, atmospheric temperature, user set parameter from the planet of origin (Kerbal is 20C at ground)
top_of_atmosphere = 70000  # meters, highest altitude that atmosphere affects flight & drag, user set parameter from the planet of origin (Kerbal is 70km)

delta_time = 0.01  # seconds, time step, Controls simulation accuracy vs computational cost. ???maybe a user set parameter??? ???maybe dynamically generated based on velocity or some other error calc???

flight_plan_altitude = [3000, 12000, 20000, 30000, 80000]  # meters, flight plan altitude points, user set values for their launch plan
flight_plan_attitude = [85, 72, 40, 30, 5, 0]  # degrees, flight plan attitude, user set values for their launch plan

current_time = 0.00  # seconds, simulation time, initialization

x_pos = 0.00  # meters, cartesian position of the vehicle relative to planet center (origin), initialization
y_pos = radius_planet  # meters, cartesian position of the vehicle relative to planet center (origin), initialization

alpha_pos = pi / 2  # radians, angular position of the vehicle relative to planet center (origin), initialization
radius_pos = radius_planet  # meters, radial position of the vehicle from planet center (origin), initialization

x_velocity = 0.00  # m/s, cartesian velocity of the vehicle, initialization
y_velocity = 0.00  # m/s, cartesian velocity of the vehicle, initialization

target_apoapsis = flight_plan_altitude[-1]


def tangential_velocity(vx, vy):
    # m/s, calculates velocity perpendicular to the position vector, assumes planar motion
    # Vx = x velocity, Vy = y velocity, alpha = angular position
    alpha = angular_position(x_pos, y_pos)
    return vx * sin(alpha) + vy * cos(alpha)


def radial_velocity(vx, vy, alpha):
    # m/s, calculates velocity co-linear to the position vector, assumes planar motion
    # Vx = x velocity, Vy = y velocity, alpha = alpha_pos
    return vx * cos(alpha) + vy * sin(alpha)


def altitude():
    # meters, calculates current altitude from the radial position
    return radius_pos - radius_planet


def attitude_deg():
    # degrees, angle, returns the attitude from the flight plan. Calculates the current altitude, reads the flight plan.
    alt = altitude()
    a = 0
    for i in flight_plan_altitude:
        if alt >= flight_plan_altitude[i]:
            a = i
    return flight_plan_attitude[a]


def attitude_rad():
    # radians, angle, calculate the attitude in radians from the flight plan in degrees
    # a = altitude()  # calculate current altitude, to pass to flight plan function
    return attitude_deg() * pi / 180


def air_pressure():
    # kPa, returns the air pressure from current altitude
    return initial_pressure * exp(-(altitude()) / scale_height)


def air_density():
    # kg/m^3, returns the air density
    return air_pressure() / (gas_const * temperature)


def angular_position(x, y):
    # radians, returns the angular position from cartesian coordinates, assumes planar position
    return atan(y / x)


def radial_position(x, y):
    # meters, returns the radial position from cartesian coordinates, assumes planar position
    return sqrt(x**2 + y**2)


def thrust():
    # Newtons, returns the engine thrust. Turns engine thrust on/off based on calculated apoapsis. Engine will burn until target apoapsis is expected, prevents overshoot.
    if apoapsis_altitude() >= target_apoapsis:
        return 0
    else:
        return engine_thrust


def velocity(vx, vy):
    return sqrt(vx**2 + vy**2)


def orbital_energy(m, v, r):
    # kg * m^2 / s^2, returns conserved orbital energy (euler? lagrangian? can't remember the technical term), requires vehicle mass, scalar velocity, and radial position
    return 1 / 2 * m * v**2 - gravity_const * mass_planet * m / r


def angular_momentum(r, vt, m):
    # kg*m^2/s, returns orbital angular momentum, requires radial position, tangential velocity, and vehicle mass
    return r * vt * m


def reduced_mass(m):
    # kg, effective mass in two body newtonian mechanics
    return mass_planet * m / (mass_planet + m)


def new_vehicle_mass():
    # kg, returns an adjusted vehicle mass for a small time step
    if thrust() > 0:
        return mass_vehicle - exhaust_flow * delta_time
    else:
        return mass_vehicle


def central_force():
    # kg * m^3 / s^2, returns the scalar force value for the planetary gravity, without direction or distance
    return gravity_const * mass_planet * mass_vehicle


def eccentricity():
    # unitless, returns the orbital eccentricity of the vehicle around the planet
    m = mass_vehicle
    mr = reduced_mass(m)
    v = velocity(x_velocity, y_velocity)
    r = radial_position(x_pos, y_pos)
    f = central_force()
    vt = tangential_velocity(x_velocity, y_velocity)
    return sqrt(1 + 2 * orbital_energy(m, v, r) * (angular_momentum(r, vt, m))**2 / (mr * f**2))


def semi_major_axis(v, r):
    # meters, returns the radius of the semi major axis of the vehicle orbit around the planet, requires scalar velocity and radial position
    return -gravity_const * mass_planet / (2 * (v**2 / 2 - gravity_const * mass_planet / r))


def apoapsis_total():
    # meters, returns the radius of the orbital apoapsis
    v = velocity(x_velocity, y_velocity)
    r = radial_position(x_pos, y_pos)
    sma = semi_major_axis(v, r)
    ecc = eccentricity()
    return sma * (1 + abs(ecc))


def apoapsis_altitude():
    # meters, returns the radius minus the planetary diameter (altitude of the apoapsis)
    return apoapsis_total() - radius_planet


def periapsis_total():
    # meters, returns the radius of the orbital periapsis
    v = velocity(x_velocity, y_velocity)
    r = radial_position(x_pos, y_pos)
    sma = semi_major_axis(v, r)
    ecc = eccentricity()
    return sma * (1 - abs(ecc))


def periapsis_altitude():
    # meters, returns the radius minus the planetary diameter (altitude of the apoapsis)
    return periapsis_total() - radius_planet


def apoapsis_velocity():
    # m/s, returns velocity at apoapsis
    v = velocity(x_velocity, y_velocity)
    r = radial_position(x_pos, y_pos)
    sma = semi_major_axis(v, r)
    ecc = eccentricity()
    return sqrt((1-ecc) * gravity_const * mass_planet / ((1+ecc) * sma))


def apo_gravity_accel():
    # m/s^2, returns the free fall acceleration due to gravity at apoapsis
    apo = apoapsis_total()
    return gravity_const * mass_vehicle / (apo**2)


def circular_orbital_apoapsis_velocity():
    # m/s, returns the required orbital velocity for a circular orbit at apoapsis
    aga = apo_gravity_accel()
    apo = apoapsis_total()
    return sqrt(aga * apo)


def velocity_deficit():
    # m/s, returns the required velocity needed to circulize orbit at apoapsis
    orbv = circular_orbital_apoapsis_velocity()
    apov = apoapsis_velocity()
    return orbv - apov


def drag_force():
    # Newtons, returns atmospheric drag force scalar
    v = velocity(x_velocity, y_velocity)
    return 1/2 * drag_coefficient * air_density() * (v**2)


def gravity_force():
    # Newtons, returns scalar force due to gravity
    r = radial_position(x_pos, y_pos)
    return gravity_const * mass_planet * mass_vehicle / (r**2)


def x_force():
    # Newtons, returns sum of forces in the cartesian x direction, assumes planar forces
    alpha = angular_position(x_pos, y_pos)
    tfx = thrust() * cos(attitude_rad() - alpha + pi/2)  # thrust in the x direction
    v = velocity(x_velocity, y_velocity)
    dfx = drag_force() * cos(acos(x_velocity / v))  # drag in the x direction
    gfx = gravity_force() * cos(alpha)  # gravity in the x direction
    return tfx - dfx - gfx


def y_force():
    # Newtons, returns sum of forces in the cartesian y direction, assumes planar forces
    alpha = angular_position(x_pos, y_pos)
    tfy = thrust() * sin(attitude_rad() - alpha + pi/2)  # thrust in the y direction
    v = velocity(x_velocity, y_velocity)
    dfy = drag_force() * sin(asin(y_velocity / v))  # drag in the y direction
    gfy = gravity_force() * sin(alpha)  # gravity in the x direction
    return tfy - dfy - gfy


def delta_v():
    # m/s, returns the "delta-v" left remaining to the vehicle
    m1 = mass_vehicle  # current vehicle mass
    m2 = mass_vehicle - initial_fuel_mass  # dry vehicle mass
    i = specific_impulse
    g = g_accel
    return math.log(m1/m2) * i * g


def write_header():
    with open("flight_log.txt", "w") as file:
        file.write("time,s  mass,kg  attitude,deg")


def write_data():
    with open("flight_log.txt", "w") as file:
        value_list = [current_time, mass_vehicle, attitude_deg()]
        for i in value_list:
            file.write(str(value_list[i]))


# record initial conditions
# attitude deg, air pressure, air density, drag force, gravity force, thrust, vx, vy, vt, vr, x, y, alpha, radial, altitude, apoapsis, periapsis, vap, remaining v, delta-v

# Begin main loop
# while apoapsis <= target_ap and altitude() <= top_of_atmosphere:
# increment time step
# increment new states
# store new states
# loop

# file.close()
