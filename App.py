# Andrew Priest, Jan/27/2020
# Ascent Profile Analyzer is a tool to determine the best flight plan for a rocket on Kerbal Space Program.
# The tool uses Kerbal Physics to calculate the relevant free body forces (thrust, gravity, drag), and complete a path integral by simple time step methods.
# The end result is an optimized flight plan that minimizes fuel consumption (or if you prefer, delta-V), when launching through an atmosphere and its aerodynamic drag.

#     Copyright (C) 2020  Andrew Priest
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>


import math
# import json
from math import *

# these hard coded constants are proper
gravity_const = 6.6743 * (10 ** -11)  # m^3/kg/s^2, universal gravitational constant
gas_const = 0.287053  # kJ/kg/K, universal gas constant
g_accel = 9.81  # m/s^2, acceleration due to gravity at the surface (should only be used for normalization)

# these parameters need a user interface front end
# the drag coefficient is especially difficult, and should have a separate tool to help guide the user to determine the value
# it is a known inaccuracy in this model, that if the vehicle attitude does not equal vehicle velocity direction, then the coefficient of drag will be different and could be significantly higher. It is possible to build a more accurate model of KSP drag,  but not a priority
# this model assumes that the vehicle aerodynamics are stable & neutral, with no flipping, no lift.
# this input is currently shown as a single stage launch, but a new set of inputs for additional stages will be necessary. This complication should use all the same math, but the main loop just breaks up which global constants to use
# specific impulse can be complicated if there are multiple rocket types mixed together in the same stage. A separate simple tool may be necessary for users who don't know how to perform the arithmetic.
engine_thrust = 10000  # N, maximum engine thrust, user set parameter from the rocket design
drag_coefficient = 0.05  # m^2, drag coefficient "Cd * Area", user set parameter from the rocket design
initial_mass_vehicle = 500  # kg, user set parameter from the rocket design
initial_fuel_mass = 410  # kg, user set parameter from the rocket design
specific_impulse = 340  # seconds, user set parameter from the rocket design
target_apoapsis = 80000  # meters, user set parameter from the mission goal !!!! Data quality requirement: must be greater than top of atmosphere!!!!!
# add exception to exit program if target apoapsis <= top of atmosphere.
# launch_heading can be a user parameter, when we are ready to take into account planetary rotation. We can reasonably assume a constant value for analysis purposes
exhaust_flow = engine_thrust / (g_accel * specific_impulse)  # kg/s, derived parameter of the rocket design

# these parameters need a user interface front end
# prefer a drop down menu for selecting all the Kerbal planets, with preloaded values from a table.
# prefer to have an additional option for custom inputs
# current model does not consider planetary rotation and its effects on adding or subtracting to the tangential velocity at launch. This will require additional modification of the drag calcs too.
radius_planet = 600000.00  # meters, user set parameter from the planet of origin (Kerbal is 600km)
mass_planet = 5.2915158 * (10 ** 22)  # kg, user set parameter from the planet of origin (Kerbal is 5.2915158e22 kg)
initial_pressure = 101.3  # kPa, atmospheric initial pressure, user set parameter from the planet of origin
scale_height = 5600  # meters, atmospheric pressure constant of KSP's atmosphere model, user set parameter from the planet of origin
temperature = 293  # Kelvin, atmospheric temperature, user set parameter from the planet of origin (Kerbal is 20C at ground)
top_of_atmosphere = 70000  # meters, highest altitude that atmosphere affects flight & drag, user set parameter from the planet of origin (Kerbal is 70km)

delta_time = 0.1  # seconds, time step, Controls simulation accuracy vs computational cost. ???maybe a user set parameter??? ???maybe dynamically generated based on velocity or some other error calc???
refinement = 3  # the number of times to execute the full optimization algorithm. Controls simulation accuracy vs computational cost. ???maybe a user set parameter??? ???maybe dynamically generated based on velocity or some other error calc???
resolution_bottom = 20 # the number of altitude points to divide the flight plan in the bottom thick part of atmosphere. The atmosphere is divided by equal increments of pressure (density).
                    # By increasing the resolution, the last element will be higher in the atmosphere; as each division will be a small density change.
resolution_top = 10 # the number of altitude points to divide the flight plan in the upper atmosphere. This starts where the bottom resolution left off, and divides the distance to top of atmosphere by these equal points.
                    # !!!! Data quality requirement: we probably want resolution_bottom + resolution_top > 5 !!!!!

# initial conditions
# these should be the only global variables modified in main loop

# mass_vehicle = initial_mass_vehicle  # kg, vehicle mass
# current_time = 0.00  # seconds, simulation time, initialization
# x_pos = 0.0  # meters, cartesian position of the vehicle relative to planet center (origin), initialization
# y_pos = radius_planet  # meters, cartesian position of the vehicle relative to planet center (origin), initialization
# x_velocity = 0.0  # m/s, cartesian velocity of the vehicle, initialization
# y_velocity = 0.0  # m/s, cartesian velocity of the vehicle, initialization


def tangential_velocity(x, y, vx, vy):
    # m/s, returns velocity perpendicular to the position vector, assumes planar motion
    v = velocity(vx, vy)
    vr = radial_velocity(x, y, vx, vy)
    return sqrt(v ** 2 - vr ** 2)


def radial_velocity(x, y, vx, vy):
    # m/s, returns velocity co-linear to the position vector, assumes planar motion
    # vx = x_velocity
    # vy = y_velocity
    alpha = angular_position(x, y)
    return vx * cos(alpha) + vy * sin(alpha)


def altitude(x, y):
    # meters, returns current altitude from the radial position
    return radial_position(x, y) - radius_planet


def attitude_deg(x, y, flight_plan):
    # degrees, angle, returns the attitude from the flight plan. Calculates the current altitude, reads the flight plan.
    alt = altitude(x, y)
    a = 0
    for j in range(len(flight_plan[0])):
        if alt >= flight_plan[0][j]:
            a = j
    return flight_plan[1][a]


def attitude_rad(x, y, flight_plan):
    # radians, angle, calculate the attitude in radians from the flight plan in degrees
    return attitude_deg(x, y, flight_plan) * pi / 180


def air_pressure(x, y):
    # kPa, returns the air pressure from current altitude
    return initial_pressure * exp(-(altitude(x, y)) / scale_height)


def air_density(x, y):
    # kg/m^3, returns the air density
    return air_pressure(x, y) / (gas_const * temperature)


def angular_position(x, y):
    # radians, returns the angular position from cartesian coordinates, assumes planar position
    # x = x_pos
    # y = y_pos
    if x == 0:
        return pi / 2
    else:
        return atan(y / x)


def radial_position(x, y):
    # meters, returns the radial position from cartesian coordinates, assumes planar position
    # x = x_pos
    # y = y_pos
    return sqrt(x ** 2 + y ** 2)


def fuel_remain(mass_vehicle):
    # kg, returns the amount of fuel left
    return initial_fuel_mass - (initial_mass_vehicle - mass_vehicle)


def thrust(x, y, vx, vy, mass_vehicle):
    # Newtons, returns the engine thrust. Turns engine thrust on/off based on calculated apoapsis. Engine will burn until target apoapsis is expected, prevents overshoot.
    if apoapsis_altitude(x, y, vx, vy, mass_vehicle) >= target_apoapsis or fuel_remain(mass_vehicle) <= 0:
        return 0
    else:
        return engine_thrust


def velocity(vx, vy):
    # m/s, returns the scalar velocity
    # vx = x_velocity
    # vy = y_velocity
    return sqrt(vx ** 2 + vy ** 2)


def orbital_energy(x, y, vx, vy, mass_vehicle):
    # kg * m^2 / s^2, returns conserved orbital energy (LaGrangian), requires vehicle mass, scalar velocity, and radial position
    m = mass_vehicle
    v = velocity(vx, vy)
    r = radial_position(x, y)
    return 1 / 2 * m * (v ** 2) - gravity_const * mass_planet * m / r


def angular_momentum(x, y, vx, vy, mass_vehicle):
    # kg*m^2/s, returns orbital angular momentum
    r = radial_position(x, y)
    vt = tangential_velocity(x, y, vx, vy)
    m = mass_vehicle
    return r * vt * m


def reduced_mass(mass_vehicle):
    # kg, effective mass in two body newtonian mechanics
    m = mass_vehicle
    return mass_planet * m / (mass_planet + m)


def new_vehicle_mass(x, y, vx, vy, mass_vehicle):
    # kg, returns an adjusted vehicle mass for a small time step
    if thrust(x, y, vx, vy, mass_vehicle) > 0:
        return mass_vehicle - exhaust_flow * delta_time
    else:
        return mass_vehicle


def central_force(mass_vehicle):
    # kg * m^3 / s^2, returns the scalar force value for the planetary gravity, without direction or distance
    return gravity_const * mass_planet * mass_vehicle


def eccentricity(x, y, vx, vy, mass_vehicle):
    # unitless, returns the orbital eccentricity of the vehicle around the planet
    mr = reduced_mass(mass_vehicle)
    f = central_force(mass_vehicle)
    return sqrt(1 + 2 * orbital_energy(x, y, vx, vy, mass_vehicle) * (angular_momentum(x, y, vx, vy, mass_vehicle) ** 2) / (mr * (f ** 2)))


def semi_major_axis(x, y, vx, vy):
    # meters, returns the radius of the semi major axis of the vehicle orbit (ellipse) around the planet. Semi major axis is measured from the center of the ellipse, whereas the apoapsis is measured from the center of the planet.
    v = velocity(vx, vy)
    r = radial_position(x, y)
    return -gravity_const * mass_planet / (2 * (v ** 2 / 2 - gravity_const * mass_planet / r))


def apoapsis(x, y, vx, vy, mass_vehicle):
    # meters, returns the radius of the orbital apoapsis
    sma = semi_major_axis(x, y, vx, vy)
    ecc = eccentricity(x, y, vx, vy, mass_vehicle)
    return sma * (1 + abs(ecc))


def apoapsis_altitude(x, y, vx, vy, mass_vehicle):
    # meters, returns the radius minus the planetary diameter (altitude of the apoapsis)
    return apoapsis(x, y, vx, vy, mass_vehicle) - radius_planet


def periapsis(x, y, vx, vy, mass_vehicle):
    # meters, returns the radius of the orbital periapsis
    sma = semi_major_axis(x, y, vx, vy)
    ecc = eccentricity(x, y, vx, vy, mass_vehicle)
    return sma * (1 - abs(ecc))


def periapsis_altitude(x, y, vx, vy, mass_vehicle):
    # meters, returns the radius minus the planetary diameter (altitude of the apoapsis)
    return periapsis(x, y, vx, vy, mass_vehicle) - radius_planet


def apoapsis_velocity(x, y, vx, vy, mass_vehicle):
    # m/s, returns velocity at apoapsis
    sma = semi_major_axis(x, y, vx, vy)
    ecc = eccentricity(x, y, vx, vy, mass_vehicle)
    return sqrt((1 - ecc) * gravity_const * mass_planet / ((1 + ecc) * sma))


def apo_gravity_accel(x, y, vx, vy, mass_vehicle):
    # m/s^2, returns the free fall acceleration due to gravity at apoapsis
    apo = apoapsis(x, y, vx, vy, mass_vehicle)
    return gravity_const * mass_planet / (apo ** 2)


def circular_orbital_apoapsis_velocity(x, y, vx, vy, mass_vehicle):
    # m/s, returns the required orbital velocity for a circular orbit at apoapsis
    aga = apo_gravity_accel(x, y, vx, vy, mass_vehicle)
    apo = apoapsis(x, y, vx, vy, mass_vehicle)
    return sqrt(aga * apo)


def tangential_velocity_error(x, y, vx, vy, mass_vehicle):
    # m/s, returns the additional velocity needed to circularize orbit at apoapsis
    orbv = circular_orbital_apoapsis_velocity(x, y, vx, vy, mass_vehicle)
    apov = apoapsis_velocity(x, y, vx, vy, mass_vehicle)
    return orbv - apov


def drag_force(x, y, vx, vy):
    # Newtons, returns atmospheric drag force scalar
    v = velocity(vx, vy)
    if altitude(x, y) >= top_of_atmosphere:
        return 0
    else:
        return 1 / 2 * drag_coefficient * air_density(x, y) * (v ** 2)


def gravity_force(x, y, mass_vehicle):
    # Newtons, returns scalar force due to gravity
    r = radial_position(x, y)
    return gravity_const * mass_planet * mass_vehicle / (r ** 2)


def x_force(x, y, vx, vy, mass_vehicle, flight_plan):
    # Newtons, returns sum of forces in the cartesian x direction, assumes planar forces
    alpha = angular_position(x, y)
    ftx = thrust(x, y, vx, vy, mass_vehicle) * cos(attitude_rad(x, y, flight_plan) - alpha + pi / 2)  # thrust in the x direction
    v = velocity(vx, vy)
    if v == 0:
        fdx = 0
    else:
        fdx = drag_force(x, y, vx, vy) * cos(acos(vx / v))  # drag in the x direction
    fgx = gravity_force(x, y, mass_vehicle) * cos(alpha)  # gravity in the x direction
    return ftx - fdx - fgx


def y_force(x, y, vx, vy, mass_vehicle, flight_plan):
    # Newtons, returns sum of forces in the cartesian y direction, assumes planar forces
    alpha = angular_position(x, y)
    fty = thrust(x, y, vx, vy, mass_vehicle) * sin(attitude_rad(x, y, flight_plan) - alpha + pi / 2)  # thrust in the y direction
    v = velocity(vx, vy)
    if v == 0:
        fdy = 0
    else:
        fdy = drag_force(x, y, vx, vy) * sin(asin(vy / v))  # drag in the y direction
    fgy = gravity_force(x, y, mass_vehicle) * sin(alpha)  # gravity in the y direction
    return fty - fdy - fgy


def delta_v(mass_vehicle):
    # m/s, returns the "delta-v" left remaining to the vehicle
    m1 = mass_vehicle  # current vehicle mass
    m2 = initial_mass_vehicle - initial_fuel_mass  # dry vehicle mass
    isp = specific_impulse
    g = g_accel
    return math.log(m1 / m2) * isp * g


def height_potential(x, y):
    # m/s, returns the instantaneous radial velocity needed to achieve the target height (orbital altitude)
    return sqrt(2 * gravity_const * mass_planet * (1/radial_position(x, y) - 1/target_apoapsis))


def radial_velocity_error(x, y, vx, vy):
    # m/s, returns the additional radial velocity needed to attain the target altitude
    return height_potential(x, y) - radial_velocity(x, y, vx, vy)


def write_header():
    with open("flight_log.txt", "w") as file:
        file.write(
            "time,s  mass,kg  attitude,deg  air_pressure,kPa  air_density,kg/m^3  drag_force,N  gravity_force,N  thrust,N  Fx,N  Fy,N  Vx,m/s  Vy,m/s  Vt,m/s  Vr,m/s  X,m  Y,m  Alpha,rad  "
            "Radial,m  Altitude,m  Apoapsis,m  Periapsis,m  Apoapsis_Velocity,m/s  Remaining_Orbital_Velocity,m/s  Delta-V,m/s  \n")


def write_data(x, y, vx, vy, mass_vehicle, current_time, flight_plan):
    with open("flight_log.txt", "a") as file:
        value_list = [current_time, mass_vehicle, attitude_deg(x, y, flight_plan), air_pressure(x, y), air_density(x, y), drag_force(x, y, vx, vy),
                  gravity_force(x, y, mass_vehicle), thrust(x, y, vx, vy, mass_vehicle), x_force(x, y, vx, vy, mass_vehicle, flight_plan), y_force(x, y, vx, vy, mass_vehicle, flight_plan), vx, vy, tangential_velocity(x, y, vx, vy),
                  radial_velocity(x, y, vx, vy), x, y, angular_position(x, y),
                  radial_position(x, y), altitude(x, y), apoapsis(x, y, vx, vy, mass_vehicle), periapsis(x, y, vx, vy, mass_vehicle), apoapsis_velocity(x, y, vx, vy, mass_vehicle), tangential_velocity_error(x, y, vx, vy, mass_vehicle),
                  delta_v(mass_vehicle)]
    for j in range(len(value_list)):
        file.write(str(value_list[j]) + ", ")
    file.write("\n")
# when development is done, we don't need this output. It can be trimmed down, reduced to every few time steps, or removed entirely if we can self generate a flight profile plan


def flight_analysis(flight_plan, state):
    # returns the fuel consumption of the simulated flight plan
    x = state[0]
    y = state[1]
    vx = state[2]
    vy = state[3]
    mass_vehicle = state[4]
    current_time = state[5]

    step_state = list(state)
    # write_header()
    # file = open("flight_log.txt", "a")
    # while (radial_velocity_error() > 5 or tangential_velocity_error() > 5 or altitude() <= top_of_atmosphere) and fuel_remain() > 0:
    for t in range(int(initial_fuel_mass / exhaust_flow / delta_time)):
        # increment time step & new states
        # write_data()
        if altitude(x, y) < flight_plan[0][altitude_points + 1]
            step_state = [x, y, vx, vy, mass_vehicle, current_time]  # saving the current state variable, so it doesn't have to be recalculated later
        current_time += delta_time
        mass_vehicle = new_vehicle_mass(x, y, vx, vy, mass_vehicle)
        vx = vx + x_force(x, y, vx, vy, mass_vehicle, flight_plan) / mass_vehicle * delta_time
        vy = vy + y_force(x, y, vx, vy, mass_vehicle, flight_plan) / mass_vehicle * delta_time
        x = x + vx * delta_time
        y = y + vy * delta_time
        if fuel_remain(mass_vehicle) <= 0 or (radial_velocity_error(x, y, vx, vy) < 5 and tangential_velocity_error(x, y, vx, vy, mass_vehicle) < 5 and altitude(x, y) >= top_of_atmosphere):
            break

    # write_data(x, y, vx, vy, mass_vehicle)  # need to record the last update.
    # file.close()
    fa = [fuel_remain(mass_vehicle), radial_velocity_error(x, y, vx, vy), tangential_velocity_error(x, y, vx, vy, mass_vehicle), step_state]
    return fa

def seed_flight_plan():
    fp = []
    for altitudes in range(resolution_bottom):
        fp[0][altitudes] = int(scale_height * math.log(1 + altitudes/count))
    last = fp[0][resolution_bottom - 1]
    for altitudes in range(resolution_bottom, resolution_bottom + extend):
        fp[0][altitudes] = top_of_atmosphere - last

    return
# this is a seed to generate the optimal flight profile as an algorithm output that the user can copy
# flight_plan_altitude = [0, 287, 590, 910, 1250, 1611, 1997, 2412, 2861, 3348, 3882, 4472, 5131, 5879, 6742, 7763, 9013, 10624, 12894, 16776, 25789, top_of_atmosphere, target_apoapsis]  # meters, flight plan altitude points
# these altitude points should be generated by a function based on the planet parameter: scale height
# altitude = int( scale_height * ln (1+i/20)) (i from 0 to 19)
# altitude = int( scale_height * ln (100)) (for element 20)
# flight_plan_attitude = [90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # degrees, flight plan attitude




# Begin main loop
# seed flight plan
# optimized_flight_plan = flight_plan_seed()
# set the state conditions
# state = initial_conditions()
initial_flight_plan = seed_flight_plan()

for altitude_points in range(len(initial_flight_plan[1])):
    for theta in range(last_angle, -10, -5):
        # analyze the flight plan
        fa1 = list(flight_analysis(flight_plan, start))
        # modify flight_plan[1][i]
        flight_plan[1][altitude_points] = theta
        # analyze the flight plan
        fa2 = list(flight_analysis(flight_plan, start))
        # if better --> modify the same element again, incrementing another angle
        # if worse/same --> revert the change, save the state (and start all analysis from there), go to next element
        if fa1[0] > 0 and 0 < fa2[0] <= fa1[0]:
            flight_plan[1][altitude_points] += 5
            state = list(fa1[3])
            break
        if fa1[0] <= 0 or fa2[0] <= 0:
            fa1_v_error = fa1[1]+fa1[2]
            fa2_v_error = fa2[1]+fa2[2]
            if fa2_v_error >= fa1_v_error:
                flight_plan[1][altitude_points] += 5
                state = list(fa1[3])
                break
        if theta == -5:
            state = list(fa1[3])