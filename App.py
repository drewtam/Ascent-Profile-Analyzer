import math
# import json
from math import *

# these hard coded constants are proper
gravity_const = 6.6743 * (10**-11)  # m^3/kg/s^2, universal gravitational constant
gas_const = 0.287053  # kJ/kg/K, universal gas constant
g_accel = 9.81  # m/s^2, acceleration due to gravity at the surface (should only be used for normalization)

# these parameters need a user interface front end
# the drag coefficient is especially difficult, and should have a separate tool to help guide the user to determine the value
# it is a known inaccuracy in this model, that if the vehicle attitude does not equal vehicle velocity direction, then the coefficient of drag will be different and could be significantly higher. It is possible to build a more accurate model of KSP drag,  but not a priority
# this model assumes that the vehicle aerodynamics are stable & neutral, with no flipping, no lift.
# this input is currently shown as a single stage launch, but a new set of inputs for additional stages will be necessary. This complication should use all the same math, but the main loop just breaks up which global constants to use
engine_thrust = 10000  # N, maximum engine thrust, user set parameter from the rocket design
drag_coefficient = 0.05  # m^2, drag coefficient "Cd * Area", user set parameter from the rocket design
initial_mass_vehicle = 500  # kg, user set parameter from the rocket design
initial_fuel_mass = 410  # kg, user set parameter from the rocket design
specific_impulse = 340  # seconds, user set parameter from the rocket design
exhaust_flow = engine_thrust / (g_accel * specific_impulse)  # kg/s, derived parameter of the rocket design

# these parameters need a user interface front end
# prefer a drop down menu for selecting all the Kerbal planets, with preloaded values from a table.
# prefer to have an additional option for custom inputs
radius_planet = 600000.00  # meters, user set parameter from the planet of origin (Kerbal is 600km)
mass_planet = 5.2915158 * (10**22)  # kg, user set parameter from the planet of origin (Kerbal is 5.2915158e22 kg)
initial_pressure = 101.3  # kPa, atmospheric initial pressure, user set parameter from the planet of origin
scale_height = 5600  # meters, atmospheric pressure constant of KSP's atmosphere model, user set parameter from the planet of origin
temperature = 293  # Kelvin, atmospheric temperature, user set parameter from the planet of origin (Kerbal is 20C at ground)
top_of_atmosphere = 70000  # meters, highest altitude that atmosphere affects flight & drag, user set parameter from the planet of origin (Kerbal is 70km)

delta_time = 0.05  # seconds, time step, Controls simulation accuracy vs computational cost. ???maybe a user set parameter??? ???maybe dynamically generated based on velocity or some other error calc???

# this is a temporary method, prefer to generate the optimal flight profile as an algorithm output that the user can copy
flight_plan_altitude = [0, 3000, 8000, 12000, 20000, 80000]  # meters, flight plan altitude points, user set values for their launch plan
flight_plan_attitude = [85, 72, 40, 20, 0, 0]  # degrees, flight plan attitude, user set values for their launch plan
target_apoapsis = flight_plan_altitude[-1]  # meters, target apoapsis is assumed to the be the last element of the flight plan
# need a check to verify flight plan is appropriately sized to each other, & the last element is above planet's top of atmosphere


# initial conditions
# these should be the only global variables modified in main loop
mass_vehicle = initial_mass_vehicle  # kg, vehicle mass
current_time = 0.00  # seconds, simulation time, initialization
x_pos = 0.0  # meters, cartesian position of the vehicle relative to planet center (origin), initialization
y_pos = radius_planet  # meters, cartesian position of the vehicle relative to planet center (origin), initialization
x_velocity = 0.0  # m/s, cartesian velocity of the vehicle, initialization
y_velocity = 0.0  # m/s, cartesian velocity of the vehicle, initialization


def tangential_velocity():
    # m/s, returns velocity perpendicular to the position vector, assumes planar motion
    v = velocity()
    vr = radial_velocity()
    return sqrt(v**2 - vr**2)


def radial_velocity():
    # m/s, returns velocity co-linear to the position vector, assumes planar motion
    vx = x_velocity
    vy = y_velocity
    alpha = angular_position()
    return vx * cos(alpha) + vy * sin(alpha)


def altitude():
    # meters, returns current altitude from the radial position
    return radial_position() - radius_planet


def attitude_deg():
    # degrees, angle, returns the attitude from the flight plan. Calculates the current altitude, reads the flight plan.
    alt = altitude()
    a = 0
    for i in range(len(flight_plan_altitude)):
        if alt >= flight_plan_altitude[i]:
            a = i
    return flight_plan_attitude[a]


def attitude_rad():
    # radians, angle, calculate the attitude in radians from the flight plan in degrees
    return attitude_deg() * pi / 180


def air_pressure():
    # kPa, returns the air pressure from current altitude
    return initial_pressure * exp(-(altitude()) / scale_height)


def air_density():
    # kg/m^3, returns the air density
    return air_pressure() / (gas_const * temperature)


def angular_position():
    # radians, returns the angular position from cartesian coordinates, assumes planar position
    x = x_pos
    y = y_pos
    if x == 0:
        return pi/2
    else:
        return atan(y / x)


def radial_position():
    # meters, returns the radial position from cartesian coordinates, assumes planar position
    x = x_pos
    y = y_pos
    return sqrt(x**2 + y**2)


def fuel_remain():
    # kg, returns the amount of fuel left
    return initial_fuel_mass - (initial_mass_vehicle - mass_vehicle)


def thrust():
    # Newtons, returns the engine thrust. Turns engine thrust on/off based on calculated apoapsis. Engine will burn until target apoapsis is expected, prevents overshoot.
    if apoapsis_altitude() >= target_apoapsis or fuel_remain() <= 0:
        return 0
    else:
        return engine_thrust


def velocity():
    # m/s, returns the current scalar velocity
    vx = x_velocity
    vy = y_velocity
    return sqrt(vx**2 + vy**2)


def orbital_energy():
    # kg * m^2 / s^2, returns conserved orbital energy (LaGrangian), requires vehicle mass, scalar velocity, and radial position
    m = mass_vehicle
    v = velocity()
    r = radial_position()
    return 1/2 * m * (v**2) - gravity_const * mass_planet * m / r


def angular_momentum():
    # kg*m^2/s, returns orbital angular momentum
    r = radial_position()
    vt = tangential_velocity()
    m = mass_vehicle
    return r * vt * m


def reduced_mass():
    # kg, effective mass in two body newtonian mechanics
    m = mass_vehicle
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
    mr = reduced_mass()
    f = central_force()
    return sqrt(1 + 2 * orbital_energy() * (angular_momentum()**2) / (mr * (f**2)))


def semi_major_axis():
    # meters, returns the radius of the semi major axis of the vehicle orbit (ellipse) around the planet. Semi major axis is measured from the center of the ellipse, whereas the apoapsis is measured from the center of the planet.
    v = velocity()
    r = radial_position()
    return -gravity_const * mass_planet / (2 * (v**2 / 2 - gravity_const * mass_planet / r))


def apoapsis():
    # meters, returns the radius of the orbital apoapsis
    sma = semi_major_axis()
    ecc = eccentricity()
    return sma * (1 + abs(ecc))


def apoapsis_altitude():
    # meters, returns the radius minus the planetary diameter (altitude of the apoapsis)
    return apoapsis() - radius_planet


def periapsis():
    # meters, returns the radius of the orbital periapsis
    sma = semi_major_axis()
    ecc = eccentricity()
    return sma * (1 - abs(ecc))


def periapsis_altitude():
    # meters, returns the radius minus the planetary diameter (altitude of the apoapsis)
    return periapsis() - radius_planet


def apoapsis_velocity():
    # m/s, returns velocity at apoapsis
    sma = semi_major_axis()
    ecc = eccentricity()
    return sqrt((1-ecc) * gravity_const * mass_planet / ((1+ecc) * sma))


def apo_gravity_accel():
    # m/s^2, returns the free fall acceleration due to gravity at apoapsis
    apo = apoapsis()
    return gravity_const * mass_planet / (apo**2)


def circular_orbital_apoapsis_velocity():
    # m/s, returns the required orbital velocity for a circular orbit at apoapsis
    aga = apo_gravity_accel()
    apo = apoapsis()
    return sqrt(aga * apo)


def velocity_deficit():
    # m/s, returns the additional velocity needed to circularize orbit at apoapsis
    orbv = circular_orbital_apoapsis_velocity()
    apov = apoapsis_velocity()
    return orbv - apov


def drag_force():
    # Newtons, returns atmospheric drag force scalar
    v = velocity()
    if altitude() >= top_of_atmosphere:
        return 0
    else:
        return 1/2 * drag_coefficient * air_density() * (v**2)


def gravity_force():
    # Newtons, returns scalar force due to gravity
    r = radial_position()
    return gravity_const * mass_planet * mass_vehicle / (r**2)


def x_force():
    # Newtons, returns sum of forces in the cartesian x direction, assumes planar forces
    alpha = angular_position()
    ftx = thrust() * cos(attitude_rad() - alpha + pi/2)  # thrust in the x direction
    v = velocity()
    if v == 0:
        fdx = 0
    else:
        fdx = drag_force() * cos(acos(x_velocity / v))  # drag in the x direction
    fgx = gravity_force() * cos(alpha)  # gravity in the x direction
    return ftx - fdx - fgx


def y_force():
    # Newtons, returns sum of forces in the cartesian y direction, assumes planar forces
    alpha = angular_position()
    fty = thrust() * sin(attitude_rad() - alpha + pi/2)  # thrust in the y direction
    v = velocity()
    if v == 0:
        fdy = 0
    else:
        fdy = drag_force() * sin(asin(y_velocity / v))  # drag in the y direction
    fgy = gravity_force() * sin(alpha)  # gravity in the y direction
    return fty - fdy - fgy


def delta_v():
    # m/s, returns the "delta-v" left remaining to the vehicle
    m1 = mass_vehicle  # current vehicle mass
    m2 = initial_mass_vehicle - initial_fuel_mass  # dry vehicle mass
    i = specific_impulse
    g = g_accel
    return math.log(m1/m2) * i * g


# calculate the needed radial velocity needed to achieve the target height
# calculate the needed tangential velocity to maintain orbit
# calculate the deficit (target - actual) for each velocity component
# calculate the change in deficit (deltaDV = sqrt(dvr^2 + dvt^2) from last step to current step
# divide by change in deltaV (maybe not necessary)
# iterate current step through several attitude angles which maximizes deficit reduction (possible sign error checking here, also needs to constrain between 90deg and 0deg)
# record optimal attitude through flight; this probably doesn't need to be calculated every time step
# average optimal attitude in 5km increments, output as separate flight plan recommendation


def write_header():
    with open("flight_log.txt", "w") as file:
        file.write("time,s  mass,kg  attitude,deg  air_pressure,kPa  air_density,kg/m^3  drag_force,N  gravity_force,N  thrust,N  Fx,N  Fy,N  Vx,m/s  Vy,m/s  Vt,m/s  Vr,m/s  X,m  Y,m  Alpha,rad  "
                   "Radial,m  Altitude,m  Apoapsis,m  Periapsis,m  Apoapsis_Velocity,m/s  Remaining_Velocity,m/s  Delta-V,m/s  \n")


def write_data():
    with open("flight_log.txt", "a") as file:
        value_list = [current_time, mass_vehicle, attitude_deg(), air_pressure(), air_density(), drag_force(), gravity_force(), thrust(), x_force(), y_force(), x_velocity, y_velocity, tangential_velocity(), radial_velocity(), x_pos, y_pos, angular_position(),
                      radial_position(), altitude(), apoapsis(), periapsis(), apoapsis_velocity(), velocity_deficit(), delta_v()]
        for i in range(len(value_list)):
            file.write(str(value_list[i]) + ", ")
        file.write("\n")
# this is opening and closing the file every line, fix it!
# when development is done, we don't need this output. It can be trimmed down, reduced to every few time steps, or removed entirely if we can self generate a flight profile plan


# record initial conditions
write_header()

# Begin main loop (you can take the boy out of C, but can't take the C out of the boy)
while (apoapsis_altitude() <= target_apoapsis or altitude() <= top_of_atmosphere) and fuel_remain() > 0:
    # increment time step & new states
    write_data()
    current_time += delta_time
    mass_vehicle = new_vehicle_mass()
    x_velocity = x_velocity + x_force() / mass_vehicle * delta_time
    y_velocity = y_velocity + y_force() / mass_vehicle * delta_time
    x_pos = x_pos + x_velocity * delta_time
    y_pos = y_pos + y_velocity * delta_time
