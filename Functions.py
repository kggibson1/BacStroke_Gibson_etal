'''
This script contains functions used throughout the other scripts contained within
the BacStroke Repository.
'''

# Imports #####################################################################

# modules
import numpy as np
import pandas as pd
import math
import random
import os
import matplotlib.pyplot as plt

###############################################################################

# def initialise_swimming_direction():
#     '''
#     Initialises an initial unit vector direction of bacteria swimming
#     by generating two random angles on the unit sphere and converting
#     to cartesian coordinates.
#     '''
    
#     # generating random angles from the circular coordinate system
#     # cos theta is generated as this is uniform unlike theta
#     costheta = np.random.uniform(0, 1)
#     phi = np.random.uniform(0, 2*np.pi)
    
#     # components of 3D unit vector on unit circle
#     # converting from circular to cartesian coordinates
#     x = np.sin(phi)*costheta
#     y = np.sin(phi)*np.sqrt(1-costheta**2)
#     z = np.cos(phi)
    
#     # random unit vector generated from the circular coordinate system
#     return np.array([x, y, z])

# def round_to_first_non_zero(x):
#     '''
#     Rounds floats to first non zero didget.
    
#     x: float
#     float to be rounded
#     '''
#     if x == 0:
#         return 0
#     # Calculate the magnitude (order of magnitude)
#     magnitude = -int(f"{x:e}".split('e')[-1])
#     # Round the number to that magnitude
#     return round(x, magnitude)
# def round_to_first_non_zero(x):
#     if x == 0:
#         return '0'
        
#     # Calculate the magnitude (order of magnitude)
#     magnitude = -int(f"{x:e}".split('e')[-1])
    
#     # Round the number to that magnitude
#     rounded_value = round(x, magnitude)
    
#     # Format the rounded value to prevent scientific notation
#     formatted_value = f"{rounded_value:.{magnitude}f}"
    
#     # Remove trailing zeros after the decimal point
#     return formatted_value.rstrip('0').rstrip('.')


# should work with numbers larger than 1 now
def round_to_first_non_zero(x):
    if x == 0:
        return '0'
    
    # Determine the magnitude
    magnitude = -int(f"{x:e}".split('e')[-1])
    
    if x >= 1:
        # For numbers >= 1, we round to the nearest whole number
        rounded_value = round(x)
        # Convert to string without scientific notation
        formatted_value = str(rounded_value)
    else:
        # For numbers < 1, we keep the same logic
        rounded_value = round(x, magnitude)
        # Format the rounded value
        formatted_value = f"{rounded_value:.{magnitude}f}"
    
    # Remove trailing zeros after the decimal point
    return formatted_value.rstrip('0').rstrip('.')



def initialise_swimming_direction():
    '''
    Initialises an initial unit vector direction of bacteria swimming
    by generating two random angles on the unit sphere and converting
    to cartesian coordinates.
    '''
    
    # generating random angles from the circular coordinate system
    # cos theta is generated as this is uniform unlike theta
    theta = np.random.uniform(0, np.pi)
    phi = np.random.uniform(0, 2*np.pi)
    
    # components of 3D unit vector on unit circle
    # converting from circular to cartesian coordinates
    x = np.sin(phi)*np.cos(theta)
    y = np.sin(phi)*np.sin(theta)
    z = np.cos(phi)
    
    # random unit vector generated from the circular coordinate system
    return np.array([x, y, z])


# def initialise_swimming_direction():
#     '''
#     Initialises an initial unit vector direction of bacteria swimming
#     by generating two random angles on the unit sphere and converting
#     to cartesian coordinates.
#     '''
    
#     # generating random angles from the circular coordinate system
#     # cos theta is generated as this is uniform unlike theta
#     theta = np.random.uniform(0, np.pi)
#     phi = np.random.uniform(0, 2*np.pi)
    
#     # components of 3D unit vector on unit circle
#     # converting from circular to cartesian coordinates
#     x = np.sin(phi)*np.cos(theta)
#     y = np.sin(phi)*np.sin(theta)
#     z = np.cos(phi)
    
#     # random unit vector generated from the circular coordinate system
#     return np.array([x, y, z])


def tangential_velocity(current_position, planar_position, current_velocity):
    '''
    Inakes the current velocity of a bacterium and returns the tangential component
    of that velocity.
    
    :param current_position: numpy array [3], xyz components of bacteriums current position
    :param planar_position
    '''
    
    # radial unit vector
    radial_unit_vec = planar_position / np.linalg.norm(planar_position)
    
    # radial component of velocity
    radial_magnitude = np.dot(radial_unit_vec, current_velocity)
    
    # 
    tangential_vel = current_velocity - (radial_magnitude*radial_unit_vec)
    
    return tangential_vel


def convert(lst):
    '''
    This function takes a list and converts it into a space separated 
    string.
    
    :param lst: list of any length
    '''
    return ' '.join(lst)


def change_swimming_speed(initial_conditions_file_path, swimming_speeds):
    '''
    This function opens and changes the swimming speed within an inital
    condition file. The swimming velocity is a supplied parameter, if there
    is multiple lines the initial conditions file then the number of elements
    in the swimming_velocities list or array must match the length of lines in 
    the initial conditions file.
    
    :param initial_conditions_file_path: string, file path to initial conditions file
    :param swimming_speeds: array or list (length equal to number of lines in
    initial conditions file) containing the swimming velocities that are to be
    appended.
    
    '''
    
    # opening initial conditions file
    with open(initial_conditions_file_path, 'r') as file:
        # read each line in configuration file
        data = file.readlines()       
    
    # list of data where each element is a line in the initial conditions file
    lists = []    
    
    # splitting each line and reassigning the desired swimming speed value
    for i in range(len(data)):
        params = data[i].split()
        params[-2] = str(swimming_speeds[i]) #+ '\n'
        
        # converting split list back into string format
        combined = convert(params)
        # putting params back into correct format for putting back into file
        lists.append(combined)
        
    # and write everything back to config file
    with open(initial_conditions_file_path, 'w') as file:
        file.writelines(lists)
    
        
def change_starting_coords(initial_conditions_file_path, starting_positions):
    '''
    This function opens and changes the starting coordinates within an inital
    condition file. The swimming velocity is a supplied parameter, if there
    is multiple lines the initial conditions file then the number of elements
    in the swimming_velocities list or array must match the length of lines in 
    the initial conditions file.
    
    :param initial_conditions_file_path: string, file path to initial conditions file
    :param starting_positions: n x 3 array where n is the number of lines in the 
    initial conditions file, positions in xyz format.
    '''
    
    # opening initial conditions file
    with open(initial_conditions_file_path, 'r') as file:
        # read each line in configuration file
        data = file.readlines()       
    
    # list of data where each element is a line in the initial conditions file
    lists = []    
    
    # splitting each line and reassigning the desired swimming speed value
    for i in range(len(data)):
        params = data[i].split()
        params[2] = str(starting_positions[0]) # x coord
        params[3] = str(starting_positions[1]) # y coord
        params[4] = str(starting_positions[2]) # z coord
        params[6] = str(params[6]) + '\n' # making sure new line is added
        #print(params)
        
        # converting split list back into string format
        combined = convert(params)
        # putting params back into correct format for putting back into file
        lists.append(combined)
        
    # and write everything back to config file
    with open(initial_conditions_file_path, 'w') as file:
        file.writelines(lists)  
        
        
def radial_velocity(planar_position, current_velocity):
    '''
    This function calculates the magnitude and direction of the radial velocity 
    component of a bacteriums current velocity. 
    
    :param current_position: numpy array [3], xyz components of bacteriums current position
    :param planar_position: numpy array [3], xyz components of bacteriums position with z component set to 0
    :param current_velocity: numpy array [3], xyz components of bacteriums velocity
    
    each of the above parameters need to be from the same timestep
    
    :returns radial_magnitude: magnitude of bacteriums velocity vector in the radial direction [same units as velocity vector]
    :returns radial_unit_vec: direction of bacterium in radial direction from its current position
    '''

    # set velocity and orientation - !!!
    # radial unit vector
    planar_radial_unit_vec = planar_position / np.linalg.norm(planar_position)
    
    # radial component of velocity
    radial_magnitude = np.dot(planar_radial_unit_vec, current_velocity)
    
    return radial_magnitude, planar_radial_unit_vec


def get_unit_vector(vector):
    '''
    Get the unit vector of any unit vector in any coordinate system.
    
    :param vector: array of any length that the unit vector is desired.
    
    :returns unit: array of length vector, unit vector of input vector 
    '''
    
    mag = np.linalg.norm(vector)
    
    unit = vector/mag
    unit_mag = np.linalg.norm(unit)

    return unit, unit_mag


def peclet(swimming_vel, terminal_vel, centripetal_vel, translational_diffusion_coefficent, bacterium_radius):
    '''
    Calculates the three peclet numbers assositated with the bacterium within
    the system
    '''
    
    # common factor present in each peclet number
    common_factor = bacterium_radius/translational_diffusion_coefficent
    
    Ps = swimming_vel*common_factor # swimming peclet number
    Pg = terminal_vel*common_factor # gravitational peclet number
    Pc = centripetal_vel*common_factor # centripetal peclet number
    
    return Ps, Pg, Pc


# def sample_hollow_cylinder(r, R, h):
#     '''
#     This function randomly samples coordinated from a hollow cylindrical system
#     and then converts those coorindates into an xyz coordinate frame.

#     Parameters
#     ----------
#     r : float
#         Inner radius of cylinder [m].
#     R : float
#         Outer radius of cylinder [m].
#     h : float
#         length of cylinder [m].

#     Returns
#     -------
#     xyz : numpy array, shape 3
#         xyz coordinates generated from a uniform distribution within a hollow
#         cylinder.

#     '''
    
#     # generate random angle and radius (from central axis) from a uniform 
#     # distribution - this gets the radius on circular face of the cylinder
#     theta = random.uniform(0, 2*np.pi)
#     radii = random.uniform(r, R) # radius is larger than r, smaller than R
    
#     # get x and y coords from theta and radii
#     x = radii*np.cos(theta)
#     y = radii*np.sin(theta)
    
#     # generate z coord from uniform distribution
#     z = random.uniform(0, h)
    
#     # turn coords into a numpy array
#     xyz = np.array([x, y, z])
    
#     return xyz

# def sample_from_hollow_cylinder(radius_inner, radius_outer, height):
    
#     # Define bounding box enclosing the hollow cylinder
#     min_x = -radius_outer
#     max_x = radius_outer
#     min_y = -radius_outer
#     max_y = radius_outer
#     min_z = 0
#     max_z = height
    
#     while True:
#         # Sample a point uniformly in the bounding box
#         x = np.random.uniform(min_x, max_x)
#         y = np.random.uniform(min_y, max_y)
#         z = np.random.uniform(min_z, max_z)
        
#         # Compute radial distance from z-axis
#         radial_dist = np.sqrt(x**2 + y**2)
        
#         # Reject points outside the annular region defined by inner and outer radii
#         if radial_dist >= radius_inner and radial_dist <= radius_outer:

#             return np.array([x, y, z])
        
def sample_from_hollow_cylinder(radius_inner, radius_outer, height):
    
    u = np.random.uniform(radius_inner**2, radius_outer**2)
    r = np.sqrt(u)

    theta = np.random.uniform(0, 2*np.pi)
    z = np.random.uniform(0, height)

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    return np.array([x, y, z])


def measure_steady_state_time(folder_path, position_file_name):
    '''
    Measure the average steady state time of a collection of datasets 
    corresponding to a single configuration.(in the y direction)
    
    This is done by measuring the standard deviaiton (sigma) of the set at the end of 
    the dataset. Then find where a line of y = <y> + 1sigma passes through the average of the
    dataset. The time at which this occurs is taken as the steady state time.
    
    Parameters
    ----------
    
    folder_path: string
        path to folder containing subfolders with positional files
        
    time_path: string
        path to time file corresponding to all positions file (one file must work
                                                               for all positions files)
        
    position_file_name
        name of position file. Must be the same for each run

    
    Returns
    -------
    ind: integer
        index corresponding to steady state time of dataset

    '''
    
    # get all of the subfolders within folder path
    runs = os.listdir(folder_path)
    no_runs = len(runs)
    
    # find out how many data points are in each file
    position_path = folder_path + '/' + runs[0] + '/positions.csv'
    data = pd.read_csv(position_path)
    positions = np.array(data)
    no_pos = len(positions)
    
    
    # compile each set of y values into large array
    y = np.zeros([no_pos, no_runs]) # rows, cols
    #print(no_pos)
    
    # read positions file of each run & add y values to y array
    for i in range(no_runs):
        
        # path to positions file
        ith_positions_path = folder_path + '/' + runs[i] + '/positions.csv'
        
        # open positions file
        positions_r = pd.read_csv(ith_positions_path) # read csv with pandas
        positions = np.array(positions_r)
        
        # store the y data
        y[:, i] = positions[:, 1]
    
    #print(y)
    
    # get mean of each row
    ymean = np.mean(y, axis = 1)
    
    p90 = int((no_pos/100)*95)
    #print(p90)
    #print(y[p90:, :])
    #print(y[p90:, :].shape)
    
    # measure mean & standard deviation of last 90% of points (take mean again to get a singular value, but all should be very similar anyway)
    mean90 = np.mean(np.mean(y[p90:, :], axis = 1))
    stddev90 = np.mean(np.std(y[p90:, :], axis = 1)) # much more like it value wise
    print(stddev90)
    
    c = mean90 + (np.abs(stddev90))
    print(c)
    arr = np.zeros(no_pos) + (c) # straight line
    
    # this could be done more accuratley using interpolation of the y data 
    # to get a function that we could equate the straight line too but this is
    # as accurate that is needed at the moment
    
    # return a boolean array stating if the tolerance line is close to the real 
    # data
    mask = np.isclose(arr, ymean, atol = 10E-5)
    
    # index of first true
    ind = np.argmax(mask)
    
    # ignore next 3 lines ### 
    # find first index where y data drops below the 1 sigma threshold (as we
    # know that the non 0 g follow rough exponential distributions)
    #crossing_index = np.argmax(ymean <= (c))
    
    # read in time data
    #time = pd.read_csv(time_path)
    #times = np.array(time)
    #print(times)
    #print(len(times), len(arr), len(ymean))
    
    # steady state time is the time corresponding to the first true index
    #steady_state_time = times[ind]
    
    #print(ind,times[ind])
    
    #plt.scatter(times, arr, s = 2, zorder = 20)
    #plt.scatter(times, ymean, s = 2)
    
    return ind
    
    
#measure_steady_state_time('Studies/Exp1/Exp1_data/Data/g=0.098,vs=0.0', 'Studies/Exp1/Exp1_data/time.csv','positions.csv')

def steady_state_time(y, no_pos, time_path):
    
    '''
    y: numpy array
        
    no_pos: integer
        number of positions recorded per run
    '''
    
    # get mean of each row
    ymean = np.mean(y, axis = 1)
    
    p90 = int((no_pos/100)*95)
    #print(p90)
    #print(y[p90:, :])
    #print(y[p90:, :].shape)
    
    # measure mean & standard deviation of last 90% of points (take mean again to get a singular value, but all should be very similar anyway)
    mean90 = np.mean(np.mean(y[p90:, :], axis = 1))
    stddev90 = np.mean(np.std(y[p90:, :], axis = 1)) # much more like it value wise
    print(stddev90)
    
    c = mean90 + (2*np.abs(stddev90))
    print(c)
    arr = np.zeros(no_pos) + (c) # straight line
    
    # this could be done more accuratley using interpolation of the y data 
    # to get a function that we could equate the straight line too but this is
    # as accurate that is needed at the moment
    
    # return a boolean array stating if the tolerance line is close to the real 
    # data
    mask = np.isclose(arr, ymean, atol = 10E-5)
    
    # index of first true
    ind = np.argmax(mask)
    
    # ignore next 3 lines ### 
    # find first index where y data drops below the 1 sigma threshold (as we
    # know that the non 0 g follow rough exponential distributions)
    #crossing_index = np.argmax(ymean <= (c))
    
    # read in time data
    #time = pd.read_csv(time_path)
    #times = np.array(time)
    #print(times)
    #print(len(times), len(arr), len(ymean))
    
    # steady state time is the time corresponding to the first true index
    #steady_state_time = times[ind]
    
    #print(ind,times[ind])
    
    #plt.scatter(times, arr, s = 2, zorder = 20)
    #plt.scatter(times, ymean, s = 2)
    
    return ind
           

def read_config(config_path):
    '''
    Generates an array of the constants contained within the configuration file. 

    Parameters
    ----------
    config_path : string
        path to configuration file

    Returns
    -------
    config : Numpy array
        constants read from config file

    '''
    # opening config file containing all file paths needed to execute code
    file = open(config_path, "r")
    
    # reading every entry from configuration file (containing constants etc)
    config = file.readlines()[1::3]
    
    # removing trailing new line (\n)
    for i in range(len(config)):
        config[i] = config[i].rstrip()
        
    return config


def read_ic(ic_path):
    '''
    Read the parameters from an initial condition file and have them returned 
    in an array format. 

    Parameters
    ----------
    ic_path : STRING
        absolute path to an initial condition file.

    Returns
    -------
    params : LIST
        List of each item in the initial condition file.

    '''
    
    # opening initial conditions file
    with open(ic_path, 'r') as file:
        # read each line in configuration file
        data = file.readlines()       
    
    # list of data where each element is a line in the initial conditions file
    lists = []    
    
    # splitting each line and reassigning the desired swimming speed value
    for i in range(len(data)):
        params = data[i].split()
    
    return params
    

def edit_config(config_path, line, new_val):
    '''
    Edit a line in a configuration file, current file is written over. 

    Parameters
    ----------
    config_path : string
        path to configuration file
    line : int
        which val should be replaced
        0 = ic_path
        1 = timestep
        2 = simulation length
        etc...
    new_val : varied
        parameter to replace value in specified line

    '''
    # opening config file containing all file paths needed to execute code
    with open(config_path, "r") as file:
        all_lines = file.readlines()
    
    # reading every entry from configuration file (containing constants etc)
    config = all_lines[1::3]
    
    # removing trailing new line (\n)
    for i in range(len(config)):
        config[i] = config[i].rstrip()
    
    # assign new value to the specified line
    if type(new_val) == str:
        config[line] = new_val
        
    else: config[line] = f'{new_val}'
    
    # Write back to config file
    value_index = 0
    with open(config_path, "w") as file:
        for i in range(len(all_lines)):
            # If this is a value line (every third line starting from index 1)
            if i % 3 == 1:
                file.write(f"{config[value_index]}\n")
                value_index += 1
            else:
                # Write comment lines as they were
                file.write(all_lines[i])
    
    return config
    
#edit_config('C:/Users/Kenzie/Documents/Education/Edinburgh_University/Summer_Masters/BacStroke/test_config.txt', -4, 200000)    
    
    
    
    
    