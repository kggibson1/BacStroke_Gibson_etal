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

def to_bool(x):
    '''
    Convert string value of True and False to bool in config files for toggle 
    centripetal force and boundaries. Catches spelling mistakes.

    Parameters
    ----------
    x : STR
        toggle indicator of parameter from config being on or off.

    Raises
    ------
    ValueError
        Spelling mistake present (input is not capital variation of True
                                  or False).

    Returns
    -------
    bool
        Boolean True or False toggle for parameters in config file.
    '''
    
    # remove whitespace and put all letters to lower case
    value = x.strip().lower()

    # return corresponding boolean
    if value == "true":
        return True
    if value == "false":
        return False
    
    # indicate possible spelling error
    raise ValueError(f"Expected True or False, got {x!r}")


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


# not used
def sample_from_hollow_cylinder(radius_inner, radius_outer, height):
    
    u = np.random.uniform(radius_inner**2, radius_outer**2)
    r = np.sqrt(u)

    theta = np.random.uniform(0, 2*np.pi)
    z = np.random.uniform(0, height)

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    return np.array([x, y, z])           


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
    
    