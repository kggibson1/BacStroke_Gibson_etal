'''
BacStroke.py
'''

import numpy as np
from numba import njit
import Functions as f
import os
import argparse

# inputs ######################################################################

# FUNCTIONS ###################################################################

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


@njit
def rotate_xy(x, y, theta):
    '''
    Apply rotation matrix:
        
        R(θ)=(cosθ ​−sinθ)
             (sinθ​  cosθ​)
             
    On r = (x
            y).
    
    Output is r' = R(θ)r

    Parameters
    ----------
    x : FLOAT
        x position in lab frame.
    y : FLOAT
        y position in lab frame.
    theta : float
        omega (rotational velocity) [rad/s] * dt [s].

    Returns
    -------
    x_new : FLOAT
        x position in rotated frame.
    y_new : FLOAT
        y position in rotated frame.

    '''
    
    # components of rotation matrix
    c = np.cos(theta)
    s = np.sin(theta)
    
    # calculate new coords from rotation matrix
    x_new = c * x - s * y # x rotated
    y_new = s * x + c * y # y rotated

    return x_new, y_new


@njit
def random_unit_vector():
    '''
    Generate random 3D unit vector on a unit sphere.
    
    x = sin(θ) cos(φ)
    y = sin(θ) sin(φ)
    z = cos(θ)

    where:
        θ ∈ [0, π]   (theta, polar angle)
        φ ∈ [0, 2π]  (phi, azimuthal angle)

    Returns
    -------
    np.array([x, y, z])
        random 3D unit vector on a unit sphere in cartesian coordinates.

    '''
    theta = np.random.uniform(0.0, np.pi)
    phi = np.random.uniform(0.0, 2.0*np.pi)

    x = np.sin(phi)*np.cos(theta)
    y = np.sin(phi)*np.sin(theta)
    z = np.cos(phi)

    return np.array([x, y, z])


@njit
def tumble_probability(dt, tumbling_rate, dps):
    '''
    Generate random number

    Parameters
    ----------
    dt : float
        Simulation timestep.
    tumbling_rate : integer
        how many times a bacterium tumble per second
    dps : integer
        how many decimal places are in the .

    Returns
    -------
    integer
        1 = bacterium tumbles, 0 = bacterium does not tumble.

    '''
    
    # generate random number 
    rand = np.random.uniform(0.0, 1.0 - dt)
    
    # round random number to same number of decimal places (dps) as dt
    rand = np.round(rand, dps)

    # threshold that allows bacterium to tumble 
    tumble_prob = 1.0 - (tumbling_rate * dt)
    
    # determine if bacterium tumbles
    # if randomly generated number greater than or equal to theshold a tumble
    # occurs
    # otherwise no tumble occurs
    return 1 if rand >= tumble_prob else 0

# MAIN SIMULATION #############################################################

@njit
def simulate_single_bacterium(
    pos, # initial position [x, y, z]
    swim_dir, # initial swimming direction [x, y, z]
    swim_speed, # swimming speed [m/s]
    bm, # bouyant mass
    radius, # radius of bacterium [m]
    dt, #timestep [s]
    num_steps, # total sim time / dt (how many steps simulation carries out)
    omega, # rotation rate of clinostat [rad/s]
    g, # gravitational acceleration [m/s^2]
    viscosity, # medium viscosity within clinostat [Pa/s]
    diffusion_coeff, # diffusion coeff of medium at room temperature [m^2/s]
    rot_diff_coeff, # [1/s]
    tumbling_rate, # how many times does a bacterium tumble [1/s]
    dps, # number of decimal places in the timestep
    R_outer, # outer radius of clinostat [m]
    R_inner, # inner radius of clinostat [m]
    H, # length of clinostat [m]
    use_centripetal, # is centripetal force present?
    use_boundaries, # are boundary conditions on?
    output_every, # data output interval [s]
    #seed
):
    #np.random.seed(seed)

    vel = np.zeros(3)
    rot_vel = np.zeros(3)
    swim_vel = np.zeros(3)
    term_vel = np.zeros(3)
    cent_vel = np.zeros(3)

    drag = 6.0*np.pi*radius

    # terminal velocity
    VTy = bm*g/(viscosity*drag)
    term_vel[1] = -VTy

    # output sizing
    num_saved = num_steps//output_every + 1
    traj = np.zeros((num_saved, 4))

    save_idx = 0 # track which sim step for saving into output array
    time = 0.0 # [s]
    
    # update all velocity components to determine new total velocity and 
    # position
    for t in range(num_steps + 1):

        # -------------------
        # rotation
        # -------------------
        
        x = pos[0]
        y = pos[1]

        # exact rotation over dt 
        theta = omega*dt
        x_rot, y_rot = rotate_xy(x, y, theta)
        
        rot_vel[0] = (x_rot-x)/dt
        rot_vel[1] = (y_rot-y)/dt
        rot_vel[2] = 0.0

        # -------------------
        # centripetal force
        # -------------------
        if use_centripetal:
            cent_vel[0] = omega*omega*x
            cent_vel[1] = omega*omega*y
            cent_vel[2] = 0.0

        else:
            cent_vel[0] = 0.0
            cent_vel[1] = 0.0
            cent_vel[2] = 0.0

        # -------------------
        # tumbling and rotational diffusion
        # -------------------
        tumble = tumble_probability(dt, tumbling_rate, dps)

        if tumble == 1:
            new_dir = random_unit_vector()
            norm = np.sqrt(np.dot(new_dir, new_dir))
            swim_dir[:] = new_dir/norm

        else:
            ex = swim_dir[0]
            ey = swim_dir[1]

            rotation = np.array([-ey*omega, ex*omega, 0.0])

            coeff = np.sqrt((2.0*rot_diff_coeff)/dt)
            noise = np.random.normal(0.0, 1.0, 3)

            outer = np.empty((3, 3))
            outer[0, :] = swim_dir[0]*swim_dir
            outer[1, :] = swim_dir[1]*swim_dir
            outer[2, :] = swim_dir[2]*swim_dir

            diffusion_term = coeff*((np.eye(3) - outer)@noise)

            dedt = rotation + diffusion_term
            new_dir = swim_dir + dedt*dt

            norm = np.sqrt(np.dot(new_dir, new_dir))
            swim_dir[:] = new_dir/norm

        swim_vel[:] = swim_speed*swim_dir

        # -------------------
        # translational diffusion
        # -------------------
        noise = np.random.normal(0.0, 1.0, 3)
        diffusion = noise * np.sqrt(2.0*diffusion_coeff/dt)

        # -------------------
        # total velocity
        # -------------------
        vel[:] = term_vel + rot_vel + swim_vel + diffusion

        # if use_centripetal:
        #     vel[0] += omega * omega * pos[0] * dt
        #     vel[1] += omega * omega * pos[1] * dt
        
        if use_centripetal:
            coeff = bm / (6.0*np.pi*viscosity*radius)
            vel[0] += coeff*pos[0]*omega**2
            vel[1] += coeff*pos[1]*omega**2

        # -------------------
        # update position from total velocity
        # -------------------
        pos[:] += vel*dt

        # -------------------
        # boundary conditions
        # -------------------
        if use_boundaries:
            r_xy = np.sqrt(pos[0]**2 + pos[1]**2)
            
            # outer radius boundary
            if r_xy > R_outer:
                scale = R_outer / r_xy
                pos[0] *= scale
                pos[1] *= scale
            
            # inner radius boundary
            elif r_xy < R_inner:
                scale = R_inner / r_xy
                pos[0] *= scale
                pos[1] *= scale
            
            # z boundaries
            if pos[2] < 0.0:
                pos[2] = 0.0
            elif pos[2] > H:
                pos[2] = H

        # -------------------
        # save current position and time
        # -------------------
        if t % output_every == 0:
            traj[save_idx, 0] = pos[0]
            traj[save_idx, 1] = pos[1]
            traj[save_idx, 2] = pos[2]
            traj[save_idx, 3] = time
            save_idx += 1

        time += dt

    return traj


# INPUT PARSE #################################################################

# parser = argparse.ArgumentParser()


# parser.add_argument("--config", 
#                     required=True, 
#                     help="Absolute path to config file")

# parser.add_argument("--output",
#                     required=True,
#                     help="Absolute path to output file")

# args = parser.parse_args()

# print("Config path:", args.config)
# print("Output folder:", args.output)

# config_path = args.config
# output_folder = args.output

#python BacStroke_numba_test_single_bacterium.py --config configs/omega=0.txt --output results/run1

# load config #################################################################
config_path = 'C:/Users/Kenzie/Documents/Education/Edinburgh_University/Summer_Masters/BacStroke/test_config.txt'

config = f.read_config(config_path)

dt = float(config[1])
total_time = float(config[2])
num_steps = int(round(total_time / dt))
dps = str(dt)[::-1].find('.')
print(type(dps))

omega = float(config[3]) * 2*np.pi/60
R = float(config[4])
r = float(config[5])
H = float(config[6])
density = float(config[7])
g = float(config[8])
rot_diff = float(config[9])
viscosity = float(config[10])
diffusion = float(config[11])
tumbling_rate = int(config[12])
use_centripetal = to_bool(config[13])
use_boundaries = to_bool(config[14])

#print(use_centripetal, type(use_centripetal), use_boundaries, type(use_boundaries))

output_dt = float(config[15])
output_every = int(round(output_dt / dt))


# --- load IC ---
with open(config[0], "r") as f_ic:
    p = f_ic.readline().split()

radius = float(p[1])
pos0 = np.array([p[2], p[3], p[4]], dtype=np.float64)
swim_speed = float(p[5])
organism_density = float(p[6])

bm = (4/3)*np.pi*density*(radius**3)*((organism_density/density)-1)

# initial direction
swim_dir0 = f.initialise_swimming_direction()

print('Starting Simulation')

# print simulation parameters
sim_params = {
    "Timestep, dt [s]": dt,
    "Simulation duration [s]": total_time,
    "Output interval [s]": output_dt,
    "Omega [rad/s]": omega,
    "Outer clinostat radius, R [m]": R,
    "Inner clinostat radius, Rₒ [m]": r,
    "Clinostat length, H [m]": H,
    "Medium density [kg/m³]": density,
    "Gravitational acc., g [m/s²]": g,
    "Rot. diff. coef., D_R [1/s]": rot_diff,
    "Viscosity, η [Pa s]": viscosity,
    "Transl. diff. coef, D_T [m²/s]": diffusion,
    "Tumbling rate, k [1/s]": tumbling_rate,
    "Centripetal force status": use_centripetal,
    "Boundary condition status": use_boundaries,
    "Organism radius, r [m]": radius,
    "Swimming speed, v_s [m/s]": swim_speed,
    "Organism density, ρ_b [kg/m³]": organism_density,
    "Buoyant mass [kg]": bm,
    "Initial position [x, y, z]": pos0,
    "Initial orientation": swim_dir0,
}

print("\n=== Simulation Parameters ===")
for key, value in sim_params.items():
    print(f"{key:30s}: {value}") # format and print each parameter
print("="*29)

traj = simulate_single_bacterium(
    pos0.copy(),
    swim_dir0.copy(),
    swim_speed,
    bm,
    radius,
    dt,
    num_steps,
    omega,
    g,
    viscosity,
    diffusion,
    rot_diff,
    tumbling_rate,
    dps,
    R,
    r,
    H,
    use_centripetal,
    use_boundaries,
    output_every,
    #seed=42
)

print("Sim finished")

# save file
def save_runs(N, traj, output_folder):
    
    if os.path.isdir(output_folder) == False:
        os.makedirs(output_folder)
    
    for i in range(N-1):
        filepath = os.path.join(output_folder, f"run_{i+1}.csv")
        
        np.savetxt(
            filepath,
            traj[:, i, :],
            delimiter=",")
        
    
def save_run(traj, output_folder):

    print("Entered save_run")

    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    filepath = os.path.join(output_folder, "run_1.csv")

    print("Saving to:", filepath)
    print("Trajectory shape:", traj.shape)

    np.savetxt(filepath, traj, delimiter=",")

    print("Save complete")
    
#np.save("numba_test_results".replace(".csv", ".npy"), traj)
       
output_folder = 'C:/Users/Kenzie/Documents/Education/Edinburgh_University/Summer_Masters/Gibson_et_al/BacStroke_output'

save_run(traj, output_folder)