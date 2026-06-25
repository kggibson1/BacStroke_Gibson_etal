# BacStroke.py

## Overview

BacStroke simulates the motion of a single motile bacterium inside a rotating clinostat. The bacterium is modelled as an active Brownian particle that can:

* Swim at a constant speed.
* Reorient through rotational diffusion.
* Undergo stochastic tumbling events.
* Sediment under gravity.
* Diffuse translationally through Brownian motion.
* Be advected by clinostat rotation.
* Experience centrifugal drift (optional).
* Interact with clinostat boundaries (optional).

The simulation outputs the bacterium trajectory as a CSV file containing position and time.

---
## Run BacStroke

The simulation is run from the command line as

```bash
python BacStroke.py --config CONFIG --output OUTPUT
```

where
* `CONFIG` is the path to a configuration file
* `OUTPUT` is the directory into which simulation results will be written.

For example,

```bash
python BacStroke.py --config test_config.txt --output Results
```

See test_config.txt for an example config file.

## Dependencies

External libraries:

NumPy and Numba are both required.

pip install numpy numba

Note: Functions.py must be in the same directory as BacStroke.py.

## Model description

Bacterial positions are computed numerically by integrating the Langevin equations of motion for bacterial position $r$ and orientation $\hat{e}_s$.

  $\dot{r}=v = v_g + v_R + v_c + v_T + v_s,$

  $\dot{\hat{e}}_s = \omega_R + \omega_T + \text{tumbles}.$

Integration is performed using a velocity Verlet algorithm with simulation timestep $\Delta t$.

The first term in $\dot{r}$ is the sedimentation velocity $v_g$, given by

  $v_g = -v_g \hat{e}_y,$

with $v_g=\Delta m g/\zeta$, $\Delta m = 4\pi a^3\Delta \rho/3$ and $\zeta = 6\pi \eta a$. The velocity $v_R$ due to advection by the clinostat rotation is

  $v_R = \omega \times r=v_R\hat{e}_\theta,$

The centrifugal velocity $v_c$ is

  $v_c = \frac{\Delta m \omega^2 r'}{\zeta} \hat{e}_{r},$

with $r'$ being the distance of an organisms centre of mass from the axis of rotation in the $x-y$ plane.

When considering diffusion and swimming, the Langevin equations are those of active Brownian particles. The organism's translational diffusivity $D_T=k_Bt/\zeta$ and rotational diffusivity $D_R=3D_T/4a^2$ are both standard expressions for the thermal diffusivity of a sphere. Diffusion is modeled as white noise, so the instantaneous thermal velocity can diverge. Hence, the thermal velocity $v_T$ is obtained by integrating over a small time interval, specifically the simulation timestep $\Delta t$. Hence, 

  $v_T=\sqrt{\frac{2D_T}{\Delta t}}\xi$

where the noise $\xi$ is uncorrelated 3D Gaussian white noise for which $\langle \xi\rangle = 0$ and $\langle \xi(t) \otimes \xi (t') \rangle = I \delta(t-t')$, where $\langle \cdot \rangle$ is the temporal average, $t$ and $t'$ are two time points, $\delta(\cdot)$ is the Kronecker-delta, $I$ is the 3D identity matrix and $\otimes$ indicates the tensor product. 

The bacterial swimming velocity $v_s$ is

  $v_s=v_s \hat{e}_s\,$

where $v_s$ is the swimming speed and $\hat{e}_s$ is the unit vector in the swimming direction, whose dynamics are given by \eq{Langevin equation - orientation}. Here, the first term, $\omega_R=\omega\times\hat{e}_s$ corresponds to the deterministic re-orientation produced by the clinostat rotating at rate $\omega=|\omega|$ around the axis $\hat{z}$. The $\omega_T$ term indicates stochastic (random) reorientation due to Brownian motion. As with the thermal velocity, the stochastic term is obtained as an average over the simulation timestep

  $\omega_T=\sqrt{\frac{2D_R}{\Delta t}}\xi\cdot(I - \hat{e}_s\otimes\hat{e}_s)$
  
where the combination $I-\hat{e}_s\otimes\hat{e}_s$ ensures that rotational diffusion does not change the magnitude of the orientation unit vector.

Tumbling is not included in the rotation vector $\dot{\hat{e}}$; instead, the `tumbles' term means that, at certain times, the bacterium will tumble, instantaneously reorienting. During tumbles, a new orientation is drawn uniformly from the unit sphere. These tumbling events occur randomly with rate $k$, i.e., within each small simulation time-step of length $\Delta t$ there is a probability $k\Delta t$ that the bacterium will tumble.

Boundary conditions are implemented such that a bacterium cannot pass through the walls of the clinostat. If the distance $r'$ of a bacterium's center of mass from the axis of rotation in the $x-y$ plane exceeds $R - a$, then its center of mass position is moved radially to lie $1.1a$ inside the wall. The radial component of the velocity is also removed such that the velocity of the bacterium becomes tangential to the clinostat's wall. An equivalent condition is applied at the inner boundary. Likewise, in the $x-z$ plane the same boundary conditions are applied when the $z$ component of the bacterium's center of mass is within 1$a$ from either end wall, i.e $z > H - a$ or $z < a$. 

## Input files

A Config and initial conditions file are required for the simulation. 

The config file contains environmental parameters such as clinostat dimensions and rotation, gravitational acceleration, tumbling frequency and more.

See `test_config.txt` for an example.

The initial conditions file contains the organism parameters, taking `Object_initial_conditions/Ecoli.txt` as an example:

1E-12 1E-6 0.5E-2 0.5E-2 0.5E-2 0 1027

Mass (kg), radius (m), $x_i$, $y_i$, $z_i$, swimming speed ($v_s$ (m/s)), organism density ( kg m(^{-3}))

where ($x_i$, $y_i$, $z_i$) is the initial position of the organism within the clinostat.


## Output

Simulation output is written as

```text
run_1.csv
```

with columns

| Column |    |
| ------ | -------------- |
| 1      | x position (m) |
| 2      | y position (m) |
| 3      | z position (m) |
| 4      | time (s)       |

Data are written at the output interval specified in the configuration file.

## Before running

All quantities should be in SI units. Lengths should be specified in metres, times in seconds, densities in kg m(^{-3}), viscosities in Pa s and diffusion coefficients in their standard SI units.

Ensure that the initial position lies inside the clinostat dimensions,

$R_{\mathrm{inner}} < \sqrt{x^2+y^2} < R_{\mathrm{outer}}$,

with

$0 \le z \le H.$

If boundary conditions are enabled and the initial position lies outside the allowed region, the particle will immediately be projected onto a boundary.

Ensure the timestep is small enough to resolve any behaviour needing to be observed e.g tumbling, thermal or ballistic motion.

For simulations including rotation we generally followed the rule:

$\Delta t \lesssim \frac{1}{10}\frac{2\pi}{\omega}$

but this may need to be smaller when tumbling and diffusion (behaviours that occur on much smaller timescales) are desired to be observed.

## Reproducing published results

The repository includes example config files (See /Configs) used to produce the results presented in [the accompanying publication (hopefully!)]. Because the model has stochastic processes, individual trajectories will differ between runs, but on averaged over multiple runs results should remain similar.

