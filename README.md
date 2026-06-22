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

Example configuration files are included with the repository.

## Model description

The bacterium position **r** evolves according to


$\frac{d\mathbf r}{dt}$ = $\mathbf v_g$ + $\mathbf v_R$ + $\mathbf v_c$ + $\mathbf v_T$ + $\mathbf v_s$

where $\mathbf v_g$ is the sedimentation velocity, $\mathbf v_R$ is the velocity induced by clinostat rotation, $\mathbf v_c$ is the centrifugal drift velocity, $\mathbf v_T$ is the translational Brownian velocity and $\mathbf v_s$ is the swimming velocity.

Swimming occurs at constant speed $v_s$,


$\mathbf v_s = v_s \hat{\mathbf e}$,


where $\hat{\mathbf e}$ is the swimming direction.

The swimming direction evolves through a combination of deterministic rotation and rotational diffusion,


$\frac{d\hat{\mathbf e}}{dt} = \boldsymbol{\omega}\times\hat{\mathbf e} + \sqrt{\frac{2D_R}{\Delta t}}(I-\hat{\mathbf e}\hat{\mathbf e}^{T})\boldsymbol{\xi}.$

The first term describes reorientation due to clinostat rotation and the second describes rotational Brownian motion. The projection operator $I-\hat{\mathbf e}\hat{\mathbf e}^{T}$ ensures that diffusion changes the orientation of the bacterium without changing the magnitude of the unit swimming vector.

In addition to rotational diffusion, the bacterium can undergo stochastic tumbling events. During a tumble a new random swimming direction is selected.

## Input files

The simulation requires a configuration file and an initial-condition file.

The configuration file specifies all simulation parameters, including the timestep, simulation duration, clinostat geometry, fluid properties and organism properties. The initial-condition file specifies the organism radius, initial position, swimming speed and density.

See `test_config.txt` for an example.

## Output

Simulation output is written as

```text
run_1.csv
```

with columns

| Column | Description    |
| ------ | -------------- |
| 1      | x position (m) |
| 2      | y position (m) |
| 3      | z position (m) |
| 4      | time (s)       |

Data are written at the output interval specified in the configuration file.

## Before running

A few simple checks can prevent most user errors.

First, ensure that all quantities are supplied in SI units. Lengths should be specified in metres, times in seconds, densities in kg m(^{-3}), viscosities in Pa s and diffusion coefficients in their standard SI units.

Second, ensure that the initial position lies inside the clinostat geometry,

$R_{\mathrm{inner}} < \sqrt{x^2+y^2} < R_{\mathrm{outer}}$,

with

$0 \le z \le H.$

If boundary conditions are enabled and the initial position lies outside the allowed region, the particle will immediately be projected onto a boundary.

Ensure the timestep is small enough to resolve any behaviour needing to be observed e.g tumbling, thermal or ballistic motion.

For simulations including rotation we generally followed the rule:

$\Delta t \lesssim \frac{1}{10}\frac{2\pi}{\omega}$

## Reproducing published results

The repository includes example config files (See /Configs) used to produce the results presented in [the accompanying publication (hopefully!)]. Because the model has stochastic processes, individual trajectories will differ between runs, but on averaged over multiple runs results should remain similar.

