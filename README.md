# Documentation
This package contains several python codes to run simulations of
auto-chemotatic droplet, in two-dimensional or three-dimensional space. These codes can compute the
chemical concentration distribution, chemical gradient, obtain the droplet velocities and position for deterministic or stochastic
dynamic simulations, using a mesh-free method.


We explain in the next sections how to use the package.

**Note:** Please refer to the note for the governing equations and detailed algorithm.

## 1. Prepare the package
The codes are implemented in python (version 3.x) and it is not necessary to compile the package to use it.

## 2. Droplet body configuration
The coordinates of the discretized points forming a 2D/3D droplet surface in the default configuration
(location (0, 0) and quaternion (1, 0, 0, 0)) are given to the codes
through `*.vertex` files. The format of these files is:

```
number_of_discretized_points_on_droplet_surface
vector_location_point_0
vector_location_point_1
.
.
.
```

For example, the file `structures/circle_R_1_N10.vertex` gives the
structure of a 2D circular particle discretized by 10 points.

We use a vector (2 numbers for 2D; 3 numbers for 3D) and a quaternion (4 numbers for 3D) to represent the
location (2D and 3D) and orientation of each body (3D only, and knowledge for quaternion please see Ref. [1](http://dx.doi.org/10.1063/1.4932062) for details).
This information is saved by the code in the `*.clones` files,
with format (2D):

```
vector_location_body_0
vector_location_body_1
.
.
.
```
3D format:
```
vector_location_body_0 quaternion_body_0
vector_location_body_1 quaternion_body_1
.
.
.
```

## 3. Run dynamic simulations
Here, we explain how to use the main
code which allows to run deterministic and stochastic simulations for droplet in 2D case.

First, create a directory to store your simulation data, like `simulation_results`;
Then inspect the input file
`test.txt`:

---

```
# Job description
job_name                    test_2d
job_type                    dynamic
domain                      2D
scheme                      history_local_compose_2d
acceleration                parallel
core                        8
numerical_method            forward_euler

# Parameters specification
mobility_alpha              12
radius                      1
intrinsic_linear_velocity   1
intrinsic_angular_velocity  1
emission_rate               1
persistence_length          1

peclet_number               120

initial_position            0 0
initial_orientation         1 0

# Numerical simulation configuration
droplet_num                 1
n_steps                     80000
dt                          0.00390625

# Output configuration
output_name                 simulation_results/run
save_clones				    one_file
initial_step                0

# Discretization info
structure                   structures/circle_R_1_N60.vertex
```

---

With this input we can run a simulation with one 2D droplet;
see structures given to the options `structure`. To run the simulation use

`
python main.py --input-file test.txt
`

Now, you can inspect the outputs, `ls simulation_results/run.*`. The output files are:

* `.config`: For each time step saved the
code saves a file with the location of the droplet. The name format is (output_name + time_step + .config)
The format of the files is the same that in the input .config files.

* `.chemforce.dat`: For each time step saved the
code saves a file with the chemical force applied to the droplet. The name format is (output_name + time_step + .chemforce.dat)
The format of the files is as following.

* `.inputfile`: a copy of the input file.

* `.time`: the wall-clock time elapsed during the simulation (in seconds).
* `.time.log`: the wall-clock time elapsed per step (in seconds).

**List of options for the input file:**
* `domain` (string). Options: `2D` and `3D`. 3D codes are not updated, try 2D first.
* `scheme` (string). Option: `history_local_compose_2d`. (No effects on codes now)
* `acceleration` (string). Options: `numba` and `parallel`. Numba acceleration is recommended for total step <= 10000;
parallel acceleration is recommended for EXTRA-long simulation and fine grid of the structure (or even 3D cases).
* `core` (int (default 1)). Number of cores used for parallel processing. Only effective for the case `acceleration` used `parallel`.
* `numerical_method` (string). Options: `forward_euler` and `adams_bashforth_2`

| Name | Solver type | Notes |
| ---- | ----------- | ----- |
| forward_euler               | Iterative    | first order accuracy            |
| adams_bashforth_2             | Iterative    | second order accuracy           |
| stochastic_first_order_RFD                | Iterative    | NOT implemented yet|


* `mobility_alpha` (float (default 1)): In the JCP paper and my note, we use notation $\Lambda$, which is a mobility parameter to determine
the magnitude of the chemical force.
* `radius` (float (default 1)): The geometric radius of the droplet. Normally we use non-dimensionlized parameter R=1.
* `intrinsic_linear_velocity` (float (default 1)): The intrinsic linear velocity of the droplet. 
Normally we use non-dimensionlized parameter v_0=1.
* `intrinsic_angular_velocity` (float (default 1)): The intrinsic angular velocity of the droplet.
Normally we use non-dimensionlized parameter \omega_0=1.
* `emission_rate` (float (default 1)): The emission rate of the chemical substance. 
Normally we use non-dimensionlized parameter Q_0=1.
* `persistence_length` (float (default 1))).
* `peclet_number`(float (default 1)): $Pe = Rv_o/D$, 
is the ratio of self-propelling rate of the droplet to diffusion rate of emitted solute.
* `initial_position`(float (vector default 0 0)): Vector format, 2D in format (x_0, y_0), 3D in format (x_0, y_0, z_0)
* `initial_orientation`(float (vector default 0 0 )): Vector format, 2D in format (R\cos\theta, R\sin\theta), 3D use quaternion format (NOT implemented yet).
* `droplet_num`(int (default 1)): Currently this code only support single particle cases.
* `n_steps`(int (default 1)): Number of simulation steps.
* `dt`(float): time step length to advance the simulation.
* `save_clones`(string (default `one_file`)) :options
`_one_file_per_step_` and `one_file`. With the option
`_one_file_per_step_` the clones configuration are saved in one file per time step. With the option
`one_file` the code saves one file with the
configurations of all the time steps.
* `initial_step`(int (default 0)): Use this option to restart a simulation.
If `initial_step > 0` the code will run from time step `initial_step` to
`n_steps`. Also, the code will try to load `.config` files with the name
(output_name + structure_name + initial_step + .config). (This restart function has NOT implemented yet).
* `structure`(string): The file path under main directory and file name of the discretized surface points `.vertex` file.

## 4. Software organization
* **body/**: it contains a class to handle a single droplet body. `body_2D.py` for 2D cases and `body_3D.py` for 
3D cases.
* **integrator/**: it has a small class to handle quaternions and
the schemes to integrate the equations of motion.
* **read_input/**: it has a small class to read and handle input information and vertex information.
* **structures/**: it stores `.vertex` files.
* **tools/**: start-up useful tools (NOT necessary).
* **main.py**: it calls, processes and advances for simulations.