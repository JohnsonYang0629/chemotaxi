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
Chemical substance distributions on the spherical droplet surface are given by the `*.chem_dist.dat` files. 
Each line in these files represents the chemical substance emitting rate for the corresponding line in the `*.vertex` files
for node positions.
```
1
1
0
0
.
.
.
```

## 3. Run dynamic simulations
Here, we explain how to use the main
code which allows to run deterministic and stochastic simulations for droplet in 2D case.

First, create a directory to store your simulation data, like `simulation_results`;
Then inspect the input file
`test_2d.txt`:

---

```
# Job description
job_name                    test_2d
job_type                    dynamic
domain                      2D
scheme                      history_local_compose_2d
acceleration                numba
core                        8
numerical_method            stochastic_first_order

# Parameters specification
mobility_alpha              4
radius                      1
intrinsic_linear_velocity   1
intrinsic_angular_velocity  1
emission_rate               1
translational_noise_gamma   500
rotational_noise_gamma      500

peclet_number               40

initial_position_2D                   0 0
initial_orientation_2D_vector         1 0

# Numerical simulation configuration
droplet_num                 1
n_steps                     2000
dt                          0.1

# Output configuration
output_name                 simulation_results/chemo_pe_40_lambda_4_noise_500
save_clones	            one_file
initial_step                0

# Discretization info
structure                   structures/circle_R_1_N60.vertex```

---
```
With this input we can run a simulation with one 2D droplet;
see structures given to the options `structure`. To run the simulation use

`
python main.py --input-file test_2d.txt
`

Now, you can inspect the outputs, `ls simulation_results/outputname.*`. The output files are:

* `.config`: For each time step saved the code saves a file with the location 
(and quaternion for 3D cases) of the droplet.

* `.velocity.dat`: For each time step saved the code saves a file with the velocity of the droplet. 

* `.chemforce.dat`: For each time step saved the code saves a file with the chemical force 
* (chemcial torque for 3D cases) applied to the droplet. 

* `.inputfile`: a copy of the input file.

* `.time`: the wall-clock time elapsed during the simulation (in seconds).
* `.time.log`: the wall-clock time elapsed per step (in seconds).
* `.random_state`: the file with the state of the random generator from current simulation.

**List of options for the input file:**
* `domain` (string). Options: `2D` and `3D`. 
* `scheme` (string). Option: `history_local_compose_2d` and `history_local_compose_3d`.
* `particle_type` (string). Option for 3D: `non_janus` and `janus`.
* `acceleration` (string). Options: `numba` and `parallel`. Numba acceleration is recommended for total step <= 10000;
Parallel acceleration is recommended for EXTRA-long simulation and fine grid of the structure (or even 3D cases).
* `core` (int (default 1)). Number of cores used for parallel processing. Only effective for the case `acceleration` used `parallel`.
* `numerical_method` (string). Options: `forward_euler`， `adams_bashforth_2` and `stochastic_first_order`.

| Name | Solver type | Notes |
| ---- | ----------- | ----- |
| forward_euler               | Iterative    | first order accuracy            |
| adams_bashforth_2             | Iterative    | second order accuracy           |
| stochastic_first_order               | Iterative    |  deterministically first order accuracy |


* `mobility_alpha` (float (default 1)): In the JCP paper and my note, we use notation $\Lambda$, which is a mobility parameter to determine
the magnitude of the chemical force.
* `radius` (float (default 1)): The geometric radius of the droplet. Normally we use non-dimensionlized parameter $R=1$.
* `intrinsic_linear_velocity` (float (default 1)): The intrinsic linear velocity of the droplet. 
Normally we use non-dimensionlized parameter $v_0=1$.
* `intrinsic_angular_velocity` (float (default 1)): The intrinsic angular velocity of the droplet.
Normally we use non-dimensionlized parameter $\omega_0=1$.
* `emission_rate` (float (default 1)): The emission rate of the chemical substance. 
Normally we use non-dimensionlized parameter $Q_0=1$.
* `rotational_noise_gamma` (float (default 500)): 2D cases: dθ/dt = Ω₀ + √(2/Γ) * ξ(t); 3D cases: Δθ = Ω₀ * τ̂ * Δt + √(2/Γ) * ξ(t) * √(Δt).
* `translational_noise_gamma` (float (default 500)): drₚ/dt = p̂ + F꜀ + √(2/Γₜ) * η(t).
* `peclet_number`(float (default 1)): $Pe = Rv_o/D$, 
is the ratio of self-propelling rate of the droplet to diffusion rate of emitted solute.
* `initial_position_2D`(float (vector default 0 0)) or `initial_position_3D`(float (vector default 0 0 0)): Vector format, 2D in format $(x_0, y_0)$, 3D in format $(x_0, y_0, z_0)$
* `initial_orientation_2D_vector`(float (vector default 0 0 )) or `initial_orientation_3D_quaternion`(float (vector default 1 0 0 0 )): Vector format, 2D in format $(R\cos\theta, R\sin\theta)$; 3D use Quaternion format.
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
* `chemical_distribution`(string): The file path under main directory and file name of the chemical substance distribution `.chem_dist.dat` file.
* `random_state` (string): name of a file with the state of the random generator from a previous simulation. It can be used to generate the same random numbers in different simulations.

## 4. Software organization
* **body/**: it contains a class to handle a single droplet body. `body_2D.py` for 2D cases and `body_3D.py` for 
3D cases.
* **integrator/**: it has a small class to handle quaternions and
the schemes to integrate the equations of motion.
* **read_input/**: it has a small class to read and handle input information and vertex information.
* **structures/**: it stores `.vertex` files.
* **tools/**: start-up useful tools (NOT necessary).
* **main.py**: it calls, processes and advances for simulations.
* **chem_functions.py**: it calculates related chemical gradient forces (to be called).

## 5. Notes
* Control MKL's Threading: To run the simulation by forcing MKL to run in single-threaded mode, 
letting Numba handle all the high-level parallelization with `prange`:
`
MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_NUM_THREADS=1 python main.py --input-file test_3d.txt
`
