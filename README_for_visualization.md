# Visualization Tools for Chemotaxi Simulation

These scripts act as post-processing tools designed to work with the output of the chemotaxi simulation package.

## Prerequisites

- Python 3.x
- Numpy
- Matplotlib (`mplot3d` toolkit included)
- Numba
- Scipy (for `.mat` colormap loading)
- The `read_input` and `read_vertex_file` modules from the original package.

Install required libraries:
```bash
pip install numpy matplotlib numba scipy
```

---

## 1. Chemical Concentration Field Visualization

The script: `visualize_concentration_allkind_acc_v6.py` calculates and visualizes the chemical concentration field at a specific moment in time (or as a continuous animation), based on the complete trajectory of the simulated particles.

### Features

- Reads simulation parameters directly from the `*.inputfile`, including support for individual parameters (like Peclet numbers) for multiple particles.
- Loads and parses multi-particle trajectory data (`*.config` file).
- Supports both **2D** and **3D** simulation domains:
  - **2D Mode**: Generates static heatmaps or full animations of the concentration evolution. Supports overlaying gradient vector fields (`-∇C`) with quiver arrows.
  - **3D Mode**: Efficiently computes and renders two key 2D concentration slices (YZ and XZ) alongside 3D trajectories, avoiding expensive full 3D volumetric calculations.
- **Particle Models Supported**:
  - **`non-janus`**: Point particle sources.
  - **`janus`**: Spherical particles with defined surface chemical and mobility distributions, properly handling translational and rotational dynamics (via quaternions).
  - Configurable 2D particle rendering: solid circle, chemical active regions, mobility regions, or both.
- Accelerated with `numba` parallel JIT compilation (`@njit(parallel=True, fastmath=True)`).
- Flexible CLI options: support for file prefixes, custom bounding ranges/centers, and various slice centering options in 3D.

### How to Use

#### For 2D Simulations (Static and Animated Outputs)

**Example 1: Generate a 2D Static Plot using File Prefix & Gradient Overlay**
```bash
python visualize_concentration_allkind_acc_v6.py \
    --file-prefix simulation_results/chemo_2d_test \
    --mode static \
    --time 100.0 \
    --plot-gradient \
    --visualize-particle all \
    --output-file static_2d_concentration
```
This will generate `static_2d_concentration.png`.

**Example 2: Generate a 2D Animation (GIF format)**
```bash
python visualize_concentration_allkind_acc_v6.py \
    --input-file simulation_results/chemo_2d_test.inputfile \
    --trajectory-file simulation_results/chemo_2d_test.config \
    --mode animation \
    --frame-interval 1.0 \
    --writer gif \
    --visualize-particle chem \
    --output-file animation_2d_chemo
```
This will generate `animation_2d_chemo.gif`.

#### For 3D Simulations (Cross-Section Slices)

1. **Ensure your `*.inputfile` is configured for Janus particles (if applicable):**
   ```ini
   domain = 3D
   particle_type = janus
   chemical_distribution_file = path/to/your/surface_sources.chem_dist.dat
   ```

2. **Ensure your `*.config` trajectory file contains 7 columns:** `x y z qw qx qy qz`.

3. **Run the script:**
   ```bash
   python visualize_concentration_allkind_acc_v6.py \
       --input-file simulation_results/janus_3d.inputfile \
       --trajectory-file simulation_results/janus_3d.config \
       --time 100.0 \
       --output-file concentration_3d_janus \
       --resolution 100 \
       --slice-center auto
   ```
   This generates a 3-panel plot containing YZ and XZ cross-sections alongside the 3D trajectory and slice planes.

### Command-Line Arguments

- `--file-prefix` *(optional)*: Path and prefix for input files. If specified, `--input-file` and `--trajectory-file` are automatically matched.
- `--input-file` *(optional)*: Path to the simulation input file (`*.inputfile`).
- `--trajectory-file` *(optional)*: Path to the trajectory file (`*.config`).
- `--output-file` *(optional)*: Base name of output file (without extension). Defaults to prefix or `concentration_field`.
- `--mode` *(optional)*: Visualization mode: `static` or `animation`. Default: `static`.
- `--time` *(required for static mode)*: Specific simulation time $t_1$ to visualize.
- `--slice-center` *(optional for 3D)*: Method to center slices: `auto`, `end`, or explicit coordinates `x,y,z`. Default: `auto`.
- `--frame-interval` *(optional for animation)*: Simulation time between frames. Default: `5.0`.
- `--resolution` *(optional)*: Grid resolution along each axis. Default: `100`.
- `--padding` *(optional)*: Padding distance around trajectories for grid boundaries. Default: `5.0`.
- `--writer` *(optional for animation)*: Animation writer: `gif` (Pillow) or `ffmpeg` (MP4). Default: `gif`.
- `--visualize-particle` *(optional)*: Visualize particles in 2D: `none`, `solid`, `chem`, `mobility`, or `all`. Default: `none`.
- `--hide-markers` *(optional)*: Hide start/end trajectory markers in static mode.
- `--center` *(optional)*: Center the 2D grid at `(X0, Y0)` keeping computed width/height.
- `--range` *(optional)*: Explicitly set 2D domain range `XMIN XMAX YMIN YMAX` (overrides `--center`).
- `--plot-gradient` *(optional)*: Overlay the concentration gradient field (`-∇C`) via quiver arrows.
- `--quiver-density` *(optional)*: Sampling step size for quiver arrows. Default: `4`.

---

## 2. Trajectory Visualization ()

The script:`visualize_trajectory_v2.py` visualizes 2D and 3D multi-particle trajectory curves from the block-formatted simulation output file.

### Features

- Handles multiple trajectories simultaneously for both **2D** and **3D** systems.
- Trajectory coloration:
  - **Time Gradient Mode**: Displays trajectory lines with smooth temporal color transitions and per-particle colorbars.
  - **Solid Color Mode**: Plots distinct solid colors for each trajectory without colorbars.
- Supports custom colormaps loaded from a MATLAB `.mat` file (`colormap.mat`), falling back to Matplotlib defaults (`viridis`, `plasma`, etc.).
- Filter trajectories by time range (`--time_period t1 t2`).
- Custom bounding box and aspect ratio controls via `--center` and `--range`.
- Marks starting positions with **green** circles and current/ending positions with **red** circles.

### How to Use

**Example 1: Default Gradient Plot**
```bash
python visualize_trajectory_v2.py simulation_results/test_2d.config
```

**Example 2: Plot a Specific Time Range**
```bash
python visualize_trajectory_v2.py simulation_results/test_2d.config --time_period 10.0 50.0
```

**Example 3: Select Specific Colormaps from a `.mat` File**
```bash
python visualize_trajectory_v2.py simulation_results/test_2d.config --colormap_file colormap.mat --cmap_indices 1 3 5
```

**Example 4: Specify Manual Center and Bounding Box**
```bash
python visualize_trajectory_v2.py simulation_results/test_2d.config --center 0 0 5 --range 20 20 10
```

**Example 5: Solid Color Trajectory Plot**
```bash
python visualize_trajectory_v2.py simulation_results/test_2d.config --solid_color
```

### Command-Line Arguments

- `trajectory_file` *(positional, required)*: Path to the trajectory file (`*.config`).
- `--dt` *(optional)*: Simulation time step between consecutive data frames. Default: `1.0`.
- `--colormap_file` *(optional)*: Path to `.mat` file containing custom colormaps. Default: `colormap.mat`.
- `--cmap_indices` *(optional)*: Space-separated list of colormap indices to extract from the `.mat` file (e.g. `2 4 8`).
- `--time_period` *(optional)*: Plot trajectory within time window `t1 t2`.
- `--center` *(optional)*: Center coordinates of the plot box (`X0 Y0 Z0`). Must be paired with `--range`.
- `--range` *(optional)*: Dimensions / total lengths of the plot box (`Lx Ly Lz`). Must be paired with `--center`.
- `--solid_color` *(optional)*: Render trajectory lines in solid colors instead of temporal gradients.