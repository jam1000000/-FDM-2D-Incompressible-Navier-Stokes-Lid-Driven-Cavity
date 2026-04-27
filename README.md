# Lid-Driven Cavity Simulation — Finite Difference Method

A 2D incompressible Navier–Stokes solver for the benchmark problem implemented from scratch in Python using finite difference discretisation and a fractional-step (projection) method.

---

## Overview

The lid-driven cavity is one of the most widely used benchmark problems in computational fluid dynamics (CFD). A square domain is filled with viscous fluid, and the top wall moves at a constant horizontal velocity while the remaining three walls are stationary (no-slip condition). Despite its simple geometry, the flow develops vortices that are sensitive to the Reynolds number.

This notebook solves the 2D incompressible Navier–Stokes equations:

$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu \nabla^2 \mathbf{u}$$
$$\nabla \cdot \mathbf{u} = 0$$

using an explicit time-stepping scheme on a uniform Cartesian grid.

---

## Method

### Spatial Discretisation
All spatial derivatives are approximated with second-order central differences on a uniform N×N grid spanning the unit square [0,1]².

### Time Integration — Fractional-Step (Projection) Method
Each time step proceeds in three stages:

1. Tentative velocity step — advance momentum with convection and diffusion terms, ignoring pressure:

$$\mathbf{u}^* = \mathbf{u}^n + \Delta t \left[ -(\mathbf{u}^n \cdot \nabla)\mathbf{u}^n + \nu \nabla^2 \mathbf{u}^n \right]$$

2. Pressure Poisson solve — enforce the incompressibility constraint by solving:

$$\nabla^2 p^{n+1} = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*$$

using Gauss–Seidel-style iterative relaxation.

3. Velocity correction — project the tentative velocity onto a divergence-free field:

$$\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla p^{n+1}$$

### Boundary Conditions
| Boundary | u (horizontal) | v (vertical) | p (pressure) |
|----------|---------------|--------------|--------------|
| Top wall (lid) | `U_top = 1.0` | 0 (no-slip) | Dirichlet: p = 0 |
| Bottom wall | 0 (no-slip) | 0 (no-slip) | Neumann: ∂p/∂n = 0 |
| Left wall | 0 (no-slip) | 0 (no-slip) | Neumann: ∂p/∂n = 0 |
| Right wall | 0 (no-slip) | 0 (no-slip) | Neumann: ∂p/∂n = 0 |

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N_POINTS` | 128 | Grid resolution (N×N) |
| `DOMAIN_SIZE` | 1.0 | Side length of the cavity |
| `N_ITERATIONS` | 40,000 | Number of time steps |
| `TIME_STEP_LENGTH` | 0.0001 | Δt |
| `KINEMATIC_VISCOSITY` | 0.0003125 | ν |
| `DENSITY` | 1.0 | ρ |
| `HORIZONTAL_VELOCITY_TOP` | 1.0 | Lid velocity |
| `N_PRESSURE_POISSON_ITERATIONS` | 128 | Inner Poisson solver iterations per step |
| Reynolds number | 3200 | Re = U·L / ν (derived) |

---

## Outputs

The simulation produces three plots at the end of the run:

1. Vertical centerline u-velocity profile — u(y) at x = 0.5. Useful for comparison against the Ghia et al. (1982) benchmark data.
2. Horizontal centerline v-velocity profile — v(x) at y = 0.5. Also benchmarked against Ghia et al.
3. Pressure contour + velocity vector field — filled contour plot of pressure overlaid with a white quiver plot of the velocity field, giving a visual picture of the recirculation structure.

---

## Requirements

```
python >= 3.8
numpy
matplotlib
tqdm
```

Install dependencies:

```bash
pip install numpy matplotlib tqdm
```

---

## Usage

### In Jupyter

Open the notebook and run all cells:

```bash
jupyter notebook Lid-Driven_Cavity_Sim_1_0_2__FDM_.ipynb
```

### As a Script

The notebook is structured with a `main()` guard and can be run directly:

```bash
jupyter nbconvert --to script Lid-Driven_Cavity_Sim_1_0_2__FDM_.ipynb
python Lid-Driven_Cavity_Sim_1_0_2__FDM_.py
```

> **Note:** At Re = 3200 and a 128×128 grid, a full 40,000-step run can take several minutes on a standard CPU. Reduce `N_ITERATIONS` or `N_POINTS` for quick exploratory runs.

---

## Stability Considerations

The explicit scheme is conditionally stable. The simulation includes overflow guards that halt the loop early if any field variable becomes `NaN` or `Inf`. If you change parameters, keep the CFL condition in mind:

$$\Delta t \lesssim \min\left(\frac{\Delta x}{U},\ \frac{\Delta x^2}{2\nu}\right)$$

The defaults satisfy this constraint. Increasing `N_POINTS` without reducing `TIME_STEP_LENGTH` proportionally will likely cause blow-up.

---

## Validation

Centerline velocity profiles can be compared against the widely cited benchmark:

> Ghia, U., Ghia, K. N., & Shin, C. T. (1982). *High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method.* Journal of Computational Physics, 48(3), 387–411.

---

## References

- Ghia et al. (1982) — canonical benchmark data for the lid-driven cavity at Re = 100–10,000
- Chorin, A. J. (1968) — original projection method for incompressible flow
- Ferziger, J. H. & Perić, M. — *Computational Methods for Fluid Dynamics*
