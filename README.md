# 🧪 Biomedical Flow Simulation Project: Lung Airflow Modeling

## 📌 Project Overview

This project was developed alongside our theoretical coursework, focusing on the **numerical simulation of airflow in the lungs** using **finite difference methods**. Specifically, we solved the **2D Stokes and Navier-Stokes equations** in an **asymmetric domain** representing lung airways, with **imposed pressure conditions** at the left and right boundaries to simulate inhalation and exhalation pressures.

Simulating airflow in the lungs helps better understand respiratory mechanics, gas exchange efficiency, and effects of airway deformations in diseases such as asthma or COPD.

---

## 🧠 Objectives

- Implement a **finite difference solver** for the Stokes equations.
- Extend the simulation to the **Navier-Stokes equations** in 2D.
- Model an **asymmetric airway domain** reflecting realistic lung geometry.
- Apply **Dirichlet boundary conditions** with prescribed pressures simulating breathing cycles.
- **Visualize** velocity vectors, pressure fields, streamlines, and compare results with known analytical solutions when available.

---

## 🧰 Code Libraries and Tools Used

- **NumPy**: for efficient numerical computations and linear algebra.
- **Matplotlib**: for plotting velocity fields, pressure distributions, and streamlines.
- **SciPy** (optional): for advanced numerical solvers.
- **Python 3.x**: chosen for its readability and scientific ecosystem.
- **Jupyter Notebook** (optional): for interactive development and visualization.

---

## 🧮 Mathematical Model

We solve the steady-state incompressible flow equations modeling airflow:

### 🔹 Stokes Equations (Low Reynolds number):
\[
-\nabla p + \mu \nabla^2 \vec{u} = 0 \\
\nabla \cdot \vec{u} = 0
\]

### 🔹 Navier-Stokes Equations (General case):
\[
\rho (\vec{u} \cdot \nabla)\vec{u} = -\nabla p + \mu \nabla^2 \vec{u} \\
\nabla \cdot \vec{u} = 0
\]

Where:
- \( \vec{u} = (u_x, u_y) \): airflow velocity vector in lung airways
- \( p \): pressure
- \( \mu \): dynamic viscosity of air
- \( \rho \): air density

---

## 🧰 Numerical Method

- **Discretization**: Finite Difference Method on a 2D grid
- **Domain**: Asymmetric geometry representing lung airways
- **Boundary Conditions**:
  - Prescribed pressures at inlet and outlet boundaries to simulate breathing
  - No-slip condition on airway walls
- **Solver**: Linear least squares solver using NumPy
- **Post-processing**: Visualization of velocity, pressure, and flow streamlines; error analysis against analytical profiles

---

## 🧑‍💻 Project Structure

```plaintext
📁 project_root/
│
├── stokes_solver_no_deformation.py
├── navier_stokes_solver_with_deformation.py
├── streamline_plotter.py
├── poiseuille_law_reference.py
├── velocity_vector_visualization.py
└── README.md
