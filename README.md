# Solid Tank Ray‑Tracing Simulation

This repository provides MATLAB functions for simulating and evaluating a solid‑tank fan‑beam optical CT scanner design. You can compute key performance metrics (effective radius ratio, image uniformity, and ray‑order preservation) and generate visualization plots for any given lens geometry and matching fluid.

---

## Features

* **`SolidTankObj`**: Compute objective metrics for a candidate design.
* **`SolidTankPlot`**: Generate ray‑tracing visualizations and save plots.

---

## Prerequisites

* MATLAB R2024a (or later)

---



## Usage

### 1. Evaluate objectives with `SolidTankObj`

```matlab
% Design variables x (1×10 vector)
x = [200, 0, 60, 80, 2.5, 300, 1, 1, 1, 1];

% Fidelity: number of rays for quick tests
testRays = 1000;  

% Matching fluid refractive index (e.g., FlexiDos3D)
gelIndex = 1.4225;  

% Rear‑lens type: 'conic', 'polynomial', or 'aspheric'
lensType = 'aspheric';

[effRad, intensityCV, tau] = SolidTankObj(x, testRays, gelIndex, lensType);

disp(['Effective Radius: ', num2str(effRad)]);
disp(['Intensity CV: ', num2str(intensityCV)]);
disp(['Kendall \tau: ', num2str(tau)]);
```

### 2. Generate visualization with `SolidTankPlot`

```matlab
% Supply same design parameters as above
x = [200, 0, 60, 80, 2.5, 300, 1, 1, 1, 1];
gelIndex = 1.4225;
lensType = 'aspheric';

% Output file path (e.g. 'plots/design1.png')
savePath = 'plots/design1.png';

SolidTankPlot(x, gelIndex, lensType, savePath);
```

---

## Function Reference

### `SolidTankObj(x, rayCount, gelIndex, rearLensType)`

* **Inputs**:

  * `x` (1×8 numeric vector for 'conic', 1×10 numeric vector for 'polynomial' and 'aspheric'): lens geometry parameters
  * `rayCount` (scalar): number of rays for fidelity control
  * `gelIndex` (scalar): refractive index of matching fluid
  * `rearLensType` (string): 'conic' | 'polynomial' | 'aspheric'
* **Outputs**:

  * `effRad` (double): effective radius ratio
  * `intensityCV` (double): coefficient of variation of detector intensities
  * `tau` (double): Absolute Value of Kendall’s tau for ray‐order preservation

### `SolidTankPlot(x, gelIndex, rearLensType, savePath)`

* **Inputs**:

  * `x` (1×8 numeric vector for 'conic', 1×10 numeric vector for 'polynomial' and 'aspheric'): lens geometry parameters
  * `gelIndex` (scalar): matching fluid refractive index
  * `rearLensType` (string): lens type
  * `savePath` (string): destination filename for the plot

---

## Contact

Zhongda Huang – \[[dominic.huang@ubc.ca](mailto:dominic.huang@ubc.ca)]

Project Link: [https://github.com/ViggleH/MO-STD](https://github.com/ViggleH/MO-STD)
