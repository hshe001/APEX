# APEX Reconstruction Toolbox

MATLAB code for simulation, demultiplexing, and reconstruction in **Axially Parallel Excitation (APEX) microscopy**.

This repository provides an implementation of the APEX forward model and reconstruction framework described in our work. The code can be used both for numerical simulations and for the reconstruction of experimentally acquired APEX data.

---

## Overview

APEX employs multiplexed axial excitation to simultaneously encode signals from multiple imaging planes into a reduced number of measurements. The original volumetric information can then be recovered through computational demultiplexing and reconstruction.

This repository demonstrates:

- Generation of a synthetic 3D ground-truth object
- Simulation of APEX multiplexed measurements
- PMT shot-noise and dark-current noise modeling
- Direct matrix inversion reconstruction
- Non-negative least-squares (NNLS) reconstruction
- Maximum-likelihood estimation (MLE) reconstruction
- Quantitative image-quality evaluation using PSNR and SSIM

---

## Repository Structure

```text
.
├── main.m
├── createGroundTruth.m
├── addPMTShotNoise3D.m
└── README.md
```

### main.m

Main script for:

- generating the ground-truth volume;
- simulating multiplexed APEX measurements;
- adding realistic PMT noise;
- reconstructing volumetric data using different algorithms;
- calculating quantitative reconstruction metrics.

### createGroundTruth.m

Generates a synthetic 3D object used for numerical simulations.

### addPMTShotNoise3D.m

Simulates PMT photon-counting statistics by incorporating:

- Poisson-distributed shot noise;
- PMT dark-current noise.

---

## Reconstruction of Experimental APEX Data

Although the example included in this repository is based on numerical simulations, the reconstruction framework is identical to that used for experimental APEX measurements.

To reconstruct experimental data:

1. Replace the simulated multiplexed measurements with experimentally acquired APEX images.
2. Specify the corresponding sampling matrix used during acquisition.
3. Run the reconstruction section in `main.m`.

The implemented NNLS and MLE reconstruction methods can be directly applied to experimental datasets without modification.

---

## Requirements

- MATLAB R2022a or later
- Optimization Toolbox
- Image Processing Toolbox
- Parallel Computing Toolbox

---

## Example Workflow

```matlab
% Generate ground truth
ground_truth = createGroundTruth(...);

% Simulate PMT noise
noise_3D = addPMTShotNoise3D(...);

% Run reconstruction
main
```

The script automatically generates:

- Ground-truth volume
- Noisy measurements
- Reconstructed volumes
- PSNR values
- SSIM values

---

## Output

The reconstruction results include:

| Method | Description |
|----------|----------|
| Direct Inversion | Matrix inversion reconstruction |
| NNLS | Non-negative least-squares reconstruction |
| MLE | Maximum-likelihood reconstruction under Poisson statistics |

---
