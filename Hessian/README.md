# Hessian2DfILTERS — Multiscale Hessian Filtering (2D)

This repository provides a **clean, well-tested implementation of multiscale Hessian-based filters** for 2D images.

The focus is on **geometric feature detection**, not classification or segmentation.

---

## Features

- Multiscale Hessian analysis (Gaussian scale space)
- Robust eigenvalue handling
- Max-over-scales aggregation (stateless)
- Clean separation of filter semantics
- Unit-tested invariants
- MATLAB-native, no toolboxes required beyond Image Processing

---

## Supported Filters

| FilterType | Detects | Typical Applications |
|-----------|--------|---------------------|
| vesselness | Tubular structures | Blood vessels, pipes |
| ridge | Line-like structures | Ridges, fibers |
| neuriteness | Curvilinear continuity | Neurites, filaments |
| blob | Isotropic blobs | Vesicles, spots |
| plate | Thick elongated regions | Membranes, bands |

---

## Algorithm Overview

For each scale σ:

1. Smooth image with Gaussian(σ)
2. Compute Hessian matrix:
   \[
   H = \begin{bmatrix} I_{xx} & I_{xy} \\ I_{xy} & I_{yy} \end{bmatrix}
   \]
3. Compute eigenvalues (λ₁, λ₂) and eigenvectors
4. Evaluate a scalar response function
5. Keep the maximum response across scales

This design ensures:
- scale invariance via max selection
- no inter-scale feedback
- predictable behavior

---

## Usage

### Basic Vesselness

```matlab
I = imread('image.png');
I = im2single(I);

[R, scale, dir] = hessian2DFilters(I, ...
    'FilterType','vesselness');
imshow(R,[]);
