# MFAT – Multiscale Fractional Anisotropy Tensor

## Overview
MFAT is a deterministic tensor-based framework for curvilinear structure
enhancement, with optional probabilistic and uncertainty-aware extensions.

## Core Components
- MFAT-λ: deterministic backbone
- MFAT-Prob: probabilistic confidence
- Entropy weighting: conservative modifier
- Fractional shaping: weak-signal enhancement

## Usage
```matlab
cfg = mfatConfig();
R = mfatLambda(I, sigmas);
P = mfatProb(I, sigmas);

