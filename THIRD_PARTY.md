# Third-party code in CurvilinearFilters

The repository as a whole is GPL-2.0-or-later (see `README.md` and `LICENSE`). The items below are based on other people's code or algorithms and keep their own terms.

| Path | Author(s) | Licence | Source |
|---|---|---|---|
| `MFAT/` | H. Alhasson, M. Alharbi, B. Obara. Refactoring and extensions by M. D. Fricker | Original code: Apache License 2.0. Modifications: GPL-2.0-or-later | https://github.com/Haifafh/MFAT |
| `PhaseCongruency/src/phaseCong3Optimised.m`, `phaseSymOptimised.m`, `utils/lowpassfilter.m` | Peter Kovesi; optimised by M. D. Fricker | MIT-style (notice in each file) | https://www.peterkovesi.com/matlabfns/ |
| `Hessian/original/FrangiFilter2D.m`, `Hessian/src/applyHessian2D.m`, `Vesselness/src/vesselnessFilter2DOptimised.m` (and Kroon-style `eig2image.m`, `applyHessian3D*.m`) | Marc Schrijver; Dirk-Jan Kroon (University of Twente) | BSD (MATLAB File Exchange terms) | MATLAB File Exchange 24409 "Hessian based Frangi Vesselness filter" |

Apache License 2.0 is compatible with GPL version 3 but not version 2. Copies that include `MFAT/` should therefore be distributed under GPLv3 terms, which the "or later" clause allows. Apache 2.0 also requires a copy of its licence text to go with the code; that copy is still to be added (see open items).

These are wrappers only; nothing is bundled:
- ERnet-v2 (GPL-3.0). The model files and model code are downloaded at setup, and `ERnet/python/ernetServer.py` loads them at run time.
- nERdy. The package is installed into its own Python environment.

## Open items

- Add the Apache 2.0 licence text as `MFAT/LICENSE` (copy `LICENSE.txt` from the upstream repository).
- Add Kroon's BSD notice to the Kroon-derived files listed above.
- `Hessian/src/legacy/HessianVectorField.m` and `NeuritenessFilter2D.m` (B. Obara) have no licence. The same applies to `EigenMatrix2x2M.m`, `Gaussian2ndDerivatives.m` and `HessianMatrix2D.m` in the same folder, which are probably from Obara's utilities. Either get permission or rewrite them.
- `SOAGK/legacy/AGlineDetector.m` and `createGaussianFilter.m` (C. Lopez-Molina) have no licence. Either get permission, or drop them and keep only the re-implementations in `SOAGK/src`.
- `BowlerHat/src/BowlerHatFilter.m` is probably derived from Sazak, Nelson and Obara's published code, with no attribution or licence. Confirm its origin.
