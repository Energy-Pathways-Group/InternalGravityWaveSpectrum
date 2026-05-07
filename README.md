# InternalGravityWaveSpectrum

`InternalGravityWaveSpectrum` is a MATLAB class for constructing internal gravity wave spectra from ocean stratification profiles using vertical modes and a logarithmically spaced horizontal wavenumber grid.

![Garrett-Munk power spectrum example](img/PowerSpectrumGM_allStations.png)

## Overview

This repository provides tools for building and diagnosing spectral models of internal gravity waves in a stratified ocean. Given a buoyancy frequency profile `N2(z)` and an ocean depth, the main class computes vertical modes, modal structure functions, dispersion relations, modal amplitudes, and energy components.

The resulting energy fields can be analyzed as horizontal kinetic energy (`HKE`), vertical kinetic energy (`VKE`), potential energy (`PE`), or total energy (`TE`). The code is intended for research workflows where the user wants to compare theoretical internal-wave spectra with model output, observations, or idealized stratification profiles.

## Features

- Computes vertical modes with `InternalModesSpectral` and `InternalModesWKBSpectral`.
- Builds a log-spaced radial horizontal wavenumber grid.
- Evaluates modal structure functions, dispersion relations, and eigendepths.
- Supports energy components `HKE`, `VKE`, `PE`, and `TE`.
- Includes Garrett-Munk, generalized, tidal, near-inertial, and frequency-localized spectra.
- Provides interpolation and diagnostic tools for vertical profiles, horizontal wavenumber spectra, frequency spectra, and mode spectra.

## Requirements

- MATLAB, with support for class definitions and argument validation blocks.
- The `src/` directory from this repository on the MATLAB path.
- External classes used by the main model:
  - `InternalModesSpectral`
  - `InternalModesWKBSpectral`

Some examples and tests may also require:

- `matlab.unittest`
- `WVTransformBoussinesq`
- `VerticalModeAtlas`
- Local `.mat` or `.nc` data files used by specific research scripts.

## Installation

Clone the repository and add the source directory to your MATLAB path:

```matlab
addpath("src")
```

Also add any external dependency directories that provide `InternalModesSpectral` and `InternalModesWKBSpectral`.

## Quick Start

```matlab
Lz = 4000;
N0 = 3 * 2*pi/3600;
Lgm = 1300;
N2 = @(z) N0^2 * exp(2*z/Lgm);

igws = InternalGravityWaveSpectrum(N2, Lz, ...
    latitude=33, ...
    nModes=64, ...
    nK=64);

S = igws.gmSpectrum();
igws = igws.assignEnergySpectrum(S);

z = linspace(-Lz, 0, 500);
TE = igws.verticalVariance('TE', zVector=z);

figure
plot(TE, z)
xlabel('Total energy variance')
ylabel('Depth [m]')
grid on
```

## Main Usage Pattern

1. Define a stratification profile as a function handle, `N2(z)`.
2. Initialize the model with depth, latitude, and grid resolution.
3. Choose a built-in spectrum or define a custom function handle `S(k,j)`.
4. Assign the spectrum with `assignEnergySpectrum`.
5. Evaluate diagnostics such as vertical variance, horizontal wavenumber spectra, frequency spectra, or energy by vertical mode.

## Available Spectra

The class includes these spectrum constructors:

- `gmSpectrum`: Garrett-Munk type spectrum normalized to the GM energy level.
- `generalSpectrum`: tunable spectrum with mode and wavenumber slopes.
- `frequencyLocalizedSpectrum`: Lorentzian spectrum localized around a target frequency and vertical mode.
- `tidalSpectrum`: M2 tidal wrapper around `frequencyLocalizedSpectrum`.
- `fSpectrum`: near-inertial wrapper centered at the Coriolis frequency.
- `fM2Spectrum`: wrapper centered at `f0 + omega_M2`.
- `randomRealization`: stochastic realization based on the assigned variance field.

You can inspect the available named spectrum constructors from MATLAB:

```matlab
igws.listAvailableSpectra();
```

## Coordinate and Energy Conventions

- Depth is negative below the surface: `z < 0`, with `z = 0` at the surface.
- Latitude should be away from the equator; the model rejects latitudes within about 5 degrees of the equator.
- Very high latitudes are also rejected by the constructor.
- Frequencies satisfy approximately `f0 < omega < sqrt(N2max)`.
- Horizontal wavenumber `K` is radial and log-spaced.
- Energy terms are requested with `'TE'`, `'HKE'`, `'VKE'`, or `'PE'`.

## Examples

Live-script examples are available in:

- `src/Example/example.mlx`
- `src/Example/QuadraturePoints.mlx`

Additional rendered figures are available in `img/`, including quadrature-point diagnostics, radius-of-deformation diagnostics, and frequency/power spectrum examples.

## Testing

Tests are located in `src/unitTest/`. A basic test suite can be run with:

```matlab
addpath("src")
result = run(matlab.unittest.TestSuite.fromClass(?TestIGWSInitialization));
table(result)
```

Some tests rely on external packages, local atlas files, or model-output data paths. If a test fails because a data file or external class is missing, first check whether the relevant dependency is available on your MATLAB path.

## Known Limitations

- The model is not valid near the equator.
- Some tests and examples currently depend on local research data paths.
- Some diagnostics currently assume scalar vertical-position queries.
- Large `.mat` and `.nc` files can make the repository heavy if committed directly to Git.

## References

Please cite the relevant internal-wave, vertical-mode, and Garrett-Munk spectrum literature used in your analysis. The class documentation also refers to Jeffrey et al. 2021 for the squared linear wave-equation energy coefficients.

## Author

Developed by Leticia Fabre de Lima for research on internal gravity wave spectra.

