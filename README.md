# InternalGravityWaveSpectrum

`InternalGravityWaveSpectrum` is a MATLAB toolbox for computing non-hydrostatic ocean internal gravity wave vertical modes and assigning wave energy from a prescribed spectrum `S(K,j)`.

## Overview

This repository provides a numerical implementation for constructing internal-wave spectral model. Given buoyancy frequency squared `N2(z)`, total depth, latitude, and a spectral function `S(K,j)`, the toolbox computes the vertical modes and the associated horizontal kinetic energy (`HKE`), vertical kinetic energy (`VKE`), and potential energy (`PE`) as functions of horizontal wavenumber, vertical mode, and depth.

The main class, `InternalGravityWaveSpectrum`, first computes the non-hydrostatic vertical modes `F(z)` and `G(z)`. The vertical-mode calculation is performed through `GLOceanKit`, which uses stretched vertical coordinates and Chebyshev-polynomial quadrature to resolve sharp vertical structure such as pycnoclines and turning-depth behavior. The quadrature points define the effective vertical grid used by the toolbox.

After the modes are computed, the class builds energy coefficients for each energy component. A user-prescribed spectrum `S(K,j)` is then normalized and integrated over horizontal-wavenumber bands and vertical modes. Multiplying the band energy by the modal energy coefficients gives `HKE`, `VKE`, and `PE` on the native `(z,K,j)` grid.

## Dependencies

Core use requires:

- MATLAB.
- This repository's `src/` directory on the MATLAB path.
- `GLOceanKit`: <https://github.com/Energy-Pathways-Group/GLOceanKit>

The examples in `src/Example/example.mlx` also use:

- `VerticalModeAtlas`: <https://github.com/JeffreyEarly/vertical-mode-atlas>

`VerticalModeAtlas` is only needed to reproduce the example workflow that uses atlas-derived stratification profiles for Ocean Station PAPA and the Agulhas region. It is not required if you provide your own `N2(z)` function handle.

## Installation

Clone this repository, add its source directory to your MATLAB path, and add `GLOceanKit` to your MATLAB path:

```matlab
addpath("src")
addpath(genpath("path/to/GLOceanKit"))
```

If you want to run `src/Example/example.mlx`, also add `VerticalModeAtlas` and make the atlas NetCDF file available in the path used by the example.

The local `data/` and `old/` directories are not part of the GitHub repository workflow and are not required for normal use.

## Quick Start

This minimal example uses an idealized exponential stratification and does not require `VerticalModeAtlas`.

```matlab
latExp = 33;                % Latitude [degrees]
LzExp  = 4000;              % Ocean depth [m]

L_gm = 1300;                % E-folding scale depth [m]
N0   = 3 * 2*pi / 3600;     % Surface buoyancy frequency [rad/s]

% Buoyancy frequency squared N^2(z), with z < 0 below the surface
N2funcExp = @(z) (N0.^2) .* exp(2*z/L_gm);

nModes = 128;
nK     = 256;

igwExp = InternalGravityWaveSpectrum( ...
    N2funcExp, LzExp, ...
    nModes=nModes, nK=nK, latitude=latExp);

SExp = igwExp.generalSpectrum( ...
    j_star=3, ...
    slope_j=1, ...
    slope_k=1, ...
    A=1);

igwExp.assignEnergySpectrum(SExp);

zvect = linspace(-LzExp, 0, 1000);

HKE = igwExp.verticalVariance('HKE', 'zVector', zvect);
VKE = igwExp.verticalVariance('VKE', 'zVector', zvect);
PE  = igwExp.verticalVariance('PE',  'zVector', zvect);

figure
plot(HKE*1e4, zvect, VKE*1e4, zvect, PE*1e4, zvect)
xlabel('Variance [cm^2 s^{-2}]')
ylabel('Depth [m]')
legend('HKE', 'VKE', 'PE')
grid on
```

## Example Workflow

The live script `src/Example/example.mlx` demonstrates the workflow used for comparing three stratification profiles:

- idealized exponential stratification,
- Ocean Station PAPA, using `N2(z)` from `VerticalModeAtlas`,
- Agulhas region, using `N2(z)` from `VerticalModeAtlas`.

The example initializes one `InternalGravityWaveSpectrum` object per stratification, applies the same generalized spectrum to each case, assigns energy to the model, and compares how stratification changes the resulting energy distributions.

The main diagnostics used in the example are:

```matlab
verticalVariance('HKE', 'zVector', zvect)
verticalVariance('VKE', 'zVector', zvect)
verticalVariance('PE',  'zVector', zvect)

energyAtHorizontalWavenumber(z, 'HKE')
energyAtHorizontalWavenumber(z, 'VKE')
energyAtHorizontalWavenumber(z, 'PE')

energyAtFrequencies(z, 'HKE', omegaVector=omegaVector)
energyAtFrequencies(z, 'VKE', omegaVector=omegaVector)
energyAtFrequencies(z, 'PE',  omegaVector=omegaVector)
```

Rendered figures from this workflow are included in `img/`.

## Spectral Functions

The toolbox assigns energy using a spectrum function handle `S(K,j)`. The default workflow uses `generalSpectrum`, which provides tunable slopes in horizontal wavenumber and vertical mode:

```matlab
S = igw.generalSpectrum(j_star=3, slope_j=1, slope_k=1, A=1);
igw.assignEnergySpectrum(S);
```

You can also provide a custom spectrum function handle, for example from a region-specific characterization of the internal-wave field. The spectrum is normalized and integrated over the model's wavenumber bands before energy is assigned to each `(K,j)` bin.

Available spectrum tools include:

- `generalSpectrum`: tunable spectrum with vertical-mode and horizontal-wavenumber slopes.
- `gmSpectrum`: Garrett-Munk type internal-wave spectrum.
- `frequencyLocalizedSpectrum`: spectrum localized around a target frequency and vertical mode.
- `tidalSpectrum`: M2 tidal wrapper around `frequencyLocalizedSpectrum`.
- `fSpectrum`: near-inertial wrapper centered at the Coriolis frequency `f0`.
- `fM2Spectrum`: wrapper centered at `f0 + omega_M2`.

You can list the named spectrum constructors from MATLAB:

```matlab
igw.listAvailableSpectra();
```

## Random Realizations

The deterministic spectrum gives the expected total energy in each `(K,j)` band. The method `randomRealization` produces a stochastic realization by treating the internal-wave field as an ensemble of independent linear waves with zero-mean Gaussian amplitudes.

Each band energy is distributed among the independent waves in that band. The summed squared amplitude then follows a chi-squared statistic, so the random realization preserves the expected band energy while allowing finite-sample variability.

## Coordinate Conventions and Limitations

- Depth is negative below the surface: `z < 0`, with `z = 0` at the surface.
- Latitude must be away from the equator; the constructor rejects latitudes within about 5 degrees of the equator.
- Very high latitudes are also rejected by the constructor.
- Horizontal wavenumber `K` is radial and log-spaced.
- Frequencies satisfy approximately `f0 < omega < sqrt(N2max)`.
- Energy components are requested with `'HKE'`, `'VKE'`, `'PE'`, or `'TE'`.
- Some diagnostics currently assume scalar vertical-position queries.

## Examples and Figures

Examples are available in:

- `src/Example/example.mlx`
- `src/Example/QuadraturePoints.mlx`

Figures in `img/` show example outputs, including vertical-structure comparisons, horizontal-wavenumber spectra, frequency spectra, radius-of-deformation diagnostics, and quadrature-point diagnostics.

## Testing Status

The files in `src/unitTest/` are currently outdated and should not be treated as the recommended user-facing validation workflow. They are kept as development history and may need updates before being used as a reliable test suite.

## References

For the numerical vertical-mode calculation, see Early (2020) and Jeffrey et al. (2021).

## Author

Developed by Leticia Fabre-Lima, Jeffrey Early and Miles Sundermeyer for research on internal gravity wave spectra.
