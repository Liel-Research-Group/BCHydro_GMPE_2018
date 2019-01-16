# BCHydro_GMPE_2018
This is a MatLab function for the 2018 Update to BC Hydro by Abrahamson et al. 2018.

## Running the Function
Call the function "BCHydro2018.m" from the matlab command line with the appropriate functions inputs to run the method and generate outputs. Output data is returned through the function handle. 

## Inputs
T = Period [sec]
M = Magnitude [Mw]
R = Rupture Distance [km]
Ztor = Depth to Slab [km]
Vs30 = Average shear wave velocity in top 30 m [m/s]
F is a inslab flag; F = 1 => Inslab, F = 0 => Interface

## Outputs
Sa = Spectral Acceleration [g]
sig = Sigma in ln space [unitless]
