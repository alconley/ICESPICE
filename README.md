# ICESPICE Geant4 Simulation

This simulation is inspired by the geant4's advanced example named `purging_magnet`.

## Overview
ICESPICE (Internal Conversion Electron SPectrometer In Coincidence Experiments) is a mini-orange spectrometer in devolpment at the John D. Fox Lab at Florida State University. ICESPICE uses permanent magnets to separate electrons from other types of particles, reducing background noise and enhancing the study of internal conversion electrons. The concept was first introduced by [van Klinken (1972)](https://doi.org/10.1016/0029-554X(72)90416-8), with modern examples including [fIREBALL](https://doi.org/10.1016/j.nima.2023.168288) and the [Mini-Orange Spectrometer at CIAE](https://iopscience.iop.org/article/10.1088/1674-1137/40/8/086002/pdf).

The goal of this simulation is to import an external magnetic field from [COMSOL](https://www.comsol.com) and simulate electron trajectories and detector responses using [Geant4](https://geant4.web.cern.ch).

## Installation and Running

**Prerequisites:**
- Geant4 version: 11.2.0
- ROOT version: 6.30/06
- CMake version: 3.29.2

**Installing Geant4**
I have included a script to install Geant4 on MacOS. A similar script exists at [geant4_setup_tools](https://github.com/eli-temanson/geant4_setup_tools) for Ubuntu 22.04.

**Building and running the simulation:**
```bash
mkdir build
cd build
cmake ..
make
./ICESPICE
```

## COMSOL Magnetic Field Configuration

In [COMSOL](https://www.comsol.com), the magnetic field is simulated for an arrangement involving five 1"x1"x1/8" N42 Neodymium magnets. These magnets are arranged around a tantalum attenuator, which is a rod with a diameter of 1/8" and a height of approximately 30mm, placed within a vacuum volume. The setup is designed to create a toroidal magnetic field, optimizing the path of electrons from a positive z position (source) towards a focused negative z position (detector).

### Magnetic Field Data Export

The magnetic field data from COMSOL is exported to a .csv file, which includes spatial and vector field data at half-millimeter intervals within a 10cm x 10cm x 14cm region centered at the origin. Given the extensive data coverage and fine resolution, the file size is typically large (>1 GB) and could not be put on github. The COMSOL output must be converted to the correct format (same as the "purging_magnet" example) for the code to read it. This script can be found in ./scripts/comsol_to_geant_table.py.

## Geometry

The geometry of ICESPICE is easily imported into Geant4 using [CADMESH](https://github.com/christopherpoole/CADMesh). From SolidWorks, the assembly geometry is exported to a .step file. This is then imported into [FreeCad](https://www.freecad.org/) and then exported as a .obj file. This is annoying but SolidWorks doesn't allow you to export an assembly as an .obj file. The .obj file can be easily read into geant4 using CADMESH. I had to change the groups ('g ') in the .obj file to objects ('o '). The name of the objects must not have spaces too.

## Scripts
Multiple scripts exist in different folders for the simulations of different purposes.

