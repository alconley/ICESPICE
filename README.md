# Mini Orange Spectrometer Simulation

This simulation is inspired by the advanced example named `purging_magnet`.

## Overview
A Mini-Orange Spectrometer uses permanent magnets to separate electrons from other types of particles, reducing background noise and enhancing the study of internal conversion electrons. The concept was first introduced by [van Klinken (1972)](https://doi.org/10.1016/0029-554X(72)90416-8), with modern examples including [fIREBALL](https://doi.org/10.1016/j.nima.2023.168288) and the [Mini-Orange Spectrometer at CIAE](https://iopscience.iop.org/article/10.1088/1674-1137/40/8/086002/pdf).

The goal of this simulation is to import an external magnetic field from [COMSOL](https://www.comsol.com) and simulate electron trajectories and detector responses using [GEANT4](https://geant4.web.cern.ch).

## Installation and Running

The field map, due to its large size (>1 GB), is not included in the repository. If required, please contact me for access via OneDrive or similar services.

**Prerequisites:**
- GEANT4 version: 11.2.0
- ROOT version: 6.30/06
- CMake version: 3.29.2

**Building and running the simulation:**
```bash
mkdir build
cd build
cmake ..
make
./ICESPICE
```

## COMSOL Magnetic Field Configuration

In [COMSOL](https://www.comsol.com), the magnetic field is simulated for an arrangement involving five 1"x1"x1/8" N42 Neodymium magnets. These magnets are arranged around a tantalum attenuator, which is a rod with a diameter of 1/8" and a height of approximately 30mm, placed within a vacuum volume. The setup is strategically designed to create a toroidal magnetic field, optimizing the path of electrons from a positive z position (source) towards a focused negative z position (detector).

The magnetic field configuration aims to ensure that when the source is positioned optimally, electrons are deflected in the negative phi direction (using cylindrical coordinates), enabling precise focusing at the detector's location.

### Magnetic Field Data Export

The magnetic field data from COMSOL is exported to a .csv file, which includes detailed spatial and vector field data at half-millimeter intervals within a 10cm x 10cm x 10cm region centered at the origin. Given the extensive data coverage and fine resolution, the file size is typically large (>1 GB). 

### Data Conversion for GEANT4

To facilitate the use of this COMSOL-generated magnetic field data in GEANT4, a Python script (`comsol_to_geant_table.py`) is provided. This script reformats the data from the COMSOL output to match the input requirements for the GEANT4 simulation, following the format used in the `purging_magnet` example. Users need to adjust the filepath in the script to point to the correct COMSOL .csv file.

For more information or to request the large field map file, please reach out via the contact methods provided in this repository.

## PIPS Detectors

This project includes simulations of [PIPS Detectors](https://www.mirion.com/products/technologies/spectroscopy-scientific-analysis/research-education-and-industrial-solutions/passivated-implanted-planar-silicon-pips-detectors/standard-pips-detectors/pips-detectors-passivated-implanted-planar-silicon-detectors).

### Detector Specifications

We utilize four different PIPS detectors, each characterized by a unique thickness:
- **100 microns**
- **300 microns**
- **500 microns**
- **1000 microns**

All detectors share an active area of 50 mm², making them suitable for detailed particle interaction studies.

### Simulation Configuration

#### Adjusting Detector Thickness

To simulate different detector behaviors based on their thickness, use the following command in the GEANT4 GUI:

```bash
/ICESPICE/Detector/Thickness value
```

Here, `value` can range from 0 to 3000 micrometers. 


#### Changing Detector Position

The position of the detector can also be modified to better understand its detection capabilities under different spatial configurations:

```bash
/ICESPICE/Detector/Position value
```

Where `value` can range from -60 to 0 millimeters and represents the value from the origin to the surface of the detector.

#### Changing Source Properties

Adjust the properties of the simulation source to match specific experimental conditions or hypotheses. The following commands allow you to configure the primary aspects of the source:

- **Energy**: Set the initial energy of the particles emitted by the source.

```bash
/gun/energy value keV
```

Replace `value` with the desired energy level in kilo-electron volts (keV).

- **Position**: Specify the starting position of the particles along the z-axis.

```bash
/gun/position 0 0 zvalue mm
```

Replace `zvalue` with the position in millimeters from the origin along the z-axis.

- **Particle Type**: Choose the type of particle to be emitted.


```bash
/gun/particle particle 

```

Replace `particle` with the type of particle you wish to simulate (e.g., `e+`, `e-`, `proton`, `neutron`, `gamma`). 

### Data Analysis

Energy deposition data for each detected event is stored in a ROOT file, which can be analyzed post-simulation to assess performance characteristics. To visualize the data, use the command:

```bash
root -x plotHisto.C
```

This command initiates a ROOT session that executes the `plotHisto.C` script, providing a graphical representation of the energy spectrum captured by the detectors.
