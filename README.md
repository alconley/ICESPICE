# ICESPICE GEANT4 Simulation

This simulation is inspired by the geant4's advanced example named `purging_magnet`.

## Overview
ICESPICE (Internal Conversion Electron SPectrometer In Coincidence Experiments) is a mini-orange spectrometer in devolpment at the John D. Fox Lab at Florida State University. ICESPICE uses permanent magnets to separate electrons from other types of particles, reducing background noise and enhancing the study of internal conversion electrons. The concept was first introduced by [van Klinken (1972)](https://doi.org/10.1016/0029-554X(72)90416-8), with modern examples including [fIREBALL](https://doi.org/10.1016/j.nima.2023.168288) and the [Mini-Orange Spectrometer at CIAE](https://iopscience.iop.org/article/10.1088/1674-1137/40/8/086002/pdf).

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

#### Adjusting Detector 

The detector geometry is build with SOLIDWORKS and exported as a .PLY ascii file. The detector can be easily modified in the ICESPICEDetectorConstruction.cc file.

#### Changing Detector Position

The position of the detector can also be modified to better understand its detection capabilities under different spatial configurations:

```bash
/ICESPICE/Detector/Position value
```

Where `value` can range from -60 to 0 millimeters and represents the value from the origin to the surface of the detector. (g)

#### Changing Source Properties

Adjust the properties of the simulation source to match specific experimental conditions or hypotheses. The following commands allow you to configure the primary aspects of the source:

- **Energy**: Set the initial energy of the particles emitted by the source.

```bash
/gun/energy value keV
```

Replace `value` with the desired energy level in kilo-electron volts (keV).

- **Position**: Specify the starting position of the particles along the z-axis. (f)

```bash
/gun/position 0 0 zvalue mm
```

Replace `zvalue` with the position in millimeters from the origin along the z-axis.

- **Particle Type**: Choose the type of particle to be emitted.


```bash
/gun/particle particle 

```

Replace `particle` with the type of particle you wish to simulate (e.g., `e+`, `e-`, `proton`, `neutron`, `gamma`). 

## Data Analysis

As I prefer not to use ROOT, I have opted to perform my data analysis using Python. The Python scripts that I utilize are stored in the 'python_scripts' directory. These scripts are primarily involved in generating macro files tailored to various simulations.

One of the key functions, transmission_histogram in the analysis.py file, allows for parsing output files using pandas. This function enables me to plot the energy deposited in the detector across different simulation energies using matplotlib. It also facilitates the calculation of the transmission probability, the probability that the full energy is deposited in the detector, and the detector's efficiency. Future enhancements will include unfolding the spectrum to correct for the detector's resolution.

Another function, transmission_probability, processes different file paths for various simulation energies under identical detector and f-position settings. This function plots the results against the transmission probability, illustrating the likelihood of full energy deposition in the detector.

Lastly, the plot_transmission_summary function mirrors the operations of transmission_probability but aggregates and plots all variations of the g-values together.