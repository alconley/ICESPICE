### Mini Orange Spectrometer Simulation

This simulation is based off the advanced example named purging_magnet. 

A mini-orange spectrometer uses permanent magnets to seperate electrons from other types of particles. This reduces background and is used to study internal conversion electrons. The first mini-orange spectrometer was by [van Klinken](https://doi.org/10.1016/0029-554X(72)90416-8). Modern day examples include [fIREBALL](https://doi.org/10.1016/j.nima.2023.168288), [Mini-orange spectrometer at CIAE](https://iopscience.iop.org/article/10.1088/1674-1137/40/8/086002/pdf), and more.

The goal is to take an external magnetic field from [COMSOL](https://www.comsol.com) and simulate electron trajectories and the detector response in [GEANT4](https://geant4.web.cern.ch). 

### COMSOL Magnetic Field

In [COMSOL](https://www.comsol.com), I simulate the magnetic field of five 1"x1"x1/8" N42 Neodymium magnets around an attenuator made of tantalum (rod with a diameter of 1/8" and a height of ~30mm) in some vacuum volume. The attenuator is placed at the origin with the height in the z-direction. The magnets are placed in the circular pattern at some radius to form an toroidal magnetic field. If the source is placed at some positive z position away from the origin (f), the magnetic field must go in the negative phi (cylindrical coordinates) direction for the electrons to be focused at some negative z position away from the origin (g). 

The magnetic field is then exported to a .csv file. In this example I exported multiple values (X, Y, Z, Bx, By, Bz, B_norm/H_norm) at every 1/2mm point in a 10cm x 10cm x 10cm region, centered about the origin. For these dimensions, the file size is quite large (>1 Gb). I then have a python script [comsol_to_geant_table.py (just change the filepath to point to the correct comsol .csv file)] that changes the format for the file to the field file in the purging_magnet example.

### PIPS Detectors

A goal is also to create a simulation for our [PIPS Detectors](https://www.mirion.com/products/technologies/spectroscopy-scientific-analysis/research-education-and-industrial-solutions/passivated-implanted-planar-silicon-pips-detectors/standard-pips-detectors/pips-detectors-passivated-implanted-planar-silicon-detectors). At FSU, we have 4 different detectors. All detectors have an active area of 50mm^2, but have thicknesses 100, 300, 500, 1000 microns. The thickness of the detectors can be changed in the MiniOrangeDetectorConstruction.cc file. The energy deposited in the detector for each event is stored in a root file. The histogram can be plotted using

`root -x plotHisto.C`

