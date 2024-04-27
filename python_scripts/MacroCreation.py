# Define the file path for the macro
macro_path = './build/ICESPICE.mac'


def macro_creation(n_particles: int, macro_path: str, thickness:int, position: int):
    # Open the file in write mode to add the initial commands and the loops
    with open(macro_path, 'w') as file:
        # Write the initial verbose and initialization commands
        file.write('/control/verbose 1\n')
        file.write('/run/verbose 0\n')
        file.write('/event/verbose 0\n')
        
        # Write the command to set the thickness
        file.write(f'/ICESPICE/Detector/Thickness {thickness}\n')
    
        # Write the command to set the position
        file.write(f'/ICESPICE/Detector/Position {position}\n')
        
        # Loop through energies from 50 keV to 2000 keV in steps of 50 keV
        for energy in range(100, 2100, 50):
            # Write the command to set the energy
            file.write(f'/gun/energy {energy} keV\n')
            
            # Write the command to set the filename of the output root file
            file.write(f'/analysis/setFileName ICESPICE_PIPS{thickness}_f50mm_g{abs(position)}mm_{energy}.root\n')
            
            file.write('/run/initialize\n')
            
            # Write the command to run the simulation
            file.write(f'/run/beamOn {n_particles}\n')



for position in range(-20,-55, -5):
    n_particles = 10000
    path = f'./build/ICESPICE_{abs(position)}mm.mac'
    thickness = 1000

    macro_creation(n_particles, path, thickness, position)

# create a bash script to run all the macros
with open('./build/run_all.sh', 'w') as file:
    file.write('#!/bin/bash\n')

    for position in range(-20,-55, -5):

        # print out what position we are running in the bash script
        file.write(f'echo Running simulation for position {abs(position)}mm\n')
        file.write(f'./ICESPICE ICESPICE_{abs(position)}mm.mac\n')