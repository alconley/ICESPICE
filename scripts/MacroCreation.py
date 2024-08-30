import glob

# def macro_creation(n_particles: int, macro_path: str, thickness:int, f_position:int, g_position: int, folded_fwhm: float = 2.0):
def macro_creation(macro_path: str, thickness:int, folded_fwhm: float = 2.0):
    
    # Open the file in write mode to add the initial commands and the loops
    with open(macro_path, 'w') as file:
        # Write the initial verbose and initialization commands

        file.write('/run/initialize\n')
        file.write('/control/verbose 1\n')
        file.write('/event/verbose 0\n')

        file.write('/tracking/storeTrajectory 0\n')
        
        for f_position in [70]:
        
            # Write the commmand to change the gun position
            file.write(f'/gun/position 0 0 {f_position} mm\n')

            # # # Write the command to set the thickness
            # file.write(f'/ICESPICE/Detector/Thickness {thickness}\n')
            
            for g_position in range(-20,-50, -5): # -20 mm to -50 mm in steps of -5 mm
                n_particles = 100000
        
                # # Write the command to set the position
                file.write(f'/ICESPICE/Detector/Position {g_position}\n')

                # Loop through energies from 100 keV to 2000 keV in steps of 100 keV
                for energy in range(100, 2100, 100):

                    # Write the command to set the energy
                    file.write(f'/gun/energy {energy} keV\n')
                    
                    # Write the command to set the filename of the output root file
                    # file.write(f'/analysis/setFileName ICESPICE_PIPS{thickness}_f{f_position}mm_g{abs(g_position)}mm_{energy}.root\n')
                    
                    # Write the command to run the simulation
                    file.write(f'/run/beamOn {n_particles}\n')
                    
                    file.write(f"/control/shell root -x 'Plots.C({energy}, {folded_fwhm}, {f_position}, {abs(g_position)}, {thickness})'\n")
                    

def all_macros():
    # Loop through the positions from -20 to -55 in steps of -5
    # for detector in [100, 300, 500, 1000]:
    for detector in [1000]:
        for f_position in [70]:
            for g_position in range(-20,-50, -5): # -20 mm to -50 mm in steps of -5 mm
                n_particles = 10000
                # path = f'./build/MACRO_ICESPICE_PIPS{detector}_f{f_position}mm_g{abs(g_position)}mm.mac'
                path = f'./build/MACRO_ICESPICE_PIPS{detector}_f{f_position}mm_g{abs(g_position)}mm.mac'

                macro_creation(n_particles, macro_path=path, thickness=detector, f_position=f_position, g_position=g_position)

# all_macros()

macro_creation('./build/MACRO_ICESPICE_PIPS1000.mac', 1000, 2.0)

# create a bash script to run all the macros
with open('./build/run_all.sh', 'w') as file:
    file.write('#!/bin/bash\n')

    # look for all files that start with MACRO_ICESPICE
    for macro_file in glob.glob('./build/MACRO_ICESPICE*.mac'):
        # get rid of the build/ part of the path
        macro_file = macro_file[8:]
        file.write(f'./ICESPICE {macro_file}\n')
    
    # write a command to remove a file called 'ICESPICE.root' if it exists
    file.write('rm ICESPICE.root\n')
    
    # hadd all the files together
    file.write('hadd ICESPICE.root *.root\n')
    
        