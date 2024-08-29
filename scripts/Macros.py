import glob

# def macro_creation(n_particles: int, macro_path: str, thickness:int, f_position:int, g_position: int, folded_fwhm: float = 2.0):
def macro_creation(macro_path: str, thickness:int, folded_fwhm: float = 2.0):
    # Open the file in write mode to add the initial commands and the loops
    with open(macro_path, 'w') as file:
        # Write the initial verbose and initialization commands
        
        n_particles = 1000000

        file.write('/run/initialize\n')
        file.write('/control/verbose 2\n')
        file.write('/run/verbose 1\n')
        
        file.write('/gps/pos/type Point\n')
        file.write('/gps/ang/type iso\n')
        file.write('/gps/particle e-\n')
        
        for f_position in [70]:
            file.write(f'/gps/pos/centre 0 0 {f_position} mm\n')
            
            for g_position in range(-20,-50, -5):
            
                file.write(f'/ICESPICE/Detector/Position {g_position}\n')

                for energy in range(100, 2100, 100):
                    file.write(f'/gps/energy {energy} keV\n')
                    file_name = f'ICESPICE_PIPS{thickness}_f{f_position}mm_g{abs(g_position)}mm_{energy}keV.root'
                    file.write(f'/analysis/setFileName {file_name}\n')
                    
                    file.write(f'/run/printProgress 10000\n')
                    file.write(f'/run/beamOn {n_particles}\n')
                    
                    file.write(f"/control/shell root -x 'Plots.C('{file_name}', {energy}, {folded_fwhm}, {f_position}, {abs(g_position)}, {thickness})'\n")
                    
                    


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
    file.write('rm *t*.root\n')
    
    # hadd all the files together
    file.write('hadd ICESPICE.root *.root\n')
    
        