
def macro_creation(macro_path: str):
    # Open the file in write mode to add the initial commands and the loops
    with open(macro_path, 'w') as file:
        # Write the initial verbose and initialization commands
        
        file.write('/run/initialize\n')
        file.write('/control/verbose 0\n')
        file.write('/run/verbose 0\n')
        file.write('/event/verbose 0\n')
        file.write('/gps/pos/type Point\n')
        file.write('/gps/ang/type iso\n')
        file.write('/gps/pos/centre 0 0 70 mm\n')
        file.write('/ICESPICE/Detector/Position -30\n')
        file.write('/gps/particle ion\n')
        file.write('/gps/ion 83 207\n')
        file.write('/gps/ene/type Mono\n')
        file.write('/gps/ene/mono 0 eV\n')
        
        total_particle_goal = int(1e9)
        particles_per_run = int(1e6)
        print_progress = int(1e5)
        
        n_runs = int(total_particle_goal / particles_per_run)
        print("There will be ", n_runs, " runs")
        
        file.write(f'/run/printProgress {print_progress}\n')
        
        for i in range(n_runs):
            
            # print out the progress, must be done with /control/echo 
            file.write(f'/control/echo "Run {i+1}/{n_runs}"\n')
            
            file_name = f'Bi207_f70mm_g30mm_{i}.root'
            file.write(f'/analysis/setFileName {file_name}\n')
            file.write(f'/run/beamOn {particles_per_run}\n')
                        
            # Copy the file to the Git repository and push it
            file.write(f'/control/shell cp {file_name} /home/alconley/git_workspace/ICESPICE_analysis/\n')
            file.write(f'/control/shell cd /home/alconley/git_workspace/ICESPICE_analysis/ && git add {file_name} && git commit -m "Add {file_name}" && git push\n')


macro_creation('./build/MACRO_ICESPICE_PIPS1000_207Bi.mac')