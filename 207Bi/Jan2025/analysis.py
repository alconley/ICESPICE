import matplotlib.pyplot as plt

import sys
sys.path.append('/home/alconley/git_workspace/ICESPICE')
from geant4analyzer import Geant4Analyzer

# source $(brew --prefix root)/bin/thisroot.sh  # for ROOT

if __name__ == "__main__":
    simulation_root_file = "./207Bi/Sept2024_LSU/geant_sim/run_98_ICESPICE_RadDecay_z83_a207_e0keV_f70mm_g30mm_n100000000_PIPS1000_AllProcesses_Si02Window50nm_Source500nmThick.root"
    simulation_histogram_name = "Esil"

    experimental_root_file = "./207Bi/Sept2024_LSU/exp_data/207Bi_ICESPICE_f70mm_g30mm_run_14_15.root"
    experimental_histogram_name = "PIPS1000EnergyCalibrated"

    analyzer = Geant4Analyzer(geant4_root_file_path=simulation_root_file, geant4_histogram_name=simulation_histogram_name, experimental_root_file_path=experimental_root_file, experimental_histogram_name=experimental_histogram_name)

    # axes = analyzer.plot(experiment=True, simulation=True)
    # axes[1].set_yscale('log')
    # plt.show()

    analyzer.gaussian_smear_simulation(fwhm=10)
    analyzer.scale_geant4_to_experiment(scaling_range=(462, 1075))
    analyzer.apply_threshold_to_geant4(threshold_range=(300, 462))
    axes = analyzer.plot(experiment=True, simulation=True, residuals=True, same_axes=True)
    
    for ax in axes:
        ax.set_xlim(200, 1200)

    plt.show()
