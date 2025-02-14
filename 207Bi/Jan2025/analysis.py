import matplotlib.pyplot as plt

import sys
sys.path.append('/Users/alconley/Projects/ICESPICE')
sys.path.append('/home/alconley/git_workspace/ICESPICE')
from geant4analyzer import Geant4Analyzer
import numpy as np

# source $(brew --prefix root)/bin/thisroot.sh  # for ROOT

if __name__ == "__main__":

    geant4_path = "./207Bi/Jan2025/run_15_Jan2025_207BiDecay_PIPS1000_f9mm_g0mm_n75000000_Source100nmThick.root"
    # geant4_path = None
    no_icespice_analyzer = Geant4Analyzer(
        experimental_root_file_path="./207Bi/Jan2025/run_80_noICESPICE_f10mm_g0mm.root", 
        experimental_histogram_name="ICESPICE/PIPS1000/PIPS1000EnergyCalibrated",
        geant4_root_file_path=geant4_path,
        geant4_histogram_name="Esil"
    )
    
    fig, axes = plt.subplots(2, 3, sharex=True, figsize=(15, 8))
    axes = axes.flatten()

    exp_axes = axes[0]
    no_icespice_analyzer.set_experimental_histogram_noise_to_zero(25)
    no_icespice_analyzer.plot_experiment(ax=exp_axes, color='black', label="PIPS1000 w/o ICESPICE (run 56)")
    no_icespice_analyzer.GaussianFit(name='Experiment: 564-K', region_markers=(460,500), peak_markers=[481], ax=exp_axes)
    no_icespice_analyzer.GaussianFit(name='Experiment: 564-LM', region_markers=(535,585), peak_markers=[553, 582], ax=exp_axes)
    no_icespice_analyzer.GaussianFit(name='Experiment: 1064-K', region_markers=(940,1000), peak_markers=[975], ax=exp_axes)
    no_icespice_analyzer.GaussianFit(name='Experiment: 1064-LM', region_markers=(1030,1080), peak_markers=[1047, 1059], ax=exp_axes)
    average_fwhm = no_icespice_analyzer.average_fwhm()

    exp_axes.text(0.95, 0.95, f"Average FWHM {average_fwhm:.4}", transform=exp_axes.transAxes, verticalalignment='top', horizontalalignment='right')
    exp_axes.set_title(f"{no_icespice_analyzer.experimental_root_file_path[:50]}\n{no_icespice_analyzer.experimental_root_file_path[50:]}", fontsize=8)

    raw_sim_axes = axes[3]
    no_icespice_analyzer.plot_simulation(ax=raw_sim_axes)
    raw_sim_axes.text(0.95, 0.95, f"N Simulation = {sum(no_icespice_analyzer.geant4_bin_content)}", transform=raw_sim_axes.transAxes, verticalalignment='top', horizontalalignment='right')
    raw_sim_axes.set_title(f"{no_icespice_analyzer.geant4_root_file_path[:50]}\n{no_icespice_analyzer.geant4_root_file_path[50:]}", fontsize=8)
    raw_sim_axes.set_yscale('log')
    raw_sim_axes.set_ylim(bottom=1)

    no_icespice_analyzer.plot_experiment(ax=axes[1], color='black', label="PIPS1000 w/o ICESPICE")

    no_icespice_analyzer.gaussian_smear_simulation(average_fwhm)
    no_icespice_analyzer.plot_simulation(ax=axes[1], label="Geant4: Gaussion Smeared", linewidth=1)

    no_icespice_analyzer.scale_geant4_to_experiment(scaling_range=(370, 1075))
    no_icespice_analyzer.plot_simulation(ax=axes[1], label=f"Geant4: Best Scale={no_icespice_analyzer.scale:.5}", color='mediumpurple')

    no_icespice_analyzer.apply_threshold_to_geant4(threshold_range=(200,370), initial_guess=(200, 10000, 200))
    no_icespice_analyzer.plot_simulation(ax=axes[1], label=f"Geant4: w/threshold {no_icespice_analyzer.threshold_range}", color='deepskyblue', plot_uncertainity=False)

    no_icespice_analyzer.plot_simulation(ax=axes[4], color='black', plot_uncertainity=False)
    no_icespice_analyzer.GaussianFit(name='Simulation: 564-K', region_markers=(460,500), peak_markers=[481], ax=axes[4], use_geant4_data=True)
    no_icespice_analyzer.GaussianFit(name='Simulation: 564-LM', region_markers=(535,585), peak_markers=[553, 582], ax=axes[4], use_geant4_data=True)
    no_icespice_analyzer.GaussianFit(name='Simulation: 1064-K', region_markers=(940,1000), peak_markers=[975], ax=axes[4], use_geant4_data=True)
    no_icespice_analyzer.GaussianFit(name='Simulation: 1064-LM', region_markers=(1030,1080), peak_markers=[1047, 1059], ax=axes[4], use_geant4_data=True)


    axes[1].set_ylim(bottom=1)
    axes[1].set_yscale('log')
    axes[1].legend()

    no_icespice_analyzer.plot_residuals(ax=axes[5])


    # get fit parameters for Experiment
    _,_,_,_, area_564K, area_564K_uncertainity = no_icespice_analyzer.fit_stats(name='Experiment: 564-K')
    _,_,_,_, area_564L, area_564L_uncertainity = no_icespice_analyzer.fit_stats(name='Experiment: 564-LM', index=0)
    _,_,_,_, area_564M, area_564M_uncertainity = no_icespice_analyzer.fit_stats(name='Experiment: 564-LM', index=1)
    _,_,_,_, area_1064K, area_1064K_uncertainity = no_icespice_analyzer.fit_stats(name='Experiment: 1064-K')
    _,_,_,_, area_1064L, area_1064L_uncertainity = no_icespice_analyzer.fit_stats(name='Experiment: 1064-LM', index=0)
    _,_,_,_, area_1064M, area_1064M_uncertainity = no_icespice_analyzer.fit_stats(name='Experiment: 1064-LM', index=1)
    
    # get fit parameters for Simulation
    _,_,_,_, sim_area_564K, sim_area_564K_uncertainity = no_icespice_analyzer.fit_stats(name='Simulation: 564-K')
    _,_,_,_, sim_area_564L, sim_area_564L_uncertainity = no_icespice_analyzer.fit_stats(name='Simulation: 564-LM', index=0)
    _,_,_,_, sim_area_564M, sim_area_564M_uncertainity = no_icespice_analyzer.fit_stats(name='Simulation: 564-LM', index=1)
    _,_,_,_, sim_area_1064K, sim_area_1064K_uncertainity = no_icespice_analyzer.fit_stats(name='Simulation: 1064-K')
    _,_,_,_, sim_area_1064L, sim_area_1064L_uncertainity = no_icespice_analyzer.fit_stats(name='Simulation: 1064-LM', index=0)
    _,_,_,_, sim_area_1064M, sim_area_1064M_uncertainity = no_icespice_analyzer.fit_stats(name='Simulation: 1064-LM', index=1)
    
    # Compare areas of peaks
    ratio_564K = sim_area_564K / area_564K
    ratio_564L = sim_area_564L / area_564L
    ratio_564M = sim_area_564M / area_564M
    ratio_1064K = sim_area_1064K / area_1064K
    ratio_1064L = sim_area_1064L / area_1064L
    ratio_1064M = sim_area_1064M / area_1064M
    
    # Compare uncertainities of peaks
    ratio_564K_uncertainity = ratio_564K * np.sqrt((sim_area_564K_uncertainity/sim_area_564K)**2 + (area_564K_uncertainity/area_564K)**2)
    ratio_564L_uncertainity = ratio_564L * np.sqrt((sim_area_564L_uncertainity/sim_area_564L)**2 + (area_564L_uncertainity/area_564L)**2)
    ratio_564M_uncertainity = ratio_564M * np.sqrt((sim_area_564M_uncertainity/sim_area_564M)**2 + (area_564M_uncertainity/area_564M)**2)
    ratio_1064K_uncertainity = ratio_1064K * np.sqrt((sim_area_1064K_uncertainity/sim_area_1064K)**2 + (area_1064K_uncertainity/area_1064K)**2)
    ratio_1064L_uncertainity = ratio_1064L * np.sqrt((sim_area_1064L_uncertainity/sim_area_1064L)**2 + (area_1064L_uncertainity/area_1064L)**2)
    ratio_1064M_uncertainity = ratio_1064M * np.sqrt((sim_area_1064M_uncertainity/sim_area_1064M)**2 + (area_1064M_uncertainity/area_1064M)**2)
    
    # plot ratio vs energy
    ratio_x = [481.6935, 553.8372, 565.8473, 975.651, 1047.795, 1059.805]
    ratio_y = [ratio_564K, ratio_564L, ratio_564M, ratio_1064K, ratio_1064L, ratio_1064M]
    ratio_y_uncertainity = [ratio_564K_uncertainity, ratio_564L_uncertainity, ratio_564M_uncertainity, ratio_1064K_uncertainity, ratio_1064L_uncertainity, ratio_1064M_uncertainity]
    
    axes[2].errorbar(ratio_x, ratio_y, yerr=ratio_y_uncertainity, fmt='o', color='black', markersize=2)
    axes[2].set_xlabel("Energy (keV)")
    axes[2].set_ylabel(r"Area$_{Geant4}$ / Area$_{Experiment}$")
    

    # no_icespice_analyzer.plot_simulation(ax=axes[3])


    no_icespice_analyzer.print_fit_parameters()

    # icespice_analyzer.plot_experiment(ax=axes[0], color='black', label="PIPS1000 w/ ICESPICE (run 55)")
    # icespice_analyzer.GaussianFit(name='564-K', region_markers=(465,495), peak_markers=[481], ax=axes[0])
    # icespice_analyzer.GaussianFit(name='564-LM', region_markers=(540,580), peak_markers=[553, 582], ax=axes[0])
    # icespice_analyzer.GaussianFit(name='1064-K', region_markers=(960,995), peak_markers=[975], ax=axes[0])
    # icespice_analyzer.GaussianFit(name='1064-LM', region_markers=(1030,1080), peak_markers=[1047, 1059], ax=axes[0])


    # icespice_analyzer.print_fit_parameters()

    # n_564K_ratio = icespice_analyzer.fits['564-K'].params['g0_area'].value / no_icespice_analyzer.fits['564-K'].params['g0_area'].value
    # n_564L_ratio = icespice_analyzer.fits['564-LM'].params['g0_area'].value / no_icespice_analyzer.fits['564-LM'].params['g0_area'].value
    # n564M_ratio = icespice_analyzer.fits['564-LM'].params['g1_area'].value / no_icespice_analyzer.fits['564-LM'].params['g1_area'].value
    # n_1064K_ratio = icespice_analyzer.fits['1064-K'].params['g0_area'].value / no_icespice_analyzer.fits['1064-K'].params['g0_area'].value
    # n_1064L_ratio = icespice_analyzer.fits['1064-LM'].params['g0_area'].value / no_icespice_analyzer.fits['1064-LM'].params['g0_area'].value
    # n1064M_ratio = icespice_analyzer.fits['1064-LM'].params['g1_area'].value / no_icespice_analyzer.fits['1064-LM'].params['g1_area'].value

    # # plot ratio vs energy
    # ratio_x = [481.6935, 553.8372, 565.8473, 975.651, 1047.795, 1059.805]
    # ratio_y = [n_564K_ratio, n_564L_ratio, n564M_ratio, n_1064K_ratio, n_1064L_ratio, n1064M_ratio]

    # axes[2].plot(ratio_x, ratio_y, 'o', color='black', markersize=2)
    # axes[2].set_xlabel("Energy (keV)")
    # axes[2].set_ylabel(r"$N_{ICESPICE}$ / $N_{Without\,ICESPICE}$")

    # for ax in axes:
    #     if ax != axes[0]:
    #         ax.set_xlim(200, 1200)

    fig.tight_layout()
    
    plt.show()
