import matplotlib.pyplot as plt

import sys
sys.path.append('/Users/alconley/Projects/ICESPICE')
sys.path.append('/home/alconley/git_workspace/ICESPICE')
from geant4analyzer import Geant4Analyzer

# source $(brew --prefix root)/bin/thisroot.sh  # for ROOT

if __name__ == "__main__":

    geant4_path = "../207Bi/Jan2025/run_7_Jan2025_207BiDecay_PIPS1000_f10mm_g0mm_n150000000_Source100nmThick.root"
    # geant4_path = None
    no_icespice_analyzer = Geant4Analyzer(
        experimental_root_file_path="../207Bi/Jan2025/run_80_noICESPICE_f10mm_g0mm.root", 
        experimental_histogram_name="ICESPICE/PIPS1000/PIPS1000EnergyCalibrated",
        geant4_root_file_path=geant4_path,
        geant4_histogram_name="Esil"
    )
    
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(10, 8))

    no_icespice_analyzer.set_experimental_histogram_noise_to_zero(25)
    no_icespice_analyzer.plot_experiment(ax=axes[0], color='black', label="PIPS1000 w/o ICESPICE (run 56)")
    no_icespice_analyzer.GaussianFit(name='564-K', region_markers=(460,500), peak_markers=[481], ax=axes[0])
    no_icespice_analyzer.GaussianFit(name='564-LM', region_markers=(535,585), peak_markers=[553, 582], ax=axes[0])
    no_icespice_analyzer.GaussianFit(name='1064-K', region_markers=(960,995), peak_markers=[975], ax=axes[0])
    no_icespice_analyzer.GaussianFit(name='1064-LM', region_markers=(1030,1080), peak_markers=[1047, 1059], ax=axes[0])

    no_icespice_analyzer.print_fit_parameters()
    average_fwhm = no_icespice_analyzer.average_fwhm()

    no_icespice_analyzer.gaussian_smear_simulation(average_fwhm)
    no_icespice_analyzer.scale_geant4_to_experiment(scaling_range=(370, 1075))
    no_icespice_analyzer.apply_threshold_to_geant4(threshold_range=(200,370))

    no_icespice_analyzer.plot_simulation(ax=axes[0])

    # no_icespice_analyzer.GaussianFit(name='564-K', region_markers=(460,500), peak_markers=[481], ax=axes[0], use_geant4_data=True)
    # no_icespice_analyzer.GaussianFit(name='564-LM', region_markers=(535,585), peak_markers=[553, 582], ax=axes[0], use_geant4_data=True)
    # no_icespice_analyzer.GaussianFit(name='1064-K', region_markers=(960,995), peak_markers=[975], ax=axes[0], use_geant4_data=True)
    # no_icespice_analyzer.GaussianFit(name='1064-LM', region_markers=(1030,1080), peak_markers=[1047, 1059], ax=axes[0], use_geant4_data=True)

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
    #     ax.set_xlim(200, 1200)

    fig.tight_layout()
    
    plt.show()
