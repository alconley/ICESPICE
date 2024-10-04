import polars as pl
import matplotlib.pyplot as plt
import numpy as np
import ROOT
from scipy.optimize import minimize
import argparse  
import os
import glob
import lmfit
import warnings

# for virtual environment
# source $(brew --prefix root)/bin/thisroot.sh  # for ROOT

# Function to parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Run analysis on a specified ROOT file")
    parser.add_argument("root_file_path", nargs='?', default=None, help="Path to the ROOT file")
    parser.add_argument("--icespice", action="store_true", default=False, help="If ICESPICE in the simultation")
    parser.add_argument("--fwhm", type=float, default=10, help="FWHM value for Gaussian smearing (default: 10.0)")
    parser.add_argument("--save-pic", action="store_true", default=False, help="Save the plot as an image if this flag is set (default: False)")
    parser.add_argument("--save-path", type=str, default="picture.png", help="Path to save the plot image (default: picture.png)")

    return parser.parse_args()

# Function to get list of .root files
def list_root_files(directory):
    return glob.glob(os.path.join(directory, "*.root"))

# Function to get the histogram data from a ROOT file
def get_root_hist_data(root_file_path: str, histogram_name: str, print=False):

    root_file = ROOT.TFile(root_file_path, "READ") # Open the ROOT file    
    histogram = root_file.Get(histogram_name) # Get the histogram by name
    
    sim_hist = []
    sim_hist_uncertainity = []
    sim_bin_edges = []
    sim_bin_centers = []
    
    
    for i in range(1, histogram.GetNbinsX() + 1):  # ROOT bins are 1-indexed
        bin_content = histogram.GetBinContent(i) # Get the bin content
        bin_center = histogram.GetBinCenter(i) * 1000  # Convert to keV
        bin_width = histogram.GetBinWidth(i) * 1000  # Convert to keV
        bin_uncertainty = histogram.GetBinError(i) # Get the bin uncertainty
        bin_edge = bin_center - 0.5 * bin_width # Get the bin edge from the center and the width
        
        if print:
            print(f"Bin {i}: Center = {bin_center}, Lower Edge = {bin_edge}, Width = {bin_width}, Content = {bin_content}, Uncertainty = {bin_uncertainty}")
        
        sim_hist.append(bin_content)
        sim_hist_uncertainity.append(bin_uncertainty)
        sim_bin_edges.append(bin_edge)
        sim_bin_centers.append(bin_center)
        
        if i == 2000:
            bin_edge = bin_center + 0.5 * bin_width # add the last bin edge
            sim_bin_edges.append(bin_edge)
            break
    
    return sim_hist, sim_bin_centers, sim_bin_edges, sim_hist_uncertainity

# Function to smear a histogram bin with a Gaussian
def gaussian_smear(bin_contents, bin_centers, bin_uncertainity, fwhm):
    """
    Apply Gaussian smearing to histogram data.
    
    Parameters:
    - bin_centers: List of bin centers.
    - bin_contents: List of bin contents (counts).
    - fwhm: Full width at half maximum (FWHM) of the Gaussian distribution.

    Returns:
    - smeared_contents: Gaussian smeared bin contents.
    """
    
    # Ensure bin_centers is a numpy array for mathematical operations
    bin_centers = np.array(bin_centers)
    
    total_counts = np.sum(bin_contents)  # Total counts before smearing
    
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to sigma
    
    smeared_bin_contents = np.zeros_like(bin_contents)
    smeared_bin_contents_low = np.zeros_like(bin_contents)
    smeared_bin_contents_high = np.zeros_like(bin_contents)
    
    # Loop over all bins and distribute each bin's content over neighboring bins
    for i in range(len(bin_centers)):
        bin_center = bin_centers[i]
        
        bin_content = bin_contents[i]
        bin_content_low = bin_contents[i] - bin_uncertainity[i]
        bin_content_high = bin_contents[i] + bin_uncertainity[i]
        
        # Generate a Gaussian centered at bin_center
        weights = np.exp(-0.5 * ((bin_centers - bin_center) / sigma) ** 2)
        weights /= np.sum(weights)  # Normalize the weights to sum to 1
        
        smeared_bin_contents += bin_content * weights  # Apply the smearing with normalized weights
        smeared_bin_contents_low += bin_content_low * weights  # Apply the smearing with normalized weights
        smeared_bin_contents_high += bin_content_high * weights  # Apply the smearing with normalized weights
        
    
    smeared_total_counts = np.sum(smeared_bin_contents)  # Total counts after smearing
    
    if total_counts != smeared_total_counts:
        warnings.warn(f"Total counts before and after smearing do not match: {total_counts} vs {smeared_total_counts}")
    
    return smeared_bin_contents, smeared_bin_contents_low, smeared_bin_contents_high

def chi_squared(params, exp_hist, exp_unc, sim_hist, sim_unc, exp_bins, sim_bins):
    """
    Chi-squared function to minimize the difference between experimental and simulated data, including uncertainties.

    Parameters:
    - params: lmfit Parameters object containing the 'scale' parameter.
    - exp_hist: Experimental histogram counts.
    - exp_unc: Experimental histogram uncertainties (e.g., sqrt(counts)).
    - sim_hist: Simulated histogram counts.
    - sim_unc: Simulated histogram uncertainties (e.g., sqrt(counts)).
    - exp_bins: Experimental bin edges.
    - sim_bins: Simulated bin edges.

    Returns:
    - chi2: The chi-squared value.
    """
    scale = params['scale'].value  # Extract the scaling factor
    
    # Scale the simulated histogram by the scale factor
    scaled_sim_hist, _ = np.histogram(sim_bins[:-1], bins=exp_bins, weights=sim_hist * scale)

    # Scale the simulated uncertainties
    scaled_sim_unc = sim_unc * scale

    # Combine the uncertainties (experimental and simulated)
    total_unc = np.sqrt(exp_unc ** 2 + scaled_sim_unc ** 2)

    # Calculate the chi-squared value
    chi2 = np.sum(((exp_hist - scaled_sim_hist) ** 2) / total_unc ** 2)

    return chi2

# Main part of your code
if __name__ == "__main__":
    
    # Parse the command-line arguments
    args = parse_args()

    # Use the parsed root_file_path from the arguments, or list files if not provided
    if args.root_file_path is None:
        root_file_dir = "../207Bi/geant_sim/"
        root_files = list_root_files(root_file_dir)
        if not root_files:
            print(f"No .root files found in {root_file_dir}")
            exit(1)
        print("Available .root files:")
        for i, file in enumerate(root_files):
            print(f"{i + 1}: {file}")
        file_index = int(input("Select a file by number: ")) - 1
        root_file_path = root_files[file_index]
    else:
        root_file_path = args.root_file_path

    df_withICESPICE = pl.read_parquet("../207Bi/exp_data/207Bi_ICESPICE_f70mm_g30mm_run_*.parquet")
    df_withoutICESPICE = pl.read_parquet("../207Bi/exp_data/207Bi_noICESPICE_f9mm_g0mm_run_13.parquet")

    # Energy calibration of m=0.5395 and b=2.5229
    df_withICESPICE = df_withICESPICE.with_columns([(pl.col("PIPS1000Energy") * 0.5395 + 2.5229).alias("PIPS1000EnergyCalibrated")])
    df_withoutICESPICE = df_withoutICESPICE.with_columns([(pl.col("PIPS1000Energy") * 0.5395 + 2.5229).alias("PIPS1000EnergyCalibrated")])

    ###############################################################################################################
    
    linewidth = 1
    fs = 10
    
    exp_fig, exp_axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True, num="Experimental Data")
    exp_axs = exp_axs.flatten()
    
    ###############################################################################################################
    # The experimental histogram with and without ICESPICE
    
    exp_bins = 1100
    exp_range = [0, 1100]

    exp_hist_withICESPICE, exp_bin_edges_withICESPICE = np.histogram(df_withICESPICE["PIPS1000EnergyCalibrated"], bins=exp_bins, range=exp_range)
    exp_hist_withoutICESPICE, exp_bin_edges_withoutICESPICE = np.histogram(df_withoutICESPICE["PIPS1000EnergyCalibrated"], bins=exp_bins, range=exp_range)
    
    exp_hist_withICESPICE_uncertainity = np.sqrt(exp_hist_withICESPICE)
    exp_hist_withoutICESPICE_uncertainity = np.sqrt(exp_hist_withoutICESPICE)
    
    exp_axs[0].stairs(values=exp_hist_withoutICESPICE, edges=exp_bin_edges_withoutICESPICE, color="#A6192E", label=r"$^{207}$Bi without ICESPICE", linewidth=linewidth)
    exp_axs[1].stairs(values=exp_hist_withICESPICE, edges=exp_bin_edges_withICESPICE, color="#5CB8B2", label=r"$^{207}$Bi with ICESPICE", linewidth=linewidth)
    
    # scale of the hist of the data with out ICESPICE so the counts match the 975 peak
    scale = 51700/65625
    exp_axs[2].stairs(values=exp_hist_withoutICESPICE*scale, edges=exp_bin_edges_withoutICESPICE, color="#A6192E", label=r"$^{207}$Bi without ICESPICE (Scaled)", linewidth=linewidth)
    exp_axs[2].stairs(values=exp_hist_withICESPICE, edges=exp_bin_edges_withICESPICE, color="#5CB8B2", label=r"$^{207}$Bi with ICESPICE", linewidth=linewidth)
    
    bbox = dict(facecolor='white', edgecolor='none')
    exp_axs[2].text(481, 3300, r"570K", horizontalalignment='center', verticalalignment='bottom', fontsize=fs, rotation=90, bbox=bbox)
    exp_axs[2].text(553, 1500, r"570L", horizontalalignment='center', verticalalignment='bottom', fontsize=fs, rotation=90, bbox=bbox)
    exp_axs[2].text(565, 1300, r"570M", horizontalalignment='left', verticalalignment='bottom', fontsize=fs, rotation=90, bbox=bbox)
    exp_axs[2].text(975, 5150, r"1064K", horizontalalignment='center', verticalalignment='bottom', fontsize=fs, rotation=90, bbox=bbox)
    exp_axs[2].text(1047, 1400, r"1064L", horizontalalignment='center', verticalalignment='bottom', fontsize=fs, rotation=90, bbox=bbox)
    exp_axs[2].text(1059, 700, r"1064M", horizontalalignment='left', verticalalignment='bottom', fontsize=fs, rotation=90, bbox=bbox)

    exp_axs[2].set_ylim(0, 6500)
    
    for ax in exp_axs:
        ax.set_xlabel(r"Energy [keV]")
        ax.set_ylabel(r"Counts/keV")
        ax.axvline(x=481.6935, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=553.8372, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=565.8473, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=975.651, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=1047.795, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=1059.805, color='green', linestyle='--', linewidth=1)
        
        ax.legend(loc='upper left', shadow=False, frameon=True, fancybox=False, edgecolor='none', facecolor='none')
        ax.set_xlim(0, 1100)   
        ax.minorticks_on()
        ax.tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=3)
        ax.tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=5)
    
    exp_fig.tight_layout()
    
    ###############################################################################################################
    geant_fig, geant_axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True, num="Geant4 Simulation")
    
    # Get the simulation data from the ROOT file and plot the raw data
    
    sim_hist, sim_bin_centers, sim_bin_edges, sim_hist_uncertainity = get_root_hist_data(root_file_path, "Esil", print=False)

    geant_axs[0].stairs(values=sim_hist, edges=sim_bin_edges, color="dodgerblue", label="Geant4 simulation (Raw)", linewidth=linewidth)

    n_sim_particles = np.sum(sim_hist)
    # in scientific notation
    geant_axs[0].text(0.05, 0.95, f"Total Counts: {n_sim_particles:.1e}", transform=geant_axs[0].transAxes, ha='left', va='top')

    fwhm = args.fwhm  # Get the FWHM value from the argument
    
    # Smear the simulation data with a Gaussian and plot the smeared data
    sim_hist_smeared, sim_hist_smeared_low, sim_hist_smeared_high = gaussian_smear(bin_contents=sim_hist, 
                                                                                   bin_centers=sim_bin_centers, 
                                                                                   bin_uncertainity=sim_hist_uncertainity, 
                                                                                   fwhm=fwhm)
    
    geant_axs[1].stairs(values=sim_hist_smeared, edges=sim_bin_edges, color="dodgerblue", label=f"Geant4 simulation (Smeared: {fwhm} keV)", linewidth=linewidth)
    geant_axs[1].fill_between(sim_bin_centers, sim_hist_smeared_low, sim_hist_smeared_high, color='dodgerblue', alpha=0.2, label="Uncertainty")
    geant_axs[1].set_ylim(top=sim_hist_smeared.max() * 1.1)
    
    for ax in geant_axs:
        ax.legend(loc='upper right', shadow=False, frameon=True, fancybox=False, edgecolor='none', facecolor='none')
        ax.minorticks_on()
        ax.tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=3)
        ax.tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=5)
        ax.set_xlabel(r"Energy [keV]")
        ax.set_ylabel(r"Counts/keV")
        ax.set_yscale("log")
        ax.set_ylim(bottom=1)
        ax.set_xlim(0, 2000)
        
        ax.axvline(x=481.6935, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=553.8372, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=565.8473, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=975.651, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=1047.795, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=1059.805, color='green', linestyle='--', linewidth=1)
        
    geant_fig.tight_layout()
    
    ###############################################################################################################
    if args.icespice:
        exp_v_geant_fig,  exp_v_geant_axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True, num="Experiment vs Simulation with ICESPICE")
        exp_v_geant_axs = exp_v_geant_axs.flatten()
        
        exp_v_geant_axs[0].stairs(values=exp_hist_withICESPICE, edges=exp_bin_edges_withICESPICE, color="#5CB8B2", label=r"$^{207}$Bi with ICESPICE", linewidth=linewidth)
    else:
        exp_v_geant_fig,  exp_v_geant_axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True, num="Experiment vs Simulation without ICESPICE")
        exp_v_geant_axs = exp_v_geant_axs.flatten()
        
        exp_v_geant_axs[0].stairs(values=exp_hist_withoutICESPICE, edges=exp_bin_edges_withoutICESPICE, color="#A6192E", label=r"$^{207}$Bi without ICESPICE", linewidth=linewidth)
        
        scale_range = [280, 380]
        
        
    for ax in exp_v_geant_axs:
        ax.set_xlabel(r"Energy [keV]")
        ax.set_ylabel(r"Counts/keV")
        ax.axvline(x=481.6935, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=553.8372, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=565.8473, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=975.651, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=1047.795, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=1059.805, color='green', linestyle='--', linewidth=1)
        
        ax.legend(loc='upper left', shadow=False, frameon=True, fancybox=False, edgecolor='none', facecolor='none')
        ax.set_xlim(0, 1100)   
        ax.minorticks_on()
        ax.tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=3)
        ax.tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=5)

    exp_v_geant_fig.tight_layout()
    
    ###############################################################################################################
    
    save_pic = args.save_pic
    save_path = args.save_path

    # If save-pic is False, ignore save-path
    if save_pic:
        print(f"Saving the plot to {save_path}")
        # Save the plot to the provided path
        plt.savefig(save_path, dpi=300)
        
    plt.show()



        
        
    # Filter the simulation data so the binning is the same
    # filtered_sim_hist = sim_hist_smeared[(sim_bin_centers >= exp_range[0]) & (sim_bin_centers <= exp_range[1])]
        
    # geant_bin_centers, geant_bin_counts, geant_bin_counts_low, geant_bin_counts_high, n_sim_particles, n_interactions = get_root_hist_and_gauss_smear(root_file_path, "Esil", fwhm=fwhm)
    # geant_count_uncertainty = (geant_bin_counts_high - geant_bin_counts_low) / 2
    
    # # convert bin_centers to left edges
    # geant_bin_edges = geant_bin_centers - 0.5 * (geant_bin_centers[1] - geant_bin_centers[0])

    # # Filter simulation data to only include energies between 200-1200 keV
    # filtered_geant_bin_centers = geant_bin_centers[(geant_bin_centers >= 400) & (geant_bin_centers <= 1200)]
    # filtered_geant_bin_counts = geant_bin_counts[(geant_bin_centers >= 400) & (geant_bin_centers <= 1200)]
    # filtered_geant_bin_counts_low = geant_bin_counts_low[(geant_bin_centers >= 400) & (geant_bin_centers <= 1200)]
    # filtered_geant_bin_counts_high = geant_bin_counts_high[(geant_bin_centers >= 400) & (geant_bin_centers <= 1200)]
    # filtered_geant_bin_count_uncertainty = (filtered_geant_bin_counts_high - filtered_geant_bin_counts_low) / 2

    # # Make sure the binning of the real and sim data matches
    # sim_hist, sim_bins = np.histogram(filtered_geant_bin_centers, bins=exp_bin_edges, weights=filtered_geant_bin_counts)

    # # Initialize lmfit parameters
    # params = lmfit.Parameters()
    # params.add('scale', value=10.0)  # Initial scale guess

    # # Perform chi-squared minimization using lmfit
    # minimizer = lmfit.Minimizer(chi_squared_func, params, fcn_args=(exp_hist, exp_hist_uncertainity, sim_hist, filtered_geant_bin_counts_low[:-1], filtered_geant_bin_counts_high[:-1]))
    # result = minimizer.minimize()
    
    

    # best_scale = result.params['scale'].value    
    # result.params.pretty_print()
    
    # print(f"Best scale: {best_scale:.3f}")
    
    # # # Apply the best scale to the simulation
    # scaled_sim_filtered_hist = sim_hist * best_scale

    # ## scale the original simulation data (geant_bin_counts) 
    # scaled_sim_og_hist = geant_bin_counts * best_scale

    # plot the histogram of the Geant4 simulation without ICESPICE
    # if args.icespice:
    #     axs[0].hist(df_withICESPICE["PIPS1000EnergyCalibrated"], bins=1200, range=[0, 1200], histtype="step", color='#5CB8B2', label="without ICESPICE", linewidth=linewidth)
    # else:
    #     axs[0].hist(df_withoutICESPICE["PIPS1000EnergyCalibrated"], bins=1200, range=[0, 1200], histtype="step", color='#5CB8B2', label="without ICESPICE", linewidth=linewidth)
    

    # axs[0].step(geant_bin_centers, geant_bin_counts * best_scale, color="dodgerblue", label="Scaled Geant4 simulation", linewidth=1, where="mid")    
    # axs[0].fill_between(geant_bin_centers, geant_bin_counts_low * best_scale, geant_bin_counts_high * best_scale, color='dodgerblue', alpha=0.2, label="Uncertainty")

    # axs[0].text(0.49, 0.95, f"Scale factor: {best_scale:.3f}\nTotal Counts: {n_sim_particles}\nSimulation Counts in Detector: {n_interactions}", transform=axs[0].transAxes, ha='left', va='top')



    # axs[0].set_ylim(1, 100000)
    # axs[0].set_yscale("log")

# \nSimulation particles: {n_sim_particles:.1e}
    ###############################################################################################################



    # axs[1].plot(real_bins[:-1] + 0.5, (real_hist - scaled_sim_filtered_hist)/real_hist, marker='o', markersize=1, color="black", linestyle='None', label="Exp-Sim/Exp")
    
    # # axs[1].errorbar(real_bins[:-1] + 0.5, (real_hist - scaled_sim_filtered_hist), yerr= np.sqrt(filtered_geant_bin_count_uncertainty[:-1]**2+real_hist_error**2), fmt='none', color='black', capsize=2, capthick=1, elinewidth=1)
    # axs[1].set_xlabel(r"Energy [keV]")
    # axs[1].set_ylabel(r"Exp-Sim/Exp")
    # # draw a horizontal line y=0
    # axs[1].axhline(y=0, color='black', linestyle='--', linewidth=1)
    
    # axs[1].fill_between(real_bins[:-1] + 0.5, -0.1, 0.1, color='green', alpha=0.2, label="10%")
    
    # axs[1].fill_between(real_bins[:-1] + 0.5, 0.1, 0.2, color='yellow', alpha=0.2, label="10-20%")
    # axs[1].fill_between(real_bins[:-1] + 0.5, -0.1, -0.2, color='yellow', alpha=0.2)
    
    # axs[1].fill_between(real_bins[:-1] + 0.5, 0.2, 0.3, color='red', alpha=0.2, label="20-30%")
    # axs[1].fill_between(real_bins[:-1] + 0.5, -0.2, -0.3, color='red', alpha=0.2)
    
    
    # axs[1].set_ylim(-1, 1)
    ###############################################################################################################

    # save the figure
    # plt.savefig("../207Bi/207Bi/207Pb_1663keV_decay_f9mm_g0mm_n10mil.png", dpi=300)