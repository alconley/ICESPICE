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

# for virtual environment on mac
# source $(brew --prefix root)/bin/thisroot.sh  # for ROOT

def parse_args():
    parser = argparse.ArgumentParser(description="Run analysis on a specified ROOT file")
    parser.add_argument("root_file_path", nargs='?', default=None, help="Path to the ROOT file")
    parser.add_argument("--icespice", action="store_true", default=False, help="If ICESPICE in the simultation")
    parser.add_argument("--fwhm", type=float, default=10, help="FWHM value for Gaussian smearing (default: 10.0)")
    parser.add_argument("--save-pic", action="store_true", default=False, help="Save the plot as an image if this flag is set (default: False)")

    return parser.parse_args()

def list_root_files(directory):
    return glob.glob(os.path.join(directory, "*.root"))

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

def experimental_results(plot=True):
    df_withICESPICE = pl.read_parquet("../207Bi/exp_data/207Bi_ICESPICE_f70mm_g30mm_run_*.parquet")
    df_withoutICESPICE = pl.read_parquet("../207Bi/exp_data/207Bi_noICESPICE_f9mm_g0mm_run_13.parquet")

    # Energy calibration of m=0.5395 and b=2.5229
    df_withICESPICE = df_withICESPICE.with_columns([(pl.col("PIPS1000Energy") * 0.5395 + 2.5229).alias("PIPS1000EnergyCalibrated")])
    df_withoutICESPICE = df_withoutICESPICE.with_columns([(pl.col("PIPS1000Energy") * 0.5395 + 2.5229).alias("PIPS1000EnergyCalibrated")])
    
    exp_bins = 1100
    exp_range = [0, 1100]

    exp_hist_withICESPICE, exp_bin_edges_withICESPICE = np.histogram(df_withICESPICE["PIPS1000EnergyCalibrated"], bins=exp_bins, range=exp_range)
    exp_hist_withoutICESPICE, exp_bin_edges_withoutICESPICE = np.histogram(df_withoutICESPICE["PIPS1000EnergyCalibrated"], bins=exp_bins, range=exp_range)
    exp_hist_withICESPICE_uncertainity = np.sqrt(exp_hist_withICESPICE)
    exp_hist_withoutICESPICE_uncertainity = np.sqrt(exp_hist_withoutICESPICE)
    
    if plot: 
        linewidth = 1
        fs = 10
        
        exp_fig, exp_axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True, num="Experimental Data")
        exp_axs = exp_axs.flatten()
        
        ###############################################################################################################
        # The experimental histogram with and without ICESPICE
        

        
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
            ax.set_xlim(200, 1100)   
            ax.minorticks_on()
            ax.tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=3)
            ax.tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=5)
        
        exp_fig.tight_layout()
    
    return exp_hist_withICESPICE, exp_bin_edges_withICESPICE, exp_hist_withICESPICE_uncertainity, exp_hist_withoutICESPICE, exp_bin_edges_withoutICESPICE, exp_hist_withoutICESPICE_uncertainity

def geant4_simulation_results(root_file_path, fwhm, plot=True):
    # Get the simulation data from the ROOT file and plot the raw data
    sim_hist, sim_bin_centers, sim_bin_edges, sim_hist_uncertainity = get_root_hist_data(root_file_path, "Esil", print=False)
    
    if plot:
        geant_fig, geant_axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True, num="Geant4 Simulation")
        
        geant_axs[0].stairs(values=sim_hist, edges=sim_bin_edges, color="dodgerblue", label="Geant4 simulation (Raw)", linewidth=linewidth)

        n_sim_particles = np.sum(sim_hist)
        # in scientific notation
        geant_axs[0].text(0.05, 0.95, f"Total Counts: {n_sim_particles:.1e}", transform=geant_axs[0].transAxes, ha='left', va='top')
        
        # Smear the simulation data with a Gaussian and plot the smeared data
        sim_hist_smeared, sim_hist_smeared_low, sim_hist_smeared_high = gaussian_smear(bin_contents=sim_hist, 
                                                                                    bin_centers=sim_bin_centers, 
                                                                                    bin_uncertainity=sim_hist_uncertainity, 
                                                                                    fwhm=fwhm)
        
        sim_hist_smeared_uncertainity = (sim_hist_smeared_high - sim_hist_smeared_low) / 2
        
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
        
    return np.array(sim_hist_smeared),  np.array(sim_bin_centers),  np.array(sim_bin_edges),  np.array(sim_hist_smeared_uncertainity)

def arctan_threshold(x, width, scale, phase):
    """
    Arctan threshold function that sets negative values to 0, mimicking a threshold.
    
    The 'start' value is where the transition begins, and 'end' is where it ends.

    Parameters:
    - x: Input data (energy values).
    - start: The energy value where the threshold starts (below this, output is 0).
    - end: The energy value where the threshold ends (above this, output is 1).
    - width: Controls the smoothness of the threshold.
    - scale: Scaling factor for fitting the experimental counts.

    Returns:
    - Threshold values applied to x, with negative values set to 0.
    """
    # Compute arctan values
    arctan_values = scale * (np.arctan(x / width - phase*np.pi) / np.pi)
    
    # set negative values to 0
    arctan_values[arctan_values < 0] = 0
    
    return arctan_values

def threshold_residuals(params, x, data):
    """
    Residuals function for fitting the experimental data.

    Parameters:
    - params: lmfit Parameters object containing 'start', 'end', 'width', and 'scale'.
    - x: Experimental bin centers (energy values).
    - data: Experimental bin contents (counts).

    Returns:
    - Residuals: Difference between the experimental data and the model.
    """
    width = params['width'].value
    scale = params['scale'].value
    phase = params['phase'].value

    model = arctan_threshold(x, width, scale, phase)
    return data - model

def scale_residuals(params, exp_hist, exp_unc, sim_hist, sim_unc):
    """
    Residuals function for fitting the simulated data to the experimental data.
    
    Parameters:
    - params: lmfit Parameters object containing 'scale'.
    - exp_hist: Experimental histogram counts.
    - exp_unc: Experimental histogram uncertainties.
    - sim_hist: Simulated histogram counts.
    - sim_unc: Simulated histogram uncertainties.
    
    Returns:
    - Residuals: The difference between the scaled simulated histogram and the experimental histogram, normalized by uncertainties.
    """
    scale = params['scale'].value  # Extract the scaling factor

    # Scale the simulated histogram
    scaled_sim_hist = sim_hist * scale

    # Combine uncertainties (experimental and simulated uncertainties)
    total_unc = np.sqrt(exp_unc ** 2 + sim_unc ** 2)

    # Compute the residuals (normalized by uncertainties)
    residuals = (exp_hist - scaled_sim_hist) / total_unc

    return residuals


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

    
    linewidth = 1
    fs = 10
    
    ###############################################################################################################

    exp_hist_withICESPICE, exp_bin_edges_withICESPICE, exp_hist_withICESPICE_uncertainity, \
    exp_hist_withoutICESPICE, exp_bin_edges_withoutICESPICE, exp_hist_withoutICESPICE_uncertainity = experimental_results(plot=True)
    
    ###############################################################################################################

    sim_hist, sim_bin_centers, sim_bin_edges, sim_hist_uncertainity = geant4_simulation_results(root_file_path, args.fwhm, plot=True)
        
    ###############################################################################################################
    if args.icespice:
        exp_v_geant_fig,  exp_v_geant_axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True, num="Experiment vs Simulation with ICESPICE")
        exp_v_geant_axs = exp_v_geant_axs.flatten()
        
        exp_v_geant_axs[0].stairs(values=exp_hist_withICESPICE, edges=exp_bin_edges_withICESPICE, color="black", label=r"$^{207}$Bi with ICESPICE", linewidth=linewidth)
        exp_bin_centers = (exp_bin_edges_withICESPICE[:-1] + exp_bin_edges_withICESPICE[1:]) / 2
        exp_hist_data = exp_hist_withICESPICE  # Replace with the appropriate experimental data
        exp_hist_uncertainity = exp_hist_withICESPICE_uncertainity
        
        threshold_ranges = [100, 450]
        scale_range = [450, 1100]
    else:
        exp_v_geant_fig,  exp_v_geant_axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True, num="Experiment vs Simulation without ICESPICE")
        exp_v_geant_axs = exp_v_geant_axs.flatten()
        
        exp_v_geant_axs[0].stairs(values=exp_hist_withoutICESPICE, edges=exp_bin_edges_withoutICESPICE, color="black", label=r"$^{207}$Bi without ICESPICE", linewidth=linewidth)
        
        exp_bin_centers = (exp_bin_edges_withoutICESPICE[:-1] + exp_bin_edges_withoutICESPICE[1:]) / 2
        exp_hist_data = exp_hist_withoutICESPICE  # Replace with the appropriate experimental data
        exp_hist_uncertainity = exp_hist_withoutICESPICE_uncertainity
        
        threshold_ranges = [300, 390]
        scale_range = [390, 1080]
            
    # Section for getting the threshold function
    threshold_mask = (exp_bin_centers >= threshold_ranges[0]) & (exp_bin_centers <= threshold_ranges[1])
    threshold_exp_bin_centers = exp_bin_centers[threshold_mask]
    threshold_exp_hist_data = exp_hist_data[threshold_mask]

    threshold_params = lmfit.Parameters()
    threshold_params.add('width', value=90.0)      # Smoothness of the threshold
    threshold_params.add('scale', value=2000.0)  # Scaling factor for the counts
    threshold_params.add('phase', value=1.0)  # Phase shift for the arctan function
    threshold_result = lmfit.minimize(threshold_residuals, threshold_params, args=(threshold_exp_bin_centers, threshold_exp_hist_data))

    print("\nFitted Parameters for Threshold Function (Arctan):")
    lmfit.report_fit(threshold_result.params)
    fitted_width = threshold_result.params['width'].value
    fitted_scale = threshold_result.params['scale'].value
    fitted_phase = threshold_result.params['phase'].value
    
    # exp_v_geant_axs[0].plot(threshold_exp_bin_centers, arctan_threshold(threshold_exp_bin_centers, fitted_width, fitted_scale, fitted_phase), color='green', label="Fitted Threshold Function")
    
    # Section for getting the best scale
    
    exp_scale_mask = (exp_bin_centers >= scale_range[0]) & (exp_bin_centers <= scale_range[1])
    scale_exp_bin_centers = exp_bin_centers[exp_scale_mask]
    scale_exp_hist_data = exp_hist_data[exp_scale_mask]
    scale_exp_hist_uncertainity = exp_hist_uncertainity[exp_scale_mask]
    
    sim_scale_mask = (sim_bin_centers >= scale_range[0]) & (sim_bin_centers <= scale_range[1])
    scale_sim_bin_centers = sim_bin_centers[sim_scale_mask]
    scale_sim_hist = sim_hist[sim_scale_mask]
    scale_sim_hist_uncertainity = sim_hist_uncertainity[sim_scale_mask]

    # Initialize lmfit parameters
    scale_params = lmfit.Parameters()
    scale_params.add('scale', value=1.0)  # Initial scale guess

    # Perform residuals minimization using lmfit
    scale_result = lmfit.minimize(scale_residuals, scale_params, args=(scale_exp_hist_data, scale_exp_hist_uncertainity, scale_sim_hist, scale_sim_hist_uncertainity))

    # Get the best scale factor
    best_scale = scale_result.params['scale'].value
    scale_uncertainty = scale_result.params['scale'].stderr
    
    print("\nFitted Parameters for Scaling:")
    lmfit.report_fit(scale_result.params)

    # Apply the scaling factor to the simulated histogram within the scale range
    scaled_sim_hist_mask = scale_sim_hist * best_scale
    scaled_sim_hist_uncertainty_mask = scaled_sim_hist_mask * np.sqrt(
        (scale_sim_hist_uncertainity / scale_sim_hist) ** 2 + (scale_uncertainty / best_scale) ** 2
    )

    # Scale the raw data and get the uncertainty
    scaled_sim_hist = sim_hist*best_scale
    scaled_sim_hist_uncertainty = sim_hist * np.sqrt(
        (sim_hist_uncertainity / sim_hist) ** 2 + (scale_uncertainty / best_scale) ** 2
    )

    # Apply the threshold function to the scaled_sim_hist for energy values below 390 keV
    below_threshold_mask = sim_bin_centers < threshold_ranges[1]  

    # Apply the arctan threshold function to the scaled_sim_hist for energies below 390 keV
    scaled_sim_hist_below_threshold = scaled_sim_hist[below_threshold_mask]
    scaled_sim_hist_threshold = arctan_threshold(
        sim_bin_centers[below_threshold_mask], fitted_width, fitted_scale, fitted_phase
    )

    # Combine the thresholded data below 390 keV with the original data above 390 keV
    scaled_sim_hist_with_threshold = np.copy(scaled_sim_hist)
    scaled_sim_hist_with_threshold[below_threshold_mask] = scaled_sim_hist_threshold

    # Plot the scaled simulated histogram with threshold
    exp_v_geant_axs[0].stairs(values=scaled_sim_hist_with_threshold, edges=sim_bin_edges, color='dodgerblue', label="Thresholded & Scaled Geant4 Simulation")
    exp_v_geant_axs[0].fill_between(
        scale_sim_bin_centers,
        scaled_sim_hist_mask - scaled_sim_hist_uncertainty_mask,
        scaled_sim_hist_mask + scaled_sim_hist_uncertainty_mask,
        color='dodgerblue',
        alpha=0.2,
        label="Uncertainty"
    )
    
    
    if args.icespice:
        activity = (np.sum(sim_hist)*best_scale)/52200
    else:   
        activity = (np.sum(sim_hist)*best_scale)/4080
        
    exp_v_geant_axs[0].text(0.5, 0.95, f"N: {np.sum(sim_hist):.1e}\nBest Scale: {best_scale:.2f} ± {scale_uncertainty:.2f}\nActivity Estimate: {activity:.0f} Bq", transform=exp_v_geant_axs[0].transAxes, ha='center', va='top')
    
    # Set the y-axis to log scale and adjust limits
    exp_v_geant_axs[0].set_ylim(0, 7000)
    exp_v_geant_axs[0].set_xlabel(r"Energy [keV]")
    exp_v_geant_axs[0].set_ylabel(r"Counts/keV")
        
    exp_minus_sim = scale_exp_hist_data - scaled_sim_hist_mask
    exp_minus_sim_uncertainity = np.sqrt(scale_exp_hist_uncertainity**2 + scaled_sim_hist_uncertainty_mask**2)
    
    percent_diff = 100 * exp_minus_sim / scale_exp_hist_data
    percent_difference_uncertainty = percent_diff * np.sqrt( (exp_minus_sim_uncertainity / exp_minus_sim)**2 + (scale_exp_hist_uncertainity / scale_exp_hist_data)**2)
        
    exp_v_geant_axs[1].plot(scale_exp_bin_centers, percent_diff, marker='o', markersize=1, color="black", linestyle='None', label="Exp-Sim/Exp [%]")
    exp_v_geant_axs[1].fill_between(scale_exp_bin_centers, 
                                    percent_diff - np.abs(percent_difference_uncertainty),
                                    percent_diff + np.abs(percent_difference_uncertainty),
                                    color='red',
                                    alpha=0.2,
                                    label="Uncertainty"
    )
    
    exp_v_geant_axs[1].set_ylim(-100, 100)
    exp_v_geant_axs[1].axhline(y=0, color='black', linestyle='--', linewidth=1)
    exp_v_geant_axs[1].set_xlabel(r"Energy [keV]")
    exp_v_geant_axs[1].set_ylabel(r"Exp-Sim/Exp [%]")
    
    
    exp_v_geant_axs[2].errorbar(scale_exp_bin_centers, exp_minus_sim, yerr=exp_minus_sim_uncertainity, ecolor='red', capsize=1)
    exp_v_geant_axs[2].axhline(y=0, color='black', linestyle='--', linewidth=1)
    exp_v_geant_axs[2].set_xlabel(r"Energy [keV]")
    exp_v_geant_axs[2].set_ylabel(r"Residuals")
                                    
    for ax in exp_v_geant_axs:
        ax.axvline(x=481.6935, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=553.8372, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=565.8473, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=975.651, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=1047.795, color='green', linestyle='--', linewidth=1)
        ax.axvline(x=1059.805, color='green', linestyle='--', linewidth=1)
        
        ax.legend(loc='upper left', shadow=False, frameon=True, fancybox=False, edgecolor='none', facecolor='none')
        ax.set_xlim(270, 1080)   
        ax.minorticks_on()
        ax.tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=3)
        ax.tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=5)

    exp_v_geant_fig.tight_layout()
    
    ###############################################################################################################
    
    save_pic = args.save_pic
    # If save-pic is False, ignore save-path
    if save_pic: 
        # take the file name without the extension
        save_path = "../207Bi/geant_sim/" + root_file_path.split(".root")[0] + ".png"
        print(f"Saving the plot to {save_path}")
        # Save the plot to the provided path
        plt.savefig(save_path, dpi=300)
        
    plt.show()