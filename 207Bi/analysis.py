import polars as pl
import matplotlib.pyplot as plt
import numpy as np
import ROOT
from scipy.optimize import minimize
import argparse  
import os
import glob

# for virtual environment
# source $(brew --prefix root)/bin/thisroot.sh  # for ROOT

# Function to parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Run analysis on a specified ROOT file")
    parser.add_argument("root_file_path", nargs='?', default=None, help="Path to the ROOT file")
    parser.add_argument("--icespice", action="store_true", default=True, help="If ICESPICE in the simultation")
    parser.add_argument("--fwhm", type=float, default=10.0, help="FWHM value for Gaussian smearing (default: 10.0)")
    parser.add_argument("--save-pic", action="store_true", default=False, help="Save the plot as an image if this flag is set (default: False)")
    parser.add_argument("--save-path", type=str, default="picture.png", help="Path to save the plot image (default: picture.png)")

    return parser.parse_args()

# Function to get list of .root files
def list_root_files(directory):
    return glob.glob(os.path.join(directory, "*.root"))

# Function to compute residuals
def residuals(scale, real_hist, sim_hist):
    """
    Calculate the residuals between real and scaled simulation histograms.
    
    Parameters:
    - scale: Scale factor to apply to the simulation histogram.
    - real_hist: The real data histogram.
    - sim_hist: The simulation data histogram.
    
    Returns:
    - The sum of squared residuals.
    """
    # Scale the simulation histogram
    scaled_sim_hist = sim_hist * scale
    
    # Compute the residuals (difference between real and scaled simulation)
    res = real_hist - scaled_sim_hist
    
    # Return the sum of squared residuals (minimize this)
    return np.sum(res**2)

# Function to smear a histogram bin with a Gaussian
def gaussian_smear(bin_centers, bin_contents, fwhm):
    """
    Apply Gaussian smearing to histogram data.
    
    Parameters:
    - bin_centers: List of bin centers.
    - bin_contents: List of bin contents (counts).
    - fwhm: Full width at half maximum (FWHM) of the Gaussian distribution.

    Returns:
    - smeared_contents: Gaussian smeared bin contents.
    """
    
    total_counts = np.sum(bin_contents)  # Total counts before smearing
    
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to sigma
    
    smeared_contents = np.zeros_like(bin_contents)
    
    # Loop over all bins and distribute each bin's content over neighboring bins
    for i in range(len(bin_centers)):
        bin_center = bin_centers[i]
        bin_content = bin_contents[i]
        
        # Generate a Gaussian centered at bin_center
        weights = np.exp(-0.5 * ((bin_centers - bin_center) / sigma) ** 2)
        weights /= np.sum(weights)  # Normalize the weights to sum to 1
        
        smeared_contents += bin_content * weights  # Apply the smearing with normalized weights
    
    smeared_total_counts = np.sum(smeared_contents)  # Total counts after smearing
    
    # print(f"Total counts before smearing: {total_counts}")
    # print(f"Total counts after smearing: {smeared_total_counts}")
    
    return smeared_contents

# Function to extract bin centers and contents from a ROOT histogram and apply Gaussian smearing
def get_root_hist_and_gauss_smear(root_file_path, histogram_name, fwhm):
    """
    Extracts bin centers and bin contents from a ROOT histogram and applies Gaussian smearing.

    Parameters:
    - root_file_path: Path to the ROOT file.
    - histogram_name: Name of the histogram in the ROOT file.
    - sigma: Standard deviation for the Gaussian smearing.

    Returns:
    - bin_centers: List of bin centers.
    - smeared_bin_contents: List of Gaussian smeared bin contents (counts).
    """
    # Open the ROOT file
    root_file = ROOT.TFile(root_file_path, "READ")
    
    # Get the histogram by name
    histogram = root_file.Get(histogram_name)
    
    if not histogram:
        raise ValueError(f"Histogram '{histogram_name}' not found in {root_file_path}")
    
    # Prepare lists to store bin centers and smeared contents
    bin_centers = []
    bin_contents = []
    
    n_sim_particles = 0
    n_interactions = 0
    # Loop over the bins of the histogram
    for i in range(1, histogram.GetNbinsX() + 1):  # ROOT bins are 1-indexed
        bin_center = histogram.GetBinCenter(i) * 1000  # Convert to keV if necessary
        if bin_center == 0.5:
            n_sim_particles = histogram.GetBinContent(i)
            bin_content = 0
        else:
            bin_content = histogram.GetBinContent(i)
            n_sim_particles += bin_content
            n_interactions += bin_content
        
        bin_centers.append(bin_center)
        bin_contents.append(bin_content)

    # Close the ROOT file
    root_file.Close()
    
    # Convert to numpy arrays for easy manipulation
    bin_centers = np.array(bin_centers)
    bin_contents = np.array(bin_contents)
    
    # Apply Gaussian smearing
    smeared_bin_contents = gaussian_smear(np.array(bin_centers), np.array(bin_contents), fwhm)
    
    return bin_centers, smeared_bin_contents, n_sim_particles, n_interactions

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
    df_withICESPICE = df_withICESPICE.with_columns([
        (pl.col("PIPS1000Energy") * 0.5395 + 2.5229).alias("PIPS1000EnergyCalibrated")
    ])

    df_withoutICESPICE = df_withoutICESPICE.with_columns([
        (pl.col("PIPS1000Energy") * 0.5395 + 2.5229).alias("PIPS1000EnergyCalibrated")
    ])

    ###############################################################################################################

    fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    axs = axs.flatten()

    linewidth = 1

    # scale of the hist of the data with out ICESPICE so the counts match the 975 peak using a gaussian fit (from gNat program)
    # data_scale = 51700/65625

    # axs[0].hist(df_withoutICESPICE["PIPS1000EnergyCalibrated"], bins=1000, range=[200, 1200], histtype="step", color='#5CB8B2', weights=[data_scale]*len(df_withoutICESPICE), label="without ICESPICE (scaled)", linewidth=linewidth)
    # axs[0].hist(df_withICESPICE["PIPS1000EnergyCalibrated"], bins=1000, range=[200, 1200], histtype="step", color="#A6192E", label="with ICESPICE", linewidth=linewidth)
    # axs[0].set_xlabel(r"Energy [keV]")
    # axs[0].set_ylabel(r"Counts/keV")

    ###############################################################################################################
    
    if args.icespice:
        real_hist, real_bins = np.histogram(df_withICESPICE["PIPS1000EnergyCalibrated"], bins=799, range=[400, 1199])
    else:
        real_hist, real_bins = np.histogram(df_withoutICESPICE["PIPS1000EnergyCalibrated"], bins=799, range=[400, 1199])

    fwhm = args.fwhm  # Get the FWHM value from the argument

    geant_bin_centers, geant_bin_counts, n_sim_particles, n_interactions = get_root_hist_and_gauss_smear(root_file_path, "Esil", fwhm=fwhm)
    
    # convert bin_centers to left edges
    geant_bin_edges = geant_bin_centers - 0.5 * (geant_bin_centers[1] - geant_bin_centers[0])

    # Filter simulation data to only include energies between 200-1200 keV
    filtered_geant_bin_centers = geant_bin_centers[(geant_bin_centers >= 400) & (geant_bin_centers <= 1200)]
    filtered_geant_bin_counts = geant_bin_counts[(geant_bin_centers >= 400) & (geant_bin_centers <= 1200)]

    # Make sure the binning of the real and sim data matches
    sim_hist, sim_bins = np.histogram(filtered_geant_bin_centers, bins=real_bins, weights=filtered_geant_bin_counts)

    # Minimize the residuals
    initial_scale = 10.0
    result = minimize(residuals, initial_scale, args=(real_hist, sim_hist))
    best_scale = result.x[0]

    # # Apply the best scale to the simulation
    scaled_sim_hist = sim_hist * best_scale

    # plot the histogram of the Geant4 simulation without ICESPICE
    if args.icespice:
        axs[0].hist(df_withICESPICE["PIPS1000EnergyCalibrated"], bins=1000, range=[200, 1200], histtype="step", color='#5CB8B2', label="without ICESPICE", linewidth=linewidth)
    else:
        axs[0].hist(df_withoutICESPICE["PIPS1000EnergyCalibrated"], bins=1000, range=[200, 1200], histtype="step", color='#5CB8B2', label="without ICESPICE", linewidth=linewidth)
    
    axs[0].step(filtered_geant_bin_centers[:-1], scaled_sim_hist, color="#425563", label="Scaled Geant4 simulation", linewidth=1, where="mid")
    axs[0].text(0.50, 0.95, f"Scale factor: {best_scale:.3f}\nSimulation Counts in Detector: {n_interactions}", transform=axs[0].transAxes, ha='center', va='top')
    axs[0].set_xlabel(r"Energy [keV]")
    axs[0].set_ylabel(r"Counts/keV")

# \nSimulation particles: {n_sim_particles:.1e}
    ###############################################################################################################

    # # Residual plot
    axs[1].step(real_bins[:-1], real_hist - scaled_sim_hist, where="mid", color="black", label="Residuals", linewidth=1)
    axs[1].set_xlabel(r"Energy [keV]")
    axs[1].set_ylabel(r"Residuals")

    ###############################################################################################################

    for ax in axs:
        ax.legend(loc='upper left', shadow=False, frameon=True, fancybox=False, edgecolor='none', facecolor='none')
        ax.set_xlim(200, 1200)
            
        ax.minorticks_on()
        ax.tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=3)
        ax.tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=5)

    fig.tight_layout()
    
    # save the figure
    # plt.savefig("../207Bi/207Bi/207Pb_1663keV_decay_f9mm_g0mm_n10mil.png", dpi=300)
    save_pic = args.save_pic
    save_path = args.save_path

    # If save-pic is False, ignore save-path
    if save_pic:
        print(f"Saving the plot to {save_path}")
        # Save the plot to the provided path
        plt.savefig(save_path, dpi=300)
        
    plt.show()
