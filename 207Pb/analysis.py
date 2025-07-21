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
            
            
        # print(f"Bin content: {bin_content}, Bin uncertainty: {bin_uncertainty}")
        
        
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
            ax.set_ylim(bottom=0.1)
            ax.set_xlim(30, 2000)
            
            ax.axvline(x=481.6935, color='green', linestyle='--', linewidth=1)
            ax.axvline(x=553.8372, color='green', linestyle='--', linewidth=1)
            ax.axvline(x=565.8473, color='green', linestyle='--', linewidth=1)
            ax.axvline(x=975.651, color='green', linestyle='--', linewidth=1)
            ax.axvline(x=1047.795, color='green', linestyle='--', linewidth=1)
            ax.axvline(x=1059.805, color='green', linestyle='--', linewidth=1)
            
        geant_fig.tight_layout()
        
    return np.array(sim_hist_smeared),  np.array(sim_bin_centers),  np.array(sim_bin_edges),  np.array(sim_hist_smeared_uncertainity)

# Main part of your code
if __name__ == "__main__":
    
    # Parse the command-line arguments
    args = parse_args()

    # Use the parsed root_file_path from the arguments, or list files if not provided
    if args.root_file_path is None:
        root_file_dir = "../207Pb/"
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

    sim_hist, sim_bin_centers, sim_bin_edges, sim_hist_uncertainity = geant4_simulation_results(root_file_path, args.fwhm, plot=True)
    
    ###############################################################################################################
    
    save_pic = args.save_pic
    # If save-pic is False, ignore save-path
    if save_pic: 
        # take the file name without the extension
        save_path = "../207Pb/" + root_file_path.split(".root")[0] + ".png"
        print(f"Saving the plot to {save_path}")
        # Save the plot to the provided path
        plt.savefig(save_path, dpi=300)
        
    plt.show()