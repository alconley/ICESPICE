import ROOT
import numpy as np
import matplotlib.pyplot as plt
import warnings
import os
from matplotlib.cm import get_cmap

class Geant4Analyzer:
    def __init__(self, file: str, histogram_name: str):
        self.file = file
        self.histogram_name = histogram_name
        self.root_file = ROOT.TFile(file, "READ")
        self.histogram = self.root_file.Get(histogram_name)
        
        if not self.histogram:
            raise ValueError(f"Histogram '{histogram_name}' not found in file '{file}'.")
        
        # Extract and store histogram data during initialization
        self.bin_content, self.bin_centers, self.bin_edges, self.bin_uncertainties = self._extract_histogram_data()
    
    def _extract_histogram_data(self):
        bin_content, bin_centers, bin_edges, bin_uncertainties = [], [], [], []

        for i in range(1, self.histogram.GetNbinsX() + 1):  # ROOT bins are 1-indexed
            content = self.histogram.GetBinContent(i)
            center = self.histogram.GetBinCenter(i) * 1000  # Convert to keV
            width = self.histogram.GetBinWidth(i) * 1000  # Convert to keV
            uncertainty = self.histogram.GetBinError(i)
            edge = center - 0.5 * width
            
            bin_content.append(content)
            bin_centers.append(center)
            bin_uncertainties.append(uncertainty)
            bin_edges.append(edge)
        
        # Add the last bin edge
        last_bin_edge = bin_centers[-1] + 0.5 * self.histogram.GetBinWidth(self.histogram.GetNbinsX()) * 1000
        bin_edges.append(last_bin_edge)
        
        return (
            np.array(bin_content),
            np.array(bin_centers),
            np.array(bin_edges),
            np.array(bin_uncertainties),
        )
    
    def gaussian_smear(self, fwhm, plot=True):
        bin_contents = self.bin_content
        bin_centers = self.bin_centers
        bin_uncertainity = self.bin_uncertainties
        
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        smeared_bin_contents = np.zeros_like(bin_contents)
        smeared_bin_contents_low = np.zeros_like(bin_contents)
        smeared_bin_contents_high = np.zeros_like(bin_contents)
        
        bin_centers = np.array(bin_centers)
        total_counts = np.sum(bin_contents)
        
        for i in range(len(bin_centers)):
            bin_center = bin_centers[i]
            
            bin_content = bin_contents[i]
            bin_content_low = bin_contents[i] - bin_uncertainity[i]
            bin_content_high = bin_contents[i] + bin_uncertainity[i]
            
            weights = np.exp(-0.5 * ((bin_centers - bin_center) / sigma) ** 2)
            weights /= np.sum(weights)
            
            smeared_bin_contents += bin_content * weights
            smeared_bin_contents_low += bin_content_low * weights
            smeared_bin_contents_high += bin_content_high * weights
            
        smeared_total_counts = np.sum(smeared_bin_contents)
        
        if total_counts != round(smeared_total_counts):
            warnings.warn(f"Total counts before and after smearing do not match: {total_counts} vs {smeared_total_counts}")
            
        if plot:
            geant_fig, geant_axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True, num="Geant4 Simulation")
            
            geant_axs[0].stairs(values=self.bin_content, edges=self.bin_edges, color="dodgerblue", label="Geant4 simulation (Raw)", linewidth=0.5)

            n_sim_particles = np.sum(self.bin_content)
            
            geant_axs[0].text(0.05, 0.95, f"Total Counts: {n_sim_particles:.1e}", transform=geant_axs[0].transAxes, ha='left', va='top')            
            geant_axs[1].stairs(values=smeared_bin_contents, edges=self.bin_edges, color="dodgerblue", label=f"Geant4 simulation (Smeared: {fwhm} keV)", linewidth=0.5)
            geant_axs[1].fill_between(self.bin_centers, smeared_bin_contents_low, smeared_bin_contents_high, color='dodgerblue', alpha=0.2, label="Uncertainty")
            geant_axs[1].set_ylim(top=smeared_bin_contents.max() * 1.1)
            
            for ax in geant_axs:
                ax.legend(loc='upper right', shadow=False, frameon=True, fancybox=False, edgecolor='none', facecolor='none')
                ax.minorticks_on()
                ax.tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=3)
                ax.tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=5)
                ax.set_xlabel(r"Energy [keV]")
                ax.set_ylabel(r"Counts/keV")
                ax.set_yscale("log")
                ax.set_ylim(bottom=1)                
            geant_fig.tight_layout()
            
            plt.show()
        
        return smeared_bin_contents, smeared_bin_contents_low, smeared_bin_contents_high

# Analyze and plot results for all `PIPS1000` files
def spectrums(files, fwhm, ax: plt.Axes):
    
    cmap = get_cmap("tab20")  # Choose a colormap (e.g., "tab20", "viridis", "plasma", etc.)
    num_files = len(files)
    
    for idx, file in enumerate(files):
        # Extract the energy from the file name
        base_name = os.path.basename(file)
        energy_label = base_name.split('_')[1]  # Assumes format "PIPS1000_100keV_SpectrumTemplete..."
        
        # Create analyzer for each file
        analyzer = Geant4Analyzer(file, "Esil")
        
        edges = analyzer.bin_edges
        smeared_content, _, _ = analyzer.gaussian_smear(fwhm, plot=False)
        
        # Extract energy bin centers and plot
        ax.stairs(
            values=smeared_content, 
            edges=edges, 
            label=f"{energy_label} keV", 
            color=cmap(idx / num_files)  # Assign a unique color based on file index
        )
        
        

pips1000_files = []
pips500_files = []
pips300_files = []
pips100_files = []

base_path = "./PIPSSpectrumTempletes/data"

detectors = ["PIPS1000", "PIPS500", "PIPS300", "PIPS100"]
    
for energy in range(100, 1600, 100):
    for detector in detectors:
        file = os.path.join(base_path, f"{detector}_{energy}keV_SpectrumTemplete_f0mm_g10mm_n100000_PointSource_Zdirection.root")
        if not os.path.exists(file):
            continue
        if detector == "PIPS1000":
            pips1000_files.append(file)
        elif detector == "PIPS500":
            pips500_files.append(file)
        elif detector == "PIPS300":
            pips300_files.append(file)
        elif detector == "PIPS100":
            pips100_files.append(file)

# Plot all results with Gaussian smearing
# spectrums(pips1000_files, fwhm=10)


fig, axs = plt.subplots(4, 1, figsize=(10, 12), num="Geant4 Simulation")
axs = axs.flatten()
spectrums(pips1000_files, fwhm=10, ax=axs[0])
spectrums(pips500_files, fwhm=10, ax=axs[1])
spectrums(pips300_files, fwhm=10, ax=axs[2])
spectrums(pips100_files, fwhm=10, ax=axs[3])

axs[0].text(0.5, 0.95, f"PIPS1000", transform=axs[0].transAxes, ha='center', va='top')
axs[1].text(0.5, 0.95, f"PIPS500", transform=axs[1].transAxes, ha='center', va='top')
axs[2].text(0.5, 0.95, f"PIPS300", transform=axs[2].transAxes, ha='center', va='top')
axs[3].text(0.5, 0.95, f"PIPS100", transform=axs[3].transAxes, ha='center', va='top')


# Collect handles and labels from all axes for a global legend
handles, labels = [], []
for ax in axs:
    for handle, label in zip(*ax.get_legend_handles_labels()):
        if label not in labels:  # Avoid duplicate labels
            handles.append(handle)
            labels.append(label)

# Add a single legend at the top
fig.legend(handles, labels, loc='upper center', ncol=6)

for ax in axs:
    ax.set_ylabel("Counts/keV")
    ax.set_yscale("log")
    ax.set_ylim(1, 1e5)
    ax.set_xlim(0, 1550)
    
axs[3].set_xlabel("Energy (keV)")

fig.subplots_adjust(top=0.916,
bottom=0.049,
left=0.07,
right=0.985,
hspace=0.19,
wspace=0.2)

plt.savefig("PIPSSpectrumTempletes/spectrum_templetes.png", dpi=300)

plt.show()
