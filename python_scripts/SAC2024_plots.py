import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from math import pi, sqrt

tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    "font.serif" : ["CMR10"],
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 10,
    "font.size": 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8
}

plt.rcParams.update(tex_fonts)
matplotlib.rcParams['axes.unicode_minus'] = False

def gaussian(x, mu, sig, area):
    normalization = area / (sig * sqrt(2 * pi))
    return normalization * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

# Implement the folding function
def fold(counts, bin_edges, fwhm_keV):
    fwhm_mev = fwhm_keV / 1000  # Convert keV to MeV
    sigma = fwhm_mev / 2.355

    folded_counts = np.zeros_like(bin_edges[:-1])
    
    for i, count in enumerate(counts):
        bin_center = (bin_edges[i] + bin_edges[i + 1]) / 2
        folded_counts += gaussian(bin_edges[:-1], bin_center, sigma, count * (bin_edges[1] - bin_edges[0]))

    return folded_counts

def transmission_histogram_plot(file_path: str, ax: plt.Axes):

    parts = file_path.split('_')
    # Extract relevant details based on position
    # Example filename: 'ICESPICE_PIPS1000_f50mm_g20mm_1000_h1_Esil.csv'
    detector = parts[1]  # e.g., 'PIPS1000'
    f = parts[2]   # e.g., 'f50mm'
    g = parts[3]            # e.g., 'g20mm'
    simulation_energy = int(parts[4])   # e.g., 1000

        # Open the file to read the bin information
    with open(file_path, 'r') as file:
        lines = file.readlines()
        bin_info = lines[3]  # This is the line with bin axis information

    # Extracting the number of bins and the range from the bin_info line
    _, _, bins, start, end = bin_info.split()
    bins = int(bins) + 3  # Adjusted bins for plt.step
    start = float(start)
    end = float(end)

    # Load the histogram data
    df = pd.read_csv(file_path, skiprows=7, names=['counts', 'Sw', 'Sw2', 'Sxw0', 'Sx2w0'])

    # Extract counts
    counts = df['counts']

    total_counts = counts.sum()
    counts_4pi = total_counts * 2

    transmission_counts = counts[2:].sum()
    transmission_probability = transmission_counts / counts_4pi * 100

    full_energy_depoisted_counts = counts[simulation_energy]
    full_energy_transmission_probability = full_energy_depoisted_counts / counts_4pi * 100
    full_energy_efficiency = full_energy_depoisted_counts / transmission_counts * 100

    textstr = fr'''
    Geant4 Simulation of ICESPICE

    Particle: Electron
        Energy: {simulation_energy} keV
        Number (4$\pi$):  {total_counts*2}
        Distance to ICESPICE [f]: 70 mm

    Detector: Passivated Implanted Planar Silicon (PIPS)
        Thickness: 1000 $\mu$m
        Active area: 50 mm$^2$
        Distance from ICESPICE [g]: 25 mm
        Transmission Counts: {transmission_counts}
        Transmission Probability: {transmission_probability:.2f}\%
        Full Energy Deposited Counts: {full_energy_depoisted_counts}
        Transmission Probability (Full Energy Deposited): {full_energy_transmission_probability:.2f}\%
        Full Energy Deposited Counts/Transmission Counts = {full_energy_efficiency:.2f}\%
        '''

    # Generate the x-axis for the bins
    # bin_edges = np.linspace(start, end, bins)  # Corrected to 'bins' for edges calculation

    # Create the step plot
    # ax.step(bin_edges[:-1], counts, where='post', linewidth=0.5, color='#782F40')  # 'post' ensures the value at each step is constant until the next edge
    
    # Generate the x-axis for the bins
    bin_edges = np.linspace(start, end, bins + 1)  # Corrected to 'bins + 1' for edges calculation

    # Fold the counts with the Gaussian
    folded_counts = fold(counts, bin_edges, 2.2)

    # sum the folded counts
    total_folded_counts = folded_counts.sum()

    # print
    print(f'Total folded counts: {total_folded_counts}\n\n\n')


    # Create the step plot
    ax.step(bin_edges[:-1], folded_counts, where='mid', linewidth=0.5, color='#782F40')  # 'mid' ensures the value at each step is centered

    ax.set_xlabel(r'e$^{-}$ Energy [MeV]')
    ax.set_ylabel(r'Counts/keV')
    ax.set_yscale('log')  # Logarithmic scale for the y-axis
    ax.set_ylim(0.1, 1e3)
    ax.set_xlim(0, 1.1)
    ax.text(0.02, 0.99, textstr, transform=ax.transAxes, verticalalignment='top', horizontalalignment='left')

    return

def transmission_histogram(file_path, plot=True):

    parts = file_path.split('_')
    # Extract relevant details based on position
    # Example filename: 'ICESPICE_PIPS1000_f50mm_g20mm_1000_h1_Esil.csv'
    detector = parts[1]  # e.g., 'PIPS1000'
    f = parts[2]   # e.g., 'f50mm'
    g = parts[3]            # e.g., 'g20mm'
    simulation_energy = int(parts[4])   # e.g., 1000

        # Open the file to read the bin information
    with open(file_path, 'r') as file:
        lines = file.readlines()
        bin_info = lines[3]  # This is the line with bin axis information

    # Extracting the number of bins and the range from the bin_info line
    _, _, bins, start, end = bin_info.split()
    bins = int(bins) + 3  # Adjusted bins for plt.step
    start = float(start)
    end = float(end)

    # Load the histogram data
    df = pd.read_csv(file_path, skiprows=7, names=['counts', 'Sw', 'Sw2', 'Sxw0', 'Sx2w0'])

    # Extract counts
    counts = df['counts']

    total_counts = counts.sum()
    counts_4pi = total_counts * 2

    transmission_counts = counts[2:].sum()
    transmission_probability = transmission_counts / counts_4pi * 100

    full_energy_depoisted_counts = counts[simulation_energy]
    full_energy_transmission_probability = full_energy_depoisted_counts / counts_4pi * 100

    full_energy_efficiency = full_energy_depoisted_counts / transmission_counts * 100

        # Text annotation in the plot
    textstr = f'''File: {file_path}
        Total counts: {total_counts}
        4π Counts: {counts_4pi}
        Transmission counts: {transmission_counts}
        4π transmission probability: {transmission_probability:.2f}%
        Full energy deposited counts: {full_energy_depoisted_counts}
        4π Probability for full energy deposition: {full_energy_transmission_probability:.2f}%
        Full energy efficiency: {full_energy_efficiency:.2f}%'''
    
    # print(textstr)

    if plot:
        # Generate the x-axis for the bins
        bin_edges = np.linspace(start, end, bins)  # Corrected to 'bins' for edges calculation

        # Create the step plot
        plt.figure(figsize=(10, 6))
        plt.step(bin_edges[:-1], counts, where='post', linewidth=0.5)  # 'post' ensures the value at each step is constant until the next edge
        plt.xlabel('Energy (MeV)')
        plt.ylabel('Entries')
        plt.yscale('log')  # Logarithmic scale for the y-axis

        plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, fontsize=8, verticalalignment='top', horizontalalignment='right')

        plt.show()

    return simulation_energy, transmission_probability, full_energy_transmission_probability, full_energy_efficiency

def get_file_paths(detector, f, g):
    import os
    files = []
    for file in os.listdir('./analysis/data/ICESPICEwithf70mm/'):
        if file.endswith('.csv'):
            if file.startswith(f'ICESPICE_{detector}_f{f}mm_g{g}mm'):
                files.append(f'./analysis/data/ICESPICEwithf70mm/{file}')
    return files

def transmission_probability(files, title=None, plot=True):
    energy = []
    transmission_prob = []
    full_energy_transmission_prob = []
    efficiency = []

    for file in files:
        simulation_energy, transmission_probability, full_energy_transmission_probability, full_energy_efficiency = transmission_histogram(file, plot=False)
        energy.append(simulation_energy)
        transmission_prob.append(transmission_probability)
        full_energy_transmission_prob.append(full_energy_transmission_probability)
        efficiency.append(full_energy_efficiency)

    # Sort data by energy
    zipped_data = zip(energy, transmission_prob, full_energy_transmission_prob, efficiency)
    sorted_data = sorted(zipped_data, key=lambda x: x[0])  # Sort by the first element in each tuple, which is energy
    energy, transmission_prob, full_energy_transmission_prob, efficiency = zip(*sorted_data)  # Unzip back into separate lists


    if plot:
        fig, axs = plt.subplots(1, 3, figsize=(12, 5))

        if title:
            fig.suptitle(title)

        axs[0].plot(energy, transmission_prob, marker='.', linestyle='-')
        axs[0].set_title("Transmission Probability")
        axs[0].set_xlabel('Energy (keV)')
        axs[0].set_ylabel('Probability (%)')

        axs[1].plot(energy, full_energy_transmission_prob, marker='.', linestyle='-')
        axs[1].set_title("Full Energy Deposited Transmission Probability")
        axs[1].set_xlabel('Energy (keV)')
        axs[1].set_ylabel('Probability (%)')

        axs[2].plot(energy, efficiency, marker='.', linestyle='-')
        axs[2].set_title("Efficiency")
        axs[2].set_xlabel('Energy (keV)')
        axs[2].set_ylabel('Efficiency (%)')

        fig.tight_layout()

        plt.show()

    return energy, transmission_prob, full_energy_transmission_prob, efficiency

def plot_transmission_summary(detector, f, g_values, ax: plt.Axes):

    data = {}
    for g in g_values:
        files = get_file_paths(detector, f, g)
        energy, transmission_prob, full_energy_transmission_prob, efficiency = transmission_probability(files, title=None, plot=False)
        
        data[f'g{g}mm'] = {
            'energy': energy,
            'transmission_prob': transmission_prob,
            'full_energy_transmission': full_energy_transmission_prob,
            'efficiency': efficiency
        }

    marker_shapes = ['o', 's', 'D', '^', 'v', '<', '>', 'p', 'P', '*', 'X', 'd', 'H', 'h', '+', '|', '_', '.', ',']

    marker_index = 0
    for key, values in data.items():
        # rename the key to remove the 'mm' from the label
        key = key.replace('mm', ' mm')
        key = key.replace('g', 'g = ')
        axs[0].plot(values['energy'], values['transmission_prob'], marker=marker_shapes[marker_index], markersize = 4, linestyle='-', label=key, linewidth=0.5)
        marker_index += 1


    ax.set_ylabel(r'Transmission Probability [\%]')
    ax.set_xlabel(r'e$^{-}$ Energy [keV]')

    ax.legend(loc='upper left', frameon=False)

fig, axs = plt.subplots(1, 2, figsize=(10, 5))
fig.subplots_adjust(wspace=0.12,
                    hspace=0.1,
                    left=0.051,
                    bottom=0.082,
                    right=0.996,
                    top=0.98,)

transmission_histogram_plot('./analysis/data/ICESPICEwithf70mm/ICESPICE_PIPS1000_f70mm_g25mm_1000_h1_Esil.csv', axs[1])
plot_transmission_summary('PIPS1000', '70', [20, 25, 30, 35, 40, 45], axs[0])

for i in range(len(axs)):

    axs[i].minorticks_on()
    axs[i].tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=2)
    axs[i].tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=4)

plt.savefig('./analysis/plots/ICESPICE_transmission_SAC2024.png', dpi=900, bbox_inches='tight')

plt.show()
