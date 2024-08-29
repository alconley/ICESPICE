import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def transmission_histogram(file_path, plot=True):

    parts = file_path.split('_')
    print(parts)
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
    
    print(textstr)

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

def get_file_paths(detector, f, g):
    import os
    files = []
    for file in os.listdir('./analysis/data/ICESPICEwithf70mm/'):
        if file.endswith('.csv'):
            if file.startswith(f'ICESPICE_{detector}_f{f}mm_g{g}mm'):
                files.append(f'./analysis/data/ICESPICEwithf70mm/{file}')
    return files

def plot_transmission_summary(detector, f, g_values):

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

    fig, axs = plt.subplots(1, 3, figsize=(12, 5))
    axs = axs.flatten()

    fig.suptitle(f'{detector} | f={f}mm')

    for key, values in data.items():
        axs[0].plot(values['energy'], values['transmission_prob'], marker='.', linestyle='-', label=key)
        axs[1].plot(values['energy'], values['full_energy_transmission'], marker='.', linestyle='-', label=key)
        axs[2].plot(values['energy'], values['efficiency'], marker='.', linestyle='-', label=key)

    axs[0].set_title("4π Transmission Probability")
    axs[0].set_ylabel('Probability (%)')

    axs[1].set_title("Full Energy Deposited Transmission Probability")
    axs[1].set_ylabel('Probability (%)')

    axs[2].set_title("Efficiency")
    axs[2].set_ylabel('Efficiency (%)')

    # Add legends
    for ax in axs:
        ax.set_xlabel('Energy (keV)')

    # Add a single legend
    handles, labels = axs[0].get_legend_handles_labels() 
    fig.legend(handles, labels, loc='lower center', ncol=len(g_values), frameon=False) 

    plt.tight_layout() 
    fig.subplots_adjust(bottom=0.17)  # Adjust the figure to make room for the legend

    plt.savefig(f'./analysis/plots/{detector}_f{f}mm_summary.png', dpi=900)

    # plt.show()
    
# Usage examples:

# plot the histogram for a single file
# transmission_histogram('./analysis/data/ICESPICEwithf70mm/ICESPICE_PIPS1000_f70mm_g30mm_1000_h1_Esil.csv')

# # plot the transmission probability for a detector, f value, and g value
# pips1000_f50_g30_files = get_file_paths(detector='PIPS1000', f='50', g='30')
# transmission_probability(pips1000_f50_g30_files, title='PIPS1000 | f=50mm | g=30mm')

pips1000_f70_g35_files = get_file_paths(detector='PIPS1000', f='70', g='35')
transmission_probability(pips1000_f70_g35_files, title='PIPS1000 | f=70mm | g=35mm')

# plot_transmission_summary(detector='PIPS1000', f='70', g_values=[20, 25, 30, 35, 40, 45])

# plot_transmission_summary(detector='PIPS1000', f='50', g_values=[20, 25, 30, 35, 40, 45, 50])
# plot_transmission_summary(detector='PIPS500', f='50', g_values=[20, 25, 30, 35, 40, 45, 50])
# plot_transmission_summary(detector='PIPS300', f='50', g_values=[20, 25, 30, 35, 40, 45, 50])
# plot_transmission_summary(detector='PIPS100', f='50', g_values=[20, 25, 30, 35, 40, 45, 50])

plt.show()
