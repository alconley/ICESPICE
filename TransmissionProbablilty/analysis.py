import uproot
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

import uproot
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

def transmission_histogram(file_path, simulation_energy, plot=True):
    # Open the ROOT file and retrieve histogram data
    with uproot.open(file_path) as file:
        histogram = file["Esil"]  # Assuming the histogram is named "Esil"
        counts, edges = histogram.to_numpy()

    # Convert edges to keV
    edges_keV = edges * 1000

    # Calculate total counts and counts above 0 bin
    total_counts = counts.sum()
    transmission_counts = counts[1:].sum()  # Counts above 0 bin
    transmission_probability = (transmission_counts / total_counts) * 100 if total_counts > 0 else 0
    transmission_error = (np.sqrt(transmission_counts) / total_counts) * 100 if total_counts > 0 else 0

    # Define the 10 keV range around simulation_energy
    lower_bound = simulation_energy - 10
    upper_bound = simulation_energy + 10

    # Find bins within the 10 keV range
    in_range = (edges_keV[:-1] >= lower_bound) & (edges_keV[1:] <= upper_bound)
    full_energy_deposited_counts = counts[in_range].sum()
    full_energy_transmission_probability = (full_energy_deposited_counts / total_counts) * 100 if total_counts > 0 else 0
    full_energy_error = (np.sqrt(full_energy_deposited_counts) / total_counts) * 100 if total_counts > 0 else 0

    return (
        simulation_energy, transmission_probability, transmission_error,
        full_energy_transmission_probability, full_energy_error
    )

def batch_transmission_probability(files_with_energies):
    energies, trans_probs, trans_errors, full_energy_trans_probs, full_energy_errors, efficiencies, efficiency_errors = [], [], [], [], [], [], []

    for file_path, simulation_energy in files_with_energies:
        sim_energy, transmission_probability, transmission_error, full_energy_transmission_probability, full_energy_error = transmission_histogram(
            file_path, simulation_energy, plot=False
        )
        energies.append(sim_energy)
        trans_probs.append(transmission_probability)
        trans_errors.append(transmission_error)
        full_energy_trans_probs.append(full_energy_transmission_probability)
        full_energy_errors.append(full_energy_error)
        
        # Calculate efficiency (FEP transmission / total transmission probability)
        if transmission_probability > 0:
            efficiency = (full_energy_transmission_probability / transmission_probability) * 100
            efficiency_error = efficiency * np.sqrt(
                (transmission_error / transmission_probability)**2 + (full_energy_error / full_energy_transmission_probability)**2
            )
        else:
            efficiency = 0
            efficiency_error = 0
        
        efficiencies.append(efficiency)
        efficiency_errors.append(efficiency_error)

    # Sort by energy for plotting
    sorted_data = sorted(zip(energies, trans_probs, trans_errors, full_energy_trans_probs, full_energy_errors, efficiencies, efficiency_errors))
    energies, trans_probs, trans_errors, full_energy_trans_probs, full_energy_errors, efficiencies, efficiency_errors = zip(*sorted_data)

    return energies, trans_probs, trans_errors, full_energy_trans_probs, full_energy_errors, efficiencies, efficiency_errors

def batch_plot_by_fg(files_by_fg):
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    colors = ["b", "g", "r", "c", "m", "y"]

    for (f, g), files_with_energies in files_by_fg.items():
        # Run batch analysis for each set of files
        energies, trans_probs, trans_errors, full_energy_trans_probs, full_energy_errors, efficiencies, efficiency_errors = batch_transmission_probability(files_with_energies)

        # Choose a color for each (f, g) pair and add labels
        color = colors.pop(0)

        # Plot Transmission Probability with error bars
        axs[0].errorbar(energies, trans_probs, yerr=trans_errors, marker="o", linestyle="-", color=color, label=f"f={f}mm, g={g}mm")
        axs[0].set_title("Transmission Probability")
        axs[0].set_xlabel("Energy (keV)")
        axs[0].set_ylabel("Probability (%)")

        # Plot Full Energy Transmission Probability with error bars
        axs[1].errorbar(energies, full_energy_trans_probs, yerr=full_energy_errors, marker="x", linestyle="-", color=color, label=f"f={f}mm, g={g}mm")
        axs[1].set_title("Full Energy Transmission Probability")
        axs[1].set_xlabel("Energy (keV)")
        axs[1].set_ylabel("Probability (%)")

        # Plot Efficiency with error bars
        axs[2].errorbar(energies, efficiencies, yerr=efficiency_errors, marker="s", linestyle="-", color=color, label=f"f={f}mm, g={g}mm")
        axs[2].set_title("Efficiency (FEP Transmission / Total Transmission)")
        axs[2].set_xlabel("Energy (keV)")
        axs[2].set_ylabel("Efficiency (%)")

    # Add legend and adjust layout
    for ax in axs:
        ax.legend(loc="upper right")
    fig.tight_layout()
    plt.show()

# energy_1000keV, trans_prob_1000keV, fep_trans_prob_1000keV = transmission_histogram("../TransmissionProbablilty/TransmissionProb_f70mm_g30mm_n1000000_PIPS1000_PointSource_1000keV.root", 1000, plot=True)

# Define fixed parameters and base path
base_path = "../TransmissionProbablilty/data"
detector = "PIPS1000"
n = 1000000

# Define the range of energies and values for f and g
energies = range(100, 2100, 100)  # From 100 keV to 2000 keV in 100 keV steps
f_values = [70]  # Example values for f
g_values = [25, 30, 35]  # Example values for g

# Generate file paths for all combinations of f, g, and energy
files_by_fg = {
    (f, g): [
        (f"{base_path}/TransmissionProb_f{f}mm_g{g}mm_n{n}_{detector}_PointSource_{energy}keV.root", energy)
        for energy in energies
    ]
    for f, g in product(f_values, g_values)
}


# Call the plotting function with the files_by_fg dictionary
batch_plot_by_fg(files_by_fg)