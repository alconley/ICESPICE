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
        axs[0].errorbar(energies, trans_probs, yerr=trans_errors, marker="o", linestyle="-", color=color, label=f"f={f}mm, g={g} mm")
        axs[0].set_title(r"Transmission Probability")
        axs[0].set_xlabel(r"Energy (keV)")
        axs[0].set_ylabel(r"Probability (%)")

        # Plot Full Energy Transmission Probability with error bars
        axs[1].errorbar(energies, full_energy_trans_probs, yerr=full_energy_errors, marker="x", linestyle="-", color=color, label=f"f={f}mm, g={g} mm")
        axs[1].set_title(r"Full Energy Transmission Probability")
        axs[1].set_xlabel(r"Energy (keV)")
        axs[1].set_ylabel(r"Probability (%)")

        # Plot Efficiency with error bars
        axs[2].errorbar(energies, efficiencies, yerr=efficiency_errors, marker="s", linestyle="-", color=color, label=f"f={f}mm, g={g} mm")
        axs[2].set_title(r"Efficiency (FEP Transmission / Total Transmission)")
        axs[2].set_xlabel(r"Energy (keV)")
        axs[2].set_ylabel(r"Efficiency (%)")

    # Add legend and adjust layout
    for ax in axs:
        ax.legend(loc="upper right")
    fig.tight_layout()
    

def demo_plot():
    import matplotlib
    tex_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
        "font.family": "serif",
        "font.serif" : ["CMR10"],
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 8,
        "font.size": 6,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 6,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6
    }


    def set_size(width, fraction=1, rows=1):
        """Set figure dimensions to avoid scaling in LaTeX.

        Parameters
        ----------
        width: float
                Document textwidth or columnwidth in pts
        fraction: float, optional
                Fraction of the width which you wish the figure to occupy

        Returns
        -------
        fig_dim: tuple
                Dimensions of figure in inches
        """
        # Width of figure (in pts)
        fig_width_pt = width * fraction

        # Convert from pt to inches
        inches_per_pt = 1 / 72.27

        # Golden ratio to set aesthetic figure height
        # https://disq.us/p/2940ij3
        golden_ratio = (5**.5 - 1) / 2

        # Figure width in inches
        fig_width_in = fig_width_pt * inches_per_pt
        # Figure height in inches
        fig_height_in = fig_width_in * golden_ratio

        fig_dim = (fig_width_in, fig_height_in*rows)

        return fig_dim


    plt.rcParams.update(tex_fonts)
    matplotlib.rcParams['axes.unicode_minus'] = False

    fig, ax = plt.subplots(2, 1, figsize=set_size(222, rows=2))
    ax = ax.flatten()
    
    n = 1000000

    # Define the range of energies and values for f and g
    energies = range(100, 2100, 100)  # From 100 keV to 2000 keV in 100 keV steps
    f_values = [70]  # Example values for f
    g_values = [20, 25, 30, 35, 40, 45]  # Example values for g

    # Generate file paths for all combinations of f, g, and energy
    files_by_fg_1000 = {
        (f, g): [
            (f"../TransmissionProbablilty/data/5N42_1x1x1_8in_PIPS1000/TransmissionProb_f{f}mm_g{g}mm_n{n}_PIPS1000_PointSource_{energy}keV.root", energy)
            for energy in energies
        ]
        for f, g in product(f_values, g_values)
    }
    
        # Generate file paths for all combinations of f, g, and energy
    files_by_fg_300 = {
        (f, g): [
            (f"../TransmissionProbablilty/data/5N42_1x1x1_16in_PIPS300/TransmissionProb_f{f}mm_g{g}mm_n{n}_PIPS300_PointSource_{energy}keV.root", energy)
            for energy in energies
        ]
        for f, g in product(f_values, g_values)
    }
    
    colors = ["b", "g", "r", "c", "m", "y"]
    line_styles = [
        (0, (1, 1)),           # Densely dotted
        (0, (5, 5)),           # Dashed
        (0, (1, 10)),          # Dotted
        (0, (3, 5, 1, 5)),     # Dash dot
        (0, (3, 5, 1, 5)),     # Dash dot (repeated)
        (0, (3, 1, 1, 1))      # Densely dash dot
    ]

    for (f, g), files_with_energies in files_by_fg_1000.items():
        # Run batch analysis for each set of files
        energies, trans_probs, trans_errors, full_energy_trans_probs, full_energy_errors, efficiencies, efficiency_errors = batch_transmission_probability(files_with_energies)

        # Choose a color for each (f, g) pair and add labels
        color = colors.pop(0)
        line_style = line_styles.pop(0)  # Select the next line style


        # Plot Transmission Probability with error bars
        ax[0].errorbar(energies, trans_probs, yerr=trans_errors, marker="o", linestyle=line_style, color=color, label=f"g={g} mm", markersize=0.5, linewidth=1)

    # ax[0].legend(loc='upper right', shadow=False, frameon=True, fancybox=False, edgecolor='none', facecolor='none')
    ax[0].minorticks_on()
    ax[0].tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=2)
    ax[0].tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=4)
    ax[0].text(0.05, 0.92, r"5-1''x1''x1/8'' N42 magnets", transform=ax[0].transAxes, verticalalignment='top', horizontalalignment='left')    
    ax[0].set_xticklabels([])

    colors = ["b", "g", "r", "c", "m", "y"]
    line_styles = [
        (0, (1, 1)),           # Densely dotted
        (0, (5, 5)),           # Dashed
        (0, (1, 10)),          # Dotted
        (0, (3, 5, 1, 5)),     # Dash dot
        (0, (3, 5, 1, 5)),     # Dash dot (repeated)
        (0, (3, 1, 1, 1))      # Densely dash dot
    ]

    for (f, g), files_with_energies in files_by_fg_300.items():
        # Run batch analysis for each set of files
        energies, trans_probs, trans_errors, full_energy_trans_probs, full_energy_errors, efficiencies, efficiency_errors = batch_transmission_probability(files_with_energies)

        # Choose a color for each (f, g) pair and add labels
        color = colors.pop(0)
        line_style = line_styles.pop(0)  # Select the next line style

        # Plot Transmission Probability with error bars
        ax[1].errorbar(energies, trans_probs, yerr=trans_errors, marker="o", linestyle=line_style, color=color, label=f"g={g} mm", markersize=0.5, linewidth=1)

    # ax[1].legend(loc='upper right', shadow=False, frameon=True, fancybox=False, edgecolor='none', facecolor='none')
    ax[1].minorticks_on()
    ax[1].tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=2)
    ax[1].tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=4)    
    ax[1].text(0.05, 0.92, r"5-1''x1''x1/16'' N42 magnets", transform=ax[1].transAxes, verticalalignment='top', horizontalalignment='left')
    ax[1].set_xlabel(r"Energy [keV]")

    # Collect all handles and labels from both subplots
    handles, labels = ax[0].get_legend_handles_labels()

    # Create a single legend at the top of the figure
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False)

    # Adjust the layout to make space for the legend
    # fig.tight_layout(rect=[0, 0, 1, 0.9])

    for ax in ax:
        ax.set_xlim(0, 2000)
        ax.set_ylim(0, 0.79)
        ax.set_ylabel(r"T(E) [\%]")

    fig.subplots_adjust(
        top=0.916,
        bottom=0.085,
        left=0.123,
        right=0.963,
        hspace=0.025,
        wspace=0.2
    )
    plt.savefig(f"../TransmissionProbablilty/ICESPICE_demo_transmission_probabilities.pdf")
    plt.show()
    
demo_plot()

# # energy_1000keV, trans_prob_1000keV, fep_trans_prob_1000keV = transmission_histogram("../TransmissionProbablilty/TransmissionProb_f70mm_g30mm_n1000000_PIPS1000_PointSource_1000keV.root", 1000, plot=True)

# # Define fixed parameters and base path
# # base_path = "../TransmissionProbablilty/data/5N42_1x1x1_8in_PIPS1000"
# # detector = "PIPS1000"

# base_path = "../TransmissionProbablilty/data/3N42_1x1x1_16in_PIPS300"
# # base_path = "../TransmissionProbablilty/data/5N42_1x1x1_16in_PIPS300"
# # base_path = "../TransmissionProbablilty/data/6N42_1x1x1_16in_PIPS300"


# detector = "PIPS300"

# n = 1000000

# # Define the range of energies and values for f and g
# energies = range(100, 2100, 100)  # From 100 keV to 2000 keV in 100 keV steps
# f_values = [70]  # Example values for f
# g_values = [20, 25, 30, 35, 40]  # Example values for g

# # Generate file paths for all combinations of f, g, and energy
# files_by_fg = {
#     (f, g): [
#         (f"{base_path}/TransmissionProb_f{f}mm_g{g}mm_n{n}_{detector}_PointSource_{energy}keV.root", energy)
#         for energy in energies
#     ]
#     for f, g in product(f_values, g_values)
# }

# # Call the plotting function with the files_by_fg dictionary
# batch_plot_by_fg(files_by_fg)

# plt.savefig(f"../TransmissionProbablilty/3N42_1x1x1_16in_{detector}_transmission_probabilities.png", dpi=300)
# plt.show()

