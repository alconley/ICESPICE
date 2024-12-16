import uproot
import numpy as np
import matplotlib.pyplot as plt
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
    "legend.fontsize": 5,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6
}


def set_size(width, fraction=1):
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

    fig_dim = (fig_width_in, fig_height_in)

    return fig_dim


plt.rcParams.update(tex_fonts)
matplotlib.rcParams['axes.unicode_minus'] = False


def histogram(file_path, simulation_energy):
    # Open the ROOT file and retrieve histogram data
    with uproot.open(file_path) as file:
        histogram = file["Esil"]  # Assuming the histogram is named "Esil"
        counts, edges = histogram.to_numpy()

    # Convert edges to keV
    edges_keV = edges * 1000

    # Calculate total counts and counts above 0 bin
    total_counts = counts.sum()

    # Define the 10 keV range around simulation_energy
    lower_bound = simulation_energy - 1
    upper_bound = simulation_energy + 1

    # Find bins within the 10 keV range
    in_range = (edges_keV[:-1] >= lower_bound) & (edges_keV[1:] <= upper_bound)
    full_energy_deposited_counts = counts[in_range].sum()

    full_energy_probability = (full_energy_deposited_counts / total_counts) * 100 if total_counts > 0 else 0
    backscattering_probability = (total_counts - full_energy_deposited_counts) / total_counts * 100 if total_counts > 0 else 0

    return counts, edges_keV, full_energy_probability, backscattering_probability

pips1000_counts, pips1000_edges, pips1000_fep, pips1000_bs = histogram("../BackScattering/BackScattering_f0mm_g10mm_n100000_PIPS1000_PointSource_Zdirection_1000keV.root", 1000)
pips500_counts, pips500_edges, pips500_fep, pips500_bs = histogram("../BackScattering/BackScattering_f0mm_g10mm_n100000_PIPS500_PointSource_Zdirection_1000keV.root", 1000)
pips300_counts, pips300_edges, pips300_fep, pips300_bs = histogram("../BackScattering/BackScattering_f0mm_g10mm_n100000_PIPS300_PointSource_Zdirection_1000keV.root", 1000)
pips100_counts, pips100_edges, pips100_fep, pips100_bs = histogram("../BackScattering/BackScattering_f0mm_g10mm_n100000_PIPS100_PointSource_Zdirection_1000keV.root", 1000)

fig, ax = plt.subplots(1,1,figsize=set_size(222))
fig.subplots_adjust(top=0.99,
    bottom=0.17,
    left=0.126,
    right=0.989,
    hspace=0.2,
    wspace=0.2)

ax.step(pips1000_edges, np.append(pips1000_counts, 0), where="post", label=r"1000 $\mu$m", linewidth=0.5, color='#782F40')
ax.step(pips500_edges, np.append(pips500_counts, 0), where="post", label=r"500 $\mu$m", linewidth=0.5, color='#CEB888')
ax.step(pips300_edges, np.append(pips300_counts, 0), where="post", label=r"300 $\mu$m", linewidth=0.5, color='#5CB8B2')
ax.step(pips100_edges, np.append(pips100_counts, 0), where="post", label=r"100 $\mu$m", linewidth=0.5, color='#425563')


# Define the range of x and y
x_start, x_end = 10, 995
y_value = 6000

# Plot the horizontal line with caps
ax.plot([x_start, x_end], [y_value, y_value], color='black', linewidth=0.5, linestyle='--', label=r'Backscattered e$^{-}$')
ax.plot(x_start, y_value, marker='|', color='black', markersize=8)  # Start cap
ax.plot(x_end, y_value, marker='|', color='black', markersize=8)    # End cap


# Sample data for the table
data = [
    [r"1000 $\mu$m", f"{pips1000_bs:.2f} \%", f"{pips1000_fep:.2f} \%"],
    [r"500 $\mu$m", f"{pips500_bs:.2f} \%", f"{pips500_fep:.2f} \%"],
    [r"300 $\mu$m", f"{pips300_bs:.2f} \%", f"{pips300_fep:.2f} \%"],
    [r"100 $\mu$m", f"{pips100_bs:.2f} \%", f"{pips100_fep:.2f} \%"]
]
columns = [r"Thickness", r"Backscattered", r"Full energy"]

# Add the table
table = ax.table(
    cellText=data,
    colLabels=columns,
    cellLoc='center',
    loc='upper right',  # Top-right corner
    bbox=[0.38, 0.7, 0.6, 0.25]  # Adjust position and size [x, y, width, height]
)

# Beautify the table
table.auto_set_font_size(False)
table.set_fontsize(5)  # Smaller font size
table.auto_set_column_width([0, 1, 2])  # Adjust column widths

# Remove cell lines
for key, cell in table.get_celld().items():
    cell.set_linewidth(0)  # No border for cells

ax.set_xlabel(r"Energy (keV)")
ax.set_ylabel(r"Counts/keV")
ax.set_yscale("log")
ax.set_xlim(0, 1050)
ax.set_ylim(1, 5e5)

ax.legend(loc='upper left', shadow=False, frameon=True, fancybox=False, edgecolor='none', facecolor='none')
ax.minorticks_on()
ax.tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=2)
ax.tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=4)

plt.savefig("../BackScattering/PIPS_backscattering.pdf")

plt.show()

