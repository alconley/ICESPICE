import uproot
import glob
import numpy as np
import matplotlib.pyplot as plt

def open_root_file(file_path, hist_name):
    """Open a ROOT file and return the bins and edges (in keV)."""
    with uproot.open(file_path) as file:
        hist = file[hist_name]
        bins, edges = hist.to_numpy()      # safer / uproot4-style
        edges = edges * 1000.0             # MeV -> keV
    return bins, edges

def plot():

    import matplotlib

    tex_fonts = {
        # Use LaTeX to write all text
        # "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["CMR10"],
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 8,
        "font.size": 6,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 6,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
    }

    plt.rcParams.update(tex_fonts)
    matplotlib.rcParams["text.usetex"] = True

    from matplotlib import cycler

    # Create a long color cycle (20 colors)
    color_list = plt.cm.tab20.colors   # 20 distinct colors
    plt.rcParams['axes.prop_cycle'] = cycler(color=color_list)

    fig, axs = plt.subplots(2, 1, figsize=(2.8919330289193304, 3.5))

    detector_areas = [25, 50, 150, 300, 450, 600, 900, 1200, 200, 500]

    # For Si(Li) detectors, we have BOTH 2 mm and 5 mm for 200 and 500 mm²
    sili_thickness_um = {
        200: [2000, 5000],  # 2 mm and 5 mm
        500: [2000, 5000],  # 2 mm and 5 mm
    }

    f = 70
    g = 30
    n = 1000000

    # Deterministic marker assignment: (group, det) → marker
    # groups: "PIPS", "SiLi2", "SiLi5"
    marker_map = {
        ("PIPS", 25): "o",
        ("PIPS", 50): "s",
        ("PIPS", 150): "D",
        ("PIPS", 300): "^",
        ("PIPS", 450): "v",
        ("PIPS", 600): "<",
        ("PIPS", 900): ">",
        ("PIPS", 1200): "P",
        ("SiLi2", 200): "X",
        ("SiLi2", 500): "*",
        ("SiLi5", 200): "h",
        ("SiLi5", 500): "H",
    }

    # Store plotted line handles for grouped legend: group → {det_area: handle}
    curve_handles = {
        "PIPS": {},
        "SiLi2": {},
        "SiLi5": {},
    }

    for det in detector_areas:
        # For 200 and 500 mm², loop over both 2 mm and 5 mm Si(Li) thicknesses.
        # For all others, we only have PIPS (1.0 mm) and don't filter by thickness.
        thicknesses = sili_thickness_um.get(det, [None])

        for thick in thicknesses:
            if thick is None:
                # PIPS case: your original working pattern
                file_pattern = f"./data/*_{det}mm2_f{f}mm_g{g}mm_n{n}_energy*.root"
            else:
                # Si(Li) case: explicitly select 2000umThick or 5000umThick
                file_pattern = (
                    f"./data/*_{thick}umThick_*_{det}mm2_f{f}mm_g{g}mm_n{n}_energy*.root"
                )

            files = glob.glob(file_pattern)

            energies = []
            total_transmission_prob = []
            fed_transmission_prob = []

            for file in files:
                energy_val = float(file.split("energy")[1].split("keV")[0])

                bins, edges = open_root_file(file, "Esil")

                total_counts = np.sum(bins)
                transmission_counts = np.sum(bins[1:])  # Exclude the first bin
                total_transmission_probability = (
                    transmission_counts / total_counts * 100 if total_counts > 0 else 0
                )

                # find the number of counts at the energy bin ± 1 bin
                energy_bin = np.argmin(np.abs(edges[:-1] - energy_val))
                counts_plus_minus_one = np.sum(
                    bins[max(0, energy_bin - 1):min(len(bins), energy_bin + 2)]
                )

                fed_transmission_probability = (
                    counts_plus_minus_one / total_counts * 100
                    if total_counts > 0
                    else 0
                )

                energies.append(energy_val)
                total_transmission_prob.append(total_transmission_probability)
                fed_transmission_prob.append(fed_transmission_probability)

            # --------------------------
            # SORT BY ENERGY SAFELY
            # --------------------------
            if energies:  # Only if files exist
                combined = sorted(
                    zip(energies, total_transmission_prob, fed_transmission_prob),
                    key=lambda x: x[0],
                )
                energies, total_transmission_prob, fed_transmission_prob = map(
                    list, zip(*combined)
                )

                # Determine group and label base
                if thick is None:
                    group = "PIPS"
                    label_base = f"{det} mm$^2$ PIPS, 1.0 mm"
                elif thick == 2000:
                    group = "SiLi2"
                    label_base = f"{det} mm$^2$ Si(Li), 2.0 mm"
                elif thick == 5000:
                    group = "SiLi5"
                    label_base = f"{det} mm$^2$ Si(Li), 5.0 mm"
                else:
                    group = None
                    label_base = f"{det} mm$^2$"

                # Pick marker deterministically for this group + det
                marker = marker_map.get((group, det), "o")

                print(marker, group, det)

                markersize = 3

                # Plot on top axis and grab handle
                line0, = axs[0].plot(
                    energies,
                    total_transmission_prob,
                    marker=marker,
                    markersize=markersize,
                    markeredgecolor='black',
                    markeredgewidth=0.1,
                    linewidth=0.5,
                    label=label_base,
                    color=None,  # use color cycle
                )

                # Match color/marker on bottom axis
                axs[1].plot(
                    energies,
                    fed_transmission_prob,
                    marker=marker,
                    markersize=markersize,
                    markeredgecolor='black',
                    markeredgewidth=0.1,
                    linewidth=0.5,
                    color=line0.get_color(),
                    label=f"{label_base} (FED)",
                )

                # Store handle for grouped legend
                if group is not None:
                    curve_handles[group][det] = line0

    axs[0].set_ylabel(r"T$_{M}(E)$ [\%]")
    axs[1].set_ylabel(r"T$_{FED}(E)$ [\%]")

    # remove x labels from top plot
    axs[0].set_xticklabels([])
    axs[1].set_xlabel(r"Energy [keV]")

    # -------------------------------------------------------
    # Build grouped legend with marker edges visible
    # -------------------------------------------------------
    from matplotlib.lines import Line2D

    legend_handles = []
    legend_labels = []

    legend_groups = [
        ("PIPS",  r"PIPS"   + "\n" + r"1 mm"),
        ("SiLi2", r"Si(Li)" + "\n" + r"2 mm"),
        ("SiLi5", r"Si(Li)" + "\n" + r"5 mm"),
    ]

    for group_key, group_title in legend_groups:
        det_map = curve_handles[group_key]
        if not det_map:
            continue

        # ---------- Group heading (no marker) ----------
        legend_handles.append(Line2D([], [], linestyle="none"))
        legend_labels.append(group_title)

        # ---------- Sorted detector entries ----------
        for det in sorted(det_map.keys()):

            # Original plotted line handle
            orig = det_map[det]

            # NEW: Legend-safe handle that preserves marker shape + edge
            legend_handles.append(
                Line2D(
                    [], [],
                    marker=orig.get_marker(),
                    markersize=markersize,
                    markerfacecolor=orig.get_color(),          # fill color
                    markeredgecolor=orig.get_markeredgecolor(), # black outline
                    markeredgewidth=orig.get_markeredgewidth(),
                    linewidth=0,  # IMPORTANT: removes the line so edge is visible
                )
            )
            legend_labels.append(rf"{det} mm$^2$")

        # ---------- Spacer line ----------
        legend_handles.append(Line2D([], [], linestyle="none"))
        legend_labels.append("")

    # ---------- Final legend placement ----------
    axs[0].legend(
        handles=legend_handles,
        labels=legend_labels,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        frameon=False,
        handlelength=1.5,
    )

    labels = [r"(a)", r"(b)"]
    for i, ax in enumerate(axs.flatten()):
        ax.set_xlim(100, 1999)
        ax.minorticks_on()
        ax.tick_params(
            axis="both",
            which="minor",
            direction="in",
            top=True,
            right=True,
            length=2,
        )
        ax.tick_params(
            axis="both",
            which="major",
            direction="in",
            top=True,
            right=True,
            length=4,
        )
        ax.text(0.05, 0.95, labels[i], ha="left", va="top", transform=ax.transAxes)

    axs[0].set_xlim(100, 5000)
    axs[1].set_xlim(100, 5000)

    plt.tight_layout()
    fig.subplots_adjust(left=0.138, right=0.736, bottom=0.1, top=0.99, hspace=0.025)
    plt.savefig("transmission_probability_vs_energy_det_comparison.png", dpi=300)
    plt.savefig("transmission_probability_vs_energy_det_comparison.pdf", dpi=300)
    plt.show()

plot()

