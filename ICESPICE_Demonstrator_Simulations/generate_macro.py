
def trans_prob(det=1000):
    n = 1000000
    with open(f"../ICESPICE_Demonstrator_Simulations/PIPS{det}_transmission_prob_macro_with_ICESPICE.mac", "w") as file:


        file.write("# Auto-generated Geant4 macro for transmission probability simulations\n")
        file.write(f"\n/control/alias n {n}\n")
        if det == 1000 or det == 500:
            file.write(f"/control/alias pathsuffix 5N42_1x1x1_8in_PIPS{det}\n")
        elif det == 300:
            file.write(f"/control/alias pathsuffix 5N42_1x1x1_16in_PIPS{det}\n")

        file.write(f"\n/control/verbose 0\n")
        file.write(f"/run/verbose 0\n")
        file.write(f"/event/verbose 0\n")
        file.write(f"/tracking/verbose 0\n\n")
        file.write(f"/run/setCut 1 um\n")
        file.write(f"/run/initialize\n")

        file.write(f"/gps/ang/type iso\n")
        file.write(f"/gps/pos/type Point\n")
        file.write(f"/gps/particle e-\n")
        file.write(f"/gps/ene/type Mono\n")

# /control/alias pathsuffix 5N42_1x1x1_8in_PIPS{detector}/

        # for f in [68, 69, 70, 71, 72]:
        #     file.write(f"\n# --- f = {f} mm ---\n")
        #     file.write(f"\n/gps/pos/centre 0 0 {f} mm\n")
        #     for g in [28, 29, 30, 31, 32]:
        #         file.write(f"\n# --- g = {g} mm ---\n")
        #         file.write(f"/ICESPICE/DetectorPosition -{g}\n")

        #         for e in range(50, 2001, 50):
        #             file.write(f"\n/control/alias File trasmission_probability_PIPS{det}_f{f}mm_g{g}mm_n{n}_energy{e}keV\n")
        #             file.write(f"/gps/ene/mono {e} keV\n")
        #             file.write(f"/analysis/setFileName {{File}}\n")
        #             file.write(f"/run/printProgress 100000\n")
        #             file.write(f"/run/beamOn {{n}}\n")
        #             file.write(f"/control/shell cp {{File}}.root ../ICESPICE_Demonstrator_Simulations/transmission_prob_data/{{pathsuffix}}\n\n")

        for f in [68, 72]:
            file.write(f"\n# --- f = {f} mm ---\n")
            file.write(f"\n/gps/pos/centre 0 0 {f} mm\n")
            for g in [28, 32]:
                file.write(f"\n# --- g = {g} mm ---\n")
                file.write(f"/ICESPICE/DetectorPosition -{g}\n")

                for e in range(50, 2001, 25):
                    file.write(f"\n/control/alias File trasmission_probability_PIPS{det}_f{f}mm_g{g}mm_n{n}_energy{e}keV\n")
                    file.write(f"/gps/ene/mono {e} keV\n")
                    file.write(f"/analysis/setFileName {{File}}\n")
                    file.write(f"/run/printProgress 100000\n")
                    file.write(f"/run/beamOn {{n}}\n")
                    file.write(f"/control/shell cp {{File}}.root ../ICESPICE_Demonstrator_Simulations/transmission_prob_data/{{pathsuffix}}\n\n")



def generate_all_electrons_macro(
    N_total=100_000_000,
    runno=8, thickness_nm=1, detector_um=1000, f_mm=70, g_mm=30, cut=1, theta_deg=0, phi_deg=0,
    out_path="./ICESPICE_Demonstrator_Simulations/self_defined_207Bi.mac",
    add_git=True,
):
    """
    Build a Geant4 macro that simulates ALL electron lines (Auger + CE K/L/M)
    using monoenergetic /gps electrons at the energies and intensities provided.
    The /run/beamOn per line is assigned ~ proportional to the line's intensity,
    normalized so the total primaries across all lines equals N.
    """

    # (label, energy_keV, intensity_percent)
    lines = [
        # ("Auger_L_7p97",        7.97,      54.4),
        # ("Auger_K_56p7",        56.7,       2.9),
        ("CEK_240p10",         240.10,      1.88e-4),
        ("CEL_312p24",         312.24,      3.2e-5),
        ("CEM_324p25",         324.25,      7.5e-6),
        ("CEK_481p6935",       481.6935,    1.537),
        ("CEL_553p8372",       553.8372,    0.442),
        ("CEM_565p8473",       565.8473,    0.111),
        ("CEK_809p77",         809.77,      0.00246),
        ("CEL_881p91",         881.91,      0.000407),
        ("CEM_893p92",         893.92,      0.000095),
        ("CEK_975p651",        975.651,     7.08),
        ("CEL_1047p795",      1047.795,     1.84),
        ("CEM_1059p805",      1059.805,     0.44),
        ("CEK_1354p20",       1354.20,      0.000355),
        ("CEL_1426p34",       1426.34,      0.0000613),
        ("CEM_1438p35",       1438.35,      0.0000144),
        ("CEK_1682p224",      1682.224,     0.0238),
        ("CEL_1754p367",      1754.367,     0.0034),
    ]

    # Calculate beamOn per line
    alloc = [round(N_total * (p / 100.0)) for _, _, p in lines]


    # Filenames (combined + per-line)
    File = f"run_{runno}_PIPS{detector_um}_f{f_mm}mm_g{g_mm}mm_n{N_total}_207BiThickness{thickness_nm}nm_Cut{cut}um_AllElectrons_theta{theta_deg}_phi{phi_deg}_noICESPICE"

    def file_for(label):
        return f"run_{runno}_PIPS{detector_um}_f{f_mm}mm_g{g_mm}mm_n{N_total}_207BiThickness{thickness_nm}nm_Cut{cut}um_{label}_theta{theta_deg}_phi{phi_deg}_noICESPICE"

    # Header / common setup
    header = f"""# Auto-generated by generate_all_electrons_macro()

# variables
/control/alias runno {runno}
/control/alias thickness {thickness_nm}
/control/alias detector {detector_um}
/control/alias f {f_mm}
/control/alias g {g_mm}
/control/alias n {N_total}
/control/alias cut {cut}
/control/alias theta_deg {theta_deg}
/control/alias phi_deg {phi_deg}

/control/alias File {File}

# Optional per-line aliases (not strictly needed, but handy for debugging)
""" + "\n".join([f"/control/alias File_{lbl} {file_for(lbl)}" for (lbl, _, _) in lines]) + """

/run/setCut {cut} um

# Initialize the run manager (prepare the simulation geometry and physics)
/run/initialize

# These must come right after the initialization
/process/em/printParameters

# Verbosity
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0

# Geometry placement
# /ICESPICE/DetectorPosition -25.4   # use alias consistently
/ICESPICE/DetectorPosition -{g}   # use alias consistently

# Primary: electrons
/gps/particle e-

# --- thin cylindrical source geometry ---
/ICESPICE/FSU207BiSourceEnable 1
/ICESPICE/FSU207BiSourcePosition {f}     # Position the 207Bi source backing 
/ICESPICE/FSU207BiSourceThickness {thickness} nanometer
/ICESPICE/FSU207BiSourcePhi {phi_deg}
/ICESPICE/FSU207BiSourceTheta {theta_deg}

#/gps/pos/type Volume
#/gps/pos/shape Cylinder
#/gps/pos/centre 0 0 {f} mm
#/gps/pos/radius 2.5 mm
#/gps/pos/halfz {thickness} nm

# Direction: 4π isotropic
#/gps/ang/type iso

# Progress
/run/printProgress 1000000
"""

    # Blocks per line
    blocks = []
    for (label, energy_keV, _), nline in zip(lines, alloc):
        fname = file_for(label)
        block = f"""
# ---- {label}  ({energy_keV:.6f} keV)
/analysis/setFileName {fname}
/gps/ene/type Mono
/gps/ene/mono {energy_keV:.6f} keV
/run/beamOn {nline}
/control/shell cp {fname}.root ../ICESPICE_Demonstrator_Simulations
"""
        blocks.append(block)

    # hadd + git
    all_parts = " ".join([f"../ICESPICE_Demonstrator_Simulations/{file_for(lbl)}.root" for (lbl, _, _) in lines])
    git_tail = (
        f'\n/control/shell git add ../ICESPICE_Demonstrator_Simulations/{File}.root '
        f'&& git commit -m "Add {File}.root" && git push'
        if add_git else ""
    )
    tail = f"""
# ---- Combine all into one ROOT file
/control/shell hadd -f ../ICESPICE_Demonstrator_Simulations/{File}.root {all_parts}{git_tail}

/control/shell rm {all_parts}
"""

    macro = header + "".join(blocks) + tail
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(macro)
    return out_path

# generate_all_electrons_macro(
#     N_total=1_000_000_000,
#     runno=56,
#     thickness_nm=1,
#     detector_um=1000,
#     f_mm=70, # mm
#     g_mm=30, # mm
#     cut=1, # um
#     theta_deg=0, # degrees
#     phi_deg=0 # degrees
# )

# for f in [68, 69, 70 , 71, 72]:
#     for g in [28, 29, 30, 31, 32]:
#         generate_all_electrons_macro(
#             N_total=100_000_000,
#             runno=55,
#             thickness_nm=1000,
#             detector_um=1000,
#             f_mm=f, # mm
#             g_mm=g, # mm
#             cut=1, # um
#             theta_deg=0, # degrees
#             phi_deg=0 # degrees
#         )

trans_prob(det=1000)