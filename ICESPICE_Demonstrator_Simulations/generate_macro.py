with open("transmission_prob_macro.mac", "w") as f:
    for e in range(50, 2001, 50):
        f.write(f"# {e} keV\n")
        f.write(f"/control/alias File trasmission_probability_PIPS{{detector}}_f{{f}}mm_g{{g}}mm_n{{n}}_energy{e}keV\n")
        f.write(f"/gps/ene/mono {e} keV\n")
        f.write(f"/analysis/setFileName {{File}}\n")
        f.write(f"/run/printProgress 1000000\n")
        f.write(f"/run/beamOn {{n}}\n")
        f.write(f"/control/shell cp {{File}}.root ../ICESPICE_Demonstrator_Simulations/transmission_prob_data/{{pathsuffix}}\n\n")