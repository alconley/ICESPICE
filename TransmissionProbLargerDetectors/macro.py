def trans_prob():
    n = 1000000
    with open(f"../TransmissionProbLargerDetectors/transmission_prob_macro.mac", "w") as file:


        file.write("# Auto-generated Geant4 macro for transmission probability simulations\n")
        file.write(f"\n/control/alias n {n}\n")
        
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

        f = 70
        file.write(f"\n# --- f = {f} mm ---\n")
        file.write(f"\n/gps/pos/centre 0 0 {f} mm\n")

        g = 30

        file.write(f"\n# --- g = {g} mm ---\n")
        file.write(f"/ICESPICE/DetectorPosition -{g}\n")

        # det = [25, 50, 150, 300, 450, 600, 900, 1200]
        det = [200, 500]

        for det in det:

            file.write(f"\n/ICESPICE/CustomActiveArea {det} mm2\n")

            for e in range(2000, 5001, 100):
                file.write(f"\n/control/alias File transmission_probability_5000umThick_200nmSiWindow_{det}mm2_f{f}mm_g{g}mm_n{n}_energy{e}keV\n")
                file.write(f"/gps/ene/mono {e} keV\n")
                file.write(f"/analysis/setFileName {{File}}\n")
                file.write(f"/run/printProgress 100000\n")
                file.write(f"/run/beamOn {{n}}\n")
                file.write(f"/control/shell cp {{File}}.root ../TransmissionProbLargerDetectors/data/\n")


trans_prob()
