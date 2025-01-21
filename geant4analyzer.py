import ROOT
import numpy as np
import matplotlib.pyplot as plt
import warnings
import os
from lmfit import minimize, Parameters
import sigfig

# for virtual environment on mac
# source $(brew --prefix root)/bin/thisroot.sh  # for ROOT

class Geant4Analyzer:
    def __init__(self, geant4_root_file_path: str = None, geant4_histogram_name: str = None, experimental_root_file_path: str = None, experimental_histogram_name: str = None):
        
        self.scale = None
        self.scale_range = None

        self.threshold_result = None
        self.threshold_range = None

        if geant4_root_file_path is None:
            self.geant4_root_file_path = None
            self.geant4_root_file = None
            self.geant4_histogram_name = None
            self.geant4_histogram = None
            self.geant4_bin_content = None
            self.geant4_bin_centers = None
            self.geant4_bin_edges = None
            self.geant4_bin_uncertainties = None
        else:
            if geant4_histogram_name is None:
                raise ValueError("geant4_histogram_name must be provided if geant4_root_file is provided.")
            
            self.geant4_root_file_path = geant4_root_file_path

            if not os.path.exists(geant4_root_file_path):
                raise ValueError(f"File '{geant4_root_file_path}' not found.")
            
            self.geant4_histogram_name = geant4_histogram_name
            self.geant4_root_file = ROOT.TFile(geant4_root_file_path, "READ")
            self.geant4_histogram = self.geant4_root_file.Get(geant4_histogram_name)
            
            if not self.geant4_histogram:
                raise ValueError(f"Histogram '{geant4_histogram_name}' not found in file '{geant4_root_file_path}'.")
            
            self.geant4_bin_content, self.geant4_bin_centers, self.geant4_bin_edges, self.geant4_bin_uncertainties = self.extract_geant4_histogram_data()

        if experimental_root_file_path is None:
            self.experimental_root_file_path = None
            self.experimental_root_file = None
            self.experimental_histogram_name = None
            self.experimental_histogram = None
            self.experimental_bin_content = None
            self.experimental_bin_centers = None
            self.experimental_bin_edges = None
            self.experimental_bin_uncertainties = None
        else:
            if experimental_histogram_name is None:
                raise ValueError("experimental_histogram_name must be provided if experimental_root_file is provided.")
            
            self.experimental_root_file_path = experimental_root_file_path

            if not os.path.exists(experimental_root_file_path):
                raise ValueError(f"File '{experimental_root_file_path}' not found.")
            
            self.experimental_histogram_name = experimental_histogram_name
            self.experimental_root_file = ROOT.TFile(experimental_root_file_path, "READ")
            self.experimental_histogram = self.experimental_root_file.Get(experimental_histogram_name)
            
            if not self.experimental_histogram:
                raise ValueError(f"Histogram '{experimental_histogram_name}' not found in file '{experimental_root_file_path}'.")
            
            self.experimental_bin_content, self.experimental_bin_centers, self.experimental_bin_edges, self.experimental_bin_uncertainties = self.extract_experimental_histogram_data()
    
    def extract_geant4_histogram_data(self):
        bin_content, bin_centers, bin_edges, bin_uncertainties = [], [], [], []

        for i in range(1, self.geant4_histogram.GetNbinsX() + 1):  # ROOT bins are 1-indexed
            content = self.geant4_histogram.GetBinContent(i)
            center = self.geant4_histogram.GetBinCenter(i) * 1000  # Convert to keV
            width = self.geant4_histogram.GetBinWidth(i) * 1000  # Convert to keV
            uncertainty = self.geant4_histogram.GetBinError(i)
            edge = center - 0.5 * width
            
            bin_content.append(content)
            bin_centers.append(center)
            bin_uncertainties.append(uncertainty)
            bin_edges.append(edge)
        
        # Add the last bin edge
        last_bin_edge = bin_centers[-1] + 0.5 * self.geant4_histogram.GetBinWidth(self.geant4_histogram.GetNbinsX()) * 1000
        bin_edges.append(last_bin_edge)
        
        return (
            np.array(bin_content),
            np.array(bin_centers),
            np.array(bin_edges),
            np.array(bin_uncertainties),
        )
    
    def extract_experimental_histogram_data(self):
        bin_content, bin_centers, bin_edges, bin_uncertainties = [], [], [], []

        for i in range(1, self.experimental_histogram.GetNbinsX() + 1):  # ROOT bins are 1-indexed
            content = self.experimental_histogram.GetBinContent(i)
            center = self.experimental_histogram.GetBinCenter(i)
            width = self.experimental_histogram.GetBinWidth(i)
            uncertainty = self.experimental_histogram.GetBinError(i)
            edge = center - 0.5 * width
            
            bin_content.append(content)
            bin_centers.append(center)
            bin_uncertainties.append(uncertainty)
            bin_edges.append(edge)
        
        # Add the last bin edge
        last_bin_edge = bin_centers[-1] + 0.5 * self.experimental_histogram.GetBinWidth(self.experimental_histogram.GetNbinsX())
        bin_edges.append(last_bin_edge)
        
        return (
            np.array(bin_content),
            np.array(bin_centers),
            np.array(bin_edges),
            np.array(bin_uncertainties),
        )
    
    def gaussian_smear_simulation(self, fwhm):
        if self.geant4_root_file is None:
            raise ValueError("Geant4 histogram not available.")
        
        bin_contents = self.geant4_bin_content
        bin_centers = self.geant4_bin_centers
        bin_uncertainity = self.geant4_bin_uncertainties
        
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        smeared_bin_contents = np.zeros_like(bin_contents)
        smeared_bin_contents_uncertainity = np.zeros_like(bin_contents)
        
        bin_centers = np.array(bin_centers)
        total_counts = np.sum(bin_contents)
        
        for i in range(len(bin_centers)):
            bin_center = bin_centers[i]
            
            bin_content = bin_contents[i]
            bin_content_uncertainity = bin_uncertainity[i]
            
            weights = np.exp(-0.5 * ((bin_centers - bin_center) / sigma) ** 2)
            weights /= np.sum(weights)
            
            smeared_bin_contents += bin_content * weights
            smeared_bin_contents_uncertainity += bin_content_uncertainity * weights
            
        smeared_total_counts = np.sum(smeared_bin_contents)
        
        if total_counts != round(smeared_total_counts):
            warnings.warn(f"Total counts before and after smearing do not match: {total_counts} vs {smeared_total_counts}")
            print("\n")
        
        print(f"Smearing {self.geant4_histogram_name} histogram with FWHM = {fwhm} keV")
        
        self.geant4_bin_content = smeared_bin_contents
        self.geant4_bin_uncertainties = smeared_bin_contents_uncertainity

    def scale_geant4_to_experiment(self, scaling_range=None):
        """
        Scales the Geant4 histogram to minimize residuals with the experimental histogram using lmfit.

        Args:
            scaling_range (tuple, optional): A tuple specifying the (min, max) range in keV to use for scaling.
                                            If None, the entire range of the histograms will be used.

        Returns:
            float: The optimal scaling factor.
        """
        # Check that both histograms exist
        if self.geant4_bin_content is None or self.experimental_bin_content is None:
            raise ValueError("Both Geant4 and experimental histograms must exist to perform scaling.")

        # Ensure the histograms have the same number of bins
        if len(self.geant4_bin_content) != len(self.experimental_bin_content):
            raise ValueError("Geant4 and experimental histograms must have the same number of bins.")

        # Apply the scaling range
        if scaling_range is not None:
            min_range, max_range = scaling_range
            range_mask = (self.geant4_bin_centers >= min_range) & (self.geant4_bin_centers <= max_range)

            # Extract the bins within the specified range
            geant4_content = self.geant4_bin_content[range_mask]
            experiment_content = self.experimental_bin_content[range_mask]
        else:
            geant4_content = self.geant4_bin_content
            experiment_content = self.experimental_bin_content

        # Define the residuals function
        def residuals(params):
            scale_factor = params['scale']
            scaled_geant4 = geant4_content * scale_factor
            return scaled_geant4 - experiment_content

        # Define parameters
        params = Parameters()
        params.add('scale', value=1.0, min=0)  # Adjust bounds as needed

        # Perform the fit
        result = minimize(residuals, params)

        # Check the success of the fit
        if not result.success:
            raise RuntimeError("Optimization failed to find a scaling factor.")

        # Extract the optimal scale factor
        optimal_scale = result.params['scale'].value
        optimal_scale_uncertainty = result.params['scale'].stderr

        self.scale = optimal_scale
        self.scale_range = scaling_range

        print(f"Optimal scale factor: {formatted_round(optimal_scale, optimal_scale_uncertainty)}")

        # Update the Geant4 histogram with the optimal scale factor
        self.geant4_bin_content *= optimal_scale
        self.geant4_bin_uncertainties *= optimal_scale

        return optimal_scale

    def apply_threshold_to_geant4(self, threshold_range):
        """
        Fits an arctan threshold function to the experimental data and applies it to the Geant4 histogram.

        Parameters:
            threshold_range (tuple): (min, max) range of energies (in keV) to use for fitting the threshold.

        Returns:
            lmfit.model.ModelResult: The result of the threshold fitting process.
        """
        if self.experimental_bin_content is None or self.geant4_bin_content is None:
            raise ValueError("Both experimental and Geant4 histograms must be available.")

        # Mask data within the threshold range
        threshold_mask = (self.experimental_bin_centers >= threshold_range[0]) & \
                        (self.experimental_bin_centers <= threshold_range[1])

        threshold_exp_bin_centers = self.experimental_bin_centers[threshold_mask]
        threshold_exp_hist_data = self.experimental_bin_content[threshold_mask]

        # Define lmfit parameters
        threshold_params = Parameters()
        threshold_params.add('width', value=90.0)   # Initial guess for threshold smoothness
        threshold_params.add('scale', value=2000.0)  # Initial guess for scaling factor
        threshold_params.add('phase', value=1.0)   # Initial guess for phase offset

        # Perform threshold fitting
        threshold_result = minimize(threshold_residuals, threshold_params,
                                    args=(threshold_exp_bin_centers, threshold_exp_hist_data))

        # Apply the fitted threshold to the Geant4 histogram
        if not threshold_result.success:
            warnings.warn("Threshold fitting was unsuccessful.")

        fitted_width = threshold_result.params['width'].value
        fitted_width_uncertainty = threshold_result.params['width'].stderr
        
        
        fitted_scale = threshold_result.params['scale'].value
        fitted_scale_uncertainty = threshold_result.params['scale'].stderr
        
        fitted_phase = threshold_result.params['phase'].value
        fitted_phase_uncertainty = threshold_result.params['phase'].stderr

        # Apply the threshold to the Geant4 histogram for energies below the threshold max
        below_threshold_mask = self.geant4_bin_centers < threshold_range[1]
        threshold_values = arctan_threshold(self.geant4_bin_centers[below_threshold_mask],
                                            fitted_width, fitted_scale, fitted_phase)

        self.geant4_bin_content[below_threshold_mask] = np.minimum(
            self.geant4_bin_content[below_threshold_mask], threshold_values)

        # Set uncertainties to zero for bins within the threshold region
        self.geant4_bin_uncertainties[below_threshold_mask] = 0

        print(f"Threshold range: {threshold_range}")
        print(f"Threshold fit parameters: width = {formatted_round(fitted_width, fitted_width_uncertainty)}, scale = {formatted_round(fitted_scale, fitted_scale_uncertainty)}, phase = {formatted_round(fitted_phase, fitted_phase_uncertainty)}")
        self.threshold_range = threshold_range
        self.threshold_result = threshold_result
        
    def plot_experiment(self, ax: plt.Axes):
        if self.experimental_bin_content is None or self.experimental_bin_edges is None:
            raise ValueError("Geant4 simulation data is not available for plotting.")
        
        ax.stairs(
            values=self.experimental_bin_content,
            edges=self.experimental_bin_edges,
            label=f"{self.experimental_histogram_name}",
            color="dodgerblue",
            linewidth=0.5,
        )

        ax.fill_between(
            self.experimental_bin_centers,
            self.experimental_bin_content - self.experimental_bin_uncertainties,
            self.experimental_bin_content + self.experimental_bin_uncertainties,
            color='dodgerblue',
            alpha=0.5,
        )

        ax.set_ylabel("Counts")
        ax.set_xlabel("Energy [keV]")

        ax.legend(loc='upper center')
        ax.minorticks_on()
        ax.tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=3)
        ax.tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=5)

    def plot_simulation(self, ax: plt.Axes):
        if self.geant4_bin_content is None or self.geant4_bin_edges is None:
            raise ValueError("Geant4 simulation data is not available for plotting.")
        
        ax.stairs(
            values=self.geant4_bin_content,
            edges=self.geant4_bin_edges,
            label=f"Geant4: {self.geant4_histogram_name}",
            color="lightcoral",
            linewidth=0.5,
        )

        ax.fill_between(
            self.geant4_bin_centers,
            self.geant4_bin_content - self.geant4_bin_uncertainties,
            self.geant4_bin_content + self.geant4_bin_uncertainties,
            color='lightcoral',
            alpha=0.5,
        )

        ax.set_ylabel("Counts")
        ax.set_xlabel("Energy [keV]")

        ax.legend(loc='upper center')
        ax.minorticks_on()
        ax.tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=3)
        ax.tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=5)

    def plot_residuals(self, ax: plt.Axes):
        if self.experimental_bin_centers is None or self.geant4_bin_centers is None:
            raise ValueError("Experimental and Geant4 bin centers must be available for plotting residuals.")

        if len(self.experimental_bin_centers) != len(self.geant4_bin_centers):
            raise ValueError("Experimental and Geant4 bin centers must have the same length.")

        residuals = self.experimental_bin_content - self.geant4_bin_content

        residuals_uncertainty = np.sqrt(self.experimental_bin_uncertainties**2 + self.geant4_bin_uncertainties**2)

        # plot points with fill between as the error bars
        ax.plot(self.experimental_bin_centers, residuals, 'o', color='salmon', markersize=0.5, label="Residuals")
        ax.fill_between(
            self.experimental_bin_centers, 
            residuals - residuals_uncertainty, 
            residuals + residuals_uncertainty, 
            color='salmon', 
            alpha=0.5
        )
        ax.axhline(0, color='black', linestyle='-', linewidth=0.5)

        ax.set_xlabel("Energy [keV]")
        ax.set_ylabel("Residuals [keV]")

        # ax.legend(loc='upper center')
        ax.minorticks_on()
        ax.tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=3)
        ax.tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=5)

    def plot_percent_difference(self, ax: plt.Axes):
        """
        Plots the percent difference between experimental and simulation histograms.

        Args:
            ax (plt.Axes): The axis on which to plot the percent difference.
        """
        if self.experimental_bin_centers is None or self.geant4_bin_centers is None:
            raise ValueError("Experimental and Geant4 bin centers must be available for plotting percent difference.")

        if len(self.experimental_bin_centers) != len(self.geant4_bin_centers):
            raise ValueError("Experimental and Geant4 bin centers must have the same length.")

        # Calculate percent difference
        percent_difference = 100 * (self.experimental_bin_content - self.geant4_bin_content) / self.experimental_bin_content

        # Calculate uncertainties for percent difference
        relative_uncertainty = np.sqrt(
            (self.experimental_bin_uncertainties / self.experimental_bin_content) ** 2 +
            (self.geant4_bin_uncertainties / self.geant4_bin_content) ** 2
        )
        percent_difference_uncertainty = np.abs(percent_difference) * relative_uncertainty

        # Plot percent difference with error bands
        ax.plot(self.experimental_bin_centers, percent_difference, 'o', color='blue', markersize=0.5, label="Percent Difference [%]")
        ax.fill_between(
            self.experimental_bin_centers,
            percent_difference - percent_difference_uncertainty,
            percent_difference + percent_difference_uncertainty,
            color='blue',
            alpha=0.5
        )

        ax.axhline(0, color='black', linestyle='-', linewidth=0.5)
        ax.set_xlabel("Energy [keV]")
        ax.set_ylabel("Percent Difference [%]")

        ax.legend(loc='upper center')
        ax.minorticks_on()
        ax.tick_params(axis='both', which='minor', direction='in', top=True, right=True, left=True, bottom=True, length=3)
        ax.tick_params(axis='both', which='major', direction='in', top=True, right=True, left=True, bottom=True, length=5)

    def plot(self, experiment=True, simulation=True, residuals=False, percent_difference=False, same_axes=False):
        """
        Combines experimental, simulation, and residuals plots based on user input.

        Args:
            experiment (bool): Whether to plot the experimental data.
            simulation (bool): Whether to plot the simulation data.
            residuals (bool): Whether to plot residuals between experiment and simulation.
            same_axes (bool): Whether to overlay experiment and simulation plots on the same axes.
        """
        try:
            # Ensure at least one dataset is enabled
            if not (experiment or simulation or residuals):
                raise ValueError("At least one of 'experiment', 'simulation', or 'residuals' must be True to generate a plot.")

            # Create subplots based on conditions
            n_rows = 1 + residuals + percent_difference
            if same_axes:
                if n_rows > 1:
                    fig, axs = plt.subplots(n_rows, 1, figsize=(10, 6), sharex=True, gridspec_kw={"height_ratios": [2] + [1] * (n_rows - 1)})
                    ax_main = axs[0]
                    ax_residuals = axs[1] if residuals else None
                    ax_percent_difference = axs[-1] if percent_difference else None
                else:
                    fig, ax_main = plt.subplots(1, 1, figsize=(10, 6))
            else:
                fig, axs = plt.subplots(n_rows + 1, 1, figsize=(10, 6), sharex=True)
                ax_experiment, ax_simulation = axs[:2]
                ax_residuals = axs[2] if residuals else None
                ax_percent_difference = axs[3] if percent_difference else None

            # Plot scale range
            if self.scale_range:
                for ax in axs:
                    ax.axvspan(self.scale_range[0], self.scale_range[1], color='lightgreen', alpha=0.3, label="Scaling Range")

            # Plot threshold range
            if self.threshold_range:
                for ax in axs:
                    ax.axvspan(self.threshold_range[0], self.threshold_range[1], color='lightsalmon', alpha=0.3, label="Threshold Range")

            # Plot experiment
            if experiment:
                if same_axes:
                    self.plot_experiment(ax_main if not residuals else ax_main)
                else:
                    self.plot_experiment(ax_experiment)

            # Plot simulation
            if simulation:
                if same_axes:
                    self.plot_simulation(ax_main if not residuals else ax_main)
                else:
                    self.plot_simulation(ax_simulation)

            # Plot residuals
            if residuals:
                if same_axes:
                    self.plot_residuals(ax_residuals)
                else:
                    self.plot_residuals(ax_residuals)

            # Plot percent difference
            if percent_difference and ax_percent_difference:
                self.plot_percent_difference(ax_percent_difference)

            # Final layout adjustments
            fig.tight_layout()

            # Show plot
            # plt.show()

        except ValueError as e:
            warnings.warn(f"ValueError encountered: {e}")
        except Exception as e:
            warnings.warn(f"An unexpected error occurred: {e}")

        return axs
    
def arctan_threshold(x, width, scale, phase):
    """Arctan threshold function."""
    arctan_values = scale * (np.arctan((x - phase) / width) / np.pi)
    arctan_values[arctan_values < 0] = 0  # Clamp negative values to 0
    return arctan_values

def threshold_residuals(params, x, data):
    """Residuals for fitting the arctan threshold function."""
    width = params['width'].value
    scale = params['scale'].value
    phase = params['phase'].value
    model = arctan_threshold(x, width, scale, phase)
    return data - model

def formatted_round(value, uncertainty):
    import sigfig
    return sigfig.round(value, uncertainty, style='Drake', sep='external_brackets', spacer='')

if __name__ == "__main__":
    simulation_root_file = "./207Bi/Sept2024_LSU/geant_sim/run_98_ICESPICE_RadDecay_z83_a207_e0keV_f70mm_g30mm_n100000000_PIPS1000_AllProcesses_Si02Window50nm_Source500nmThick.root"
    simulation_histogram_name = "Esil"

    experimental_root_file = "./207Bi/Sept2024_LSU/exp_data/207Bi_ICESPICE_f70mm_g30mm_run_14_15.root"
    experimental_histogram_name = "PIPS1000EnergyCalibrated"

    analyzer = Geant4Analyzer(geant4_root_file_path=simulation_root_file, geant4_histogram_name=simulation_histogram_name, experimental_root_file_path=experimental_root_file, experimental_histogram_name=experimental_histogram_name)

    # axes = analyzer.plot(experiment=True, simulation=True)
    # axes[1].set_yscale('log')
    # plt.show()

    analyzer.gaussian_smear_simulation(fwhm=10)
    analyzer.scale_geant4_to_experiment(scaling_range=(462, 1075))
    analyzer.apply_threshold_to_geant4(threshold_range=(300, 462))
    axes = analyzer.plot(experiment=True, simulation=True, residuals=True, same_axes=True)
    
    for ax in axes:
        ax.set_xlim(200, 1200)

    plt.show()
