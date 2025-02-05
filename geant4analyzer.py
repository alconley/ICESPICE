import ROOT
import numpy as np
import matplotlib.pyplot as plt
import warnings
import os
from lmfit import minimize, Parameters
import lmfit
from tabulate import tabulate

# for virtual environment on mac
# source $(brew --prefix root)/bin/thisroot.sh  # for ROOT

class Geant4Analyzer:
    def __init__(self, 
                 geant4_root_file_path: str = None, 
                 geant4_histogram_name: str = None, 
                 experimental_root_file_path: str = None, 
                 experimental_histogram_name: str = None
                 ):
        
        self.scale = None
        self.scale_range = None
        self.rchi2 = None

        self.threshold_result = None
        self.threshold_range = None

        self.fits = {}

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
    
    def set_experimental_histogram_noise_to_zero(self, noise_value: float):
        # Set any bin content below the noise value that corresponds to the bin center to 0

        # check to see if the experimental histogram is available
        if self.experimental_bin_centers is None:
            raise ValueError("Experimental histogram is not available.")
        
        noise_mask = self.experimental_bin_centers < noise_value
        self.experimental_bin_content[noise_mask] = 0

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
            geant4_uncertainties = self.geant4_bin_uncertainties[range_mask]
            experiment_content = self.experimental_bin_content[range_mask]
            experiment_uncertainties = self.experimental_bin_uncertainties[range_mask]
        else:
            geant4_content = self.geant4_bin_content
            geant4_uncertainties = self.geant4_bin_uncertainties
            experiment_content = self.experimental_bin_content
            experiment_uncertainties = self.experimental_bin_uncertainties

        # Define the residuals function
        def residuals(params):
            scale_factor = params['scale']
            scaled_geant4 = geant4_content * scale_factor
            residual = (scaled_geant4 - experiment_content) / np.sqrt(geant4_uncertainties**2 + experiment_uncertainties**2)
            return residual

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

        # Calculate chi-squared
        scaled_geant4 = geant4_content * optimal_scale
        total_uncertainties = np.sqrt(geant4_uncertainties**2 + experiment_uncertainties**2)
        chi_squared = np.sum(((scaled_geant4 - experiment_content) / total_uncertainties) ** 2)

        # Degrees of freedom
        degrees_of_freedom = len(geant4_content) - len(result.var_names)  # Number of bins - number of fit parameters
        reduced_chi_squared = chi_squared / degrees_of_freedom
        self.rchi2 = reduced_chi_squared
        print(f"Reduced chi-squared: {reduced_chi_squared:.2f} (χ² = {chi_squared:.2f}, DoF = {degrees_of_freedom})")

        return optimal_scale

    def apply_threshold_to_geant4(self, threshold_range, initial_guess = (300.0, 2000.0, 300.0)):
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
        threshold_params.add('width', value=initial_guess[0])   # Initial guess for threshold smoothness
        threshold_params.add('scale', value=initial_guess[1])  # Initial guess for scaling factor
        threshold_params.add('phase', value=initial_guess[2])   # Initial guess for phase offset

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

        self.threshold_range = threshold_range
        self.threshold_result = threshold_result

        print(f"Threshold range: {threshold_range}")
        print(f"Threshold fit parameters: width = {formatted_round(fitted_width, fitted_width_uncertainty)}, scale = {formatted_round(fitted_scale, fitted_scale_uncertainty)}, phase = {formatted_round(fitted_phase, fitted_phase_uncertainty)}")

    def GaussianFit(self, name: str, region_markers: tuple, peak_markers: list, 
                            equal_sigma: bool = True, free_position: bool = True,
                            background_params: dict = None, ax: plt.Axes = None, 
                            use_geant4_data: bool = False):
        """
        Multiple Gaussian fit function with background model support.
        
        Parameters:
        - x_data, y_data: Lists of data points.
        - peak_markers: List of peak positions for the Gaussians.
        - equal_sigma: Whether to constrain all Gaussians to have the same sigma.
        - free_position: Whether to allow the positions of Gaussians to vary.
        - background_params: Dictionary containing background model type and parameters.
        """

        if use_geant4_data is False:
            if self.experimental_bin_centers is None or self.experimental_bin_content is None:
                raise ValueError("Experimental data is not available for fitting.")

            x_data = self.experimental_bin_centers
            y_data = self.experimental_bin_content
            bin_width = x_data[1] - x_data[0]
        else:
            if self.geant4_bin_centers is None or self.geant4_bin_content is None:
                raise ValueError("Geant4 data is not available for fitting.")
            
            x_data = self.geant4_bin_centers
            y_data = self.geant4_bin_content
            bin_width = x_data[1] - x_data[0]
        
        # reduce the data to the region of interest
        region_mask = (x_data >= region_markers[0]) & (x_data <= region_markers[1])
        x_data = x_data[region_mask]
        y_data = y_data[region_mask]
        
        # Default background params if none are provided
        if background_params is None:
            background_params = {
                'bg_type': 'linear',
                'slope': ("slope", -np.inf, np.inf, 0.0, True),
                'intercept': ("intercept", -np.inf, np.inf, 0.0, True),
                'a': ("a", -np.inf, np.inf, 0.0, True),
                'b': ("b", -np.inf, np.inf, 0.0, True),
                'c': ("c", -np.inf, np.inf, 0.0, True),
                'exponent': ("exponent", -np.inf, np.inf, 0.0, True),
                'amplitude': ("amplitude", -np.inf, np.inf, 0.0, True),
                'decay': ("decay", -np.inf, np.inf, 0.0, True),
            }
        
        bg_type = background_params.get('bg_type', 'linear')
        slope = background_params.get('slope')
        intercept = background_params.get('intercept')
        a = background_params.get('a')
        b = background_params.get('b')
        c = background_params.get('c')
        amplitude = background_params.get('amplitude')
        exponent = background_params.get('exponent')
        decay = background_params.get('decay')

        # Initialize the model with or without a background based on bg_type
        if bg_type == 'linear': 
            model = lmfit.models.LinearModel(prefix='bg_')
            params = model.make_params(slope=slope[3], intercept=intercept[3])
            params['bg_slope'].set(min=slope[1], max=slope[2], value=slope[3], vary=slope[4])
            params['bg_intercept'].set(min=intercept[1], max=intercept[2], value=intercept[3], vary=intercept[4])
        elif bg_type == 'quadratic':
            model = lmfit.models.QuadraticModel(prefix='bg_')
            params = model.make_params(a=a[3], b=b[3], c=c[3])
            params['bg_a'].set(min=a[1], max=a[2], value=a[3], vary=a[4])
            params['bg_b'].set(min=b[1], max=b[2], value=b[3], vary=b[4])
            params['bg_c'].set(min=c[1], max=c[2], value=c[3], vary=c[4])
        elif bg_type == 'exponential':
            model = lmfit.models.ExponentialModel(prefix='bg_')
            params = model.make_params(amplitude=amplitude[3], decay=decay[3])
            params['bg_amplitude'].set(min=amplitude[1], max=amplitude[2], value=amplitude[3], vary=amplitude[4])
            params['bg_decay'].set(min=decay[1], max=decay[2], value=decay[3], vary=decay[4])
        elif bg_type == 'powerlaw':
            model = lmfit.models.PowerLawModel(prefix='bg_')
            params = model.make_params(amplitude=amplitude[3], exponent=exponent[3])
            params['bg_amplitude'].set(min=amplitude[1], max=amplitude[2], value=amplitude[3], vary=amplitude[4])
            params['bg_exponent'].set(min=exponent[1], max=exponent[2], value=exponent[3], vary=exponent[4])
        elif bg_type == 'none':
            model = None
            params = lmfit.Parameters()
        else:
            raise ValueError("Unsupported background model")

        first_gaussian = lmfit.Model(gaussian, prefix=f'g0_')

        if model is None:
            model = first_gaussian
        else:
            model += first_gaussian
        
        if len(peak_markers) == 0:
            peak_markers = [x_data[np.argmax(y_data)]]


        peak_markers = sorted(peak_markers)  # sort the peak markers in ascending order

        estimated_amplitude = 1000
        estimated_sigma = 10

        params.update(first_gaussian.make_params(amplitude=estimated_amplitude, mean=peak_markers[0], sigma=estimated_sigma))
        params['g0_sigma'].set(min=0)  # Initial constraint for the first Gaussian's sigma
        params[f"g0_amplitude"].set(min=0)

        params.add(f'g0_fwhm', expr=f'2.35482 * g0_sigma')  # FWHM = 2 * sqrt(2 * ln(2)) * sigma
        params[f"g0_fwhm"].set(min=0)

        params.add(f'g0_area', expr=f'g0_amplitude * sqrt(2 * pi) * g0_sigma / {bin_width}')  # Area under the Gaussian
        params[f"g0_area"].set(min=0)

        if not free_position:
            params['g0_mean'].set(vary=False)

        params['g0_mean'].set(min=x_data[0], max=peak_markers[1] if len(peak_markers) > 1 else x_data[-1])

        # Add additional Gaussians
        for i, peak in enumerate(peak_markers[1:], start=1):
            g = lmfit.Model(gaussian, prefix=f'g{i}_')
            model += g

            estimated_amplitude = 1000
            params.update(g.make_params(amplitude=estimated_amplitude, mean=peak, sigma=10))

            min_mean = peak_markers[i-1]
            max_mean = peak_markers[i+1] if i + 1 < len(peak_markers) else x_data[-1]
            params[f'g{i}_mean'].set(min=min_mean, max=max_mean)

            params.add(f'g{i}_fwhm', expr=f'2.35482 * g{i}_sigma')
            params[f"g{i}_fwhm"].set(min=0)

            params.add(f'g{i}_area', expr=f'g{i}_amplitude * sqrt(2 * pi) * g{i}_sigma / {bin_width}')
            params[f"g{i}_area"].set(min=0)

            if equal_sigma:
                params[f'g{i}_sigma'].set(expr='g0_sigma')
            else:
                params[f'g{i}_sigma'].set(min=0)

            params[f'g{i}_amplitude'].set(min=0)

            if not free_position:
                params[f'g{i}_mean'].set(vary=False)

        # Fit the model to the data
        result = model.fit(y_data, params, x=x_data)

        print("\nInitial Parameter Guesses:")
        params.pretty_print()

        print("\nFit Report:")
        print(result.fit_report())

        # Store the fit result
        self.fits[name] = result

        # Plot the fit if axes are provided
        if ax is not None:
            # Plot the total fit

            ax.plot(x_data, result.best_fit, '-', color='blue', linewidth=1)

            # Plot the background model if it exists
            if bg_type != 'none':
                bg_prefix = 'bg_'
                background_fit = result.eval_components(x=x_data)[bg_prefix]
                ax.plot(x_data, background_fit, '-', color='green', linewidth=1)

            # Plot individual Gaussian components
            for i, peak in enumerate(peak_markers):
                gaussian_prefix = f'g{i}_'

                mean = result.params[f'{gaussian_prefix}mean'].value
                sigma = result.params[f'{gaussian_prefix}sigma'].value

                x = np.linspace(mean - 5 * sigma, mean + 5 * sigma, 1000)

                gaussian_fit = result.eval_components(x=x)[gaussian_prefix]
                ax.plot(x, gaussian_fit, '-', color='purple', linewidth=1)

    def print_fit_parameters(self):
        """
        Print the fit parameters (mean, FWHM, and area) for all fits stored in self.fits
        using the `tabulate` library for better formatting. Starts Fit # from 0 and 
        avoids repeating Fit # for multiple entries under the same fit.
        """
        if not self.fits:
            print("No fits available to display.")
            return
            
        print("\Gaussian Fit Parameters")

        # Table headers and data rows
        headers = ["Name", "Fit #", "Index", "Mean", "FWHM", "Area"]
        rows = []

        # Iterate over all fits

        for fit_num, (fit_name, result) in enumerate(self.fits.items()):
            # Determine the number of Gaussian components in this fit
            num_gaussians = sum(1 for key in result.params.keys() if key.startswith("g") and "_mean" in key)

            for i in range(num_gaussians):
                # Extract Gaussian parameters
                name = fit_name
                mean = result.params[f'g{i}_mean'].value
                mean_uncertainty = result.params[f'g{i}_mean'].stderr or 0.0
                fwhm = result.params[f'g{i}_fwhm'].value
                fwhm_uncertainty = result.params[f'g{i}_fwhm'].stderr or 0.0
                area = result.params[f'g{i}_area'].value
                area_uncertainty = result.params[f'g{i}_area'].stderr or 0.0

                # Format the values with uncertainties
                mean_formatted = formatted_round(mean, mean_uncertainty)
                fwhm_formatted = formatted_round(fwhm, fwhm_uncertainty)
                area_formatted = formatted_round(area, area_uncertainty)

                # Add a row to the table
                # Display the Fit # only for the first row of each fit
                # rows.append([f"Fit {fit_num}" if i == 0 else "", i, mean_formatted, fwhm_formatted, area_formatted])

                rows.append([f"{name}" if i == 0 else "", f"Fit {fit_num}" if i == 0 else "", i, mean_formatted, fwhm_formatted, area_formatted])

        # Print the table using `tabulate`
        print(tabulate(rows, headers=headers))

    def average_fwhm(self):
        """
        Calculate the average FWHM of the Gaussian fits stored in self.fits.
        """

        if not self.fits:
            print("No fits available to calculate average FWHM.")
            return None

        # Initialize variables to store the sum of FWHMs and the number of Gaussian components
        total_fwhm = 0
        num_gaussians = 0

        # Iterate over all fits
        for fit_name, result in self.fits.items():
            for param_name in result.params.keys():
                if param_name.endswith("_fwhm"):
                    total_fwhm += result.params[param_name].value
                    num_gaussians += 1

        if num_gaussians == 0:
            print("No Gaussian components found to calculate average FWHM.")
            return None

        # Calculate the average FWHM
        average_fwhm = total_fwhm / num_gaussians

        print(f"\nAverage FWHM: {average_fwhm:.2f} keV")
        return average_fwhm

    def fit_stats(self, name: str, index: int = 0):
        """
        Print the fit statistics for a given fit stored in self.fits.

        Parameters:
        - name: The name of the fit to display statistics for.
        """

        if name not in self.fits:
            print(f"Fit '{name}' not found.")
            return

        result = self.fits[name]

        mean = result.params[f'g{index}_mean'].value
        mean_uncertainty = result.params[f'g{index}_mean'].stderr or 0.0
        fwhm = result.params[f'g{index}_fwhm'].value
        fwhm_uncertainty = result.params[f'g{index}_fwhm'].stderr or 0.0
        area = result.params[f'g{index}_area'].value
        area_uncertainty = result.params[f'g{index}_area'].stderr or 0.0

        return mean, mean_uncertainty, fwhm, fwhm_uncertainty, area, area_uncertainty

    def plot_experiment(self, ax: plt.Axes, color="dodgerblue", label=None, linewidth=0.5, plot_uncertainity=False):
        if self.experimental_bin_content is None or self.experimental_bin_edges is None:
            raise ValueError("Experimental data is not available for plotting.")
        
        name = self.experimental_histogram_name if label is None else label

        ax.stairs(
            values=self.experimental_bin_content,
            edges=self.experimental_bin_edges,
            label=f"{name}",
            color=color,
            linewidth=linewidth,
        )

        if plot_uncertainity:
            ax.fill_between(
                self.experimental_bin_centers,
                self.experimental_bin_content - self.experimental_bin_uncertainties,
                self.experimental_bin_content + self.experimental_bin_uncertainties,
                color=color,
                alpha=0.5,
            )

        ax.set_ylabel("Counts")
        ax.set_xlabel("Energy [keV]")

        # ax.legend(loc='upper left')
        ax.minorticks_on()
        ax.tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=3)
        ax.tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=5)

    def plot_simulation(self, ax: plt.Axes, color="lightcoral", label=None, linewidth=0.5, plot_uncertainity=True):
        if self.geant4_bin_content is None or self.geant4_bin_edges is None:
            raise ValueError("Geant4 simulation data is not available for plotting.")
        
        if label is None:
            label = self.geant4_histogram_name
        else: 
            label = label

        ax.stairs(
            values=self.geant4_bin_content,
            edges=self.geant4_bin_edges,
            label=label,
            color=color,
            linewidth=linewidth,
        )

        if plot_uncertainity:
            ax.fill_between(
                self.geant4_bin_centers,
                self.geant4_bin_content - self.geant4_bin_uncertainties,
                self.geant4_bin_content + self.geant4_bin_uncertainties,
                color=color,
                alpha=0.5,
            )

        ax.set_ylabel("Counts")
        ax.set_xlabel("Energy [keV]")

        # ax.set_ylim(bottom=0.1)
        # ax.set_yscale('log')

        # ax.legend(loc='upper left')
        ax.minorticks_on()
        ax.tick_params(axis='both',which='minor',direction='in',top=True,right=True,left=True,bottom=True,length=3)
        ax.tick_params(axis='both',which='major',direction='in',top=True,right=True,left=True,bottom=True,length=5)

    def plot_residuals(self, ax: plt.Axes):
        """
        Plots residuals between the experimental and Geant4 histograms.

        Args:
            ax (plt.Axes): The axis on which to plot the residuals.
        """
        if self.experimental_bin_centers is None or self.geant4_bin_centers is None:
            raise ValueError("Experimental and Geant4 bin centers must be available for plotting residuals.")

        if len(self.experimental_bin_centers) != len(self.geant4_bin_centers):
            raise ValueError("Experimental and Geant4 bin centers must have the same length.")

        residuals = self.experimental_bin_content - self.geant4_bin_content
        residuals_uncertainty = np.sqrt(self.experimental_bin_uncertainties**2 + self.geant4_bin_uncertainties**2)

        # Plot points with fill between as the error bars
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
        ax.set_ylabel("Residuals")

        # Add reduced chi-squared text to the plot
        if self.rchi2 is not None:
            ax.text(
                0.95, 0.95, 
                f"Reduced χ² = {self.rchi2:.2f}", 
                transform=ax.transAxes, 
                ha='right',
                va='top',
            )

        ax.minorticks_on()
        ax.tick_params(axis='both', which='minor', direction='in', top=True, right=True, left=True, bottom=True, length=3)
        ax.tick_params(axis='both', which='major', direction='in', top=True, right=True, left=True, bottom=True, length=5)
   
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

        ax.legend(loc='upper left')
        ax.minorticks_on()
        ax.tick_params(axis='both', which='minor', direction='in', top=True, right=True, left=True, bottom=True, length=3)
        ax.tick_params(axis='both', which='major', direction='in', top=True, right=True, left=True, bottom=True, length=5)

    def plot(self, experiment=True, simulation=True, residuals=False, percent_difference=False, same_axes=False, save=None):
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
                    ax.axvspan(self.scale_range[0], self.scale_range[1], color='lightgreen', alpha=0.2, label="Scaling Range")

            # Plot threshold range
            if self.threshold_range:
                for ax in axs:
                    ax.axvspan(self.threshold_range[0], self.threshold_range[1], color='lightsalmon', alpha=0.2, label="Threshold Range")

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

            if save is not None:
                plt.savefig(save, dpi=300)
                print(f"Plot saved to {save}")

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

def gaussian(x, amplitude, mean, sigma):
    return amplitude * np.exp(-(x - mean)**2 / (2 * sigma**2))

def formatted_round(value, uncertainty):
    import sigfig
    if uncertainty is None:
        return sigfig.round(value, style='Drake', sep='external_brackets', spacer='')
    else:
        return sigfig.round(value, uncertainty, style='Drake', sep='external_brackets', spacer='')
