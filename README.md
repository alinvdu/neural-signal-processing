# Neural Signal Processing
The content of this course is based on Mike X Cohen: Complete neural signal processing and analysis: Zero to hero
https://www.udemy.com/course/solved-challenges-ants/

<img width="800" alt="NeuralSignalCourse" src="https://github.com/bobergsatoko/neural-signal-processing/assets/16021447/81838493-25e9-4565-8fc8-44268a3a72da">

**Time Domain Analysis**
1. erp-time-peaks.m -> Perform ERP (Event Related Potential).ERP is the average electrical activity of the brain in response to a specific stimulus or event.
2. voltage-fluctuations.m -> Voltage Fluctuations. In this project I analysed voltage fluctuations in Event Related Potential and in the end draw a visualization of zones before and after the onset of stimulus.

**Time Frequency Analysis**
1. SlowFastFourierTransform.m -> Perform manual fourier transform in a for loop. Plot the amplitude spectrum.
Visualization.

Visualizations:
<img width="800" alt="ex2" src="https://github.com/bobergsatoko/neural-signal-processing/assets/16021447/5eedf8e3-ba50-46ed-b197-bdffaa8dc12d">


3. TimeFrequency_plot.m -> Perform time frequency analysis. Apply convolution to EEG signal in some specifics channels and plot a time frequency in order to analyse brain activations. In this specific examples I used cycles but there are other examples with full width at half maximum.
4. InverseFTManual.m -> Perform the inverse operation of Fourier Transform manually, in a for loop.
5. ManualTimeFilter.m -> Perform a manual filter in the frequency domain after applying Welch method, this is a simulation of a line filter at a specific frequency.
6. WelchMethodManual.m -> Welch method is used to estimate power spectral density, manually implement the welch method and then use Matlab built-in tools.
7. SpectralSeparation.m -> Try to separate a simulated neural signal in two dipoles.
8. DetrendedFluctuationAnalysis.m -> Detect fractal-like scaling properties and correlations in non-stationary time series data. This is used in disease detection such as Alzheimer, Parkinson, Epilepsy, Schizophrenia and others.
9. 

**Synchronization Analysis**
1. InterSitePhaseClustering.m -> Measures connections between two electrodes in terms of phase on a specific trial.
2. LaplacianSpaceFilter.m -> Surface Laplacian Spatial Filtering - estimates the second spatial derivative in two dimensions across scalp features. This method is aimed at enhancing the spatial resolution of the signals and reducing the effects of volume conduction. This uses a custom implementation of Laplacian Filter provided by professor Mike Cohen and visualizes the resulting data.
3. GrangerCausalityBivariateRegression.m -> Analyzes directed synchrony between EEG channels using Granger Causality. It loads sample EEG data, identifies two channels (fcz and o1), and computes Granger predictions for both the entire dataset and time windows. The results are visualized in a plot.
4. ConnectivityHubs.m -> This Matlab script analyzes brain connectivity using phase-locking value (PLI) from EEG data. It loads EEG data, applies a wavelet convolution for frequency-specific activity, calculates all-to-all PLI between channels, visualizes connectivity matrices, determines a threshold for significant connectivity, binarizes the matrix, and plots hub regions.

**Permutation Statistics**
1. SimualtedFreqAndStats.m -> Simulated EEG data of two conditions is analyzed using wavelet convolution. Differences in power spectra are statistically tested via permutation testing, followed by cluster correction for significance.

**Multivariate**
1. MultivariateMultipleComponentAnalysis.m -> Performs a Multivariate Component Analysis on EEG data, computing and visualizing covariance matrices. It extracts and analyzes the top components, computes phase synchronization, and extends the analysis to eight components, visualizing their power spectra, correlations, and topographies.
