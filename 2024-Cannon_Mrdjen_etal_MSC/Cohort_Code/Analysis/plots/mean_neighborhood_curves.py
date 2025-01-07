import pandas as pd
from scipy.stats import zscore
import matplotlib.pyplot as plt
import statsmodels.api as sm
import numpy as np

# Set the backend to 'Agg'
import matplotlib

matplotlib.use('Agg')

# Load the CSV file
file_path = '/Users/bryjc/Downloads/mglia_with_neighbors.csv'  # Update with your actual file path
data = pd.read_csv(file_path)

# List of biomarkers
biomarkers = [
    "Beta.amyloid.1.40", "Beta.amyloid.1.42", "Beta.amyloid.17.24", "PHF.tau", "PHF.tau.AT8",
    "Clusterin", "CD31_CD105", "Fe", "GAD65", "GFAP", "Glutamine.Synthetase",
    "MAG_MCNPase", "MBP", "NEFL_MAP2_NEFH", "Parvalbumin", "PSD95",
    "SMA", "Synaptophysin", "VGAT", "VGlut1"
]

# Colors for disease status
disease_status_colors = {
    "Healthy": "#A6E002",
    "AD": "#172869"
}

# Filter data for required columns and z-score normalization
filtered_data = data[['Region', 'Disease_Status', 'PSEUDOTIME_NORMALIZED'] + biomarkers]
zscored_data = filtered_data.copy()
zscored_data[biomarkers] = zscore(filtered_data[biomarkers])

# Separate data by region
data_ca1 = zscored_data[zscored_data['Region'] == 'CA1']
data_caudate = zscored_data[zscored_data['Region'] == 'Caudate']


# Function to plot LOESS smoothed line with smoothed error lines
def plot_loess_with_smoothed_error(data, biomarker, region):
    plt.figure(figsize=(10, 6))

    # Add a light gray background grid
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='lightgray')

    # Add horizontal dark gray dotted lines at 0.33 and 0.66
    plt.axvline(0.33, color='darkgray', linestyle='dotted', linewidth=2)
    plt.axvline(0.66, color='darkgray', linestyle='dotted', linewidth=2)

    # Plot LOESS smoothed lines with smoothed error lines
    for status in ['AD', 'Healthy']:
        color = disease_status_colors[status]
        subset = data[data['Disease_Status'] == status]
        x = subset['PSEUDOTIME_NORMALIZED']
        y = subset[biomarker]

        # Fit LOESS model
        loess_sm = sm.nonparametric.lowess(y, x, frac=0.3)

        # Calculate standard deviation for each bin
        bins = np.linspace(x.min(), x.max(), 50)
        bin_centers = 0.5 * (bins[1:] + bins[:-1])
        bin_means = np.zeros(len(bin_centers))
        bin_stds = np.zeros(len(bin_centers))

        for i in range(len(bins) - 1):
            bin_data = y[(x >= bins[i]) & (x < bins[i + 1])]
            if len(bin_data) > 0:
                bin_means[i] = bin_data.mean()
                bin_stds[i] = bin_data.std()

        # Fit LOESS model to the upper and lower bounds
        loess_upper = sm.nonparametric.lowess(bin_means + bin_stds, bin_centers, frac=0.3)
        loess_lower = sm.nonparametric.lowess(bin_means - bin_stds, bin_centers, frac=0.3)

        # Plot LOESS smoothed line
        plt.plot(loess_sm[:, 0], loess_sm[:, 1], label=status, color=color, linewidth=2)

        # Add smoothed dashed edge lines for deviation
        plt.plot(loess_upper[:, 0], loess_upper[:, 1], color=color, linestyle=(0, (5, 10)), alpha=0.7, linewidth=2)
        plt.plot(loess_lower[:, 0], loess_lower[:, 1], color=color, linestyle=(0, (5, 10)), alpha=0.7, linewidth=2)

    # Set x-axis limits to 0 to 1
    plt.xlim(0, 1)

    plt.xlabel('PSEUDOTIME_NORMALIZED')
    plt.ylabel(f'Z-scored {biomarker}')
    plt.title(f'LOESS Smoothed Plot for {biomarker} in {region}')
    #plt.legend()
    plt.savefig(f'{biomarker}_{region}.tiff')  # Save the plot as a PNG file
    plt.close()


# Function to remove outliers beyond three standard deviations
def remove_outliers(data, column):
    mean_val = data[column].mean()
    std_val = data[column].std()
    threshold = 3
    return data[(data[column] > (mean_val - threshold * std_val)) & (data[column] < (mean_val + threshold * std_val))]


# Create plots for each biomarker
for biomarker in biomarkers:
    # CA1 region
    data_ca1_filtered = remove_outliers(data_ca1, biomarker)
    plot_loess_with_smoothed_error(data_ca1_filtered, biomarker, 'CA1')

    # Caudate region
    data_caudate_filtered = remove_outliers(data_caudate, biomarker)
    plot_loess_with_smoothed_error(data_caudate_filtered, biomarker, 'Caudate')
