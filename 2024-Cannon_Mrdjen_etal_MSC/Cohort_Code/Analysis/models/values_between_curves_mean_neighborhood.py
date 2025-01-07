import pandas as pd
from scipy.stats import zscore
import statsmodels.api as sm
import numpy as np

# Load the CSV file
file_path = '/Users/bryjc/Downloads/mglia_with_neighbors.csv' # Update with your actual file path
data = pd.read_csv(file_path)

# List of biomarkers
biomarkers = [
    "Beta.amyloid.1.40", "Beta.amyloid.1.42", "Beta.amyloid.17.24", "PHF.tau", "PHF.tau.AT8",
    "Clusterin", "CD31_CD105", "Fe", "GAD65", "GFAP", "Glutamine.Synthetase",
    "MAG_MCNPase", "MBP", "NEFL_MAP2_NEFH", "Parvalbumin", "PSD95",
    "SMA", "Synaptophysin", "VGAT", "VGlut1"
]

# Filter data for required columns and z-score normalization
filtered_data = data[['Region', 'Disease_Status', 'PSEUDOTIME_NORMALIZED'] + biomarkers]
zscored_data = filtered_data.copy()
zscored_data[biomarkers] = zscore(filtered_data[biomarkers])

# Separate data by region
data_ca1 = zscored_data[zscored_data['Region'] == 'CA1']
data_caudate = zscored_data[zscored_data['Region'] == 'Caudate']


# Function to remove outliers beyond three standard deviations
def remove_outliers(data, column):
    mean_val = data[column].mean()
    std_val = data[column].std()
    threshold = 3
    return data[(data[column] > (mean_val - threshold * std_val)) & (data[column] < (mean_val + threshold * std_val))]


# Function to get LOESS smoothed values at specified x values
def get_loess_values_at_x(data, biomarker, status, x_values):
    subset = data[data['Disease_Status'] == status]
    x = subset['PSEUDOTIME_NORMALIZED']
    y = subset[biomarker]

    # Fit LOESS model
    loess_sm = sm.nonparametric.lowess(y, x, frac=0.3)

    # Interpolate the LOESS smoothed values at the specified x values
    loess_interpolated = np.interp(x_values, loess_sm[:, 0], loess_sm[:, 1])

    return loess_interpolated


# Specify x values
x_values = np.arange(0, 1.01, 0.01)

# Initialize a dictionary to store the results
results = {"Region": [], "Biomarker": [], "Status": [], "x": [], "y": []}

# Create the table for each biomarker and each region
for biomarker in biomarkers:
    for region, region_data in zip(["CA1", "Caudate"], [data_ca1, data_caudate]):
        for status in ['AD', 'Healthy']:
            # Get LOESS values
            loess_values = get_loess_values_at_x(region_data, biomarker, status, x_values)

            # Append the results to the dictionary
            results["Region"].extend([region] * len(x_values))
            results["Biomarker"].extend([biomarker] * len(x_values))
            results["Status"].extend([status] * len(x_values))
            results["x"].extend(x_values)
            results["y"].extend(loess_values)

# Convert the results dictionary to a DataFrame
results_df = pd.DataFrame(results)

# Save the results to a CSV file
output_file_path = '/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data/loess_values.csv'  # Update with your desired file path
 # Update with your desired file path
results_df.to_csv(output_file_path, index=False)
