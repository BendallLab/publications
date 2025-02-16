import pandas as pd

# Load the CSV file generated by the previous script
input_file_path = '/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data/loess_values.csv'  # Update with your actual file path
data = pd.read_csv(input_file_path)

# Initialize a dictionary to store the results
difference_results = {"Region": [], "Biomarker": [], "x": [], "Difference_AD_vs_Healthy": []}

# Get the unique biomarkers and regions
biomarkers = data['Biomarker'].unique()
regions = data['Region'].unique()
x_values = data['x'].unique()

# Calculate the differences for each biomarker and region
for biomarker in biomarkers:
    for region in regions:
        for x in x_values:
            # Get the y values for AD and Healthy at the current x value
            y_ad = data[(data['Biomarker'] == biomarker) & (data['Region'] == region) & (data['Status'] == 'AD') & (
                        data['x'] == x)]['y'].values[0]
            y_healthy = data[
                (data['Biomarker'] == biomarker) & (data['Region'] == region) & (data['Status'] == 'Healthy') & (
                            data['x'] == x)]['y'].values[0]

            # Calculate the difference
            difference = y_ad - y_healthy

            # Append the results to the dictionary
            difference_results["Region"].append(region)
            difference_results["Biomarker"].append(biomarker)
            difference_results["x"].append(x)
            difference_results["Difference_AD_vs_Healthy"].append(difference)

# Convert the results dictionary to a DataFrame
difference_df = pd.DataFrame(difference_results)

# Save the results to a CSV file
output_file_path = '/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data/differences_ad_vs_healthy.csv'  # Update with your desired file path
difference_df.to_csv(output_file_path, index=False)

print(f"Differences saved to {output_file_path}")
