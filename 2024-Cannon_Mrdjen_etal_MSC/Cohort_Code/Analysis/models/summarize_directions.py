import pandas as pd

# Load the CSV file generated by the previous script
input_file_path = '/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data/differences_ad_vs_healthy.csv' # Update with your actual file path
data = pd.read_csv(input_file_path)

# Initialize a dictionary to store the summarized results
summary_results = {"Region": [], "Biomarker": [], "Group": [], "Total_Difference": [], "Direction": []}

# Define the x-value ranges for each group
groups = {
    "low": (0, 0.33),
    "medium": (0.34, 0.66),
    "high": (0.67, 1.0)
}

# Get the unique biomarkers and regions
biomarkers = data['Biomarker'].unique()
regions = data['Region'].unique()

# Summarize the differences for each biomarker and region
for biomarker in biomarkers:
    for region in regions:
        for group, (x_min, x_max) in groups.items():
            # Filter the data for the current biomarker, region, and x-value range
            filtered_data = data[
                (data['Biomarker'] == biomarker) & (data['Region'] == region) & (data['x'] >= x_min) & (
                            data['x'] <= x_max)]

            # Calculate the total difference for the current group
            total_difference = filtered_data['Difference_AD_vs_Healthy'].sum()

            # Determine the direction of the difference
            direction = 1 if total_difference > 0 else -1

            # Append the results to the dictionary
            summary_results["Region"].append(region)
            summary_results["Biomarker"].append(biomarker)
            summary_results["Group"].append(group)
            summary_results["Total_Difference"].append(total_difference)
            summary_results["Direction"].append(direction)

# Convert the summary results dictionary to a DataFrame
summary_df = pd.DataFrame(summary_results)

# Save the summarized results to a CSV file
output_file_path = '/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data/summarized_differences_with_direction.csv'  # Update with your desired file path
summary_df.to_csv(output_file_path, index=False)

print(f"Summarized differences with direction saved to {output_file_path}")
