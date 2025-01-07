import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your dataset
df = pd.read_csv('/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data/density_mglia_large_cohort.csv')

# Group by 'GW', 'region', and 'point_id' and count the total number of rows per point
grouped_df = df.groupby(['Region', 'Disease_Status', 'patient_label', 'point_id']).size().reset_index(name='total_rows_per_patient_point')

# Apply transformation to y-axis values (each FOV 0.16 mm2)
grouped_df['transformed_total_rows_per_patient_point'] = grouped_df['total_rows_per_patient_point'] / 0.16

# Define the order of disease status for the x-axis
disease_status_order = ['Healthy', 'AD']

# Set the color palette to ensure common regions have the same color
# Define custom colors for the regions
palette = {'Healthy': '#A6E002', 'AD': '#172869'}


# Plot the boxplot with narrower and taller layout, angled labels, and faceted by GW
plt.figure(figsize=(6, 12))  # Narrower and taller dimensions
sns.set(style="whitegrid")

# Use seaborn's FacetGrid to create faceted plots by 'GW' with the specified region order
g = sns.FacetGrid(grouped_df, col='Region', height=10, aspect=0.5)
g.map(sns.boxplot, 'Disease_Status', 'transformed_total_rows_per_patient_point', order=disease_status_order, palette=palette, linewidth=4)

# Adjust the layout: rename y-axis, center x-axis label, and rotate x-axis labels by 45 degrees
g.set_axis_labels('Disease Status', 'No. cells/mm2')
g.set_titles(col_template="Region: {col_name}")
g.set_xlabels('Disease Status', labelpad=10)

# Apply 45-degree rotation for x-axis labels
for ax in g.axes.flat:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

plt.tight_layout()

# Save the plot as a PDF file
plt.savefig('/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/plots/large_cohort_density_plot.pdf', format='pdf')

# Save the plot as a TIFF file
plt.savefig('/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/plots/large_cohort_density_plot.tif', format='tiff', dpi=300)

# Show the plot
plt.show()
