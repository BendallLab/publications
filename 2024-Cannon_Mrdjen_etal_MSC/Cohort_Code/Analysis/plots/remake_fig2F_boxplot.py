import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your dataset
df = pd.read_csv('/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data/anat_object_connect.csv')

# Group by 'GW', 'region', and 'point_id' and count the total number of rows per point
grouped_df = df.groupby(['GW', 'region', 'point_id']).size().reset_index(name='total_rows_per_point')

# Remove 'HIP_AD' region
filtered_df = grouped_df[grouped_df['region'] != 'HIP_AD']

# Apply transformation to y-axis values
filtered_df['transformed_total_rows_per_point'] = filtered_df['total_rows_per_point'] / 0.49

# Define the order of regions for the x-axis
region_order = ['HIP', 'SN', 'MFG', 'Caudate', 'Cerebellum']

# Set the color palette to ensure common regions have the same color
# Define custom colors for the regions
palette = {'HIP': '#a6e000', 'SN': '#6c6c9d', 'MFG': '#1bb6af', 'Caudate': '#c70e7b', 'Cerebellum': '#fc6882'}


# Plot the boxplot with narrower and taller layout, angled labels, and faceted by GW
plt.figure(figsize=(6, 12))  # Narrower and taller dimensions
sns.set(style="whitegrid")

# Use seaborn's FacetGrid to create faceted plots by 'GW' with the specified region order
g = sns.FacetGrid(filtered_df, col='GW', height=10, aspect=0.5)
g.map(sns.boxplot, 'region', 'transformed_total_rows_per_point', order=region_order, palette=palette, linewidth=4)

# Adjust the layout: rename y-axis, center x-axis label, and rotate x-axis labels by 45 degrees
g.set_axis_labels('Region', 'No. cells/mm2')
g.set_titles(col_template="GW: {col_name}")
g.set_xlabels('Region', labelpad=10)

# Apply 45-degree rotation for x-axis labels
for ax in g.axes.flat:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

plt.tight_layout()

# Save the plot as a PDF file
plt.savefig('/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/plots/Fig2F_microglial_density_plot.pdf', format='pdf')

# Save the plot as a TIFF file
plt.savefig('/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/plots/Fig2F_microglial_density_plot.tif', format='tiff', dpi=300)

# Show the plot
plt.show()
