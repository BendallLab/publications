# Load the dataset
latest_file_path = "/mnt/data/microglia-DM_fixedV3_wAnucleated_table_size_transformed.csv"
latest_df = pd.read_csv(latest_file_path)

# Filter for rows where Disease_Status == "Healthy"
healthy_df = latest_df[latest_df["Disease_Status"] == "Healthy"]

# Subsample 10% of rows from each patient group
sampled_df = healthy_df.groupby('patient_label').apply(lambda x: x.sample(frac=0.1, random_state=42)).reset_index(drop=True)

# Select features for clustering
features_clustering = ["ApoE (pan)", "Biotin", "C1q", "CD44", "CD45", "CD68", "HLADR", "Iba1", "P2RY12", "TMEM119", "UbiquitinK48"]
subset_clustering = sampled_df[features_clustering]

# Perform hierarchical clustering on the subsampled data
linked_sampled = linkage(subset_clustering, method='ward')

# Plot the dendrogram for the clustered subsampled data
plt.figure(figsize=(10, 7))
dendrogram(linked_sampled,
           orientation='top',
           labels=subset_clustering.index,
           distance_sort='descending',
           show_leaf_counts=True)
plt.title('Hierarchical Clustering Dendrogram of Healthy Patients')
plt.xlabel('Sample Index')
plt.ylabel('Distance')
plt.show()


# Extract cluster centroids from the sampled healthy data hierarchical clustering
cluster_centers_healthy_hier = sampled_healthy_df.groupby('Meta_Cluster').mean()

# Using the centroids, predict the closest clusters for AD samples
# Define a function to find the closest cluster center
import numpy as np

def find_closest_cluster(centers, row):
    distances = np.linalg.norm(centers - row, axis=1)
    return np.argmin(distances) + 1  # +1 to match cluster labeling starting at 1

# Apply the function to each row in the AD subset
subset_ad_df['Meta_Cluster'] = subset_ad_df.apply(lambda row: find_closest_cluster(cluster_centers_healthy_hier, row), axis=1)

subset_ad_df.head()


# Concatenate
# Add the meta-cluster labels to the original dataframes for both Healthy and AD subsets
healthy_df['Meta_Cluster'] = sampled_healthy_df['Meta_Cluster']
ad_df['Meta_Cluster'] = subset_ad_df['Meta_Cluster']

# Combine the rows for "Healthy" and "AD" into a single DataFrame with their meta-cluster labels
combined_meta_df = pd.concat([healthy_df, ad_df], axis=0)

# Prepare for CSV export
output_combined_meta_path = "/mnt/data/combined_healthy_ad_meta_clusters.csv"
combined_meta_df.to_csv(output_combined_meta_path, index=False)

output_combined_meta_path


# Re-attempt the clustering on the sampled Healthy dataset

# Sample 10% of rows for each patient_label group again
sampled_healthy_df = healthy_only_df.groupby('patient_label').apply(lambda x: x.sample(frac=0.1, random_state=42)).reset_index(drop=True)

# Select the features for clustering
features_healthy_sampled = sampled_healthy_df[features_kmeans]

# Perform hierarchical clustering on this sample
linked_healthy_sampled = linkage(features_healthy_sampled, method='ward')

# Create 10 meta-clusters from the hierarchical clustering
meta_clusters_healthy_sampled = fcluster(linked_healthy_sampled, t=10, criterion='maxclust')

# Assign meta cluster labels to the sampled data
sampled_healthy_df['Meta_Cluster'] = meta_clusters_healthy_sampled

sampled_healthy_df.head()
