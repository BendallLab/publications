# Re-import necessary libraries for clustering
from sklearn.cluster import KMeans

# Reload the dataset to ensure consistency
new_uploaded_file_path = "/mnt/data/microglia-DM_fixedV3_wAnucleated_table_size_transformed.csv"
new_uploaded_df = pd.read_csv(new_uploaded_file_path)

# Filter for rows where Disease_Status == "Healthy"
healthy_df = new_uploaded_df[new_uploaded_df['Disease_Status'] == "Healthy"]

# Select the features for K-means clustering
features_kmeans = ["ApoE (pan)", "Biotin", "C1q", "CD44", "CD45", "CD68", "HLADR", "Iba1", "P2RY12", "TMEM119", "UbiquitinK48"]
subset_healthy_df = healthy_df[features_kmeans]

# Perform K-means clustering on the filtered data
kmeans_healthy = KMeans(n_clusters=10, random_state=42)
kmeans_healthy.fit(subset_healthy_df)

# Get the cluster labels
cluster_labels_healthy = kmeans_healthy.labels_
cluster_centers_healthy = kmeans_healthy.cluster_centers_

# Add the cluster labels to the filtered DataFrame
healthy_df["Cluster"] = cluster_labels_healthy

healthy_df.head(), cluster_centers_healthy

# Filter for rows where Disease_Status == "AD"
ad_df = new_uploaded_df[new_uploaded_df['Disease_Status'] == "AD"]

# Select the same features for rows with AD
subset_ad_df = ad_df[features_kmeans]

# Use the cluster centers from the Healthy dataset to assign the AD samples to the closest clusters
# We will use the KMeans model that was trained on the healthy dataset, but predict on the AD subset
cluster_labels_ad = kmeans_healthy.predict(subset_ad_df)

# Add the cluster labels to the AD DataFrame
ad_df["Cluster"] = cluster_labels_ad

ad_df.head()

# Combine the rows for "Healthy" and "AD" into a single DataFrame with their cluster labels
combined_df = pd.concat([healthy_df, ad_df], axis=0)

# Prepare for CSV export
output_combined_path = "/mnt/data/combined_healthy_ad_clusters.csv"
combined_df.to_csv(output_combined_path, index=False)

output_combined_path
