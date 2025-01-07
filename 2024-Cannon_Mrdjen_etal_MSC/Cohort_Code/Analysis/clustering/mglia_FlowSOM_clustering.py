#### Step 1: Install Required Library ####
#pip install flowsom

#### Step 2: Import Libraries ####
import pandas as pd
from flowsom import FlowSOM
from sklearn.preprocessing import MinMaxScaler

#### Step 3: Load Data ####
# Load data
df = pd.read_csv('path_to_your_data.csv')

# Filter data
healthy_df = df[df['Disease_Status'] == 'Healthy']
ad_df = df[df['Disease_Status'] == 'AD']

# Select features
features = ["ApoE (pan)", "Biotin", "C1q", "CD44", "CD45", "CD68", "HLADR", "Iba1", "P2RY12", "TMEM119", "UbiquitinK48"]
healthy_features = healthy_df[features]

#### Step 4: Normalize Data ####

scaler = MinMaxScaler()
healthy_features_scaled = scaler.fit_transform(healthy_features)

#### Step 5: Initialize and Train FlowSOM ####
# Initialize FlowSOM
fsom = FlowSOM(healthy_features_scaled, n_clusters=10, map_type=1)  # map_type=1 for toroid grid

# Train FlowSOM
fsom.train()

#### Step 6: Map Clusters to Meta Clusters ####
# Map clusters to meta clusters
meta_clusters = fsom.meta_clustering(n_clusters=10)
healthy_df['Meta_Cluster'] = meta_clusters


#### Step 7: Map AD Data ####
ad_features = ad_df[features]
ad_features_scaled = scaler.transform(ad_features)  # Use the same scaler
ad_meta_clusters = fsom.map_data(ad_features_scaled)
ad_df['Meta_Cluster'] = ad_meta_clusters


#### Step 8: Combine and Save Data ####
combined_df = pd.concat([healthy_df, ad_df], axis=0)
combined_df.to_csv('combined_clustered_data.csv', index=False)


