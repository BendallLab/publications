BiocManager::install("FlowSOM")
BiocManager::install("flowCore")
BiocManager::install("pheatmap")


# Load the necessary libraries
library(FlowSOM)
library(flowCore)
library(pheatmap)

#set seed
seed <- 123 
set.seed(seed)

channels_MSC = c("ApoE (pan)", "Biotin", "C1q", "CD44", "CD45", "CD68", "HLADR", "Iba1", "P2RY12", "TMEM119", "UbiquitinK48")

#### Clustering ####
#### FlowSOM ########################################################################################################################################################################

# Load your data
fsom_data <- as.data.frame(mglia_to_save)

# Normalize each feature to the 99th percentile & Select only the relevant features
selected_features <- c("ApoE (pan)", "Biotin", "C1q", "CD44", "CD45", "CD68", "HLADR", "Iba1", "P2RY12", "TMEM119", "UbiquitinK48")

normalize_99 <- function(x) {
  x / quantile(x, 0.99, na.rm = TRUE)
}

fsom_data[selected_features] <- apply(fsom_data[selected_features], 2, normalize_99)

feature_data = fsom_data %>% filter(Disease_Status == "Healthy") %>%
  group_by(patient_label) %>%
  sample_frac(0.5) %>% ungroup() %>%
  select(selected_features)

# Convert the data frame to a FlowFrame
flow_frame <- flowCore::flowFrame(as.matrix(feature_data))

# Initialize a FlowSOM object with the FlowFrame
flowSOM_obj <- FlowSOM(flow_frame)


# Map clusters onto rest of data
flowSOM_obj_all <- NewData(flowSOM_obj, input = flowCore::flowFrame(as.matrix(fsom_data[, selected_features])))

# Apply metaclustering
metacl <- MetaClustering(flowSOM_obj_all$map$codes,
                         method = "metaClustering_kmeans",
                         max=10)

# Get metaclustering per cell
flowSOM.clustering <- metacl[flowSOM_obj_all$map$mapping[,1]]
fsom_data$Cluster = as.factor(flowSOM.clustering)
fsom_data$patient_label = droplevels(as.factor(fsom_data$patient_label))



#### Seurat - Leiden TBD#########################################################################################################################################################################

# Install necessary packages
BiocManager::install("Seurat")

# Load Seurat library
library
library(Seurat)

# Load the data
data <- as.data.frame(mglia_to_save)

# Select the features for clustering
selected_features <- c("ApoE (pan)", "Biotin", "C1q", "CD44", "CD45", "CD68", "HLADR", "Iba1", "P2RY12", "TMEM119", "UbiquitinK48")

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = t(data[selected_features]))

# Normalize the data (optional, depending on the data)
seurat_obj <- NormalizeData(seurat_obj)

# Scale the data
seurat_obj <- ScaleData(seurat_obj, features = selected_features)

# Perform PCA (Principal Component Analysis)
seurat_obj <- RunPCA(seurat_obj, features = selected_features)

# Choose the number of principal components to use
num_pcs <- 10  # Adjust based on your data and variance explained

# Create a neighborhood graph for clustering using PCA
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:num_pcs)

# Perform Leiden clustering
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, algorithm = 4)  # 4 is the Leiden algorithm

# Add cluster information back to the data
cluster_ids <- Idents(seurat_obj)
data$Cluster <- cluster_ids

# Print the data with cluster assignments
print(data)

# Optional: visualize the clusters using PCA, t-SNE, or UMAP
# PCA plot
DimPlot(seurat_obj, reduction = "pca")

# t-SNE plot
seurat_obj <- RunTSNE(seurat_obj, dims = 1:num_pcs)
TSNEPlot(seurat_obj)

# UMAP plot
seurat_obj <- RunUMAP(seurat_obj, dims = 1:num_pcs)
UMAPPlot(seurat_obj)


#### kMeans ########################################################################################################################################################################

data_path = "/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data"
name_kMeans_file = "mglia_clusters_kMeans.csv"
mglia_kMeans_clusters = data.table::fread(file=paste(data_path, name_kMeans_file, sep = "/"))
mglia_kMeans_clusters$Cluster = as.factor(mglia_kMeans_clusters$Cluster)
mglia_kMeans_clusters$patient_label = droplevels(as.factor(mglia_kMeans_clusters$patient_label))


#### hClust ########################################################################################################################################################################
data_path = "/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data"
name_hClust_file = "mglia_clusters_hClust.csv"
mglia_hClust_clusters = data.table::fread(file=paste(data_path, name_hClust_file, sep = "/"))
mglia_hClust_clusters$Cluster = as.factor(mglia_hClust_clusters$Cluster)
mglia_hClust_clusters$patient_label = droplevels(as.factor(mglia_hClust_clusters$patient_label))


#### Plotting ####

#### Patient Cluster distro ####

#### FlowSOM ####

  # Cluster distribution plot - per patient
  ggplot(fsom_data) +
    aes(x = patient_label, fill = Cluster) +
    geom_bar(position = "fill") +
    scale_fill_hue(direction = 1) +
    scale_color_hue(direction = 1) +
    theme_minimal() +
    facet_wrap(vars(Disease_Status, Region), scales = "free")
  
  # Cluster distribution plot - summed
  ggplot(fsom_data) +
    aes(x = Disease_Status, fill = Cluster) +
    geom_bar(position = "fill") +
    scale_fill_hue(direction = 1) +
    scale_color_hue(direction = 1) +
    theme_minimal() +
    facet_wrap(vars(Region), scales = "free")

#### kMeans ####

  # Cluster distribution plot - per patient
  ggplot(mglia_kMeans_clusters) +
    aes(x = patient_label, fill = Cluster) +
    geom_bar(position = "fill") +
    scale_fill_hue(direction = 1) +
    scale_color_hue(direction = 1) +
    theme_minimal() +
    facet_wrap(vars(Disease_Status, Region), scales = "free")

  # Cluster distribution plot - summed
  ggplot(mglia_kMeans_clusters) +
    aes(x = Disease_Status, fill = Cluster) +
    geom_bar(position = "fill") +
    scale_fill_hue(direction = 1) +
    scale_color_hue(direction = 1) +
    theme_minimal() +
    facet_wrap(vars(Region), scales = "free")


#### hClust ####

  # Cluster distribution plot - per patient
  ggplot(mglia_hClust_clusters) +
    aes(x = patient_label, fill = Cluster) +
    geom_bar(position = "fill") +
    scale_fill_hue(direction = 1) +
    scale_color_hue(direction = 1) +
    theme_minimal() +
    facet_wrap(vars(Disease_Status, Region), scales = "free")
  
  # Cluster distribution plot - summed
  ggplot(mglia_hClust_clusters) +
    aes(x = Disease_Status, fill = Cluster) +
    geom_bar(position = "fill") +
    scale_fill_hue(direction = 1) +
    scale_color_hue(direction = 1) +
    theme_minimal() +
    facet_wrap(vars(Region), scales = "free")

  
#### Cluster heatmaps ####

#### FlowSOM ####
  
  # Heatmap of clusters
  fsom_data_pheat = fsom_data %>% 
    select(c(channels_to_check, Cluster, patient_label, Disease_Status, Region)) %>% 
    filter(Disease_Status == "Healthy") %>%
    group_by(Cluster) %>%
    summarise(across(all_of(channels_MSC), median))
  # Assign row and column names
  rownames(fsom_data_pheat) <- fsom_data_pheat$Cluster # Gene names from A to E
  
  # Normalizing rows
  # Data already normalized
  
  # Creating the heatmap
  pheatmap(fsom_data_pheat[-1],
           color = colorRampPalette(viridis::magma(256))(256),
           show_rownames = TRUE,
           show_colnames = TRUE,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           scale = "column", # Data is already normalized
           main = "Heatmap of Normalized Gene Expression by Cluster")
  

#### kMeans ####
  
  # Heatmap of clusters
  mglia_kMeans_clusters_pheat = mglia_kMeans_clusters %>% 
    select(c(channels_to_check, Cluster, patient_label, Disease_Status, Region)) %>% 
    filter(Disease_Status == "Healthy") %>%
    group_by(Cluster) %>%
    summarise(across(all_of(channels_MSC), median))
  # Assign row and column names
  rownames(mglia_kMeans_clusters_pheat) <- mglia_kMeans_clusters_pheat$Cluster # Gene names from A to E
  
  # Normalizing rows
  # Data already normalized
  
  # Creating the heatmap
  pheatmap(mglia_kMeans_clusters_pheat[-1],
           color = colorRampPalette(viridis::magma(256))(256),
           show_rownames = TRUE,
           show_colnames = TRUE,
           labels_row =  rownames(mglia_kMeans_clusters_pheat),
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           scale = "column", # Data is already normalized
           main = "Heatmap of Normalized Gene Expression by Cluster")
  

#### hClust ####

  # Heatmap of clusters
  mglia_kMeans_hClust_pheat = mglia_hClust_clusters %>% 
    select(c(channels_to_check, Cluster, patient_label, Disease_Status, Region)) %>% 
    filter(Disease_Status == "Healthy") %>%
    group_by(Cluster) %>%
    summarise(across(all_of(channels_MSC), median))
  # Assign row and column names
  rownames(mglia_kMeans_hClust_pheat) <- mglia_kMeans_hClust_pheat$Cluster # Gene names from A to E
  
  # Normalizing rows
  # Data already normalized
  
  # Creating the heatmap
  pheatmap(mglia_kMeans_hClust_pheat[-1],
           color = colorRampPalette(viridis::magma(256))(256),
           show_rownames = TRUE,
           show_colnames = TRUE,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           scale = "column", # Data is already normalized
           main = "Heatmap of Normalized Gene Expression by Cluster")
  
#### General biaxials ####

  ggplot(mglia_to_save) +
   aes(x = Iba1, y = HLADR, colour = Region) +
   geom_point(shape = "circle", size = 1.5) +
   scale_x_continuous(trans = "asn") +
   scale_y_continuous(trans = "asn") +
   scale_color_hue(direction = 1) +
   theme_minimal() +
   facet_wrap(vars(Disease_Status, Region))
  
