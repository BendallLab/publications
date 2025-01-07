# Filtering files based on their names and grouping them into three categories
total_intensity_files = [f for f in extracted_files if f.endswith("total_intensity_stats.csv")]
nonzero_mean_files = [f for f in extracted_files if f.endswith("nonzero_mean_stats.csv")]
percentile_99_9_files = [f for f in extracted_files if f.endswith("percentile_99_9_stats.csv")]

# Function to read and append TMA identifier to each file
def read_and_label_files(file_list, folder_path):
    df_list = []
    for file in file_list:
        tma_name = file.split('_')[0]  # Extracting TMA name from the file name
        df = pd.read_csv(os.path.join(folder_path, file))
        df['TMA'] = tma_name  # Adding the TMA identifier as a new column
        df_list.append(df)
    return pd.concat(df_list, ignore_index=True)

# Reading and combining files in each category
total_intensity_combined = read_and_label_files(total_intensity_files, extracted_folder_path)
nonzero_mean_combined = read_and_label_files(nonzero_mean_files, extracted_folder_path)
percentile_99_9_combined = read_and_label_files(percentile_99_9_files, extracted_folder_path)

# Saving the combined dataframes to new CSV files
total_intensity_combined_path = '/mnt/data/combined_total_intensity_stats.csv'
nonzero_mean_combined_path = '/mnt/data/combined_nonzero_mean_stats.csv'
percentile_99_9_combined_path = '/mnt/data/combined_percentile_99_9_stats.csv'

total_intensity_combined.to_csv(total_intensity_combined_path, index=False)
nonzero_mean_combined.to_csv(nonzero_mean_combined_path, index=False)
percentile_99_9_combined.to_csv(percentile_99_9_combined_path, index=False)

(total_intensity_combined_path, nonzero_mean_combined_path, percentile_99_9_combined_path)
