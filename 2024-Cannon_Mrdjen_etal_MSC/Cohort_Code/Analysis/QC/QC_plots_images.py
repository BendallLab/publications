import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import zscore

def process_csv(file_path):
    # Read in the CSV file
    df = pd.read_csv(file_path)

    # Filter out rows where TMA = 0
    df_filtered = df[df['TMA'] != 0]
    df_filtered = df[df['Region'] != "HC-CA1" | "AD-CA1"]

    # Grouping the DataFrame by 'channel'
    grouped = df_filtered.groupby('channel')

    # Calculating z-score for the non-zero mean intensity for each channel
    df_filtered['z_score'] = grouped['Non-zero mean intensity'].transform(lambda x: zscore(x, ddof=1))

    # Combining TMA and Region into a single column
    df_filtered['TMA_Region_Core'] = df_filtered['TMA'].astype(str) + '_' + df_filtered['Region'].astype(str) + '_' + df_filtered['Core'].astype(str)

    # Averaging the z-scores for each TMA+Region combination
    pivot_zscore = df_filtered.pivot_table(index='channel', columns='TMA_Region_Core', values='z_score', aggfunc='mean', fill_value=0)

    # Plotting the heatmap for averaged z-scores
    #matplotlib.use('TkAgg')
    matplotlib.use('Qt5Agg')
    plt.figure(figsize=(15, 10))
    sns.heatmap(pivot_zscore, cmap='viridis', center=0, xticklabels=True, yticklabels=True)
    plt.title('Heatmap of Averaged Z-Scores per Channel and TMA+Region')
    plt.xlabel('TMA_Region')
    plt.ylabel('Channel')
    plt.xticks(rotation=90)
    plt.show()

# Replace 'your_file_path.csv' with the path to your CSV file
file_path = '/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/qc_tma_metrics/combined_csv_files/combined_nonzero_mean_stats.csv'
process_csv(file_path)
