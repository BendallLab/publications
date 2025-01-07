import os
import pandas as pd
from skimage import io
import numpy as np


def normalize_images_by_controls_tma_channels(csv_path, tma_num_list, input_dir, output_dir):
    # Read the CSV file
    df = pd.read_csv(csv_path)

    # Check and create the output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Grab the fov dirs associated with the current TMA
    for tma_num in tma_num_list:
        tma_dirs = [d for d in os.listdir(input_dir) if ''.join(['TMA', str(tma_num), '_'])]
        # For each TMA, check to the path
        for fov in tma_dirs:
            fov_path = os.path.join(input_dir, fov)
            if os.path.exists(fov_path):
                # Iterate through each image in the fov
                for img_file in os.listdir(fov_path):
                    # Split the base name and the extension and return only the name
                    filename_parts = os.path.splitext(img_file)
                    if filename_parts[1] != '.tiff':
                        continue
                    channel = filename_parts[0]
                    # Load the image
                    image_path = os.path.join(fov_path, img_file)
                    image = io.imread(image_path)

                    # Normalize the image using the z-score from the associated Control TMA & Channel
                    filtered_df = df[(df['TMA'] == tma_num) & (df['channel'] == channel)]
                    z_score = filtered_df['avg_z_score'].values.tolist()[0]

                    if z_score < 0:
                        normalized_image = image * abs(z_score)
                    elif z_score > 0:
                        normalized_image = image * 1 / z_score
                    else:
                        normalized_image = image

                    # Save the normalized image
                    # Check and create the output fov directory
                    if not os.path.exists(os.path.join(output_dir, fov)):
                        os.makedirs(os.path.join(output_dir, fov))
                    output_image_path = os.path.join(output_dir, fov, img_file)
                    io.imsave(output_image_path, normalized_image.astype(np.float32))


def normalize_images_by_controls_tma(csv_path, tma_num_list, input_dir, output_dir):
    # Read the CSV file
    df = pd.read_csv(csv_path)

    # Check and create the output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Grab the fov dirs associated with the current TMA
    for tma_num in tma_num_list:
        tma_dirs = [d for d in os.listdir(input_dir) if ''.join(['TMA', str(tma_num), '_']) in d]
        # For each TMA, check to the path
        for fov in tma_dirs:
            fov_path = os.path.join(input_dir, fov)
            if os.path.exists(fov_path):
                # Iterate through each image in the fov
                for img_file in os.listdir(fov_path):
                    # Split the base name and the extension and return only the name
                    filename_parts = os.path.splitext(img_file)
                    if filename_parts[1] != '.tiff':
                        continue
                    channel = filename_parts[0]
                    # Load the image
                    image_path = os.path.join(fov_path, img_file)
                    image = io.imread(image_path)

                    # Normalize the image using the z-score from the associated Control TMA & Channel
                    filtered_df = df[(df['TMA'] == tma_num)]
                    normalize_factor = filtered_df['normalize_factor'].values.tolist()[0]
                    normalized_image =  image / normalize_factor

                    # Save the normalized image
                    # Check and create the output fov directory
                    if not os.path.exists(os.path.join(output_dir, fov)):
                        os.makedirs(os.path.join(output_dir, fov))
                    output_image_path = os.path.join(output_dir, fov, img_file)
                    io.imsave(output_image_path, normalized_image.astype(np.float32))


def normalize_tma_seven(tma_num_list, input_dir, output_dir):
    # Check and create the output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Grab the fov dirs associated with the current TMA
    for tma_num in tma_num_list:
        tma_dirs = [d for d in os.listdir(input_dir) if ''.join(['TMA', str(tma_num), '_']) in d]
        # For each TMA, check to the path
        for fov in tma_dirs:
            fov_path = os.path.join(input_dir, fov)
            if os.path.exists(fov_path):
                # Iterate through each image in the fov
                for img_file in os.listdir(fov_path):
                    # Split the base name and the extension and return only the name
                    filename_parts = os.path.splitext(img_file)
                    if filename_parts[1] != '.tiff':
                        continue
                    # Load the image
                    image_path = os.path.join(fov_path, img_file)
                    image = io.imread(image_path)

                    # Normalize the image using the z-score from the associated Control TMA & Channel
                    normalize_factor = 0.8
                    normalized_image = normalize_factor * image

                    # Save the normalized image
                    # Check and create the output fov directory
                    if not os.path.exists(os.path.join(output_dir, fov)):
                        os.makedirs(os.path.join(output_dir, fov))
                    output_image_path = os.path.join(output_dir, fov, img_file)
                    io.imsave(output_image_path, normalized_image.astype(np.float32))


# Renamed images usage
process_run = False
if process_run:
    input_img_path = '/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/image_data'
    output_img_path = '/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/image_data_hc_median_norm'
    path_to_csv = '/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My ' \
                  'Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data/control_medians.csv'
    tma_list = [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12]
    #tma_list = [9, 10, 11, 12]
    #tma_list = [7]

    #normalize_images_by_controls_tma_channels(path_to_csv, tma_list, input_img_path, output_img_path)
    normalize_images_by_controls_tma(path_to_csv, tma_list, input_img_path, output_img_path)
    #normalize_tma_seven(tma_list, input_img_path, output_img_path)


