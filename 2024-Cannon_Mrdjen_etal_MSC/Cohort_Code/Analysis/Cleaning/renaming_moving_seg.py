import os
import shutil
import re


def organize_files_by_subfolder(base_folder):
    """
    Organize files in a folder by creating subfolders based on file name parts
    and moving the files into the subfolders.
    """
    # Compile the regex to match the desired pattern
    pattern = re.compile(r"^(.*?)_(microglia-DM\.tiff)$")

    # Create base folder if it doesn't exist
    if not os.path.exists(base_folder):
        os.makedirs(base_folder)

    # Iterate over all files in the base folder
    for filename in os.listdir(base_folder):
        # Only process files (ignore subdirectories)
        if os.path.isfile(os.path.join(base_folder, filename)):
            match = pattern.match(filename)
            if match:
                part1, part2 = match.groups()

                # Create the subfolder if it doesn't exist
                subfolder_path = os.path.join(base_folder, part1)
                if not os.path.exists(subfolder_path):
                    os.makedirs(subfolder_path)

                # Move the file into the subfolder
                source_file = os.path.join(base_folder, filename)
                destination_file = os.path.join(subfolder_path, part2)
                shutil.move(source_file, destination_file)
                print(f"Moved {source_file} to {destination_file}")


# Example usage
if __name__ == "__main__":
    base_folder_path = "/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/segmentation/final_mask_dir_DM_2"  # Replace with the actual folder path
    organize_files_by_subfolder(base_folder_path)

