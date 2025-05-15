import os
import pandas as pd

def import_alignment_summaries(path):
    # Initialize an empty dictionary to store the dataframes
    data_dict = {}

    # Loop through all files in the directory
    for filename in os.listdir(path):
        # Check if the file matches the format we're looking for
        if filename.endswith('_alignment_summary.txt'):
            # Extract the name from the filename
            name = filename.replace('_alignment_summary.txt', '')
            # Read the file into a dataframe
            filepath = os.path.join(path, filename)
            df = pd.read_csv(filepath, sep='\t')
            # Add the dataframe to the dictionary with the name as the key
            data_dict[name] = df
    return data_dict

def import_dists():
    dists = pd.read_hdf('/tf/primo/data/targets/query_target_dists_new.h5')
    return dists