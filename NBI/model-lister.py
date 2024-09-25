import os
import json
import pandas as pd

# Function to parse ranking_debug.json file and retrieve necessary information
def parse_ranking_debug(folder_path):
    folder_name = os.path.basename(folder_path)
    json_file = os.path.join(folder_path, "ranking_debug.json")

    if not os.path.exists(json_file):
        return None

    with open(json_file, "r") as f:
        data = json.load(f)

    if not data.get("plddts") or not data.get("order"):
        return None

    # Get model with the highest PLDDT score
    highest_model = data["order"][0]
    highest_score = data["plddts"][highest_model]

    # Find corresponding PDB file
    pdb_file = f"relaxed_{highest_model}.pdb"
    pdb_path = os.path.join(folder_path, pdb_file)

    return folder_name, highest_model, highest_score, pdb_path

# Main function to search folders and create table
def create_table(root_dir):
    table_data = []

    # Iterate through folders in root directory
    for folder_name in os.listdir(root_dir):
        folder_path = os.path.join(root_dir, folder_name)
        if os.path.isdir(folder_path):
            result = parse_ranking_debug(folder_path)
            if result:
                table_data.append(result)

    # Create DataFrame
    df = pd.DataFrame(table_data, columns=["Folder Name", "Model Name", "PLDDT Score", "PDB Filename"])

    # Adjust display settings to show complete contents of each column
    pd.set_option('display.max_colwidth', None)
    pd.set_option('display.max_rows', None)

    return df

# Example usage
root_directory = "./"
result_table = create_table(root_directory)
print(result_table)
