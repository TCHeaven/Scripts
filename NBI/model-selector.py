import os
import json

def find_highest_scoring_models(root_folder):
    highest_scoring_models = {}

    for subdir, _, _ in os.walk(root_folder):
        json_file = os.path.join(subdir, 'ranking_debug.json')
        if os.path.isfile(json_file):
            with open(json_file) as f:
                data = json.load(f)
            
            highest_score_model = max(data['plddts'], key=data['plddts'].get)
            highest_scoring_models[subdir] = highest_score_model

    return highest_scoring_models

if __name__ == "__main__":
    root_folder = "/Volumes/Saskia-Hogenhout/sam_mugford/alphafold"  # Change this to the path of your root folder
    highest_scoring_models = find_highest_scoring_models(root_folder)
    for subdir, highest_score_model in highest_scoring_models.items():
        print(f"Subdirectory: {subdir}, Highest Scoring Model: {highest_score_model}")
