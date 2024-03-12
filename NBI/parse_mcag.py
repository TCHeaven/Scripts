import pandas as pd
import sys

def merge_notes_xlsx(notes_file, xlsx_file, output_file):
    # Read the .xlsx file
  xlsx_df = pd.read_excel(xlsx_file, header=None, names=['ID', 'Name', 'Accession', 'Description', 'Species'])
    # Read the .notes file
  notes_df = pd.read_csv(notes_file, sep='\t', header=None, names=['ID', 'Gene_IDs'])
    # Split the 'Gene_IDs' column in the .notes file by '|' and explode it into multiple rows
  notes_df['Gene_IDs'] = notes_df['Gene_IDs'].str.split('|')
  notes_df = notes_df.explode('Gene_IDs')
    # Merge the dataframes based on the 'ID' column in the .xlsx file and the 'Gene_IDs' column in the .notes file
  merged_df = pd.merge(notes_df, xlsx_df, left_on='Gene_IDs', right_on='ID', how='left')
    # Concatenate the 'Description' column from the .xlsx file to the matched rows in the .notes file
  merged_df['Notes_with_Description'] = merged_df['ID'] + '\t' + merged_df['Gene_IDs'] + '\t' + merged_df['Description'] + '\t' + merged_df['Species']
    # Group the concatenated data back into a single row for each ID in the .notes file
  final_df = merged_df.groupby('ID')['Notes_with_Description'].apply(lambda x: '\t'.join(x)).reset_index()
    # Save the final result to the output file
  final_df.to_csv(output_file, sep='\t', header=False, index=False)

if __name__ == "__main__":
  if len(sys.argv) != 4:
    print("Usage: python script.py notes_file.xlsx xlsx_file.xlsx output_file.notes")
    sys.exit(1)
  notes_file = sys.argv[1]
  xlsx_file = sys.argv[2]
  output_file = sys.argv[3]

