import pandas as pd
import sys

def convert_xlsx_to_tsv(input_file, output_file):
  try:
    # Read the .xlsx file
    df = pd.read_excel(input_file)
    # Write the DataFrame to a .tsv file
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Conversion successful. TSV file saved as {output_file}")
  except Exception as e:
    print(f"Error converting file: {e}")

if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Usage: python convert_xlsx_to_tsv.py input_file.xlsx output_file.tsv")
  else:
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    convert_xlsx_to_tsv(input_file, output_file)
