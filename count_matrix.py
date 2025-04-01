# import pandas as pd
# import sys

# # Check for correct usage
# if len(sys.argv) != 3:
#     print("Usage: python make_count_matrix.py input_file.txt output_file.csv")
#     sys.exit(1)

# input_file = sys.argv[1]
# output_file = sys.argv[2]

# # Read file, skipping lines starting with '#'
# df = pd.read_csv(input_file, sep="\t", comment="#")

# # Keep only Geneid and the last column (sample counts)
# counts = df.iloc[:, [0, -1]]

# # Rename columns
# counts.columns = ['Geneid', 'Sample']

# # Set Geneid as index
# counts.set_index('Geneid', inplace=True)

# # Save to CSV
# counts.to_csv(output_file)

# print(f"Count matrix saved to: {output_file}")

import pandas as pd
import sys
from pathlib import Path

# Check for correct usage
if len(sys.argv) < 3:
    print("Usage: python count_matrix.py input1.txt input2.txt ... output_matrix.txt")
    sys.exit(1)

input_files = sys.argv[1:-1]
output_file = sys.argv[-1]

# Initialize list to store data frames
count_dfs = []

for file_path in input_files:
    df = pd.read_csv(file_path, sep="\t", comment="#")
    gene_col = df.iloc[:, 0]  # Geneid
    count_col = df.iloc[:, -1]  # Count column (last)

    # Use the filename stem (e.g., "Collibri_A") as sample name
    sample_name = Path(file_path).stem.replace("Collibri_", "")
    
    # Create a temporary DataFrame with counts
    sample_df = pd.DataFrame({
        "Geneid": gene_col,
        sample_name: count_col
    })

    count_dfs.append(sample_df)

# Merge all on Geneid
merged_df = count_dfs[0]
for df in count_dfs[1:]:
    merged_df = pd.merge(merged_df, df, on="Geneid", how="outer")

# Set Geneid as index
merged_df.set_index("Geneid", inplace=True)

# Save to output file
merged_df.to_csv(output_file, sep="\t")
print(f"Count matrix written to: {output_file}")
