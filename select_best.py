import sys
import pandas as pd
from pathlib import Path

# Input files (featureCounts outputs)
input_files = [Path(f) for f in sys.argv[1:-1]]
output_file = sys.argv[-1]

best_file = None
best_assigned = -1

for fc_file in input_files:
    summary_file = fc_file.with_suffix(fc_file.suffix + ".summary")  # add .summary after .txt

    if not summary_file.exists():
        print(f"Warning: Summary file not found for {fc_file}")
        continue

    # Read the summary file
    summary_df = pd.read_csv(summary_file, sep="\t", index_col=0, header=0)

    # Get 'Assigned' count
    assigned = summary_df.loc["Assigned"].values[0]

    print(f"{fc_file.name}: Assigned = {assigned}")

    if assigned > best_assigned:
        best_assigned = assigned
        best_file = fc_file

# Final decision
if best_file:
    best_df = pd.read_csv(best_file, sep="\t", comment="#")
    best_df.to_csv(output_file, sep="\t", index=False)
    print(f"\nSelected best file: {best_file.name} (Assigned: {best_assigned})")
    print(f"Saved to: {output_file}")
else:
    print("No valid summary files found. Nothing was written.")
