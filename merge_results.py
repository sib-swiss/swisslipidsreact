import pandas as pd
import glob
import os

# Step 1: Find all enumerated reaction TSV files
file_pattern = os.path.join("results_*", "*_enumerated_reactions.tsv")
files = glob.glob(file_pattern)

# Step 2: Read and concatenate all TSVs
dfs = []
for file in files:
    try:
        df = pd.read_csv(file, sep='\t')
        dfs.append(df)
    except Exception as e:
        print(f"Error reading {file}: {e}")

# Step 3: Merge into one DataFrame
if dfs:
    merged_df = pd.concat(dfs, ignore_index=True)

    # Step 4: Drop rows with missing Web-RInChIKey
    merged_df = merged_df.dropna(subset=['Web-RInChIKey'])

    # Step 5: Drop duplicates based on MASTER_ID and Web-RInChIKey
    merged_df = merged_df.drop_duplicates(subset=['MASTER_ID', 'Web-RInChIKey'])

    # Step 6: Save to file
    merged_df.to_csv("merged_enumerated_reactions.tsv", sep='\t', index=False)
    print(f"Saved merged file with {len(merged_df)} rows.")
else:
    print("No matching files found.")