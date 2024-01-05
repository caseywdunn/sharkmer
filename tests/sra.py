import argparse
import pysradb
import pandas as pd

# Set up the argument parser
parser = argparse.ArgumentParser(description="Run SRA query with specified organism.")
parser.add_argument("--organism", help="Specify the organism for the query", required=True)
parser.add_argument("--bases_per_spot_min", help="Minimum bases per spot", default=280, type=int)
parser.add_argument("--bases_per_spot_max", help="Maximum bases per spot", default=320, type=int)
parser.add_argument("--spots_min", help="Minimum", default=500000, type=int)
parser.add_argument("--spots_max", help="Maximum spots", default=10000000, type=int)
args = parser.parse_args()

# Create a SRAdb object
db = pysradb.SRAweb()

# Define the query with the organism specified at the command line
query = f'((((((("illumina"[Platform]) AND "wgs"[Strategy]) AND {args.organism}[Organism])) AND "paired"[Layout]) AND "genomic"[Source])) AND "random"[library selection]'

# Run the query and get the result as a pandas DataFrame
result_df = db.sra_metadata(query)

# Make sure run_total_bases and run_total_spots are numeric
result_df['run_total_bases'] = pd.to_numeric(result_df['run_total_bases'])
result_df['run_total_spots'] = pd.to_numeric(result_df['run_total_spots'])

result_df['run_bases_per_spot'] = result_df['run_total_bases'] / result_df['run_total_spots']

# Filter the results based on the bases per spot
result_df = result_df[(result_df['run_bases_per_spot'] >= args.bases_per_spot_min) & (result_df['run_bases_per_spot'] <= args.bases_per_spot_max)]

# Filter the results based on the number of spots
result_df = result_df[(result_df['run_total_spots'] >= args.spots_min) & (result_df['run_total_spots'] <= args.spots_max)]

# Sort the results by run_total_spots
result_df = result_df.sort_values(by=['run_total_spots'], ascending=True)


# Select a subset of columns to print
result_df_subset = result_df[['run_accession', 'organism_name', 'run_total_spots', 'study_title']]

# Print the result
print(result_df_subset)

# Write the results to a CSV file
result_df.to_csv(f"sra_{args.organism}.csv", index=False)