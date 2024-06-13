import os
import pandas as pd
import matplotlib.pyplot as plt

output_directory = './output'
data_directory = './data'
benchmark_directory = './benchmarks'




# Get read stats
fastq_files = [f for f in os.listdir(data_directory) if f.endswith('.fastq')]
reads_count = {}
for f in fastq_files:
    file_path = os.path.join(data_directory, f)
    with open(file_path, 'r') as file:
        line_count = sum(1 for line in file)  # Count lines in the file
        reads_count[f.split('.')[0]] = line_count // 4  # Divide by 4 to get the number of reads

for species, reads in reads_count.items():
    print(f'{species}: {reads}')

# Get output stats
# files have naming convention:
# Xenia_sp_16M_cnidaria_16s.fasta
# Where Xenia_sp is the species name, 16M is the number of reads, cnidaria is the clade and 16s is the gene
# loop over the files in the output directory and populate these fields in a panda dataframe
output_files = [f for f in os.listdir(output_directory) if f.endswith('.fasta')]
output_stats = []
for f in output_files:
    parts = f.split('_')
    species = parts[0] + '_' + parts[1]  # Concatenate the first two parts to get the species name
    million_reads = int(parts[2][:-1])  # Remove the 'M' and convert to int
    clade = parts[3]
    gene = parts[4].split('.')[0]

    with open(os.path.join(output_directory, f), 'r') as file:
        seq_count = sum(1 for line in file if line.startswith('>'))

    output_stats.append([f, clade, gene, million_reads, species, seq_count])

columns = ['filename', 'clade', 'gene', 'million_reads', 'species',  'sequences']
output_df = pd.DataFrame(output_stats, columns=columns)
output_df = output_df.apply(pd.to_numeric, errors='ignore')
output_df.sort_values(['clade', 'gene', 'million_reads', 'species'], inplace=True)

# For each entry in fastq_files, remove any row from data where the million_reads*1000000 is greater than the number of record in the fastq file
for species, reads in reads_count.items():
    output_df = output_df[~((output_df['species'] == species) & (output_df['million_reads']*1000000 > reads))]


# save the DataFrame to a CSV file
output_df.to_csv('output_stats.csv', index=False)

print("Output stats")
# Print the full DataFrame
print(output_df)

# Ingest benchmarks
benchmarks = []
for filename in os.listdir(benchmark_directory):
    if filename.endswith(".benchmark.txt"):
        # Split the filename to extract the required parts
        parts = filename.split('.')
        species = parts[1]
        million_reads = int(parts[2][:-1])  # Remove the 'M' and convert to int

        # Read the content of the file
        with open(os.path.join(benchmark_directory, filename), 'r') as file:
            lines = file.readlines()
            # Parse the benchmarks from the second line (assuming first line is header)
            values = lines[1].split()
            
            # Append the benchmarks
            benchmarks.append([filename, species, million_reads] + values)

# Define the column names
columns = ['filename', 'species', 'million_reads', 's', 'h:m:s', 'max_rss', 'max_vms', 'max_uss', 'max_pss', 'io_in', 'io_out', 'mean_load', 'cpu_time']

# Create a pandas DataFrame
benchmark_df = pd.DataFrame(benchmarks, columns=columns)

# convert all columns to numeric
benchmark_df = benchmark_df.apply(pd.to_numeric, errors='ignore')

# Sort the DataFrame by sample and then by million_reads
benchmark_df.sort_values(['species', 'million_reads'], inplace=True)

# For each entry in fastq_files, remove any row from data where the million_reads*1000000 is greater than the number of record in the fastq file
for species, reads in reads_count.items():
    benchmark_df = benchmark_df[~((benchmark_df['species'] == species) & (benchmark_df['million_reads']*1000000 > reads))]

# save the DataFrame to a CSV file
benchmark_df.to_csv('benchmarks.csv', index=False)

# Plotting CPU time vs million reads
plt.figure(figsize=(12, 6))
for sample, group in benchmark_df.groupby('sample'):
    plt.plot(group['million_reads'], group['cpu_time'], marker='o', label=sample)
plt.xlabel('Million Reads')
plt.ylabel('CPU Time (s)')
plt.title('CPU Time vs Million Reads')
plt.legend()
plt.grid(True)
plt.savefig('cpu_time_vs_million_reads.png')

# Plotting max RSS vs million reads
plt.figure(figsize=(12, 6))
for sample, group in benchmark_df.groupby('sample'):
    plt.plot(group['million_reads'], group['max_rss'], marker='o', label=sample)
plt.xlabel('Million Reads')
plt.ylabel('Max RSS (MB)')
plt.title('Max RSS vs Million Reads')
plt.legend()
plt.grid(True)
plt.savefig('max_rss_vs_million_reads.png')
