import os
import pandas as pd
import matplotlib.pyplot as plt

data_directory = './data'
benchmark_directory = './benchmarks'

fastq_files = [f for f in os.listdir(data_directory) if f.endswith('.fastq')]

reads_count = {}
for f in fastq_files:
    file_path = os.path.join(data_directory, f)
    with open(file_path, 'r') as file:
        line_count = sum(1 for line in file)  # Count lines in the file
        reads_count[f.split('.')[0]] = line_count // 4  # Divide by 4 to get the number of reads

print(reads_count)

# List to store the data
data = []

# Iterate over the files in the directory
for filename in os.listdir(benchmark_directory):
    if filename.endswith(".benchmark.txt"):
        # Split the filename to extract the required parts
        parts = filename.split('.')
        sample = parts[1]
        million_reads = int(parts[2][:-1])  # Remove the 'M' and convert to int

        # Read the content of the file
        with open(os.path.join(directory, filename), 'r') as file:
            lines = file.readlines()
            # Parse the data from the second line (assuming first line is header)
            values = lines[1].split()
            
            # Append the data
            data.append([filename, sample, million_reads] + values)

# Define the column names
columns = ['filename', 'sample', 'million_reads', 's', 'h:m:s', 'max_rss', 'max_vms', 'max_uss', 'max_pss', 'io_in', 'io_out', 'mean_load', 'cpu_time']

# Create a pandas DataFrame
df = pd.DataFrame(data, columns=columns)

# convert all columns to numeric
df = df.apply(pd.to_numeric, errors='ignore')

# Sort the DataFrame by sample and then by million_reads
df.sort_values(['sample', 'million_reads'], inplace=True)

# For each entry in fastq_files, remove any row from data where the million_reads*1000000 is greater than the number of record in the fastq file
for species, reads in reads_count.items():
    df = df[~((df['sample'] == species) & (df['million_reads']*1000000 > reads))]

# save the DataFrame to a CSV file
df.to_csv('benchmarks.csv', index=False)

# Plotting CPU time vs million reads
plt.figure(figsize=(12, 6))
for sample, group in df.groupby('sample'):
    plt.plot(group['million_reads'], group['cpu_time'], marker='o', label=sample)
plt.xlabel('Million Reads')
plt.ylabel('CPU Time (s)')
plt.title('CPU Time vs Million Reads')
plt.legend()
plt.grid(True)
plt.savefig('cpu_time_vs_million_reads.png')

# Plotting max RSS vs million reads
plt.figure(figsize=(12, 6))
for sample, group in df.groupby('sample'):
    plt.plot(group['million_reads'], group['max_rss'], marker='o', label=sample)
plt.xlabel('Million Reads')
plt.ylabel('Max RSS (MB)')
plt.title('Max RSS vs Million Reads')
plt.legend()
plt.grid(True)
plt.savefig('max_rss_vs_million_reads.png')
