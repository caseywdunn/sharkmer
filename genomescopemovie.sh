#!/bin/bash

# Parse the arguments
# -i is the input file
# -t is the number of threads to use, default is 1
# -o is the output directory, default is the input filename without the extension
# -k is the kmer size, default is 21
# -h prints the usage

# Parse arguments and print usage if invalid
while getopts ":i:t:o:k:h" opt; do
  case $opt in
    i) input_file="$OPTARG"
    ;;
    t) thread_count="$OPTARG"
    ;;
    o) output_dir="$OPTARG"
    ;;
    k) kmer_size="$OPTARG"
    ;;
    h) echo "Usage: genomescopemovie.sh -i input_file [-t thread_count] [-o output_dir] [-k kmer_size]"
       exit 0
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Input file does not exist"
    exit 1
fi

mkdir -p $output_dir

# Get the number of columns in the input file
num_columns=$(awk -F'\t' '{print NF}' "$input_file" | head -n1)
num_y_values=$((num_columns - 1))

echo "Number of y values: $num_y_values"

# Loop through each y value and create the sample_i.histo files
for ((i = 1; i <= num_y_values; i++))
do
    # Pad the number with zeros if needed to get to three digits
    i_padded=$(printf "%04d" "$i")
    histo_file="$output_dir/sample_$i_padded.histo"
    cut -f1,$((i + 1)) "$input_file" > "$histo_file"
    
    # Get run_name from the input filename by removing the extension
    run_name=$(basename "$histo_file" | sed 's/\.[^.]*$//')

    # If bash version is 4.0 or higher, use the built-in wait -n command, else just run genomescope2
    # This is to avoid the "wait -n" command not found error on older versions of bash, such as that on Mac OS
    if (( BASH_VERSINFO[0] >= 4 )); then
        # Run genomescope on the sample histo file
        genomescope2 -i "$histo_file" -o "$output_dir" -k $kmer_size -n "$run_name" &

        # if we have reached the maximum number of threads, wait for some processes to finish before starting new ones
        if (( $(jobs -p | wc -l) >= thread_count )); then
        wait -n
        fi
    else
        # Run genomescope on the sample histo file
        echo "genomescope2 -i $histo_file -o $output_dir -k $kmer_size -n $run_name"
        genomescope2 -i "$histo_file" -o "$output_dir" -k $kmer_size -n "$run_name"
    fi

done

# Wait for all background processes to finish before moving on to the next command
wait

# Summarize results in tsv

# Name of the output file
tsv_out="$output_dir/${run_name}_genomescope_stats.tsv"

# Create the output file and write the header to it
echo -e "Index\tFilename\tHomozygous (aa) min\tHomozygous (aa) max\tHeterozygous (ab) min\tHeterozygous (ab) max\tGenome Haploid Length min\tGenome Haploid Length max\tGenome Repeat Length min\tGenome Repeat Length max\tGenome Unique Length min\tGenome Unique Length max\tModel Fit min\tModel Fit max\tRead Error Rate min\tRead Error Rate max" > $tsv_out

# Go over all .txt files
for file_path in $output_dir/sample_*summary.txt
do

  echo "File path: $file_path"

  # Extract filename from full file path
  file=$(basename $file_path)

  # Extract the index from the filename
  index=${file/sample_/}
  index=${index/_summary.txt/}
  echo $index

  # Extract the values from the table using awk and write them to the output file
  awk -v filename="$file" -v idx="$index" '
  BEGIN{OFS="\t"}
  /Homozygous \(aa\)/{aamin=$3; aamax=$4}
  /Heterozygous \(ab\)/{bmin=$3; bmax=$4}
  /Genome Haploid Length/{hmin=$4; hmax=$6}
  /Genome Repeat Length/{rmin=$4; rmax=$6}
  /Genome Unique Length/{umin=$4; umax=$6}
  /Model Fit/{fmin=$3; fmax=$4}
  /Read Error Rate/{emin=$4; emax=$5}
  END{print idx, filename, aamin, aamax, bmin, bmax, hmin, hmax, rmin, rmax, umin, umax, fmin, fmax, emin, emax}' $file_path >> $tsv_out
done

sed -i 's/%//g; s/,//g' $tsv_out

# Make some movies
movie_base=$(basename "$input_file" | sed 's/\.[^.]*$//')

ffmpeg -framerate 20 -i \
  "$output_dir/sample_%04d_linear_plot.png" \
  -c:v libx264 -r 30 -pix_fmt yuv420p \
  "$output_dir/${movie_base}_linear_plot.mp4"

ffmpeg -framerate 20 -i \
  "$output_dir/sample_%04d_log_plot.png" \
  -c:v libx264 -r 30 -pix_fmt yuv420p \
  "$output_dir/${movie_base}_log_plot.mp4"
