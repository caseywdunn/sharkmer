#!/bin/bash

# if there are zero arguments or more than two arguments, print usage and exit
if [ $# -eq 0 ] || [ $# -gt 2 ]; then
    echo "Usage: genomovie.sh <input_file> [output_dir]"
    exit 1
fi

input_file="$1"  # Input filename as an argument

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Input file does not exist"
    exit 1
fi

# If there is a second argument, use it as the output directory
if [ -n "$2" ]; then
    output_dir="$2"
else
    output_dir=$(basename "$input_file" | sed 's/\.[^.]*$//')
    output_dir="$output_dir.output"
fi

mkdir -p $output_dir

# Get the number of columns in the input file
num_columns=$(awk -F'\t' '{print NF}' "$input_file" | head -n1)
num_y_values=$((num_columns - 1))

# Loop through each y value and create the sample_i.histo files
for ((i = 1; i <= num_y_values; i++))
do
    # Pad the number with zeros if needed to get to three digits
    i_padded=$(printf "%04d" "$i")
    histo_file="$output_dir/sample_$i_padded.histo"
    cut -f1,$((i + 1)) "$input_file" > "$histo_file"
    
    # Get run_name from the input filename by removing the extension
    run_name=$(basename "$histo_file" | sed 's/\.[^.]*$//')

    # Run genomescope on the sample histo file
    genomescope2 -i "$histo_file" -o "$output_dir" -k 21 -n "$run_name"
    
done

# Create an animation using ffmpeg from the linear_plot.png files
movie_file = $(basename "$input_file" | sed 's/\.[^.]*$//')
movie_file = "$movie_file.mp4"
ffmpeg -framerate 20 -i "$output_dir"/sample_%04d_linear_plot.png -c:v libx264 -r 30 -pix_fmt yuv420p "$movie_file"