z#!/bin/bash
# This script converts an XYZ file to a GRD file using GMT, and convert to cf format
# How to use the script: 
    #./xyz2grd_jagurs.sh input_xyz output_grd
#
first_line=$(head -n 1 $1)
last_line=$(tail -n 1 $1)
# Extract the first column values from the first and last lines
lonmin=$(echo "$first_line" | awk '{print $1}')
lonmax=$(echo "$last_line" | awk '{print $1}')
latmin=$(echo "$first_line" | awk '{print $2}')
latmax=$(echo "$last_line" | awk '{print $2}')

# Calculate the difference between the first column values
# Calculate the increment as the absolute difference between the first two values of the first column
increment=$(awk 'NR==1{a=$1;next} NR==2{print $1-a;exit}' "$1")
echo "Resolution grid: $increment degrees"
# Use GMT to create the GRD file with the -I and -R options
gmt xyz2grd "$1" -G"$2"=cf -I"$increment" -R"$lonmin/$lonmax/$latmin/$latmax"
# Check if the GRD file was created successfully
if [ -f "$2" ]; then
    echo "GRD file created successfully: $2"
else
    echo "Failed to create GRD file"
fi