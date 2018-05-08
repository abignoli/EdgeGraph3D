#!/bin/bash
set -e

if [ "$#" -ne 1 ]; then
    printf "Illegal number of parameters. Format is:\n\t$0 <.../build/>\n\n"
    exit 1
fi

printf "\n\nThis script assumes that:\n\t./imgs/ is the images folder, containing already the input images\n\t./edges/ contains the edge images\n\t./working_folder/ is the directory where EdgeGraph3D can write its output\n\t./input.json is the input SfM data\n\t./output.json is the output json (to be overwritten if present)\n\t./output_transformed.json is the coordinate-system-transformed output json (to be overwritten if present)\n\t./output_transformed.ply is the output point cloud (to be overwritten if present)\n\n"

read -p "Do you wish to continue? [Y/n] " -n 1 -r
printf "\n"
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    exit 1
fi

BUILD_DIR=$1
EDGEGRAPH3D_BIN="$BUILD_DIR/EdgeGraph3D"
CST_BIN="$BUILD_DIR/coordinate_system_transform"
JTP_BIN="$BUILD_DIR/json_to_ply"

mkdir -p working_folder

# Run EdgeGraph3D, with debug images enabled
$EDGEGRAPH3D_BIN imgs/ edges/ working_folder/ input.json output.json

# Transform the coordinate system of generated JSON
# Note: the order of this and the previous operation is interchangeable
$CST_BIN output.json target_camera_poses.txt output_transformed.json

# Generate a colored point cloud representing the result
$JTP_BIN -i imgs/ output_transformed.json output_transformed.ply



