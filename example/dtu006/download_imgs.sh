#!/bin/bash
set -e

if [ "$#" -ne 1 ]; then
    printf "Illegal number of parameters. Format is:\n\t$0 <.../imgs/>\n\n"
    exit 1
fi

IMGS_FOLDER=$1
DOWNLOAD_FOLDER="."
PATH_INSIDE_ZIP="SampleSet/MVS Data/Rectified/scan6"
FILENAME="SampleSet.zip"
FILEURL="http://roboimagedata2.compute.dtu.dk/data/MVS/SampleSet.zip"
DOWLOADED_FILE="$DOWNLOAD_FOLDER/$FILENAME"

printf "Downloading dtu006 imgs in $1\n"

# Download dataset

echo "Ready to download $FILENAME from $FILEURL in $DOWNLOAD_FOLDER. The download will require ~6.3GB disk space."

read -p "Do you wish to continue? [Y/n] " -n 1 -r
printf "\n"
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    exit 1
fi

mkdir -p "$IMGS_FOLDER"
mkdir -p "$DOWNLOAD_FOLDER"

wget -O "$DOWLOADED_FILE" "$FILEURL"

printf "\nExtracting dataset...\n"

 #!/bin/bash 
COUNTER=1
while [  $COUNTER -lt 50 ]; do
 IMGNAME="rect_$(printf %03d $COUNTER)_3_r5000.png"
 unzip -j "$DOWLOADED_FILE" "$PATH_INSIDE_ZIP/$IMGNAME" -d "$IMGS_FOLDER"
 # echo "$COUNTER : The filename is $IMGNAME"
 let COUNTER=COUNTER+2
done

# Renaming images

rename -e 's/rect_(\d{3})_3_r5000.png/sprintf("%04d.png",$1-1)/e' $IMGS_FOLDER/*.png

printf "\nDone.\n"

read -p "Do you wish to remove $DOWLOADED_FILE? [Y/n] " -n 1 -r
printf "\n"
if [[ $REPLY =~ ^[Yy]$ ]]
then
    rm "$DOWLOADED_FILE"
fi

