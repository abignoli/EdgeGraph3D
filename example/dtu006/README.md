How to use?
===========

This example has been provided to show how to use EdgeGraph3D and related utilities on a standard dataset in literature.

Dataset 006 of MVS Data Set – 2014
==================================

The source data for this example had been gathered from dataset 006 of MVS Data Set – 2014, here refenced as dtu006, of the Technical University of Denmark (DTU). All credits for the images of dtu006 go to the authors:

*Rasmus Jensen, Anders Dahl, George Vogiatzis, Engin Tola, Henrik Aanaes, “Large Scale Multi-view Stereopsis Evaluation“, CVPR, 2014*

Data
====

imgs/
-----

Due to limitations imposed by the authors of the DTU MVS Data Set, I could not include here the images of dtu006. On Linux, you can download its images using the provided download_imgs.sh script as:

    <.../download_imgs.sh> <.../imgs/>

The script downloads the sample set (~6.3 GB), extracts the necessary images in <.../imgs/> and executes the required renaming. Finally it asks the user whether to keep the downloaded dataset that contained the images of interest, or remove it (it won't be needed for the rest of this example).

Manually, you can download the smaller reduced [sample set](http://roboimagedata2.compute.dtu.dk/data/MVS/SampleSet.zip) from http://roboimagedata.compute.dtu.dk/, extract dtu006 data, and place all odd numbered images (with light level 3) in the imgs/ folder, renaming them to follow the same format and zero-padding of the images provided in edges/. 

Although all images of dtu006 could be used, I chose to restrict computation on a reduced sample of them, to present the functionality of EdgeGraph3D in a setting more representative of sparse 3D reconstruction.

edges/
------

The images provided in the edges/ folder have been produced using the Edge Detection and Image Segmentation (EDISON) System (http://coewww.rutgers.edu/riul/research/code/EDISON/doc/overview.html) implementing the algorithm presented in: P. Meer, B. Georgescu: "Edge detection with embedded confidence." IEEE Trans. Pattern Anal. Machine Intell., 23, 1351-1365, December 2001. Any edge detection algorithm can be used to produce the images to be placed in the edges/ folder, with preference for methods capable of producing accurate edges. On the same note, segment based algorithms such as LSD are strongly discouraged, as they usually represent curved edges very poorly, or not at all.

The color of the pixels in the images in edges/ should be black (not-edge) or white (edge), and non-maxima suppression should be enabled (although EdgeGraph3D has been tested to be resilient to this).

NOTE: The names of the images in the imgs/ folder, the edges/ folder and the names referenced in the SfM data JSON file should match.

SfM Data
--------

EdgeGraph3D receives as input SfM data in a JSON file, following the format used by [OpenMVG](https://github.com/openMVG/openMVG). For this example, a SfMData file, "input.json", has already been provided. This file was computed using [OpenMVG](https://github.com/openMVG/openMVG) (*Pierre Moulon and Pascal Monasse and Renaud Marlet and Others, https://github.com/openMVG/openMVG*).

target_camera_poses.txt
-----------------------

A file containing the (x,y,z) position of the each camera in some target coordinate system. This is used by the utility coordinate_system_transform to transform the coordinate system in an input SfM data (which potentially contain already the edge-points produced by EdgeGraph3D) so that the camera positions in the output file match the target camera positions.

Sample usage
============

An automated sample usage is also provided in sample_usage.sh, that can be executed as:

    <.../sample_usage.sh> <.../build/>
    
where <.../build/> contains the binaries of EdgeGraph3D, coordinate_system_transform and json_to_ply. The script runs EdgeGraph3D, executes coordinate_system_transform, and produces a colored point cloud with json_to_ply.

Individual binaries
===================

EdgeGraph3D
-----------

You can run EdgeGraph3D with the following command:

    <EdgeGraph3D-binary> [-i] <.../imgs/> <.../edges/> <.../working_folder/> <.../input.json> <.../output.json>
    
For example, assuming the shell is in this directory:

    mkdir working_folder
    <EdgeGraph3D-binary> -i imgs/ edges/ working_folder/ input.json output.json
    
The [-i] option enables the output of images that allow to observe in further detail PolyLine Graphs, polyline matching, output and others.

filter
------

EdgeGraph3D runs some outliers filtering on reference points and computed edge points based on some heuristics. It is possible to re-run the filtering on any JSON SfM data using the filter utility provided with EdgeGraph3D. Also note that by default EdgeGraph3D saves the intermediate SfM data before filtering in <.../working_folder/>before_filtering.json. This allows the user to change the filtering to be applied to the point cloud with custom parameters.

For instance, on the dtu006 dataset, a user may want to re-run the filtering process on the pre-filtering SfMData produced by EdgeGraph3D, for example, allowing a maximum value of 3 for Gauss-Newton reprojection mean squared error, and a minimum observation count of 5 for the edge-points produced by EdgeGraph3D. Considering there are 6268 reference points originally in input.json, and that, in this case, we'd like to run observation count filtering just on the generated edge-points, the command to run would be:

    <filter-binary> -e 3 -s 6268 -f 4 <.../working_folder/before_filtering.json> <.../output.json>
    
coordinate_system_transform
---------------------------

To compare different point clouds representing the same object it is necessary to put them in the same coordinate system. This is usually needed to compare a point cloud with the ground truth, if provided, associated to any given dataset. The coordinate_system_transform utility serves this purpose, through the target camera positions file, as explained earlier in this document.

The utility can be launched through the command:

    <coordinate_system_transform-binary> <.../input.json> <.../target_camera_poses.txt> <.../output.json>
    
json_to_ply
-----------

The json_to_ply utility is used to convert a SfM data JSON file to a point cloud in PLY format. The utility can produce a colored point cloud, if source images are provided. To run the utility the following command can be executed:

    <json_to_ply-binary> [-i <.../imgs/>] <.../input.json> <.../output.ply>
