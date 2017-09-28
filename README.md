EdgeGraph3D
===========

EdgeGraph3D is a system for multi-view stereo reconstruction of 3D edges. The system receives as input:

* a set of image observations of a scene or object
* a set of estimated camera poses and an initial cloud of 3D points, such as the one recovered by a Structure-from-Motion pipeline. EdgeGraph3D supports as input a JSON in the format of the output of the OpenMVG pipeline (https://github.com/openMVG/openMVG/)
* a set of edge images, one per input view, in which 2D edges detected on the corresponding image are represented by white pixels over a black background. This allows the user to integrate any 2D edge detection algorithm with EdgeGraph3D

and generate a set of 3D edges reconstructed in the observed scene. EdgeGraph3D can also generate a sampling of the recovered edges and integrate the input JSON with a cloud of 3D edge-points. This allows the user to use the output of our system in standard reconstruction algorithms that accept as input a point cloud, to produce more accurate 3D models.


Release
=======

The C++ code of EdgeGraph3D will be released in the next months.
