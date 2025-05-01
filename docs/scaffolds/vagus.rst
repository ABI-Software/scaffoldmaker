Vagus Scaffold
================

The current subject-specific vagus scaffold is ``3D Nerve 1`` built from ``class MeshType_3d_nerve1``.
The vagus scaffold needs an input file with segmentations to built from it.

Input file requirements
-----------------------

The input file should have annotated segmentations for the vagus trunk and the vagal branches.
Vagus trunk should be annotated as either 'left vagus nerve' or 'right vagus nerve'.
All vagal branches should all have distinctive names, e.g. there should not be several branches with the same name.

The input file should also have included some of the anatomical landmarks to be able to locate the vagus in the body
and to estimate vagus nerve real length. The current minimum is two included landmarks as marker datapoints.
The list of approved landmarks is as follows.

* left/right level of superior border of jugular foramen on the vagus nerve
* left/right level of inferior border of jugular foramen on the vagus nerve
* left/right level of inferior border of jugular foramen on the vagus nerve
* left/right level of inferior border of cranium on the vagus nerve
* left/right level of C1 transverse process on the vagus nerve
* left/right level of angle of the mandible on the vagus nerve
* left/right level of greater horn of hyoid on the vagus nerve
* left/right level of carotid bifurcation on the vagus nerve
* left/right level of laryngeal prominence on the vagus nerve
* left/right level of superior border of the clavicle on the vagus nerve
* left/right level of jugular notch on the vagus nerve
* left/right level of sternal angle on the vagus nerve
* left/right level of 1 cm superior to start of esophageal plexus on the vagus nerve
* left/right level of esophageal hiatus on the vagus nerve
* left/right level of aortic hiatus on the vagus nerve

Variants
--------

The vagus scaffold is provided with the following parameter sets:

* Human left vagus
* Human right vagus

Coordinates
-----------

The vagus box scaffold defines the geometric and material coordinates.

The geometric ``coordinates`` field is built using the coordinates supplied by the input file.

The material ``coordinates`` field defines a coordinate system to give permanent locations for landmarks on the vagus.
It is defined as a straight line as vagus trunk, with origin as the top of the vagus trunk and 1 as the bottom end
of the vagus trunk. The branches are defined as lines with the same orientation as when they are leaving the trunk.


Annotations
-----------

Important anatomical features of the vagus nerve are defined by groups of elements (or faces, edges and nodes/points) and
annotated with standard term names and identifiers from a controlled vocabulary.

Each subject-specific vagus box scaffold has the following defined groups of 3-D elements:

* vagus centroid
* vagus epineureum
* vagus anterior line
