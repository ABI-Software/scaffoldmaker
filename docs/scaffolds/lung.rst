Lung Scaffold
=============

The current recommended lung scaffold is ``3D Lung 2`` built from ``class MeshType_3d_Lung2``;
the human variant is shown in :numref:`fig-scaffoldmaker-human-lung`.

.. _fig-scaffoldmaker-human-lung:

.. figure:: ../_images/scaffoldmaker_human_lung.jpg
   :align: center

   Human lung scaffold.

The lung scaffold is a 3-D volumetric model of the lungs representing left lung, right lung and its lobes.
The lobes are lower lobe, middle lobe (only right lung), upper lobe and accessory lobe (only right lung), but
doesn't have representations of the pulmonary airway, blood vessels or alveoli.

.. note::

   The lung scaffold contains two (or more) independent meshes for the left lung, right lung, and their lobes.

Variants
--------

The lung scaffold is provided with parameter sets for the following four species, which differ in shape, and in particular
have different numbers of lobes:

* Human (2 lobes in the left, 3 lobes in the right lung)
* Pig (1 lobe in the left, 4 lobes in the right lung)
* Rat (1 lobe in the left, 4 lobes in the right lung)
* Mouse (1 lobe in the left, 4 lobes in the right lung)

These variants' geometry and annotations are best viewed in the **Scaffold Creator** tool in the ABI Mapping Tools.
On the web, the latest published generic lung scaffold variants can be viewed on the `SPARC Portal <https://sparc.science/>`_
by searching for ``lung``, filtering for models, selecting a variant and viewing the scaffold in its Gallery tab.

The lung scaffold script generates the scaffold mesh and geometry from half ellipsoid and triangular prism functions with
many parameters controlling the shape. The parameters were carefully tuned for each species, and it is not recommended that these be edited.

An advanced optional feature is to check *Open fissures* (set parameter to ``true``) which separates the lobes into the independent
meshes allowing *sliding elements* along the horizontal and oblique fissures.

Coordinates
-----------

The lung scaffold defines both geometric and material coordinates (field ``coordinates``) which give *modified* geometry at approximately unit
scale. Note that the scaffold has a *Unit scale* parameter (default value ``1.0``) which scales the entire scaffold efficiently.

A material coordinates field is provided, so to perform embedding at this time, it is necessary to use the generic
``material coordinates`` field, built with the standard, *unmodified* parameter set (incl. *Unit scale* ``1.0``) for the species.

The lung scaffold supports limited refinement/resampling by checking *Refine* (set parameter to ``true``) with chosen
*Refine number of elements~* parameters. Be aware that only the ``coordinates`` field is currently defined on the refined mesh
(but annotations are transferred), and the refined whole lungs is only conformant if all *Refine number of elements~* parameters have the same value.

Annotations
-----------

Important anatomical regions of the lungs are defined by groups of elements (or faces, edges and nodes/points) and
annotated with standard term names and identifiers from a controlled vocabulary.

Annotated 3-dimensional volume regions are defined by groups of 3-D elements including (using only one of the items separated by slash /):

* left/right lung
* lower/upper lobe of left/right lung
* lung
* middle lobe of right lung
* right lung accessory lobe

**Terms for volume regions such as the above are not to be used for digitized contours!** They are used for applying
different material properties in models and the strain/curvature penalty (stiffness) parameters in fitting.

Annotated 2-dimensional surface regions are defined for matching annotated contours digitized from medical images including
(where ``surface`` is the outside boundary on the meshes and only one of the items separated by slash /):

* base of left/right lung surface
* base of lower lobe of left/right lung surface
* base of middle lobe of right lung surface
* base of right lung accessory lobe surface
* base of upper lobe of left lung surface
* horizontal fissure of right lung
* horizontal fissure of middle/upper lobe of right lung
* lateral/medial surface of left/right lung
* lateral/medial surface of lower/upper lobe of left/right lung
* medial surface of middle lobe of right lung
* left/right lung surface
* lower/upper lobe of left/right lung surface
* middle lobe of right lung surface
* oblique fissure of left/right lung
* oblique fissure of lower/upper lobe of left/right lung
* oblique fissure of middle lobe of right lung
* right lung accessory lobe surface

Note that once the subgroups are used in fitting, the digitized medical images should be excluded from its supergroups.
For example, if ``lower/upper lobe of left lung`` are annotated in the digitized medical images, the
digitized medical images shouldn't be annotated in ``left lung``. This helps to avoid an overlapping effect in between
groups

Annotated 1-dimensional line regions are defined for matching annotated contours digitized from medical images including
(only one of the items separated by slash /):

* anterior border of left/right lung

Several fiducial marker points are defined on the lung scaffold, of which the followings are potentially usable when digitizing:

* apex of left/right lung
* dorsal base of right lung accessory lobe
* laterodorsal tip of middle lobe of right lung
* medial/ventral base of left/right lung
* medial/ventral base of right lung accessory lobe

At present these are defined on the outer surface of the lungs
