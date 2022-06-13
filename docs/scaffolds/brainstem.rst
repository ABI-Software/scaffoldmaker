Brainstem Scaffold
==================

The current recommended brainstem scaffold is ``3D brainstem 1`` built from ``class MeshType_3d_brainstem1``;
the human variant is shown in :numref:`fig-scaffoldmaker-human-brainstem`.

.. _fig-scaffoldmaker-human-brainstem:

.. figure:: ../_images/scaffoldmaker_human_brainstem.jpg
   :align: center

   Human brainstem scaffold.

The brainstem scaffold is a 3-D volumetric model of the brainstem representing all three parts: midbrain, pons, and
medulla oblongata.

.. note::

  The parts of the brainstem are connected to each other and they are inseparable.

Variants
--------

The brainstem scaffold is provided with parameter sets for the following six species, which differ in shape and size.

* Cat
* Human
* Mouse
* Pig
* Rat
* Sheep

These variants' geometry and annotations are best viewed in the **Scaffold Creator** tool in the ABI Mapping Tools.
On the web, the latest published generic brainstem scaffold variants can be viewed on the
`SPARC Portal <https://sparc.science/>`_ by searching for ``brainstem``, filtering for models, selecting a variant and
viewing the scaffold in its Gallery tab.

The brainstem scaffold script generates the scaffold mesh and geometry from a cylinder function with many parameters
controlling the shape and size. The parameters were carefully tuned for each species, and it is not recommended that
these be edited.

To customise the shape and size, Clicking *'Central path'* navigates a user from a 3D volumetric model to a 1D line model,
representing the 1D central path of the brainstem scaffold. If you would like to change its shape, holding *S* on your keyboard
and left-clicking nodes on the GUI to move the nodes (Make sure *Node points* in Display is checked). Similarly, if you would
like to customise the derivatives, left-clicking the nodes and their derivatives (Make sure *Node derivatives* in Display is checked).
After finishing 1D path customisation, Clicking *'<< Back'* in Central path will navigate back to a 3D volumetric model.
Additionally, Displaying ``Annotations`` is recommended along with shape and size adjustment.

Coordinates
-----------

The brainstem scaffold defines both geometric and material coordinates (field ``coordinates``) which give *modified* geometry at approximately
unit scale. Note that the scaffold has a *Unit scale* parameter (default value ``1.0``) which scales the entire
scaffold efficiently.

A material coordinates field is provided, so to perform embedding at this time, it is necessary to use the
``material coordinates`` field, built with the standard, *unmodified* parameter set (incl. *Unit scale* ``1.0``) for the species.

The brainstem scaffold supports limited refinement/resampling by checking *Refine* (set parameter to ``true``) with
chosen *Refine number of elements~* parameters. Be aware that only the ``coordinates`` field is currently defined
on the refined mesh (but annotations are transferred), and the refined brainstem is only conformant if all *Refine
number of elements~* parameters have the same value.

Annotations
-----------

Important anatomical regions of the brainstem are defined by groups of elements (or faces, edges and nodes/points) and
annotated with standard term names and identifiers from a controlled vocabulary.

Annotated 3-dimensional volume regions are defined by groups of 3-D elements:

* brainstem
* medulla oblongata
* midbrain
* pons

**Terms for volume regions such as the above are not to be used for digitized contours!** They are used for applying
different material properties in models and the strain/curvature penalty (stiffness) parameters in fitting.

Annotated 2-dimensional surface regions are defined for matching annotated contours digitized from medical images
including (where ``interface`` means the outside boundary of the brainstem connecting to different organs such as
thalamus and spinal cord, ``exterior`` is the outside boundary of the brainstem):

* brainstem exterior
* brainstem-spinal cord interface
* medulla oblongata exterior
* midbrain exterior
* pons exterior
* thalamus-brainstem interface

Annotated 1-dimensional line regions are defined for displaying a 3D annotation groups on the 1D line model. This
helps a user during shape and size customisation.

* medulla oblongata
* midbrain
* pons

Several fiducial marker points are defined on the brainstem scaffold, of which the following eight are potentially
usable when digitizing:

* brainstem dorsal midline caudal point
* brainstem ventral midline caudal point
* brainstem dorsal midline pons-medulla junction
* brainstem ventral midline pons-medulla junction
* brainstem dorsal midline midbrain-pons junction
* brainstem ventral midline midbrain-pons junction
* brainstem dorsal midline cranial point
* brainstem ventral midline cranial point

At present these are defined on the outer surface of the brainstem in each section of the brainstem.
