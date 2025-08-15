ScaffoldMaker library
=====================

The *scaffoldmaker library* contains scripts for programmatically generating anatomical scaffolds for a range of organs and other structures, and includes utility code for building common geometric shapes, annotation, and refinement. The scripts generate the scaffold model in the data structures of the underlying *Zinc library*, and through the Zinc API clients can interrogate the model or export it for subsequent use.

Most users will make scaffolds using the ABI Mapping Tools' **Scaffold Creator** user interface for this library, and its documentation gives a good introduction to anatomical scaffolds and some common features of them.

This documentation is intended to be a resource describing details about using individual scaffolds, as well as developing new ones.

.. toctree::

   install
   scaffolds/bladder
   scaffolds/brainstem
   scaffolds/colon
   scaffolds/esophagus
   scaffolds/heart
   scaffolds/lung
   scaffolds/smallintestine
   scaffolds/stomach
   scaffolds/uterus
   scaffolds/vagus
