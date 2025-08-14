import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.muscle_terms import get_muscle_term, muscle_terms
from scaffoldmaker.meshtypes.meshtype_3d_musclefusiform1 import MeshType_3d_musclefusiform1
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace

from testutils import check_annotation_term_ids


class MuscleScaffoldTestCase(unittest.TestCase):

    def test_muscle_annotations(self):
        """
        Test nomenclature of the muscle terms. 
        """
        for term_ids in muscle_terms:
            self.assertTrue(check_annotation_term_ids(term_ids), "Invalid primary term id or order not UBERON < ILX < FMA for muscle annotation term ids " + str(term_ids)) 

if __name__ == "__main__":
    unittest.main()
