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

from testutils import assertAlmostEqualList


class MuscleScaffoldTestCase(unittest.TestCase):

    def test_muscle_annotations(self):
        """
        Test nomenclature of the muscle terms. 
        """
        for term in muscle_terms:
            upper_id = term[1].upper()
            self.assertTrue(("UBERON" in upper_id) or ("ILX" in upper_id) or (upper_id == ""), "Invalid muscle term" + str(term))
            if len(term) > 2:
                uberon_index = next((i for i, v in enumerate(term) if 'UBERON' in v), -1)
                ilx_index = next((i for i, v in enumerate(term) if 'ILX' in v), -1)
                fma_index = next((i for i, v in enumerate(term) if 'FMA' in v), 10)
                self.assertTrue(fma_index > uberon_index and fma_index > ilx_index, 'FMA term should be written after UBERON and ILX terms on muscle term ' + str(term))

if __name__ == "__main__":
    unittest.main()
