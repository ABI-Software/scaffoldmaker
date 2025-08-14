import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.brainstem_terms import get_brainstem_term, brainstem_terms
from scaffoldmaker.meshtypes.meshtype_3d_brainstem import MeshType_3d_brainstem1
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace

from testutils import assertAlmostEqualList


class BrainstemScaffoldTestCase(unittest.TestCase):

    def test_brainstem_annotations(self):
        """
        Test nomenclature of the brainstem terms.  
        """
        for term in brainstem_terms:
            upper_id = term[1].upper()
            self.assertTrue(("UBERON" in upper_id) or ("ILX" in upper_id) or (upper_id == ""), "Invalid brainstem term" + str(term))
            if len(term) > 2:
                uberon_index = next((i for i, v in enumerate(term) if 'UBERON' in v), -1)
                ilx_index = next((i for i, v in enumerate(term) if 'ILX' in v), -1)
                fma_index = next((i for i, v in enumerate(term) if 'FMA' in v), 10)
                self.assertTrue(fma_index > uberon_index and fma_index > ilx_index, 'FMA term should be written after UBERON and ILX terms on brainstem term ' + str(term))

if __name__ == "__main__":
    unittest.main()
