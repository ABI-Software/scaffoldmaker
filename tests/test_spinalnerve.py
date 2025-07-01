import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.spinal_nerve_terms import get_spinal_nerve_term, spinal_nerve_terms
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace

from testutils import assertAlmostEqualList


class spinal_nerveScaffoldTestCase(unittest.TestCase):

    def test_spinal_nerve_annotations(self):
        """
        Test that all spinal_nerve terms are UBERON or ILX.
        """
        for term in spinal_nerve_terms:
            upper_id = term[1].upper()
            self.assertTrue(("UBERON" in upper_id) or ("ILX" in upper_id) or (upper_id == ""), "Invalid heart term" + str(term))

if __name__ == "__main__":
    unittest.main()
