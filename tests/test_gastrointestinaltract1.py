import unittest
from testutils import assertAlmostEqualList

from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.cecum_terms import get_cecum_term
from scaffoldmaker.annotation.colon_terms import get_colon_term
from scaffoldmaker.annotation.esophagus_terms import get_esophagus_term
from scaffoldmaker.annotation.stomach_terms import get_stomach_term
from scaffoldmaker.annotation.smallintestine_terms import get_smallintestine_term
from scaffoldmaker.meshtypes.meshtype_3d_gastrointestinaltract1 import MeshType_3d_gastrointestinaltract1
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace


class GastrointestinalTractScaffoldTestCase(unittest.TestCase):

    def test_gastrointestinaltract1(self):
        """
        Test creation of gastrointestinal tract scaffold.
        """
        scaffold = MeshType_3d_gastrointestinaltract1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ['Default', 'Human 1'])
        options = scaffold.getDefaultOptions("Default")
        self.assertEqual(10, len(options))
        networkLayoutOptions = options['Network layout']
        networkLayoutSettings = networkLayoutOptions.getScaffoldSettings()
        self.assertEqual("1-2-3-4-5-6-7.2, 8-9-10-11-7-12-13-14-15-16-17-18-19-20-21-22-23-24-25-26-27-28-29-30-31-32-"
                         "33-34-35-36-37-38-39-40-41-42-43-44-45-46-47-48-49-50-51-52-53-54-55-56-57-58-59-60-61-62-63-"
                         "64-65-66-67-68-69-70-71-72-73-74-75-76-77-78-79-80-81-82-83-84-85-86-87-88-89-90-91-92-93-94-"
                         "95-96-97-98-99-100-101-102-103-104-105-106-107-108-109-110-111-112-113-114-115-116-117-118-"
                         "119-120-121-122-123-124-125-126-127-128-129-130-131-132-133-134-135-136-137-138-139-140-141-"
                         "142-143-144-145-146-147-148-149-150-151-152-153-154-155-156-157-158-159-160-161-162-163-164-"
                         "165-166-167-168-169-170-171-172.2, 173-172-174-175-176-177-178-179-180-181-182-183-184-185-"
                         "186-187-188-189-190-191-192-193-194-195-196-197-198-199-200-201-202-203-204-205-206-207-208-"
                         "209-210-211-212-213-214-215-216-217-218-219-220-221-222",
                         networkLayoutSettings.get("Structure"))
        esoOptions = options['Esophagus']
        esoSettings = esoOptions.getScaffoldSettings()
        self.assertEqual(8, esoSettings['Number of elements around'])
        self.assertEqual(1, esoSettings['Number of elements through wall'])
        stomachOptions = options['Stomach']
        stomachSettings = stomachOptions.getScaffoldSettings()
        self.assertEqual(8, stomachSettings['Number of elements around esophagus'])
        self.assertEqual(12, stomachSettings['Number of elements around duodenum'])
        self.assertEqual(14, stomachSettings['Number of elements along'])
        self.assertEqual(1, stomachSettings['Number of elements through wall'])
        self.assertEqual(3.0, stomachSettings['Wall thickness'])
        smallIntestineOptions = options['Small intestine']
        smallIntestineSettings = smallIntestineOptions.getScaffoldSettings()
        self.assertEqual(80, smallIntestineSettings['Number of segments'])
        self.assertEqual(12, smallIntestineSettings['Number of elements around'])
        self.assertEqual(3, smallIntestineSettings['Number of elements along segment'])
        self.assertEqual(1, smallIntestineSettings['Number of elements through wall'])
        self.assertEqual(3, smallIntestineSettings['Wall thickness'])
        cecumOptions = options['Cecum']
        cecumSettings = cecumOptions.getScaffoldSettings()
        self.assertEqual(1, cecumSettings['Number of segments'])
        self.assertEqual(8, cecumSettings['Number of elements around haustrum'])
        self.assertEqual(12, cecumSettings['Number of elements along segment'])
        self.assertEqual(1, cecumSettings['Number of elements through wall'])
        self.assertEqual(1.6, cecumSettings['Wall thickness'])
        colonOptions = options['Colon']
        colonSettings = colonOptions.getScaffoldSettings()
        self.assertEqual(33, colonSettings['Number of segments'])
        self.assertEqual(True, options['Use linear through wall'])
        self.assertEqual(False, options['Refine'])
        self.assertEqual(4, options['Refine number of elements surface'])
        self.assertEqual(1, options['Refine number of elements through wall'])

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateMesh(region, options)[0]
        self.assertEqual(30, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(8448, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(32494, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(39642, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(15600, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check coordinates range
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-123.4289535529459, -190.90949289431342, 799.4386433354306], 1.0E-6)
        assertAlmostEqualList(self, maximums, [124.21567892707846, -33.094490024089325, 1403.9055789269714], 1.0E-6)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            faceMeshGroup = createFaceMeshGroupExteriorOnFace(fieldmodule, Element.FACE_TYPE_XI3_1)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, faceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 280350.6462177311, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 598505.4701647585, delta=1.0E-3)

        # check some annotationGroups:
        expectedSizes3d = {
            "esophagus": 184,
            "stomach": 168,
            "small intestine": 2916,
            "caecum": 460,
            "colon": 4752
        }
        for name in expectedSizes3d:
            term = None
            if name == "esophagus":
                term = get_esophagus_term(name)
            elif name == "stomach":
                term = get_stomach_term(name)
            elif name == "small intestine":
                term = get_smallintestine_term(name)
            elif name == "caecum":
                term = get_cecum_term(name)
            elif name == "colon":
                term = get_colon_term(name)
            group = getAnnotationGroupForTerm(annotationGroups, term)
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)

        # refine 8x8x8 and check result
        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine number of elements surface'] = 2
        options['Refine number of elements through wall'] = 1
        meshrefinement = MeshRefinement(region, refineRegion, [])
        scaffold.refineMesh(meshrefinement, options)

        refineFieldmodule.defineAllFaces()
        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(33792, mesh3d.getSize())


if __name__ == "__main__":
    unittest.main()
