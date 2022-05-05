"""
Generates a 3-D unit box mesh with variable numbers of elements in 3 directions.
"""

from __future__ import division

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldFiniteElement
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from opencmiss.zinc.element import Element, Elementbasis


class MeshType_3d_vasculature1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Vasculature 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements 1' : 1,
            'Number of elements 2' : 1,
            'Number of elements 3' : 1,
            'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements 1' : 1,
            'Refine number of elements 2' : 1,
            'Refine number of elements 3' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements 1',
            'Number of elements 2',
            'Number of elements 3',
            'Use cross derivatives',
            'Refine',
            'Refine number of elements 1',
            'Refine number of elements 2',
            'Refine number of elements 3'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements 1',
            'Number of elements 2',
            'Number of elements 3',
            'Refine number of elements 1',
            'Refine number of elements 2',
            'Refine number of elements 3']:
            if options[key] < 1:
                options[key] = 1

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """
        # elementsCount1 = options['Number of elements 1']
        # elementsCount2 = options['Number of elements 2']
        # elementsCount3 = options['Number of elements 3']
        # useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        radius = findOrCreateFieldFiniteElement(fm, "radius", components_count=1, managed=True)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.defineField(radius)
        nodetemplate.setValueNumberOfVersions(radius, -1, Node.VALUE_LABEL_VALUE, 1)

        nodetemplate1 = nodes.createNodetemplate()
        nodetemplate1.defineField(coordinates)
        nodetemplate1.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate1.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 3)
        nodetemplate1.defineField(radius)
        nodetemplate1.setValueNumberOfVersions(radius, -1, Node.VALUE_LABEL_VALUE, 1)
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        # if useCrossDerivatives:
        #     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        # if useCrossDerivatives:
        #     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
        #     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
        #     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        mesh = fm.findMeshByDimension(1)
        cubicHermiteBasis = fm.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        linearBasis = fm.createElementbasis(1, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
        eft1.setTermNodeParameter(2, 1, 1, Node.VALUE_LABEL_D_DS1, 1)
        eft2 = mesh.createElementfieldtemplate(cubicHermiteBasis)
        eft2.setTermNodeParameter(2, 1, 1, Node.VALUE_LABEL_D_DS1, 2)
        eft3 = mesh.createElementfieldtemplate(cubicHermiteBasis)
        eft3.setTermNodeParameter(2, 1, 1, Node.VALUE_LABEL_D_DS1, 3)
        eft32 = mesh.createElementfieldtemplate(cubicHermiteBasis)
        eft32.setTermNodeParameter(2, 1, 2, Node.VALUE_LABEL_D_DS1, 3)
        eftRadius = mesh.createElementfieldtemplate(linearBasis)

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate1.defineField(coordinates, -1, eft1)
        elementtemplate1.defineField(radius, -1, eftRadius)

        elementtemplate2 = mesh.createElementtemplate()
        elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate2.defineField(coordinates, -1, eft2)
        elementtemplate3 = mesh.createElementtemplate()
        elementtemplate3.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate3.defineField(coordinates, -1, eft3)
        elementtemplate2.defineField(radius, -1, eftRadius)
        elementtemplate3.defineField(radius, -1, eftRadius)



        cache = fm.createFieldcache()

        vasculature = False
        vagus = True
        if vasculature:
            x_list = [
            [20.9901, -55.2915, 1237.12],
            [21.7547, -56.1276, 1268.6],
            [24.838, -74.4412, 1299.32],
            [21.3563, -91.2, 1307.97],
            [12.9058, -111.015, 1302.16],
            [25.4696, -95.415, 1344.72],
            [39.5394, -109.4, 1435.59],
            [-10.5303, -102.14, 1331.9],
            [-27.2225, -91.4847, 1430.7],
            [-27.5775, -97.4023, 1345.7],
            [5.157545857793905, -73.7144264439337, 1067.611422937166],
            [0.7123922184414442, -76.12139177722604, 949.3444291866757],
            [1.480653991761836, -73.62566826677576, 815.1125925698087],
            [-19.50824359035452, -72.86872429654028, 929.2051634762531],
            [-20.28614848561935, -73.51853769979229, 899.0024026991548],
            [-21.61102738650197, -74.8636341463615, 841.6568853077786],
            [-18.32152657377993, -80.08693448724975, 975.6309361438332],
            [82.33929143728056, -75.38569430920019, 965.0794458269447],
            [-99.29494952824054, -84.197838408155, 963.3384151864503],
            [-14.63217490185608, -80.76093647541154, 1026.190525855144],
            [-100.6218523253879, -81.37468570139932, 1070.890294941285],
            [17.07417461634516, -120.3755886862945, 1249.284411678494]
                ]

            d1_list = [
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
            ]

            nodes_w_ver = [5,8,4,3,11,20,17]
        elif vagus:
            x_list = [
            [-3.067919428038835, -58.88831896656854, 1485.94100956215],
            [-0.6281586078592581, -62.69296504830071, 1481.347190265109],
            [19.14045007878925, -71.17249683579628, 1471.653806089609],
            [35.50375246183219, -78.35854544883568, 1463.384930602326],
            [36.38036257645483, -80.58804084795787, 1457.154198693773],
            [35.70333609041023, -82.09528266416706, 1448.653619505636],
            [34.97165431464403, -83.39378105979046, 1441.867216455742],
            [34.26489569026836, -84.45087019096185, 1435.143450297613],
            [32.29029472400944, -92.73679568270747, 1396.092960178755],
            [28.92067653338183, -109.6287737218008, 1316.75069003233],
            [28.32112681266167, -112.2390790222721, 1305.054388871628],
            [6.058014106179172, -96.19324278756504, 1168.221198705996],
            [4.49003139871628, -120.3327063642019, 982.41199971608],
            [-4.484831774677777, -76.72363699970431, 919.0809110751771],
            [-9.096316816507812, -18.3837646869268, 858.4382217409927]
                ]

            d1_list = [
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                ]

            nodes_w_ver = []
        # create nodes
        nodeIdentifier = 1
        for n3 in range(len(x_list)):
            if n3+1 in nodes_w_ver:
                nt = nodetemplate1
            else:
                nt = nodetemplate
            node = nodes.createNode(nodeIdentifier, nt)
            cache.setNode(node)

            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x_list[n3])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1_list[n3])
            if n3+1 in nodes_w_ver:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, d1_list[n3])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 3, d1_list[n3])

            radius.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 6)
            nodeIdentifier = nodeIdentifier + 1

        if vasculature:
            elem_list = [
                [2, 1],
                [3, 2],
                [4, 3],
                [5, 4],
                [3, 6],
                [4, 7],
                [5, 8],
                [8, 9],
                [8, 10],
                [1, 11],
                [11, 20],
                [20, 21],
                [20, 17],
                [17, 18],
                [17, 19],
                [11, 12],
                [12, 13],
                [17, 14],
                [14, 15],
                [15, 16],
                [22, 5],
            ]
        elif vagus:
            elem_list = [
                [1, 2],
                [2, 3],
                [3, 4],
                [4, 5],
                [5, 6],
                [6, 7],
                [7, 8],
                [8, 9],
                [9, 10],
                [10, 11],
                [11, 12],
                [12, 13],
                [13, 14],
                [14, 15],
                ]

        # create elements
        elementIdentifier = 1
        for e3 in range(len(elem_list)):
            if vasculature:
                if elementIdentifier in [7, 9, 6, 5, 16, 13, 14]:
                    elementtemplate = elementtemplate2
                    eft = eft2
                elif elementIdentifier in [8, 2, 11, 12, 15]:
                    elementtemplate = elementtemplate3
                    eft = eft3
                else:
                    elementtemplate = elementtemplate1
                    eft = eft1
            elif vagus:
                elementtemplate = elementtemplate1
                eft = eft1
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result = element.setNodesByIdentifier(eft, elem_list[e3])
            result = element.setNodesByIdentifier(eftRadius, elem_list[e3])
            elementIdentifier = elementIdentifier + 1




        fm.endChange()
        return []

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        refineElementsCount1 = options['Refine number of elements 1']
        refineElementsCount2 = options['Refine number of elements 2']
        refineElementsCount3 = options['Refine number of elements 3']
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCount1, refineElementsCount2, refineElementsCount3)
