"""
Scaffold abstract base class.
Describes method each scaffold must or may override.
"""

class Scaffold_base:
    '''
    Base class for all 
    '''

    @staticmethod
    def getName():
        return 'Unique type name for display in user interface'

    @staticmethod
    def getParameterSetNames():
        '''
        Override to return additional default parameter set names supported by getDefaultOptions().
        Always have first set name 'Default'. Do not use name 'Custom' as clients may use internally.
        '''
        return ['Default']
 
    @staticmethod
    def getDefaultOptions(parameterSetName):
        return {
            'Integer option' : 1,
            'Real option' : 1.0,
            'Boolean option' : True
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'List of parameters',
            'In order for editing in user interface'
        ]

    @staticmethod
    def checkOptions(options):
        '''
        Override to keep options within limits to prevent nonsense or errors.
        '''
        if options['Integer option'] < 1:
            options['Integer option'] = 1
        if options['Real option'] < 0.0:
            options['Real option'] = 0.0

    @classmethod
    def generateMesh(cls, region, options):
        '''
        Override to generate scaffold mesh in region using Zinc API with options.
        Scaffolds supporting refinement may switch to call other functions:
        @classmethod
        def generateBaseMesh(cls, region, options):
        # generates base high order scaffold
        @classmethod
        def refineMesh(cls, meshrefinement, options):
        # performs refinement
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: List of AnnotationGroup or None if not supported.
        '''
        return None
