"""
Scaffold abstract base class.
Describes method each scaffold must or may override.
"""

class Scaffold_base:
    '''
    Base class for scaffolds / mesh generator scripts.
    Not intended to be instantiated. Most methods must be overridden by actual scaffolds.
    '''

    @staticmethod
    def getName():
        '''
        Must override.
        :return: Unique type name for scaffold, for display in user interface.
        '''
        return None

    @staticmethod
    def getParameterSetNames():
        '''
        Optionally override to return additional default parameter set names supported by getDefaultOptions().
        Always have first set name 'Default'. Do not use name 'Custom' as clients may use internally.
        '''
        return ['Default']
 
    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        '''
        Must override to get valid initial default or other named parameter set. Must support 'Default' parameter set.
        :param parameterSetName: Name of parameter set to get, from list returned by getParameterSetNames().
        :return: Dictionary of parameter name value pairs, of integer, real or boolean type.
        '''
        return {
            'Integer option' : 1,
            'Real option' : 1.0,
            'Boolean option' : True
        }

    @staticmethod
    def getOrderedOptionNames():
        '''
        Must override to get list of parameters in order for display and editing in user interface.
        Note can omit parameter names to remove from interface.
        :return: List of parameter names in display order.
        '''
        return []

    @staticmethod
    def checkOptions(options):
        '''
        Must override to keep options within limits to prevent nonsense or errors.
        '''
        if options['Integer option'] < 1:
            options['Integer option'] = 1
        if options['Real option'] < 0.0:
            options['Real option'] = 0.0

    @classmethod
    def generateMesh(cls, region, options):
        '''
        Must override to generate scaffold mesh in region using Zinc API with options.
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
