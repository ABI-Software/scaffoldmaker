"""
Utility functions for generating a central path used to generate cylinders and tubes.
"""

import math

class CentralPath:
    """
    Modifies the initial central path. Samples elements along the central path and smooth derivatives.
    """
    def __init__(self, centralPath):
        """
        :param centralPath: 1D scaffold created using meshtype_1d_path1
        """
        self._centralPath = centralPath

    def getCentralPathNodes(self, region):
        # Central path
        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx, cd1, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)
        # for i in range(len(cx)):
        #     print(i, '[', cx[i], ',', cd1[i], ',', cd2[i],',', cd12[i], '],')
        del tmpRegion
        return