"""
Abstract base class for objects returned as construction objects by scaffold scripts.
"""
from abc import ABC, abstractmethod


class ConstructionObject:
    """
    Abstract base class for objects returned as construction objects by scaffold scripts.
    Mainly presents method for getting metadata dict.
    """

    @abstractmethod
    def getMetadata(self) -> dict:
        """
        Override to get scaffold-specific metadata.
        """
        return {}
