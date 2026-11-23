"""Finite-field conic billiards: geometric group law and elliptic-curve correspondence (research code, not a cryptographic library)."""
from .core import Conic, UndefinedReflection
from .search import CycleStats, analyze

__all__ = ['Conic', 'UndefinedReflection', 'CycleStats', 'analyze']
