"""hsv1-coords: a coordinate and annotation layer for herpes simplex virus 1."""
from .resolver import (AmbiguousName, Dictionary, Entry, RepeatPairedName, normalise,
                       repeat_pair_of)

__all__ = ["AmbiguousName", "Dictionary", "Entry", "RepeatPairedName", "normalise",
           "repeat_pair_of"]
